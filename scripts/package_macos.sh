#!/bin/sh
# Build a signed macOS installer package containing the svht_denoise executable.
#
#   package_macos.sh <version>
#
# Environment:
#   MACOS_SIGN_IDENTITY       "Developer ID Application: ..." or its SHA-1, or
#                             "-" for an ad-hoc signature. Signs the executable.
#   MACOS_INSTALLER_IDENTITY  "Developer ID Installer: ..." or its SHA-1. Signs
#                             the .pkg itself. This is a DIFFERENT certificate
#                             from the Application one, and a package signed
#                             with the Application certificate is rejected by
#                             the installer. Leave unset for an unsigned
#                             package: local testing only, cannot be notarised.
#   MACOSX_DEPLOYMENT_TARGET  oldest macOS supported (default 11.0).
#
# Releases are arm64 only. This is a compute-bound tool and the Intel Macs a
# universal build would have served are the slowest machines that could run it,
# on a platform macOS 26 is the last release to support. Intel users build from
# source with a plain `make`, which is unaffected and needs nothing from here.
#
# The package installs to /usr/local/bin, which is on the default PATH via
# /etc/paths, so the tool works immediately without the user editing a shell
# profile. That is the whole reason to prefer a package over a disk image for a
# command-line tool.
#
# Adapted from the same three-script arrangement in brainchopC (same author).
# Simpler here: one executable, no embedded weights, no OpenMP runtime to build
# and stage, and nothing linked beyond libSystem and libz -- both of which ship
# with macOS. So there is no lib/ directory, no install_name_tool rewriting and
# no second signature.
set -eu

if [ "$(uname -s)" != Darwin ]; then
	echo "package_macos.sh requires macOS" >&2
	exit 2
fi
# Releasing needs an Apple Silicon host, because this script RUNS what it builds:
# the test suite and the -version preflight both execute the arm64 binary, which
# an Intel Mac cannot do.  It already failed safely there -- non-zero, nothing
# signed -- but only after a full `make clean` and rebuild, ending in a bare
# "Bad CPU type in executable" that never names the cause.  Say it up front.
if [ "$(uname -m)" != arm64 ]; then
	echo "package_macos.sh: releases are arm64-only and this is $(uname -m)." >&2
	echo "  The script runs the test suite against the arm64 binary it builds," >&2
	echo "  which this machine cannot execute. Build the release on Apple Silicon." >&2
	echo "  To just build svht_denoise for THIS machine: make -C src" >&2
	exit 2
fi

sv_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
sv_src="$sv_root/src"
sv_version=${1:?usage: package_macos.sh VERSION}
sv_identity=${MACOS_SIGN_IDENTITY--}
sv_installer_identity=${MACOS_INSTALLER_IDENTITY:-}
sv_target=${MACOSX_DEPLOYMENT_TARGET:-11.0}
sv_arch=arm64
sv_identifier=io.github.rordenlab.svht-denoise
sv_dist="$sv_root/dist"
sv_pkg="$sv_dist/svht_denoise-$sv_version-macos-$sv_arch.pkg"

case "$sv_version" in
	*[!A-Za-z0-9._-]*|'')
		echo "package_macos.sh: VERSION may contain only letters, digits, dot, underscore, and hyphen" >&2
		exit 2
		;;
esac
case "$sv_target" in
	*[!0-9.]*|'')
		echo "package_macos.sh: deployment target must contain only digits and dots" >&2
		exit 2
		;;
esac

for sv_tool in codesign pkgbuild productbuild pkgutil otool lipo xattr security; do
	command -v "$sv_tool" >/dev/null 2>&1 || {
		echo "package_macos.sh: missing required tool: $sv_tool" >&2
		exit 2
	}
done

# A release must be THE PINNED zlib-ng build, and each of the checks below
# closes a different way of not being that.
#
# 1. The environment must not redirect the build.  `make` reads ZLIBNG_ROOT via
# `ifdef`, so it is picked up from the environment, and this script does not
# sanitize what it passes to make: with ZLIBNG_ROOT exported, a package built
# from an arbitrary libz.a outside the repo passed BOTH gates below -- the
# submodule was never even built.  The no-dynamic-libz gate cannot see it,
# because any STATIC libz.a removes the dynamic dependency just as well.
for sv_var in ZLIBNG ZLIBNG_ROOT; do
	eval "sv_val=\${$sv_var-}"
	if [ -n "$sv_val" ]; then
		echo "package_macos.sh: $sv_var is set in the environment ('$sv_val')." >&2
		echo "  A release is built from the pinned third_party/zlib-ng submodule" >&2
		echo "  and nothing else. Unset it and retry." >&2
		exit 2
	fi
done

# 2. The submodule must be present.  Otherwise `make` falls back to the system
# zlib -- right for a contributor's working copy, wrong for a package, and
# silent either way, since the fallback is ~3x slower and otherwise identical.
if [ ! -f "$sv_root/third_party/zlib-ng/configure" ]; then
	echo "package_macos.sh: the zlib-ng submodule is not checked out." >&2
	echo "  Releases link it statically; without it the build would silently" >&2
	echo "  fall back to the system zlib and ship a slower binary." >&2
	echo "  Fix with (from the repository root):" >&2
	echo "    git submodule update --init third_party/zlib-ng" >&2
	exit 2
fi

# 3. The submodule must be AT the pinned revision, with no local edits.  `git
# submodule status` prefixes '+' for "checked out at a different commit than the
# superproject records" and 'U' for conflicts; a leading space means clean.  The
# Makefile's revision stamp tracks HEAD only, so an edited working tree rebuilds
# to something that is not the pinned source and says nothing -- verified by
# appending to deflate.c, which changed the shipped library silently.
sv_sub=$(cd "$sv_root" && git submodule status third_party/zlib-ng 2>/dev/null || true)
case "$sv_sub" in
	' '*) : ;;
	'')
		echo "package_macos.sh: cannot determine the zlib-ng submodule revision." >&2
		echo "  Releases must be built from a git checkout, not an exported tree." >&2
		exit 2 ;;
	*)
		echo "package_macos.sh: the zlib-ng submodule is not at its pinned revision:" >&2
		echo "    $sv_sub" >&2
		echo "  Restore it with:  git submodule update --init --force third_party/zlib-ng" >&2
		exit 2 ;;
esac
# A DIRTY working tree is a separate question, and `git submodule status` does
# not answer it: it reports the checked-out commit, which an edited file does
# not change, so the check above passes on a tampered tree.  Measured -- an
# appended line in deflate.c built and packaged silently.
sv_dirty=$(git -C "$sv_root/third_party/zlib-ng" status --porcelain 2>/dev/null || true)
if [ -n "$sv_dirty" ]; then
	echo "package_macos.sh: the zlib-ng submodule has local modifications:" >&2
	echo "$sv_dirty" | sed 's/^/    /' >&2
	echo "  The shipped compression library would not be the pinned source." >&2
	echo "  Discard them with:  git -C third_party/zlib-ng checkout -- ." >&2
	exit 2
fi

# `${VAR--}` not `${VAR:--}`: the colon form maps an EMPTY value to "-" as well
# as an unset one, so `MACOS_SIGN_IDENTITY= make macos-pkg` silently produced an
# ad-hoc executable -- no hardened runtime, no timestamp -- inside a properly
# Developer-ID-signed package, with nothing printed.  An empty value is a
# mistake, not a request for ad-hoc; say so.
if [ -z "$sv_identity" ]; then
	echo "package_macos.sh: MACOS_SIGN_IDENTITY is set but empty." >&2
	echo "  Pass a Developer ID Application identity, or \"-\" for an ad-hoc signature." >&2
	exit 2
fi
if [ "$sv_identity" != - ]; then
	security find-identity -v -p codesigning | grep -F -- "$sv_identity" >/dev/null 2>&1 || {
		echo "package_macos.sh: signing identity not found in the current Keychain: $sv_identity" >&2
		exit 2
	}
fi
# Installer certificates are not codesigning identities, so they are absent from
# `find-identity -p codesigning`. Look at the full list instead; getting this
# wrong is easy and productbuild's own error is unhelpful.
if [ -n "$sv_installer_identity" ]; then
	security find-identity -v | grep -F -- "$sv_installer_identity" >/dev/null 2>&1 || {
		echo "package_macos.sh: installer identity not found in the current Keychain: $sv_installer_identity" >&2
		echo "  This must be a 'Developer ID Installer' certificate, which is separate" >&2
		echo "  from the 'Developer ID Application' certificate used above." >&2
		exit 2
	}
fi

sv_opt="-O3 -arch $sv_arch -mmacosx-version-min=$sv_target"

# Build AND run the test suite in ONE make invocation, so the binary that gets
# packaged is literally the one the tests passed against. Two invocations would
# not do that: the flag stamp would no longer match on the second, everything
# would be rebuilt at the default OPT, and -mmacosx-version-min would vanish --
# yielding a binary that targets the build machine's macOS rather than the 11.0
# this package advertises. (Same trap the ubsan-test targets exist to avoid --
# see the Makefile.) The minos check below would catch that particular case, but
# it is the only one it would catch; the invariant worth keeping is that the
# binary packaged is the binary tested.
# Clear the previous artifact BEFORE anything can fail.  From here to the
# productbuild at the end, every step can exit non-zero, and a stale package
# left at this exact path is one `make macos-notarize` would ship on faith.
mkdir -p "$sv_dist"
rm -f "$sv_pkg"

make -C "$sv_src" clean >/dev/null
env MACOSX_DEPLOYMENT_TARGET="$sv_target" make -C "$sv_src" OPT="$sv_opt" test

sv_exe=svht_denoise
[ -f "$sv_src/$sv_exe" ] || { echo "package_macos.sh: $sv_exe was not built" >&2; exit 1; }

# VERSION names the package and goes into its receipt; DN_VERSION is compiled
# into the binary.  Nothing links them, so `make macos-release VERSION=9.9.9`
# would happily publish a 9.9.9 package containing a binary that reports some
# entirely different version.
# Ask the binary rather than trusting the variable, and fail BEFORE anything is
# signed or submitted.
sv_reported=$("$sv_src/$sv_exe" -version 2>/dev/null | awk 'NR == 1 { print $2 }')
if [ "$sv_reported" != "$sv_version" ]; then
	echo "package_macos.sh: version mismatch." >&2
	echo "  requested VERSION=$sv_version but the binary reports '${sv_reported:-nothing}'." >&2
	echo "  Update DN_VERSION in src/svht_denoise.c to match, or pass the version it has." >&2
	exit 1
fi

# A non-system dependency would mean the installed tool breaks on a machine that
# lacks it.  verify_macos_pkg.sh checks this too, and deliberately: that one
# inspects what is INSIDE the finished package, this one the binary about to go
# in, so a failure is reported before anything is signed rather than after.
# awk on the leading tab, not a name filter: otool indents dependency lines and
# leaves the header naming the file itself flush left.  A `grep -v "$sv_exe"`
# dropped that header but also any real dependency whose path happened to
# contain the program name.
sv_foreign=$(otool -L "$sv_src/$sv_exe" | awk '/^\t/ { print $1 }' |
	grep -v '^/usr/lib/' | grep -v '^/System/' || true)
if [ -n "$sv_foreign" ]; then
	echo "package_macos.sh: $sv_exe has non-system dependencies:" >&2
	echo "$sv_foreign" | sed 's/^/  /' >&2
	exit 1
fi
# The check above cannot see this one: libz.1.dylib lives in /usr/lib, so a
# system-zlib build passes it cleanly.  Linking zlib-ng statically REMOVES that
# dependency, which makes its absence a reliable proof that the fast library is
# the one actually in the binary -- stronger than the submodule preflight, which
# only shows the source was available.
if otool -L "$sv_src/$sv_exe" | awk '/^\t/ { print $1 }' | grep -q 'libz\.'; then
	echo "package_macos.sh: $sv_exe links libz dynamically, so it was built" >&2
	echo "  against the system zlib rather than the static zlib-ng." >&2
	echo "  The build says which library it used; re-run 'make -C src' and read" >&2
	echo "  that line. A build path containing an unsafe character also falls back." >&2
	exit 1
fi
# ...and the absence of a dynamic libz only proves SOME static zlib went in.
# This proves it was the submodule's: the archive the Makefile builds is right
# here, and `make clean` above means it was produced by this run.
if [ ! -f "$sv_src/.zlibng/libz.a" ]; then
	echo "package_macos.sh: the bundled zlib-ng archive was never built, so the" >&2
	echo "  statically linked zlib came from somewhere else." >&2
	exit 1
fi
# The binary must claim the deployment target we advertise.  Collected and
# compared as one whole string rather than looped over, and that is the part to
# keep: a loop only inspects the lines it finds, so NO minos load command at all
# would have passed silently.  sort -u collapses the (now single) slice and
# would still catch two that disagreed.
# Both spellings are matched.  LC_BUILD_VERSION prints "minos" and is what an
# arm64 binary always emits, so the LC_VERSION_MIN_MACOSX arm is unreachable
# today -- it is the spelling a pre-10.14 slice used.  Left in place because it
# costs one line and fails closed; matching only "minos" once yielded zero lines
# for such a slice, which a loop read as success.
# Scoped to the load command it belongs to: a bare `$1 == "version"` also
# matches LC_SOURCE_VERSION and the linker tool version inside LC_BUILD_VERSION,
# which made this reject a perfectly good binary with "[0.0 11.0 1267.0]".
sv_minos=$(otool -l "$sv_src/$sv_exe" | awk '
	$1 == "cmd" { lc = $2 }
	lc == "LC_BUILD_VERSION"     && $1 == "minos"   { print $2 }
	lc == "LC_VERSION_MIN_MACOSX" && $1 == "version" { print $2 }' | sort -u)
[ "$sv_minos" = "$sv_target" ] || {
	echo "package_macos.sh: $sv_exe slices target macOS [$sv_minos], expected $sv_target" >&2
	exit 1
}
# Exact match, not "contains arm64": this must also reject a binary that carries
# EXTRA slices, since the package name promises arm64 and nothing else.
[ "$(lipo -archs "$sv_src/$sv_exe")" = "$sv_arch" ] || {
	echo "package_macos.sh: $sv_exe is [$(lipo -archs "$sv_src/$sv_exe")], expected $sv_arch" >&2
	exit 1
}

sv_payload=$(mktemp -d "${TMPDIR:-/tmp}/svht_denoise-pkg.XXXXXX")
trap 'rm -rf "$sv_payload"' EXIT HUP INT TERM
mkdir -p "$sv_payload/root/usr/local/bin"

cp "$sv_src/$sv_exe" "$sv_payload/root/usr/local/bin/$sv_exe"
chmod 755 "$sv_payload/root/usr/local/bin/$sv_exe"
# Clears quarantine and anything else a download or editor attached.  It cannot
# clear `com.apple.provenance`, which macOS owns and re-applies; that is why
# `pkgutil --payload-files` lists AppleDouble "._" entries beside every path.
# Those are metadata transport, not junk to be installed -- `pkgutil
# --expand-full` reassembles them into the single file, verified.
xattr -cr "$sv_payload/root/usr/local/bin/$sv_exe"
# Hardened runtime and a secure timestamp are both required for notarisation.
if [ "$sv_identity" = - ]; then
	codesign --force --sign - "$sv_payload/root/usr/local/bin/$sv_exe"
else
	codesign --force --options runtime --timestamp \
		--sign "$sv_identity" "$sv_payload/root/usr/local/bin/$sv_exe"
fi
codesign --verify --strict --verbose=2 "$sv_payload/root/usr/local/bin/$sv_exe"

pkgbuild --root "$sv_payload/root" --identifier "$sv_identifier" \
	--version "$sv_version" --install-location / \
	"$sv_payload/component.pkg"

# `productbuild --package` cannot declare a Read Me pane, so synthesize the
# distribution and add one.  Without this packaging/README-macOS.txt ships
# nowhere and no installing user ever sees it.
productbuild --synthesize --package "$sv_payload/component.pkg" "$sv_payload/dist.xml.in"

# `--synthesize` always writes hostArchitectures="x86_64,arm64".  It does NOT
# read the payload to decide that -- synthesizing from an x86_64-only component
# emits the identical line -- so with an arm64-only payload the package would
# tell macOS Installer it is installable on Intel, where it would install
# happily and then die with "Bad CPU type in executable" on every invocation.
# This attribute is the ONLY thing that refuses that install; the Read Me and
# the README can ask, and are not attached on the adhoc path anyway.  Rewritten
# for BOTH paths, before the Read Me edit, because an adhoc package is the one
# most likely to be handed to someone informally.
sed 's|hostArchitectures="[^"]*"|hostArchitectures="'"$sv_arch"'"|' \
	"$sv_payload/dist.xml.in" > "$sv_payload/dist.xml.arch"
# Fail closed if that did not take: a future productbuild that renames or drops
# the attribute would otherwise silently restore the Intel-installable package,
# which is invisible until an Intel user reports a broken binary.
grep -q "hostArchitectures=\"$sv_arch\"" "$sv_payload/dist.xml.arch" || {
	echo "package_macos.sh: could not set hostArchitectures=\"$sv_arch\" in the" >&2
	echo "  synthesized distribution; refusing to build a package that would" >&2
	echo "  install on the wrong architecture." >&2
	exit 1
}

# The BSD sed line continuation below is deliberate: a GNU-style \n in the
# replacement is not portable here.
if [ -n "$sv_installer_identity" ]; then
	sed 's|</installer-gui-script>|    <readme file="README-macOS.txt"/>\
</installer-gui-script>|' "$sv_payload/dist.xml.arch" > "$sv_payload/dist.xml"
else
	cp "$sv_payload/dist.xml.arch" "$sv_payload/dist.xml"
fi

set -- --distribution "$sv_payload/dist.xml" --package-path "$sv_payload"
if [ -n "$sv_installer_identity" ]; then
	# Only a real release gets the Read Me: it states the package is signed and
	# notarized, which is false of an ad-hoc build.
	set -- "$@" --resources "$sv_root/packaging" --sign "$sv_installer_identity"
else
	echo "package_macos.sh: MACOS_INSTALLER_IDENTITY is unset, so the package will" >&2
	echo "  be unsigned. It cannot be notarised and Gatekeeper will refuse it on" >&2
	echo "  another machine. Local testing only." >&2
fi
productbuild "$@" "$sv_pkg"

echo
echo "package_macos.sh: wrote $sv_pkg"
pkgutil --check-signature "$sv_pkg" 2>&1 | sed 's/^/  /' || true
echo "  architecture: $sv_arch"
echo "  installs /usr/local/bin/$sv_exe"
