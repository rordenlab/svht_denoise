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
#   MACOS_ARCHS               architectures (default "arm64 x86_64").
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

sv_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
sv_src="$sv_root/src"
sv_version=${1:?usage: package_macos.sh VERSION}
sv_identity=${MACOS_SIGN_IDENTITY--}
sv_installer_identity=${MACOS_INSTALLER_IDENTITY:-}
sv_target=${MACOSX_DEPLOYMENT_TARGET:-11.0}
sv_archs=${MACOS_ARCHS:-arm64 x86_64}
sv_identifier=io.github.rordenlab.svht-denoise
sv_dist="$sv_root/dist"
sv_pkg="$sv_dist/svht_denoise-$sv_version-macos-universal.pkg"

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
case "$sv_archs" in
	*[!a-z0-9_\ ]*|'')
		echo "package_macos.sh: MACOS_ARCHS must be space-separated architecture names" >&2
		exit 2
		;;
esac

for sv_tool in codesign pkgbuild productbuild pkgutil otool lipo xattr security; do
	command -v "$sv_tool" >/dev/null 2>&1 || {
		echo "package_macos.sh: missing required tool: $sv_tool" >&2
		exit 2
	}
done

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

sv_arch_flags=
for sv_arch in $sv_archs; do
	sv_arch_flags="$sv_arch_flags -arch $sv_arch"
done
sv_opt="-O3$sv_arch_flags -mmacosx-version-min=$sv_target"

# Build AND run the test suite in ONE make invocation, so the binary that gets
# packaged is literally the one the tests passed against. Two invocations would
# not do that: the flag stamp would no longer match on the second, everything
# would be rebuilt at the default OPT, and the universal binary would be
# silently replaced by a host-only one. (Same trap the ubsan-test targets exist
# to avoid -- see the Makefile.)
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
# would happily publish a 9.9.9 package containing a binary that reports 1.1.0.
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
# leaves the per-architecture "...(architecture x86_64):" headers flush left.  A
# `grep -v "$sv_exe"` dropped those headers but also any real dependency whose
# path happened to contain the program name.
sv_foreign=$(otool -L "$sv_src/$sv_exe" | awk '/^\t/ { print $1 }' |
	grep -v '^/usr/lib/' | grep -v '^/System/' || true)
if [ -n "$sv_foreign" ]; then
	echo "package_macos.sh: $sv_exe has non-system dependencies:" >&2
	echo "$sv_foreign" | sed 's/^/  /' >&2
	exit 1
fi
# Every slice must claim the deployment target we advertise.  Collected and
# compared as one string rather than looped over: a loop only inspects the lines
# it finds, so NO minos load command at all would have passed silently.  sort -u
# also catches slices that disagree with each other.
# Both spellings: LC_BUILD_VERSION prints "minos", while the older
# LC_VERSION_MIN_MACOSX that a pre-10.14 slice still uses prints "version".
# Matching only "minos" yielded zero lines for such a slice, which a loop read
# as success.
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
for sv_arch in $sv_archs; do
	lipo -archs "$sv_src/$sv_exe" | tr ' ' '\n' | grep -qx "$sv_arch" || {
		echo "package_macos.sh: $sv_exe is missing the $sv_arch slice" >&2
		exit 1
	}
done

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
# nowhere and no installing user ever sees it.  The BSD sed line continuation
# below is deliberate: a GNU-style \n in the replacement is not portable here.
productbuild --synthesize --package "$sv_payload/component.pkg" "$sv_payload/dist.xml.in"
if [ -n "$sv_installer_identity" ]; then
	sed 's|</installer-gui-script>|    <readme file="README-macOS.txt"/>\
</installer-gui-script>|' "$sv_payload/dist.xml.in" > "$sv_payload/dist.xml"
else
	cp "$sv_payload/dist.xml.in" "$sv_payload/dist.xml"
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
echo "  architectures: $(lipo -archs "$sv_src/$sv_exe")"
echo "  installs /usr/local/bin/$sv_exe"
