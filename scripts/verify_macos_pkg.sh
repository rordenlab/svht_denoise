#!/bin/sh
# Verify a built installer package WITHOUT installing it.
#
#   verify_macos_pkg.sh FILE.pkg
#
# NOT safe on a package from an untrusted source: it executes the payload.
#
# Expands the package into a scratch directory and checks the executable it
# would install: that it is correctly signed, carries no dependency outside
# macOS itself, has no build-machine path baked in, and actually runs.  It
# REPORTS the architectures rather than asserting them; the arm64-only gate is
# in package_macos.sh, before signing, so that this stays usable on any package
# including an older universal one. The last point matters most -- a package can be
# perfectly signed and still install a binary that dies looking for a library
# that is not there.
set -eu

# README.md hands this script to users, so it is the one most likely to be tried
# on a machine that cannot run it.
if [ "$(uname -s)" != Darwin ]; then
	echo "verify_macos_pkg.sh requires macOS" >&2
	exit 2
fi

sv_pkg=${1:?usage: verify_macos_pkg.sh FILE.pkg}
[ -f "$sv_pkg" ] || {
	echo "verify_macos_pkg.sh: file not found: $sv_pkg" >&2
	exit 2
}

# pkgutil exits nonzero for an unsigned package, which under `set -e` would end
# this script before it checked anything that matters. An ad-hoc build is
# legitimately unsigned, so distinguish that from a signature that is present
# and broken: warn for the first, fail for the second. notarize_macos.sh is
# where an unsigned package is actually refused.
if sv_signature=$(pkgutil --check-signature "$sv_pkg" 2>&1); then
	echo "$sv_signature" | sed 's/^/  /'
else
	case "$sv_signature" in
	*"no signature"*)
		echo "verify_macos_pkg.sh: WARNING - package is unsigned; this build" >&2
		echo "  cannot be notarised and Gatekeeper will refuse it elsewhere." >&2
		;;
	*)
		echo "$sv_signature" >&2
		echo "verify_macos_pkg.sh: package signature is invalid" >&2
		exit 1
		;;
	esac
fi

sv_work=$(mktemp -d "${TMPDIR:-/tmp}/svht_denoise-verify.XXXXXX")
trap 'rm -rf "$sv_work"' EXIT HUP INT TERM

# pkgutil --expand-full unpacks the payload as well as the metadata.
pkgutil --expand-full "$sv_pkg" "$sv_work/expanded"

sv_bin=$(find "$sv_work/expanded" -type d -path '*/usr/local/bin' -print -quit)
[ -n "$sv_bin" ] || {
	echo "verify_macos_pkg.sh: package does not install into /usr/local/bin" >&2
	exit 1
}

sv_count=0
for sv_exe in "$sv_bin"/*; do
	[ -f "$sv_exe" ] || continue
	sv_name=$(basename -- "$sv_exe")
	sv_count=$((sv_count + 1))

	codesign --verify --strict --verbose=2 "$sv_exe"

	# Anything outside /usr/lib and /System is a library the target machine is
	# not guaranteed to have.
	# See package_macos.sh: match on otool's leading tab rather than filtering by
	# name, which also dropped real dependencies under a like-named directory.
	sv_foreign=$(otool -L "$sv_exe" | awk '/^\t/ { print $1 }' |
		grep -v '^/usr/lib/' | grep -v '^/System/' || true)
	if [ -n "$sv_foreign" ]; then
		echo "verify_macos_pkg.sh: $sv_name has non-system dependencies:" >&2
		echo "$sv_foreign" | sed 's/^/  /' >&2
		exit 1
	fi
	if otool -L "$sv_exe" | grep -Eq '/opt/homebrew|/usr/local/lib'; then
		echo "verify_macos_pkg.sh: $sv_name contains a build-machine library path" >&2
		exit 1
	fi

	chmod +x "$sv_exe"   # --expand-full does not preserve the mode
	# -help rather than a denoising run: it needs no input file, exercises the
	# dynamic loader and the whole link, and is the one invocation guaranteed to
	# exist in every build.
	"$sv_exe" -help >/dev/null
	echo "  ok $sv_name ($(lipo -archs "$sv_exe"))"
done

[ "$sv_count" -gt 0 ] || {
	echo "verify_macos_pkg.sh: package installs no executables" >&2
	exit 1
}

[ "$sv_count" -eq 1 ] && sv_plural="executable" || sv_plural="executables"
echo "Verified $sv_pkg ($sv_count $sv_plural)"
