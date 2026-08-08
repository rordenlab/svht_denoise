#!/bin/sh
# Notarise and staple a signed installer package.
#
#   notarize_macos.sh FILE.pkg [KEYCHAIN_PROFILE]
set -eu

sv_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
sv_pkg=${1:?usage: notarize_macos.sh FILE.pkg [KEYCHAIN_PROFILE]}
sv_profile=${2:-svht-denoise-notary}

[ -f "$sv_pkg" ] || {
	echo "notarize_macos.sh: file not found: $sv_pkg" >&2
	exit 2
}

# An unsigned package reaches notarytool and fails there with a much less
# obvious message, so reject it here. Packages are signed with a Developer ID
# Installer certificate, which pkgutil checks -- codesign does not verify
# package signatures.
pkgutil --check-signature "$sv_pkg" >/dev/null 2>&1 || {
	echo "notarize_macos.sh: $sv_pkg is unsigned or its signature is invalid." >&2
	echo "  Set MACOS_INSTALLER_IDENTITY and rebuild with 'make macos-pkg'." >&2
	exit 2
}

xcrun notarytool submit "$sv_pkg" --keychain-profile "$sv_profile" --wait
xcrun stapler staple "$sv_pkg"
xcrun stapler validate "$sv_pkg"
"$sv_root/scripts/verify_macos_pkg.sh" "$sv_pkg"
# Installer packages are assessed under the "install" policy; "open" is for
# applications and would report a spurious rejection here.
spctl --assess --type install --verbose=4 "$sv_pkg"
shasum -a 256 "$sv_pkg"
echo "Notarized, stapled, and Gatekeeper-validated $sv_pkg"
