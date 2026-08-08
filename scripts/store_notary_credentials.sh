#!/bin/sh
# Store notarisation credentials in the Keychain, once per machine.
#
#   store_notary_credentials.sh PROFILE APPLE_ID TEAM_ID
set -eu

sv_profile=${1:?usage: store_notary_credentials.sh PROFILE APPLE_ID TEAM_ID}
sv_apple_id=${2:?usage: store_notary_credentials.sh PROFILE APPLE_ID TEAM_ID}
sv_team_id=${3:?usage: store_notary_credentials.sh PROFILE APPLE_ID TEAM_ID}

echo "Enter the app-specific password at Apple's secure prompt. It will be stored in Keychain."
echo "Do not pass the password as a command-line argument or Make variable."
xcrun notarytool store-credentials "$sv_profile" \
	--apple-id "$sv_apple_id" \
	--team-id "$sv_team_id"
