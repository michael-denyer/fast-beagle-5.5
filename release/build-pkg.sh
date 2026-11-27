#!/bin/bash
# Builds the macOS installer fast-beagle-<version>-macos-arm64.pkg from the
# unpacked release package of release/build-static.sh. The installer puts
# fast-beagle in /usr/local/bin and the licenses in
# /usr/local/share/doc/fast-beagle. It signs the binary and the installer with
# the Developer ID identities in the keychain search list, has Apple notarise
# the installer, and staples the ticket to it, so Gatekeeper accepts it offline.
#
# Use --unsigned to build a preview without accessing Apple credentials.
# Usage: release/build-pkg.sh [--unsigned] <version> <unpacked package directory> <output directory>
# Signing needs NOTARY_KEY (the path of the App Store Connect API key .p8),
# NOTARY_KEY_ID and NOTARY_ISSUER_ID.
set -euo pipefail
unsigned=false
if [[ ${1:-} == --unsigned ]]; then
  unsigned=true
  shift
fi
if [[ $# -ne 3 ]]; then
  echo "usage: $0 [--unsigned] <version> <package directory> <output directory>" >&2
  exit 2
fi
version=$1 package=$2 out=$3
pkg="$out/fast-beagle-$version-macos-arm64.pkg"

root=$(mktemp -d)
trap 'rm -rf "$root"' EXIT
mkdir -p "$root/usr/local/bin" "$root/usr/local/share/doc/fast-beagle" "$out"
cp "$package/fast-beagle" "$root/usr/local/bin/fast-beagle"
cp "$package"/LICENSE* "$package/NOTICE.md" "$root/usr/local/share/doc/fast-beagle/"

pkg_args=(--root "$root" --install-location / --identifier io.github.michael-denyer.fast-beagle --version "$version")
if [[ $unsigned == true ]]; then
  preview="${pkg%.pkg}-unsigned.pkg"
  pkgbuild "${pkg_args[@]}" "$preview"
  # pkgutil exits nonzero for an unsigned package. Require its specific verdict.
  signature=$(LC_ALL=C pkgutil --check-signature "$preview" || true)
  grep -Eq '^[[:space:]]*Status: no signature$' <<< "$signature"
  echo "$preview"
  exit 0
fi

# identity <kind>: the SHA-1 of the one "Developer ID <kind>" identity.
identity() {
  local found
  found=$(security find-identity -v | awk -v kind="\"Developer ID $1: " 'index($0, kind) { print $2 }')
  [ "$(echo "$found" | grep -c .)" -eq 1 ] \
    || { echo "build-pkg.sh: need exactly one Developer ID $1 identity, found: ${found:-none}" >&2; exit 1; }
  echo "$found"
}
application=$(identity Application)
installer=$(identity Installer)

# Notarisation requires the hardened runtime and a secure timestamp.
codesign --force --timestamp --options runtime --sign "$application" "$root/usr/local/bin/fast-beagle"
codesign --verify --strict --verbose=2 "$root/usr/local/bin/fast-beagle"

pkgbuild "${pkg_args[@]}" --sign "$installer" --timestamp "$pkg"
pkgutil --check-signature "$pkg"

notary=(--key "$NOTARY_KEY" --key-id "$NOTARY_KEY_ID" --issuer "$NOTARY_ISSUER_ID")
# Lists the earlier submissions, which fails at once on a key Apple rejects,
# before the upload and the wait.
xcrun notarytool history "${notary[@]}"
id=$(xcrun notarytool submit "$pkg" "${notary[@]}" --output-format json | sed -n 's/.*"id" *: *"\([^"]*\)".*/\1/p')
echo "notary submission $id"
# The wait can time out and notarytool exits 0 for a rejected submission, so
# read the status afterwards.
xcrun notarytool wait "$id" "${notary[@]}" --timeout 45m || true
status=$(xcrun notarytool info "$id" "${notary[@]}" --output-format json | sed -n 's/.*"status" *: *"\([^"]*\)".*/\1/p')
if [ "$status" != Accepted ]; then
  xcrun notarytool log "$id" "${notary[@]}" >&2 || true
  echo "build-pkg.sh: notary submission $id is $status, not Accepted" >&2
  exit 1
fi
xcrun stapler staple "$pkg"
xcrun stapler validate "$pkg"
spctl --assess --type install --verbose=2 "$pkg"
echo "$pkg"
