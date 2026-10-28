#!/bin/bash
set -euo pipefail

: "${MACOS_SIGNING_KEYCHAIN:?Set MACOS_SIGNING_KEYCHAIN to the signing keychain}"
: "${MACOS_SIGNING_IDENTITY:?Set MACOS_SIGNING_IDENTITY to the certificate name}"

certificate_dir=$(mktemp -d)
trap 'rm -rf "$certificate_dir"' EXIT

# Preserve system trust settings. Add only the missing intermediate certificates
# to the temporary keychain created by InstallAppleCertificate.
installed_certificates=$(security find-certificate -a -Z "$MACOS_SIGNING_KEYCHAIN")
for certificate in DeveloperIDCA DeveloperIDG2CA; do
    certificate_path=$certificate_dir/$certificate.cer
    curl --fail --location --silent --show-error --retry 3 \
        "https://www.apple.com/certificateauthority/$certificate.cer" \
        -o "$certificate_path"
    fingerprint=$(openssl x509 -inform DER -in "$certificate_path" -noout -fingerprint -sha256)
    fingerprint=${fingerprint#*=}
    fingerprint=${fingerprint//:/}
    if [[ "$installed_certificates" != *"$fingerprint"* ]]; then
        security import "$certificate_path" -t cert -k "$MACOS_SIGNING_KEYCHAIN"
    fi
done

# Report the issuer and expiry, then validate the actual signing certificate
# against the imported intermediates and macOS's existing trusted roots.
security find-certificate -c "$MACOS_SIGNING_IDENTITY" -p "$MACOS_SIGNING_KEYCHAIN" \
    > "$certificate_dir/signer.pem"
openssl x509 -in "$certificate_dir/signer.pem" -noout -issuer -dates
security verify-cert -p codeSign -L -c "$certificate_dir/signer.pem" \
    -k "$MACOS_SIGNING_KEYCHAIN"
