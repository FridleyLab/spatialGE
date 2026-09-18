## Shared environment for the cross-version comparison harness.
## Source this, don't execute it:  . dev/xver/env.sh
##
## Everything the harness needs that is not repo state lives here, so the
## individual scripts stay free of machine-specific paths.

XVER_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export XVER_REPO

# Where the two installed libraries, extracted sources, payloads and reports go.
# Deliberately outside the repo: these are large, machine-specific build
# artifacts, and nothing here should ever be a candidate for a commit.
export XVER_ROOT="${XVER_ROOT:-$HOME/.spatialge-dev/xver}"

export XVER_SRC_MAIN="$XVER_ROOT/src-main"
export XVER_SRC_ALEX="$XVER_ROOT/src-alex"
export XVER_LIB_MAIN="$XVER_ROOT/lib-main"
export XVER_LIB_ALEX="$XVER_ROOT/lib-alex"
export XVER_OUT="$XVER_ROOT/out"

# The two refs under comparison. Pinned by name, asserted by SHA in 00_gate.sh.
export XVER_REF_MAIN="${XVER_REF_MAIN:-main}"
export XVER_REF_ALEX="${XVER_REF_ALEX:-alex-refactor}"

# This network intercepts TLS with Netskope and conda's libcurl does not trust
# the intercepting root, so every download.file() from conda R fails until these
# point at a bundle containing it. See dev/environment/README.md.
XVER_CA="$HOME/.spatialge-dev/combined_cacert.pem"
if [ -f "$XVER_CA" ]; then
  export CURL_CA_BUNDLE="$XVER_CA"
  export SSL_CERT_FILE="$XVER_CA"
fi

# _R_CHECK_LIMIT_CORES_ makes count_cores() clamp on alex-refactor but not on
# main (fidelity rule 12), so leaving it set would compare two different fork
# counts and call it a version difference.
unset _R_CHECK_LIMIT_CORES_

export R_MAKEVARS_USER="$XVER_REPO/dev/environment/Makevars"

# Keep a stray .Rprofile / renv autoloader from rewriting .libPaths() out from
# under a producer and silently pointing it at the wrong spatialGE.
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
export R_PROFILE_USER=/dev/null

mkdir -p "$XVER_ROOT" "$XVER_OUT"
