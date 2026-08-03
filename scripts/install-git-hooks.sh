#!/usr/bin/env bash
# Installs the repo's git hooks. Idempotent.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
install -m 0755 "$root/scripts/pre-commit" "$root/.git/hooks/pre-commit"
echo "installed .git/hooks/pre-commit"
