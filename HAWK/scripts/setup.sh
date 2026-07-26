#!/usr/bin/env bash
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

# One-time setup: build the external deps with Nix (pinned, hash-verified —
# see nix/) and link them where the attack code expects them.
# Usage: ./scripts/setup.sh [--no-test-deps]
set -euo pipefail
cd "$(dirname "$0")/.."

NIXDIR=nix
REF_IMPL=hawk_recover/ref_impl

command -v nix-build > /dev/null 2>&1 || {
    echo "FATAL: nix-build not found. Install Nix first: https://nixos.org/download" >&2
    exit 2
}

# result-* links are GC roots; every ref_impl symlink points through them.

echo "[setup] svp_tool (BGJ sieve) — first build compiles NTL + ~25 C++ files…"
nix-build "$NIXDIR" -A svp-tool -o "$NIXDIR/result-svp-tool"
APPDIR=$REF_IMPL/bgj_sieve_src/BGJ-Sieve-AMX-main/app
mkdir -p "$APPDIR"
ln -sf "$(pwd)/$NIXDIR/result-svp-tool/bin/svp_tool" "$APPDIR/svp_tool"
ln -sf "$(pwd)/$NIXDIR/result-svp-tool/bin/run_local.sh" "$APPDIR/run_local.sh"

# Patched sieve source, linked for auditing (the attack needs only the binary).
nix-build "$NIXDIR" -A sieve-src -o "$NIXDIR/result-sieve-src"
ln -sfn "$(pwd)/$NIXDIR/result-sieve-src" "$REF_IMPL/bgj_sieve_src/source"

echo "[setup] hawk256 reference (commit-pinned) + sign/verify shim…"
nix-build "$NIXDIR" -A hawk256-ref -o "$NIXDIR/result-hawk256-ref"
nix-build "$NIXDIR" -A hawk-shim -o "$NIXDIR/result-hawk-shim"
ln -sfn "$(pwd)/$NIXDIR/result-hawk256-ref" "$REF_IMPL/hawk256_ref"
ln -sf "$(pwd)/$NIXDIR/result-hawk-shim/bin/hawk_shim" "$REF_IMPL/hawk_shim"

# pytest etc. for ./run_tests.sh — hash-pinned from PyPI into .deps/.
if [ "${1:-}" != "--no-test-deps" ] && [ ! -d .deps/pytest ]; then
    echo "[setup] test deps → .deps/ (pinned + hash-verified)…"
    nix-shell "$NIXDIR/shell.nix" --run \
        "sage -pip install --target .deps --require-hashes -r scripts/test-requirements.txt"
fi

echo "[setup] done."
echo "  attack : nix-shell nix/shell.nix --run 'sage -python -m hawk_recover <pk.hex> --out <sk.bin>'"
echo "  tests  : ./run_tests.sh"
