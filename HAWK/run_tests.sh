#!/usr/bin/env bash
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail
cd "$(dirname "$0")"
export PYTHONPATH="$PWD:$PWD/.deps${PYTHONPATH:+:$PYTHONPATH}"
if [ $# -eq 0 ]; then
    set -- hawk_recover/tests
fi
# printf %q re-quotes each argument across the nix-shell --run boundary.
exec nix-shell "$PWD/nix/shell.nix" --run "sage -python -m pytest $(printf '%q ' "$@")"
