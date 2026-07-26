# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

# Run/test environment: SageMath from the pinned nixpkgs. Test deps live in
# .deps/ (installed by setup.sh, put on PYTHONPATH by run_tests.sh).
{
  pkgs ? import ./nixpkgs.nix { },
}:

pkgs.mkShell {
  packages = [ pkgs.sage ];
}
