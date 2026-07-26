# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

# External dependencies as pinned derivations; scripts/setup.sh builds them
# and symlinks the results into hawk_recover/ref_impl/.
{
  pkgs ? import ./nixpkgs.nix { },
}:

{
  sieve-src = import ./sieve-src.nix { inherit pkgs; };
  svp-tool = import ./svp-tool.nix { inherit pkgs; };
  svp-tool-native = import ./svp-tool.nix {
    inherit pkgs;
    marchNative = true;
  };
  hawk256-ref = import ./hawk256-ref.nix { inherit pkgs; };
  hawk-shim = import ./hawk-shim.nix { inherit pkgs; };
}
