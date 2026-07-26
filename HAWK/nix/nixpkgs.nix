# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

# Pinned nixpkgs: the paper's reference environment (Sage 10.7, PARI 2.17.3,
# fpylll 0.6.4, LLVM 18.1.8). A public revision on purpose — the artifact
# must build outside any particular checkout. Override via `pkgs`.
import (
  builtins.fetchTarball {
    name = "nixpkgs-549bd84d";
    url = "https://github.com/NixOS/nixpkgs/archive/549bd84d6279f9852cae6225e372cc67fb91a4c1.tar.gz";
    sha256 = "sha256-hGdgeU2Nk87RAuZyYjyDjFL6LK7dAZN5RE9+hrDTkDU=";
  }
)
