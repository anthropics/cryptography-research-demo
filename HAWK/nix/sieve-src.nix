# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

# BGJ-Sieve-AMX source (Zhao–Ding–Yang, "Sieving with Streaming Memory
# Access", GPL-2.0) — the `bgj` backend's sieve. The only source change:
# the compiled-in thread cap is raised 112 → 192 (bgj.py clamps to the
# same limit). Upstream's install_gmp.sh / install_ntl.sh — which download
# their deps without integrity checks — are never invoked; the Nix build
# provides GMP and NTL as pinned inputs (see svp-tool.nix).
{
  pkgs ? import ./nixpkgs.nix { },
}:

pkgs.stdenv.mkDerivation {
  pname = "bgj-sieve-amx-src";
  version = "unstable-2025-01-11";

  src = pkgs.fetchFromGitHub {
    owner = "zhaoziyu0008";
    repo = "BGJ-Sieve-AMX";
    rev = "62f86a530166ab1eff908c10f84f01772353a910";
    hash = "sha256-ZIfZdwaOAGrT1mqnyPChoAiWFMjQ7Sa3wIKQilRoKMM=";
  };

  # A build-time substitution rather than a shipped patch file, so no
  # upstream (GPL) text lives in this repository. --replace-fail makes an
  # upstream rename of the constant a build error, not a silent no-op.
  postPatch = ''
    substituteInPlace include/bgj_epi8.h \
      --replace-fail '#define MAX_NTHREADS 112' '#define MAX_NTHREADS 192'
    substituteInPlace include/bgj_amx.h \
      --replace-fail '#define AMX_MAX_NTHREADS 112' '#define AMX_MAX_NTHREADS 192'
  '';

  # Unpack without chmodding directories — some sandboxed builders forbid it.
  # sourceRoot must stay unset: setting it triggers a chmod in the phase runner.
  unpackPhase = ''
    runHook preUnpack
    mkdir source
    (cd $src && tar cf - --mode=u+rwX .) | (cd source && tar xf - --no-same-owner --no-same-permissions)
    cd source
    runHook postUnpack
  '';

  # Source-only output: built by svp-tool.nix, linked for audit by setup.sh.
  dontBuild = true;
  installPhase = ''
    runHook preInstall
    cp -r . $out
    runHook postInstall
  '';

  meta = with pkgs.lib; {
    description = "BGJ lattice sieve (Zhao-Ding-Yang) source, thread cap raised for hawk_recover";
    license = licenses.gpl2Only;
    platforms = platforms.linux;
  };
}
