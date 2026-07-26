# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

# hawk_shim — the sign/verify oracle: hawk_shim.c compiled against the NIST
# reference sources. infra/shim.py runs it on every recovered key.
{
  pkgs ? import ./nixpkgs.nix { },
  hawk256Ref ? import ./hawk256-ref.nix { inherit pkgs; },
}:

pkgs.stdenv.mkDerivation {
  pname = "hawk-shim";
  version = "0.1";

  src = ../hawk_recover/ref_impl/hawk_shim.c;
  dontUnpack = true;

  buildPhase = ''
    runHook preBuild
    # Same compilation the README documents.
    cc -O2 -I${hawk256Ref} -o hawk_shim $src \
      ${hawk256Ref}/hawk_sign.c ${hawk256Ref}/hawk_vrfy.c ${hawk256Ref}/hawk_kgen.c \
      ${hawk256Ref}/ng_fxp.c ${hawk256Ref}/ng_hawk.c ${hawk256Ref}/ng_mp31.c \
      ${hawk256Ref}/ng_ntru.c ${hawk256Ref}/ng_poly.c ${hawk256Ref}/ng_zint31.c \
      ${hawk256Ref}/sha3.c -lm
    runHook postBuild
  '';

  installPhase = ''
    runHook preInstall
    install -Dm755 hawk_shim $out/bin/hawk_shim
    runHook postInstall
  '';

  meta = with pkgs.lib; {
    description = "HAWK-256 sign/verify oracle for hawk_recover";
    platforms = platforms.linux;
  };
}
