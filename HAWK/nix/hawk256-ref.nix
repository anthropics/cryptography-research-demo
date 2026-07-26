# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

# HAWK reference implementation from the HAWK team's development
# repository (MIT-licensed): the decode ground truth + the sign/verify
# oracle's sources. Pinned to an exact commit.
{
  pkgs ? import ./nixpkgs.nix { },
}:

pkgs.stdenv.mkDerivation {
  pname = "hawk256-ref";
  version = "unstable-2025-07-04";

  src = pkgs.fetchFromGitHub {
    owner = "hawk-sign";
    repo = "dev";
    rev = "1b9fef52559273fe7b40fe3e22968eaedd3a4c2a";
    hash = "sha256-UPGZuVa4WHR0z8Jgmh7/1THjfdWWXprtZ4m4bxv0fdI=";
  };

  # The portable C sources live in src/ (Reference_Implementation/ and
  # Optimized_Implementation/ are generated from it by the upstream build).
  # Read-only: work from the store path, no copy into the build dir.
  dontUnpack = true;
  dontBuild = true;

  installPhase = ''
    runHook preInstall

    ref=$src/src
    for f in hawk.h hawk_inner.h hawk_kgen.c hawk_sign.c hawk_vrfy.c \
             ng_fxp.c ng_hawk.c ng_mp31.c ng_ntru.c ng_poly.c ng_zint31.c \
             sha3.c; do
      if [ ! -f "$ref/$f" ]; then
        echo "error: expected source $f missing from the pinned HAWK repository" >&2
        ls "$ref" >&2
        exit 1
      fi
    done

    cp -r --no-preserve=mode,ownership "$ref" $out
    install -m644 "$src/LICENSE.txt" $out/LICENSE.txt
    runHook postInstall
  '';

  meta = with pkgs.lib; {
    description = "HAWK reference implementation (hawk-sign/dev), commit-pinned";
    homepage = "https://github.com/hawk-sign/dev";
    license = licenses.mit;
    platforms = platforms.all;
  };
}
