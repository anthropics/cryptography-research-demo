# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

# svp_tool — the BGJ sieve binary, clang/libc++/OpenMP with Nix GMP/NTL (the
# research runs' toolchain); bgj.py invokes the run_local.sh emitted next to
# it. cascadelake ISA by default for reproducibility; marchNative overrides.
{
  pkgs ? import ./nixpkgs.nix { },
  marchNative ? false,
  sieveSrc ? import ./sieve-src.nix { inherit pkgs; },
}:

let
  llvm = pkgs.llvmPackages_18;
  # cascadelake = oldest ISA with the required AVX512-VNNI. AMX kernels are
  # compiled in but runtime-gated (-amx), so one binary serves all machines.
  march = if marchNative then "-march=native" else "-march=cascadelake -mamx-tile -mamx-int8";
  # NTL must share the sieve's C++ stdlib ABI (libc++); stock nixpkgs ntl is
  # gcc/libstdc++, so rebuild it. GMP is C — fine as-is. unpack/sourceRoot:
  # same no-chmod treatment as below.
  ntl = (pkgs.ntl.override { stdenv = llvm.libcxxStdenv; }).overrideAttrs (old: {
    sourceRoot = null;
    unpackPhase = ''
      runHook preUnpack
      tar xf $src --no-same-owner --no-same-permissions
      cd ntl-*/src
      runHook postUnpack
    '';
  });
in
# One libc++ stdenv for everything keeps every object and the link on one ABI.
llvm.libcxxStdenv.mkDerivation {
  pname = "svp-tool";
  version = "unstable-2025-01-11";

  src = sieveSrc;

  buildInputs = [
    llvm.openmp
    pkgs.gmp
    ntl
  ];

  # Same no-directory-chmod unpack as sieve-src.nix (src is its store dir).
  unpackPhase = ''
    runHook preUnpack
    mkdir source
    (cd $src && tar cf - --mode=u+rwX .) | (cd source && tar xf - --no-same-owner --no-same-permissions)
    cd source
    runHook postUnpack
  '';

  # Upstream Makefiles hardcode ../dep/{gmp,ntl}; these variants point at the
  # Nix-provided ones. Objects and flags otherwise match upstream.
  postPatch = ''
    cat > src/Makefile.nix << MKEOF
    OBJ = config.o pool_epi8.o bucket_epi8.o bgj_epi8.o 3Sieve_parallel.o naivedh_epi8.o naivesh_epi8.o fcSieve.o bgj1.o bgjf32.o dual.o
    OBJ += vec.o quad.o utils.o lattice.o lll.o UidHashTable.o pool.o 3Sieve.o svp.o bkz.o sampler.o lsfdh_epi8.o lsfsh_epi8.o
    OBJ += bgj_amx.o bucket_amx.o booster_amx.o
    CXX = clang++
    CFLAGS = -O3 -fopenmp=libomp ${march} -pthread -stdlib=libc++ -Wno-deprecated-declarations
    INC_DIR = -I../include/ -I${ntl}/include -I${pkgs.gmp.dev}/include
    all: libllib.a
    libllib.a: \$(OBJ)
    	ar rcs libllib.a \$(OBJ)
    %.o: %.cpp
    	\$(CXX) \$(CFLAGS) \$^ \$(INC_DIR) -c
    MKEOF

    cat > app/Makefile.nix << MKEOF
    CXX = clang++
    CFLAGS = -O3 -funroll-loops -ftree-vectorize -fopenmp=libomp ${march} -stdlib=libc++ -pthread -Wno-deprecated-declarations
    LIBS = -L../src -lllib -L${ntl}/lib -lntl -L${pkgs.gmp}/lib -lgmp -L${llvm.openmp}/lib -lomp -lm
    INC_DIR = -I../include/ -I${ntl}/include -I${pkgs.gmp.dev}/include
    all: svp_tool
    %: bin_%.cpp ../src/libllib.a
    	\$(CXX) \$(CFLAGS) \$^ \$(INC_DIR) \$(LIBS) -o \$*
    MKEOF
  '';

  buildPhase = ''
    runHook preBuild
    make -C src -f Makefile.nix -j$NIX_BUILD_CORES
    make -C app -f Makefile.nix svp_tool
    runHook postBuild
  '';

  # bgj.py runs run_local.sh with a scrubbed PATH, so the wrapper uses the
  # store bash shebang and absolute paths; the binary's RPATH covers its libs.
  installPhase = ''
    runHook preInstall
    install -Dm755 app/svp_tool $out/bin/svp_tool
    cat > $out/bin/run_local.sh << WRAPEOF
    #!${pkgs.runtimeShell}
    exec "$out/bin/svp_tool" "\$@"
    WRAPEOF
    chmod 755 $out/bin/run_local.sh
    runHook postInstall
  '';

  meta = with pkgs.lib; {
    description = "BGJ lattice sieve binary (svp_tool) for hawk_recover";
    license = licenses.gpl2Only;
    platforms = [ "x86_64-linux" ];
  };
}
