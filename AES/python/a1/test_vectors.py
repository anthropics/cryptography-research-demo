# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""Golden test vectors for the A1 attack pipeline.
Run this to regenerate; diff your implementation's output against VECTORS."""

from aes7r import key_schedule, trace
from index_maps import X5_POS
from eseq_exact import e_sequence_exact
from chi_canon import chi_canon
from gf256 import inv

K = bytes.fromhex("2b7e151628aed2a6abf7158809cf4f3c")   # FIPS-197 key
P = bytes.fromhex("3243f6a8885a308d313198a2e0370734")   # FIPS-197 plaintext
rk = key_schedule(K)
x = trace(P, rk)

params = dict(
    x2_0 = x[2][0],
    x3   = [x[3][i] for i in range(4)],
    x4   = list(x[4]),
    x5   = [x[5][p] for p in X5_POS],
    x6_0 = x[6][0],
)

e = e_sequence_exact(params['x2_0'], params['x3'], params['x4'], params['x5'])
d = [0]+[inv(ew) for ew in e[1:]]
fp = chi_canon(d)

VECTORS = {
    "K": K.hex(),
    "P": P.hex(),
    "x2_0": params['x2_0'],            # 0xa4
    "x3_col0": params['x3'],           # [0xaa,0x8f,0x5f,0x03]
    "x5_diag0": params['x5'],          # [0xe0,0x63,0x35,0x01]
    "x6_0": params['x6_0'],            # 0xf1
    "e_1to8": e[1:9],
    "chi_canon_d": fp.hex(),
    "K7_adiag0": [rk[7][p] for p in (0,13,10,7)],
}

if __name__ == "__main__":
    import json
    for k,v in VECTORS.items():
        print(f"  {k:14s} = {v}")
    with open("test_vectors.json","w") as f:
        json.dump({k:(v if not isinstance(v,bytes) else v.hex()) for k,v in VECTORS.items()},
                  f, indent=2, default=lambda x: list(x) if hasattr(x,'__iter__') else x)
    print("\n→ test_vectors.json written")
