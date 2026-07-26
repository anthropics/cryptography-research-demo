# Python reference implementation of the attack layers

Pedagogical, correctness-first Python for the attack's building blocks
(standard library only).  The performance-oriented and black-box
implementations are in `../../experiments/`.

| file | what it shows |
|---|---|
| `gf256.py` | GF(256) mul/inv/pow, S-box, `M^{-1}`. Self-tested. |
| `aes7r.py` | 7-round AES-128 encrypt + `trace()`, checked against the FIPS-197 vector. |
| `index_maps.py` | SR/MC index wiring used by the e-sequence. |
| `eseq_exact.py` | the exact e-sequence from the internal state, with the full index maps. |
| `prop2.py` | naive e-sequence, Gray walk, U-table S-box cache. |
| `chi_canon.py` | the AGL-invariant χ fingerprint + `verify_agl_invariance()`. |
| `test_vectors.py` / `test_vectors.json` | known-answer test for the internal state trace. |

```bash
python3 gf256.py         # self-test
python3 aes7r.py         # FIPS vector check
python3 chi_canon.py     # AGL-invariance check
python3 test_vectors.py  # regenerate the golden vectors
```

Caveat: `prop2.py`'s diagonal index map uses a simplified closed form valid
for the diag-0 backbone used in the demo; the exact SR∘MC wiring for all 16
positions is in `eseq_exact.py` and in
`../../experiments/aes128-measured/common/prop2.h`.
