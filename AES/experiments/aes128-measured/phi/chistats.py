#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""chistats.py -- chi-square uniformity / independence report from a phi stats.txt.

Input format (written by phi.c stats_write):
  # n=... relabel=... degenerate=... zero_in_window=... a_zero=...
  hist <k> <256 counts>          (k = 0..11, Phi byte positions)
  joint <a> <b> <65536 counts>   (selected byte-position pairs)
"""
import sys, math
import numpy as np
from scipy.stats import chi2, kstwo

def main(path):
    hists, joints, meta = {}, {}, ""
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                meta = line.strip(); continue
            parts = line.split()
            if parts[0] == 'hist':
                hists[int(parts[1])] = np.array([int(x) for x in parts[2:]], dtype=np.float64)
            elif parts[0] == 'joint':
                joints[(int(parts[1]), int(parts[2]))] = np.array([int(x) for x in parts[3:]], dtype=np.float64)
    n = hists[0].sum()
    print(meta)
    print("N = %d (2^%.3f) fingerprints\n" % (n, math.log2(n)))
    print("Per-byte marginal uniformity (chi-square, 256 cells, df=255):")
    print("  pos  chi2        p-value    max|dev|/N    minbin   maxbin")
    pvals = []
    for k in sorted(hists):
        h = hists[k]
        exp = n / 256.0
        x2 = float(((h - exp) ** 2 / exp).sum())
        p = float(chi2.sf(x2, 255))
        pvals.append(p)
        dev = float(np.abs(h - exp).max() / n)
        print("  %2d   %9.2f   %8.5f   %.3e   %d  %d" % (k, x2, p, dev, int(h.min()), int(h.max())))
    # Byte 0 against the EXACT ideal-Phi marginal: P(0) = 511/65536 (zero whenever
    # D_3 = D_1, or D_2 = D_1 which triggers the relabelling rule), other values equal.
    h0 = hists[0]
    pz = 511.0 / 65536.0
    exp0 = np.full(256, n * (1 - pz) / 255.0); exp0[0] = n * pz
    x2c = float(((h0 - exp0) ** 2 / exp0).sum())
    print("  byte 0 vs ideal-Phi marginal (P(0)=511/65536): chi2=%.2f p=%.5f ; observed P(0)=%.6f vs %.6f"
          % (x2c, float(chi2.sf(x2c, 255)), h0[0] / n, pz))
    for k in (1, 2):
        print("  byte %d observed P(0)=%.6f (ideal-Phi: %.6f)" % (k, hists[k][0] / n, 1.0/256 + (255.0/256) * 256.0 ** -(k + 1)))
    # combined: KS test of the 12 p-values against uniform
    ks = kstwo.sf if False else None
    print("  (12 p-values; min=%.4f max=%.4f; under H0 ~ Uniform(0,1))" % (min(pvals), max(pvals)))
    print("\nPairwise 2-byte joint uniformity (chi-square, 65536 cells, df=65535):")
    print("  pos-pair   chi2          p-value")
    for (a, b) in sorted(joints):
        j = joints[(a, b)]
        exp = n / 65536.0
        x2 = float(((j - exp) ** 2 / exp).sum())
        p = float(chi2.sf(x2, 65535))
        print("  (%2d,%2d)   %11.2f   %8.5f" % (a, b, x2, p))
    print()

if __name__ == '__main__':
    for p in sys.argv[1:]:
        print("==== %s ====" % p)
        main(p)
