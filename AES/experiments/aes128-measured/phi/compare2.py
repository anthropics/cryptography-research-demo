#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""compare2.py A.stats B.stats -- two-sample chi-square homogeneity tests:
is the real Phi byte distribution (A) consistent with the ideal-model Phi (B)?
Per byte position (df=255) and per stored 2-byte joint (df=65535)."""
import sys, math
import numpy as np
from scipy.stats import chi2

def load(path):
    hists, joints = {}, {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split()
            if parts[0] == 'hist':
                hists[int(parts[1])] = np.array([int(x) for x in parts[2:]], dtype=np.float64)
            elif parts[0] == 'joint':
                joints[(int(parts[1]), int(parts[2]))] = np.array([int(x) for x in parts[3:]], dtype=np.float64)
    return hists, joints

def twosample(a, b):
    na, nb = a.sum(), b.sum()
    t = a + b
    mask = t > 0
    stat = float((((math.sqrt(nb/na) * a - math.sqrt(na/nb) * b) ** 2)[mask] / t[mask]).sum())
    return stat

def main(pa, pb):
    ha, ja = load(pa); hb, jb = load(pb)
    print("A = %s  (N=%d)\nB = %s  (N=%d)\n" % (pa, ha[0].sum(), pb, hb[0].sum()))
    print("Per-byte two-sample homogeneity (df=255):")
    for k in sorted(ha):
        s = twosample(ha[k], hb[k])
        print("  pos %2d: chi2=%9.2f  p=%.5f" % (k, s, chi2.sf(s, 255)))
    print("\nJoint two-sample homogeneity (df=65535):")
    for key in sorted(ja):
        if key not in jb: continue
        s = twosample(ja[key], jb[key])
        print("  pos %s: chi2=%11.2f  p=%.5f" % (key, s, chi2.sf(s, 65535)))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
