#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""
Plot the prefix collision-entropy curve H2_hat(k) vs k for the I_{m,n}
fingerprint, from clump_results.tsv produced by ../clump.c.

Output: fig_clump_entropy.pdf  (Type-3-font free; fonttype 42 = TrueType).
"""

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
# --- ensure no Type-3 fonts in the PDF/PS output -------------------------
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["text.usetex"]  = False
matplotlib.rcParams["font.family"]  = "DejaVu Sans"
import matplotlib.pyplot as plt


def load(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, parts)))
    return rows


def main():
    tsv = sys.argv[1] if len(sys.argv) > 1 else "../clump_results.tsv"
    rows = load(tsv)

    # Use the run with the largest M (most samples).
    Ms = sorted({int(r["M"]) for r in rows})
    M = Ms[-1]
    rows = [r for r in rows if int(r["M"]) == M and int(r["k"]) <= 7]
    rows.sort(key=lambda r: int(r["k"]))

    k     = np.array([int(r["k"]) for r in rows])
    C     = np.array([int(r["collisions"]) for r in rows], dtype=np.float64)
    # recompute H2 and ideal at full precision from exact M, C (TSV floats
    # are rounded to 6 dp, which is coarser than sigma at small k)
    logM2 = np.log2(float(M)) + np.log2(float(M - 1))
    with np.errstate(divide="ignore"):
        H2  = logM2 - 1.0 - np.log2(C)
        err = 1.0 / (np.sqrt(C) * np.log(2))
    ideal = k * np.log2(255.0)
    finite = C > 0

    resid = H2 - ideal
    z = resid / err

    # --- figure 1: absolute H2 vs ideal line ------------------------------
    fig1, ax = plt.subplots(figsize=(4, 3), constrained_layout=True)
    kk = np.linspace(0, k.max() + 0.5, 50)
    ax.plot(kk, kk * np.log2(255), "-", color="0.6", lw=1.2,
            label=r"ideal: $k\cdot\log_2 255$")
    ax.errorbar(k[finite], H2[finite], yerr=err[finite],
                fmt="o", ms=4, color="C0", capsize=2, lw=1,
                label=r"measured $\hat H_2(k)$")
    ax.set_xlabel(r"prefix length $k$ (bytes)")
    ax.set_ylabel(r"collision entropy $\hat H_2$ (bits)")
    ax.set_xlim(0, k.max() + 0.6)
    ax.set_ylim(0, (k.max() + 0.6) * np.log2(255))
    ax.set_xticks(range(1, k.max() + 1))
    ax.legend(frameon=False, loc="upper left", fontsize=8)
    ax.text(0.98, 0.04,
            rf"$N={M:,}$".replace(",", r"\,"),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7)
    fig1.savefig("fig_clump_entropy.pdf")
    fig1.savefig("fig_clump_entropy.png", dpi=200)

    # --- figure 2: residuals (z-scores) -----------------------------------
    fig2, axr = plt.subplots(figsize=(4, 3), constrained_layout=True)
    axr.axhline(0, color="0.6", lw=1)
    axr.axhspan(-1, 1, color="0.85", lw=0, zorder=0)
    axr.errorbar(k[finite], z[finite], yerr=np.ones_like(z[finite]),
                 fmt="o", ms=4, color="C0", capsize=2, lw=1)
    axr.set_xlabel(r"prefix length $k$ (bytes)")
    axr.set_ylabel(r"$(\hat H_2 - k\log_2 255)\,/\,\sigma_k$")
    axr.set_xlim(0, k.max() + 0.6)
    axr.set_xticks(range(1, k.max() + 1))
    axr.set_ylim(-3, 3)
    fig2.savefig("fig_clump_residuals.pdf")
    fig2.savefig("fig_clump_residuals.png", dpi=200)

    print(f"[*] wrote fig_clump_entropy.pdf, fig_clump_residuals.pdf "
          f"(and .png) — M={M:,}")
    print("    k  C                 H2_hat        resid        sigma       z")
    for i in range(len(k)):
        print(f"    {k[i]}  {int(C[i]):<17d} {H2[i]:.9f}  "
              f"{resid[i]:+.3e}  {err[i]:.3e}  {z[i]:+.2f}")


if __name__ == "__main__":
    main()
