#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""
Plot the wrong-key Hamming-distance histogram from wrongkey.c output.

Reads ../wrongkey_out.txt, extracts the hd=0..12 histogram, and plots
observed counts (bars) against the Binomial(12, 254/255) prediction (line)
on a log-y axis.  Output: fig_wrongkey_hamming.pdf (Type-3-free).
"""

import re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["font.family"]  = "DejaVu Sans"
import matplotlib.pyplot as plt
from math import comb


def parse(path):
    txt = open(path).read()
    m = re.search(r"histogram.*?:\s*\n\s*([0-9 ]+)\n", txt, re.S)
    hist = np.array([int(x) for x in m.group(1).split()], dtype=np.float64)
    M = hist.sum()
    gm = re.search(r"ghost matches.*?:\s*(\d+)", txt)
    return hist, int(M), int(gm.group(1))


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "../wrongkey_out.txt"
    hist, M, ghost = parse(path)
    h = np.arange(13)
    p = 254.0 / 255.0
    expect = np.array([comb(12, k) * p**k * (1 - p)**(12 - k) for k in h]) * M

    # --- figure 1: histogram vs binomial ---------------------------------
    fig, ax = plt.subplots(figsize=(4, 3), constrained_layout=True)
    ax.bar(h, np.maximum(hist, 0.3), width=0.7, color="C0",
           label="observed", zorder=2)
    ax.plot(h, np.maximum(expect, 0.3), "o-", color="k", ms=3.5, lw=1,
            mfc="w", label=r"Binom$(12,\frac{254}{255})$", zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(0.3, M * 2)
    ax.set_xlabel(r"Hamming distance to $I_{\mathrm{off}}$ (bytes)")
    ax.set_ylabel("count")
    ax.set_xticks(range(0, 13))
    ax.legend(frameon=False, loc="upper left", fontsize=8)
    ax.annotate("correct $k_6$", xy=(0, 1), xytext=(2.5, 1.5),
                fontsize=7, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.6))
    fig.savefig("fig_wrongkey_hamming.pdf")
    fig.savefig("fig_wrongkey_hamming.png", dpi=200)

    # --- figure 2: standardized residuals (obs-exp)/sqrt(exp) ------------
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (hist - expect) / np.sqrt(expect)
    mask = expect >= 1.0        # only bins with meaningful expectation
    fig2, axr = plt.subplots(figsize=(4, 3), constrained_layout=True)
    axr.axhline(0, color="0.6", lw=1)
    axr.axhspan(-1, 1, color="0.85", lw=0, zorder=0)
    axr.errorbar(h[mask], z[mask], yerr=np.ones(mask.sum()),
                 fmt="o", ms=4, color="C0", capsize=2, lw=1)
    axr.set_xlabel(r"Hamming distance to $I_{\mathrm{off}}$ (bytes)")
    axr.set_ylabel(r"$(\mathrm{obs}-\mathrm{exp})/\sqrt{\mathrm{exp}}$")
    axr.set_xticks(range(0, 13))
    axr.set_xlim(-0.6, 12.6)
    axr.set_ylim(-3, 3)
    fig2.savefig("fig_wrongkey_residuals.pdf")
    fig2.savefig("fig_wrongkey_residuals.png", dpi=200)

    print(f"[*] wrote fig_wrongkey_hamming.pdf, fig_wrongkey_residuals.pdf "
          f"— M={M:,}, ghost={ghost}")
    print("    h   observed        expected        z")
    for k in h:
        zs = f"{z[k]:+.2f}" if expect[k] >= 1 else "   —"
        print(f"    {k:2d}  {int(hist[k]):>14,}  {expect[k]:14.2f}   {zs}")


if __name__ == "__main__":
    main()
