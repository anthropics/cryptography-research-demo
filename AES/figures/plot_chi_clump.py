#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""
Figure for the χ* collision-entropy appendix (chi_clump.c output).

Two panels:
  (a) Ĥ₂(k) vs k, against the naive 8k line and the bit-bias prediction
      8k·h₁ where h₁ = -log₂(p²+(1-p)²), p = (1-(127/128)^255)/2.
  (b) residuals (Ĥ₂(k) - 8k·h₁)/σ_k, ±1σ band.
"""
import numpy as np, math
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["font.family"]  = "DejaVu Sans"
import matplotlib.pyplot as plt

d = np.genfromtxt("../chi_clump_results.tsv", names=True, dtype=None)
k     = d["k"].astype(int)
C     = d["collisions"].astype(float)
H2    = d["H2_hat"]
M     = float(d["M"][0])
sigma = 1.0/(np.sqrt(C)*np.log(2))

# bit-bias model (i.i.d. Bernoulli(p) bits — first order)
p  = (1 - (127/128)**255)/2
h1 = -math.log2(p*p + (1-p)*(1-p))        # ≈ 0.97375
model = 8*h1*k

# exact multinomial-parity model for nb fixed bits of χ (closed form):
#   P(match) = 2^{-nb} Σ_j C(nb,j) ((256-2j)/256)^{510}
from math import comb
def exact_H2(nb):
    s = sum(comb(nb, j) * ((256 - 2*j)/256)**510 for j in range(nb+1))
    return nb - math.log2(s)
exact = np.array([exact_H2(8*ki) for ki in k])

mask  = C >= 10
rate  = H2/k
resid = H2 - exact            # residual vs the EXACT multinomial model
z     = resid/sigma

fig,(a,b)=plt.subplots(1,2,figsize=(8,3),constrained_layout=True)

# (a) Ĥ₂(k) vs k, against the exact multinomial-parity model
a.plot(k[mask], exact[mask], "-", color="k", lw=1,
       label="multinomial-parity model")
a.errorbar(k[mask], H2[mask], yerr=sigma[mask], fmt="o", ms=4, color="C0",
           capsize=2, lw=1, label=r"measured $\hat H_2(k)$")
a.set_xlabel(r"prefix length $k$ (bytes)")
a.set_ylabel(r"$\hat H_2(k)$  (bits)")
a.set_xticks(k[mask])
a.legend(frameon=False, fontsize=8, loc="upper left")
a.text(0.98,0.04,
       fr"model: $H_2(96)={exact_H2(96):.2f}$, $H_2(104)={exact_H2(104):.2f}$",
       transform=a.transAxes,ha="right",va="bottom",fontsize=7)

# (b) residual vs the exact model, in BITS (no /σ_k)
b.axhline(0, color="k", lw=1)
b.errorbar(k[mask], resid[mask], yerr=sigma[mask],
           fmt="o", ms=4, color="C0", capsize=2, lw=1)
b.set_xlabel(r"prefix length $k$ (bytes)")
b.set_ylabel(r"$\hat H_2(k)-H_2^{\mathrm{model}}(k)$  (bits)")
b.set_xticks(k[mask])

fig.savefig("fig_chi_clump.pdf")
fig.savefig("fig_chi_clump.png",dpi=200)

print(f"[*] M={M:.0f}  p={p:.6f}  h1={h1:.6f}  8h1={8*h1:.4f}")
print(" k         C        H2_hat    iid      exact   H2-exact      z")
for i in range(len(k)):
    z = (H2[i]-exact[i])/sigma[i]
    print(f" {k[i]}  {C[i]:12.4e}  {H2[i]:8.4f}  {model[i]:7.4f}  {exact[i]:7.4f}  "
          f"{H2[i]-exact[i]:+.6f}  {z:+8.1f}")
print(f"    exact_H2(96) = {exact_H2(96):.3f}  (=> 12 bytes)")
print(f"    exact_H2(104)= {exact_H2(104):.3f}  (=> 13 bytes)")
print("    wrote fig_chi_clump.pdf/.png")
