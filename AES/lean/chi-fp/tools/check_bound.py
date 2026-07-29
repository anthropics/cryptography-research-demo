#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""Numeric sanity check for the closed inequality used in the proof of `chi_fp`.

Bsum := sum_{k=0}^{256} C(256, k) * |256 - 2k|^255

Claim (Step (iv)):  256 * 255 * Bsum * 2^200 <= 256^255 * 2^256.
"""
from math import comb, log2

n, q = 255, 256
Bsum = sum(comb(q, k) * abs(q - 2 * k) ** n for k in range(q + 1))
lhs = 256 * 255 * Bsum * 2 ** 200
rhs = 256 ** 255 * 2 ** 256
print(f"log2(Bsum)       = {log2(Bsum):.4f}")
print(f"Bsum             = {Bsum}")
print(f"inequality holds = {lhs <= rhs}")
print(f"slack (bits)     = {log2(rhs) - log2(lhs):.4f}")
assert lhs <= rhs
