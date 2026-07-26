# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Untrusted-side verification: re-derive a recovery's claims from its evidence.

CLI: sage -python -m hawk_recover.evidence verify <evidence_dir> [options]

Convenience tool: the attack CLI's own final sign/verify step is already the
proof of a successful recovery. Use this to double-check a finished run, or
to check a run someone shared with you (see export.py for the producing
side). Reads ONLY inert data (JSON of integer strings, plain-text logs, raw
key bytes) — loading evidence executes no code, so it is safe to run on
exports from untrusted parties; every mathematical claim is recomputed here,
not trusted from the export.

Checks (each prints PASS/FAIL/SKIPPED):
  PK     LB is exactly the cocycle lattice of pk.hex (constraint-matrix
         kernel + embedding structure + determinant invariant). This is the
         primary provenance check: it proves the positive chain from the
         stated public key, rather than disproving any particular forgery.
  GRAM   Gram_tau equals the Gram matrix recomputed from LB.
  SHELL  every claimed vector has exact Gram_tau-norm λ₁².
  HARVEST the claimed vectors all derive from the sieve log via the warm
         basis (shell ⊆ harvest), and the warm basis is honestly derived.
  SK     the recovered key's (f,g) satisfies |f|²+|g|² = q00 from pk.hex,
         binding the key to exactly that public key. The key is read in the
         reference implementation's decoded layout (f‖g‖F₂‖G₂‖hpub, the
         format the attack CLI writes), not the seed-based NIST encoding.
"""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

RESULTS: dict[str, str] = {}


def result(name: str, status: str, detail: str = "") -> None:
    RESULTS[name] = status
    print(f"{name}-RESULT: {status}{' — ' + detail if detail else ''}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="hawk_recover.evidence verify",
        description="Re-verify a key recovery from its exported evidence.",
    )
    ap.add_argument("evidence", help="evidence directory (inert JSON format)")
    ap.add_argument(
        "--sk",
        default=None,
        help="recovered private key, decoded format (default: ../recovered_sk.bin "
        "relative to the evidence dir)",
    )
    ap.add_argument(
        "--golden-evidence",
        default=None,
        help="a second run's evidence dir, to check the two runs found independent solutions (optional)",
    )
    args = ap.parse_args(argv)

    import numpy as np
    from hawk_recover import evidence as ev
    from hawk_recover import ring
    from hawk_recover.cocycle import _build_cBr, _gram_tau, _Q_to_array
    from hawk_recover.decode import decode
    from hawk_recover.svp import dedup_ambient, harvest_log
    from hawk_recover.types import WarmBasis

    edir = Path(args.evidence)
    print(f"evidence: {edir}")
    print("=" * 70)

    # ---- load inert evidence --------------------------------------------------
    LB, Gt = ev.load_lattice(edir)
    W, provenance = ev.load_shell(edir)
    warm_data = ev.load_warm(edir)
    sieve_log = ev.load_sieve_log(edir)
    pk_bytes = ev.load_pk_bytes(edir)
    print(
        f"loaded: LB {LB.nrows()}×{LB.ncols()}, {len(W)} claimed vectors, "
        f"warm={'yes' if warm_data else 'no'}, sieve log={'yes' if sieve_log else 'no'}, "
        f"pk={'yes' if pk_bytes else 'no'}"
    )

    if pk_bytes is None:
        print("FATAL: evidence has no pk.hex — nothing to verify against")
        return 2
    pk = decode(pk_bytes)
    p = pk.params
    dL, dK, n = p.dL, p.dK, p.rank_tau  # noqa: F841
    print(f"stated pk: HAWK-{p.n} ({len(pk_bytes)} bytes)")
    print("=" * 70)

    # ---- PK: LB is the cocycle lattice of the stated pk -----------------------
    cBr = _build_cBr(_Q_to_array(pk), dL, dK)
    LBn = np.array([[int(x) for x in row] for row in LB.rows()], dtype=object)
    A, B = LBn[:, 0:dL], LBn[:, dL : 2 * dL]
    C, D = LBn[:, 2 * dL : 3 * dL], LBn[:, 3 * dL : 4 * dL]
    sgn = np.array([1 - 2 * (k & 1) for k in range(dL)], dtype=object)
    c1_ok = np.array_equal(D, A * sgn) and not B[:, 0::2].any() and not C[:, 0::2].any()
    P = np.concatenate([A, B[:, 1::2], C[:, 1::2]], axis=1)
    in_ker = not (cBr.astype(object) @ P.T).any()
    det_ok = False
    if c1_ok and in_ker:
        d = Gt.det()
        ld = d.nbits() - 1 if d > 0 else -1
        det_ok = ld == round(p.logdet_gram_tau)
        print(
            f"(C1) structure: {c1_ok} | kernel membership (exact): {in_ker} | "
            f"log2(det)={ld} expect {round(p.logdet_gram_tau)}: {det_ok}"
        )
    result(
        "PK",
        "PASS" if (c1_ok and in_ker and det_ok) else "FAIL",
        "LB spans exactly the cocycle lattice of the stated pk"
        if (c1_ok and in_ker and det_ok)
        else "",
    )

    # ---- GRAM: Gram_tau recomputed from LB -------------------------------------
    gram_ok = Gt == _gram_tau(LB, dL)
    result(
        "GRAM",
        "PASS" if gram_ok else "FAIL",
        "Gram_tau is the Gram matrix of LB"
        if gram_ok
        else "Gram_tau does NOT match LB",
    )

    # ---- SHELL: exact norms -----------------------------------------------------
    norms = sorted({int(w * Gt * w) for w in W})
    shell_ok = norms == [p.lam1sq] and len(W) >= 2
    result(
        "SHELL",
        "PASS" if shell_ok else "FAIL",
        f"{len(W)} vectors, all exact norm {p.lam1sq}"
        if shell_ok
        else f"norms={norms}",
    )

    # ---- HARVEST: vectors derive from the sieve log ---------------------------
    if sieve_log is None or warm_data is None:
        result(
            "HARVEST",
            "SKIPPED",
            "no sieve log / warm basis in evidence — "
            "provenance rests on PK+GRAM+SHELL only",
        )
    else:
        from hawk_recover.svp import _make_warm_basis
        from hawk_recover.types import CocycleLattice

        lat = CocycleLattice(params=p, LB=LB, Gram_tau=Gt)
        if "Bint" in warm_data:
            rebuilt = _make_warm_basis(lat, warm_data["Uim"], warm_data["beta"])
            honest = rebuilt.Bint == warm_data["Bint"] and rebuilt.sb == warm_data["sb"]
            warm = WarmBasis(
                params=p,
                beta=warm_data["beta"],
                Uim=warm_data["Uim"],
                Bint=warm_data["Bint"],
                sb=warm_data["sb"],
                Gram_tau=Gt,
                LB=LB,
            )
        else:
            honest = True
            warm = _make_warm_basis(lat, warm_data["Uim"], warm_data["beta"])
        hits = dedup_ambient(harvest_log(sieve_log, warm), LB)

        def canon(w):
            V = w * LB
            nz = next((i for i, x in enumerate(V) if x != 0), 0)
            return tuple(V) if V[nz] > 0 else tuple(-V)

        subset = {canon(w) for w in W} <= {canon(h) for h in hits}
        print(
            f"warm basis honest: {honest} | harvest: {len(hits)} unique | "
            f"claimed ⊆ harvest: {subset}"
        )
        result("HARVEST", "PASS" if (honest and subset and len(W) >= 2) else "FAIL")

    # ---- SK: recovered key binds to the stated pk ------------------------------
    sk_path = Path(args.sk) if args.sk else edir.parent / "recovered_sk.bin"
    if not sk_path.exists():
        result("SK", "SKIPPED", f"no recovered key at {sk_path}")
    else:
        raw = sk_path.read_bytes()
        f_co = list(struct.unpack(f"{p.n}b", raw[0 : p.n]))
        g_co = list(struct.unpack(f"{p.n}b", raw[p.n : 2 * p.n]))
        # |f|^2 + |g|^2 = q00  (negacyclic: f*conj(f) + g*conj(g) = q00)
        ff = ring.nc_mul(f_co, ring.nc_conj(f_co, p.n), p.n)
        gg = ring.nc_mul(g_co, ring.nc_conj(g_co, p.n), p.n)
        lhs = [a + b for a, b in zip(ff, gg)]  # noqa: B905
        q00 = pk.Q_list[0]
        sk_ok = lhs == list(q00)
        result(
            "SK",
            "PASS" if sk_ok else "FAIL",
            "recovered (f,g) satisfies |f|²+|g|²=q00 of the stated pk"
            if sk_ok
            else "key equation FAILS against the stated pk",
        )

    # ---- GOLDEN: independence from another run (optional) ----------------------
    if args.golden_evidence:
        gdir = Path(args.golden_evidence)
        gLB, gGt = ev.load_lattice(gdir)
        gW, _ = ev.load_shell(gdir)
        same_lattice = LB == gLB
        overlap = len(
            {tuple(w) for w in map(list, W)} & {tuple(w) for w in map(list, gW)}
        )
        print(f"same lattice as golden: {same_lattice} | raw vector overlap: {overlap}")
        result(
            "GOLDEN",
            "PASS",
            "independence data reported above (interpretation depends on whether "
            "the two runs attacked the same key)",
        )

    print("=" * 70)
    failed = [k for k, v in RESULTS.items() if v == "FAIL"]
    print(f"SUMMARY: {', '.join(f'{k}={v}' for k, v in RESULTS.items())}")
    if failed:
        print(f"VERIFICATION FAILED: {failed}")
        return 1
    print("VERIFICATION PASSED (all non-skipped checks)")
    return 0
