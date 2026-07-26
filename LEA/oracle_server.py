#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Local chosen-plaintext encryption oracle for round-reduced LEA-128.

Serves   POST /encrypt   {"pt": "<hex>"}  ->  {"ct": "<hex>"}
where the plaintext is any multiple of 16 bytes (ECB over 16-byte blocks).

Usage:
    python3 oracle_server.py [--key <32 hex chars>] [--rounds 13] [--port 8666]

If --key is omitted, a random key is generated and printed, so the attack's
output can be checked against it afterwards. Requires numpy (the attack
requests about 2^28.6 blocks in total; a pure-Python implementation would be
far too slow to serve them).
"""

import argparse
import json
import os
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import numpy as np

MASK = 0xFFFFFFFF
DELTA = (0xC3EFE9DB, 0x44626B02, 0x79E27C8A, 0x78DF30EC)


def _rol(x, r):
    r %= 32
    return ((x << r) | (x >> (32 - r))) & MASK if r else x


def key_schedule(key16: bytes, rounds: int):
    t = [int.from_bytes(key16[4 * i : 4 * i + 4], "little") for i in range(4)]
    rk = []
    for i in range(rounds):
        d = DELTA[i & 3]
        t[0] = _rol((t[0] + _rol(d, i)) & MASK, 1)
        t[1] = _rol((t[1] + _rol(d, i + 1)) & MASK, 3)
        t[2] = _rol((t[2] + _rol(d, i + 2)) & MASK, 6)
        t[3] = _rol((t[3] + _rol(d, i + 3)) & MASK, 11)
        rk.append(tuple(t))
    return rk


def _rol_np(x, r):
    r %= 32
    if r == 0:
        return x
    return (x << np.uint32(r)) | (x >> np.uint32(32 - r))


def encrypt(pt: bytes, rk, rounds: int) -> bytes:
    """ECB-encrypt len(pt)/16 blocks, vectorized across blocks."""
    n = len(pt) // 16
    w = np.frombuffer(pt, dtype="<u4").reshape(n, 4)
    x0, x1, x2, x3 = (w[:, i].copy() for i in range(4))
    for r in range(rounds):
        t0, t1, t2, t3 = (np.uint32(v) for v in rk[r])
        y0 = _rol_np((x0 ^ t0) + (x1 ^ t1), 9)
        y1 = _rol_np((x1 ^ t2) + (x2 ^ t1), 32 - 5)
        y2 = _rol_np((x2 ^ t3) + (x3 ^ t1), 32 - 3)
        x0, x1, x2, x3 = y0, y1, y2, x0
    out = np.empty((n, 4), dtype="<u4")
    out[:, 0], out[:, 1], out[:, 2], out[:, 3] = x0, x1, x2, x3
    return out.tobytes()


def main():
    ap = argparse.ArgumentParser(description="Local LEA-128 encryption oracle")
    ap.add_argument("--key", help="32 hex chars (16-byte key); random if omitted")
    ap.add_argument("--rounds", type=int, default=13)
    ap.add_argument("--port", type=int, default=8666)
    args = ap.parse_args()

    if args.key:
        try:
            key = bytes.fromhex(args.key)
        except ValueError:
            sys.exit("--key must be exactly 32 hex characters")
        if len(key) != 16:
            sys.exit("--key must be exactly 32 hex characters")
    else:
        key = os.urandom(16)
    rk = key_schedule(key, args.rounds)
    stats = {"queries": 0, "blocks": 0}
    stats_lock = threading.Lock()

    print(f"[oracle] key    = {key.hex()}  (secret; for checking the attack's output)")
    print(f"[oracle] rounds = {args.rounds}")
    print(f"[oracle] serving on 127.0.0.1:{args.port}", flush=True)

    class Handler(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"

        def log_message(self, *a):
            pass

        def _send(self, code, obj):
            body = json.dumps(obj).encode()
            self.send_response(code)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_POST(self):
            if self.path != "/encrypt":
                self._send(404, {"error": "not found"})
                return
            try:
                n = int(self.headers.get("Content-Length", 0))
                if n > 2 * 1024 * 1024:
                    raise ValueError("request too large")
                pt = bytes.fromhex(json.loads(self.rfile.read(n))["pt"])
                if len(pt) == 0 or len(pt) % 16:
                    raise ValueError(
                        "plaintext must be a positive multiple of 16 bytes"
                    )
            except Exception as e:
                self._send(400, {"error": str(e)})
                return
            ct = encrypt(pt, rk, args.rounds)
            with stats_lock:
                stats["queries"] += 1
                stats["blocks"] += len(pt) // 16
            self._send(200, {"ct": ct.hex()})

    srv = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print(f"\n[oracle] served {stats['queries']} queries, {stats['blocks']} blocks")


if __name__ == "__main__":
    main()
