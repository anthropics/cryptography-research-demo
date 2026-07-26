# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Infrastructure: everything that is not the attack's mathematics.

Sieve backend processes, checkpointing, and the sign/verify oracle (shim).
The attack modules (decode → cocycle → svp → intertwiner → tower → encode)
live in the package root; the evidence export/verify tools in `evidence/`.
"""
