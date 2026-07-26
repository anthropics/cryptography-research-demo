# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Entry point: sage -python -m hawk_recover.evidence {export,verify} …"""

from __future__ import annotations

import sys


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    match argv:
        case ["export", *rest]:
            from hawk_recover.evidence.export import main as export_main

            return export_main(rest)
        case ["verify", *rest]:
            from hawk_recover.evidence.verify import main as verify_main

            return verify_main(rest)
        case _:
            print(
                "usage: sage -python -m hawk_recover.evidence {export,verify} …\n"
                "  export <workdir> <out_dir> --pk <pk> [--sieve-log <log>]"
                "   run checkpoints → inert evidence\n"
                "  verify <evidence_dir>"
                "   re-derive the claims present in the evidence",
                file=sys.stderr,
            )
            return 2


if __name__ == "__main__":
    sys.exit(main())
