# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Optional per-stage persistence. JSON for pure-Python stages, .sobj for Sage
matrix stages. Checkpoints are a debugging/resumption aid, not an IPC channel.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def save_json(path: str | Path, obj: dict[str, Any]) -> None:
    Path(path).write_text(json.dumps(obj, indent=None, separators=(",", ":")))


def load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text())


def save_sobj(path: str | Path, obj: Any) -> None:
    from sage.all import save

    save(obj, str(path))


def load_sobj(path: str | Path) -> Any:
    from sage.all import load

    return load(str(path))


class CheckpointDir:
    """Lightweight per-run checkpoint directory. Writes are no-ops if path is None."""

    def __init__(self, path: str | Path | None) -> None:
        self.path = Path(path) if path is not None else None
        if self.path is not None:
            self.path.mkdir(parents=True, exist_ok=True)

    def _p(self, name: str) -> Path | None:
        return self.path / name if self.path is not None else None

    def save(self, name: str, obj: Any) -> None:
        p = self._p(name)
        if p is None:
            return
        if name.endswith(".json"):
            save_json(p, obj)
        else:
            save_sobj(p, obj)

    def load(self, name: str) -> Any | None:
        p = self._p(name)
        if p is None or not p.exists():
            return None
        return load_json(p) if name.endswith(".json") else load_sobj(p)
