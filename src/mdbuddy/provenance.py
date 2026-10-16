"""Provenance: every download, command and output is recorded in mdbuddy.lock.json."""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import platform
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Optional, Sequence

from . import __version__


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


class CommandError(RuntimeError):
    def __init__(self, cmd: Sequence[str], returncode: int, log_path: Path, tail: str):
        self.cmd, self.returncode, self.log_path = list(cmd), returncode, log_path
        super().__init__(
            f"command failed (exit {returncode}): {shlex.join(cmd)}\n"
            f"full log: {log_path}\n--- last lines ---\n{tail}"
        )


class Lock:
    """Accumulates provenance for one build and writes it as JSON."""

    def __init__(self, path: Path):
        self.path = path
        self.data: dict = {
            "mdbuddy_version": __version__,
            "created": now(),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "tools": {},
            "inputs": [],
            "forcefield_files": [],
            "stages": [],
        }

    def stage(self, name: str) -> dict:
        entry = {"stage": name, "started": now(), "commands": [], "outputs": {}, "notes": []}
        self.data["stages"].append(entry)
        self.save()
        return entry

    def finish(self, entry: dict, outputs: Sequence[Path] = (), base: Optional[Path] = None) -> None:
        for p in outputs:
            p = Path(p)
            key = str(p.relative_to(base)) if base else str(p)
            entry["outputs"][key] = sha256(p)
        entry["finished"] = now()
        self.save()

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.data, indent=2) + "\n")

    def run(self, entry: dict, cmd: Sequence[str], cwd: Path, log_name: str,
            env: Optional[dict] = None) -> subprocess.CompletedProcess:
        """Run a command, tee its output to ``cwd/log_name`` and record it."""
        cmd = [str(c) for c in cmd]
        log_path = cwd / log_name
        record = {"cmd": shlex.join(cmd), "cwd": str(cwd), "log": log_name, "started": now()}
        entry["commands"].append(record)
        proc = subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True)
        log_path.write_text(proc.stdout or "")
        record["returncode"] = proc.returncode
        self.save()
        if proc.returncode != 0:
            tail = "\n".join((proc.stdout or "").splitlines()[-30:])
            raise CommandError(cmd, proc.returncode, log_path, tail)
        return proc
