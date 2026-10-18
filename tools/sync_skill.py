#!/usr/bin/env python3
"""Mirror skills/mdbuddy (canonical; Claude Code plugin layout) to .agents/skills/mdbuddy.

.agents/skills/ is the cross-agent location read by Codex, Gemini CLI and other Agent Skills
clients when they work inside this repository. A test fails if the two copies differ.
"""
import shutil
from pathlib import Path

root = Path(__file__).resolve().parents[1]
src, dst = root / "skills" / "mdbuddy", root / ".agents" / "skills" / "mdbuddy"
if dst.exists():
    shutil.rmtree(dst)
shutil.copytree(src, dst)
print(f"synced {src} -> {dst}")
