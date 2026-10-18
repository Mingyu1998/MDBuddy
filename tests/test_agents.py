"""Agent portability: MCP protocol, skill mirror, Agent Skills frontmatter."""

import filecmp
import json
import re
from pathlib import Path

import yaml

from mdbuddy import mcp_server

ROOT = Path(__file__).resolve().parents[1]


def call(method, params=None, mid=1):
    return mcp_server.handle({"jsonrpc": "2.0", "id": mid, "method": method, "params": params or {}})


def test_mcp_handshake_and_tools():
    init = call("initialize", {"protocolVersion": "2025-03-26", "capabilities": {}})["result"]
    assert init["protocolVersion"] == "2025-03-26" and init["capabilities"]["tools"]
    assert call("initialize", {"protocolVersion": "1999-01-01"})["result"]["protocolVersion"] == "2025-06-18"
    assert mcp_server.handle({"jsonrpc": "2.0", "method": "notifications/initialized"}) is None
    tools = call("tools/list")["result"]["tools"]
    names = {t["name"] for t in tools}
    assert {"doctor", "build", "write_spec", "lipids_convert", "orient_guess", "selfassembly_check"} <= names
    assert all(t["inputSchema"]["type"] == "object" for t in tools)
    assert call("nope")["error"]["code"] == -32601


def test_mcp_tool_calls(tmp_path):
    res = call("tools/call", {"name": "lipids_convert", "arguments": {"names": ["CDL2"], "source": "m2", "target": "m3"}})
    data = json.loads(res["result"]["content"][0]["text"])
    assert data["CDL2"]["candidates"][0]["name"] == "TOCL" and res["result"]["isError"] is False
    bad = call("tools/call", {"name": "write_spec", "arguments": {"spec_path": str(tmp_path / "s.yaml"),
                                                                   "yaml_text": "name: x\nprotein: {source: opm}\n"}})
    assert "protein.id" in bad["result"]["content"][0]["text"] and not (tmp_path / "s.yaml").exists()
    good = call("tools/call", {"name": "write_spec", "arguments": {"spec_path": str(tmp_path / "s.yaml"),
                                                                    "yaml_text": "name: x\nprotein: {source: opm, id: 1bl8}\n"}})
    assert json.loads(good["result"]["content"][0]["text"])["written"] and (tmp_path / "s.yaml").exists()
    wrong = call("tools/call", {"name": "validate_spec", "arguments": {"nope": 1}})
    assert wrong["result"]["isError"] is True


def test_skill_mirror_is_identical():
    a, b = ROOT / "skills" / "mdbuddy", ROOT / ".agents" / "skills" / "mdbuddy"
    files_a = sorted(p.relative_to(a) for p in a.rglob("*") if p.is_file())
    files_b = sorted(p.relative_to(b) for p in b.rglob("*") if p.is_file())
    assert files_a == files_b, "run: python tools/sync_skill.py"
    for f in files_a:
        assert filecmp.cmp(a / f, b / f, shallow=False), f"{f} differs: run python tools/sync_skill.py"


def test_skill_frontmatter_follows_agent_skills_spec():
    text = (ROOT / "skills" / "mdbuddy" / "SKILL.md").read_text()
    front = yaml.safe_load(text.split("---")[1])
    assert front["name"] == "mdbuddy" and re.fullmatch(r"[a-z0-9]+(-[a-z0-9]+)*", front["name"])
    assert 1 <= len(front["description"]) <= 1024
