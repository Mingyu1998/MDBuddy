# Using MDBuddy from different AI agents

MDBuddy works with any agent in three ways. Use whichever your agent supports, or several.

| Way | What the agent gets | Works with |
|---|---|---|
| **Skill** (`SKILL.md`, [Agent Skills](https://agentskills.io) format) | the full procedure: what to ask, how to read each build stage, the rules | Claude Code, OpenAI Codex, Gemini CLI and other Agent Skills clients |
| **MCP server** (`mdbuddy mcp`) | 10 typed tools: doctor, write/validate spec, build, lipid lookup/convert/search, orientation guess/viewer, self-assembly check | any MCP client: Claude Code/Desktop, Codex, Gemini CLI, Cursor, VS Code/Copilot, ... |
| **CLI + AGENTS.md** | shell commands with `--json` output, and repository instructions | any agent that can run a shell |

Install MDBuddy first, from a clone (the skill files live in the repository):

```bash
git clone https://github.com/Mingyu1998/mdbuddy && cd mdbuddy
pip install -e ".[build]"      # also installs martinize2, COBY, mdtraj, PDBFixer, matplotlib
```

## 1. Skill

```bash
mdbuddy skill install --agent claude     # ~/.claude/skills/mdbuddy
mdbuddy skill install --agent codex      # ~/.agents/skills/mdbuddy
mdbuddy skill install --agent gemini     # ~/.gemini/skills/mdbuddy
mdbuddy skill install --agent all        # all of the above
mdbuddy skill install --scope project --project /path/to/project   # .claude/, .agents/, .gemini/ skills there
mdbuddy skill install --dest /path/to/your/agents/skills           # any other agent
```

Inside this repository, Codex and Gemini CLI also pick the skill up from `.agents/skills/mdbuddy/`.
Claude Code users can install it as a plugin instead:
`/plugin marketplace add Mingyu1998/MDBuddy` then `/plugin install mdbuddy@mdbuddy`.

Skill locations differ between agents and change between versions. Check your agent's documentation
if it does not find the skill (Codex: `.agents/skills` in the project, user-level skills per the
Codex docs; Gemini CLI: `.gemini/skills` or `.agents/skills`; Claude Code: `.claude/skills`).

## 2. MCP server

`mdbuddy mcp-config <agent>` prints the snippet with the full path of your `mdbuddy`:

| Agent | Configuration |
|---|---|
| Claude Code | `claude mcp add mdbuddy -- mdbuddy mcp` |
| OpenAI Codex | `codex mcp add mdbuddy -- mdbuddy mcp`, or in `~/.codex/config.toml`: `[mcp_servers.mdbuddy]` `command = "mdbuddy"` `args = ["mcp"]` |
| Gemini CLI | `~/.gemini/settings.json`: `{"mcpServers": {"mdbuddy": {"command": "mdbuddy", "args": ["mcp"]}}}` |
| Cursor | `~/.cursor/mcp.json`: `{"mcpServers": {"mdbuddy": {"command": "mdbuddy", "args": ["mcp"]}}}` |
| VS Code / Copilot | `.vscode/mcp.json`: `{"servers": {"mdbuddy": {"type": "stdio", "command": "mdbuddy", "args": ["mcp"]}}}` |
| Claude Desktop | `claude_desktop_config.json`: `{"mcpServers": {"mdbuddy": {"command": "mdbuddy", "args": ["mcp"]}}}` |

The server speaks MCP over stdio (protocol versions 2025-06-18, 2025-03-26, 2024-11-05) using only
the Python standard library. It was tested with the official MCP Python SDK client. Builds can take
a minute or more, so allow a generous tool timeout if your client has one (Codex: `tool_timeout_sec`).

Tools: `doctor`, `write_spec`, `validate_spec`, `build`, `lipids_show`, `lipids_convert`,
`lipids_search`, `orient_guess`, `orient_view`, `selfassembly_check`. An agent that only has MCP still
benefits from the skill text. Point it at `skills/mdbuddy/SKILL.md`, or paste it into its instructions.

## 3. CLI + AGENTS.md

Agents that read `AGENTS.md` (Codex, Cursor, Jules, Gemini CLI with `contextFileName`, ...) get pointed
to the skill when they work in this repository. Everywhere else the CLI is enough. Useful
machine-readable commands:

```bash
mdbuddy doctor --json
mdbuddy validate spec.yaml --json
mdbuddy orient guess protein.pdb --json --html view.html
mdbuddy selfassembly check build_dir --json
mdbuddy lipids convert CDL2 POP2 --from m2 --to m3
```
