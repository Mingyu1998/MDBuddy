"""Interactive membrane-placement viewer: one self-contained HTML file (3Dmol.js from a CDN).

The page shows the protein (input frame) transformed by the current placement, the hydrophobic
slab, and sliders to tilt/rotate/shift it. It prints the spec lines for ``orientation: manual``
(membrane normal and centre in the input frame), which ``orient.manual`` applies exactly.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from .structure import Atom

THREEDMOL_URL = "https://cdn.jsdelivr.net/npm/3dmol@2.5.4/build/3Dmol-min.js"

PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Membrane Placement Viewer</title>
<style>
  :root { --bg:#f7f7f5; --panel:#ffffff; --text:#1d1d1b; --muted:#6b6b66; --line:#d9d9d4; --accent:#c2410c; }
  @media (prefers-color-scheme: dark) { :root { --bg:#151514; --panel:#1f1f1d; --text:#ecece8; --muted:#a3a39c; --line:#34342f; --accent:#fb923c; } }
  * { box-sizing: border-box; }
  body { margin:0; background:var(--bg); color:var(--text); font:14px/1.45 system-ui, sans-serif; }
  header { padding:12px 16px; border-bottom:1px solid var(--line); }
  header h1 { font-size:16px; margin:0; } header p { margin:2px 0 0; color:var(--muted); }
  main { display:grid; grid-template-columns: minmax(0,1fr) 340px; min-height: calc(100vh - 60px); }
  #view { position:relative; min-height:420px; }
  aside { background:var(--panel); border-left:1px solid var(--line); padding:14px 16px; overflow:auto; }
  label { display:block; margin:10px 0 2px; font-weight:600; }
  .row { display:flex; gap:8px; align-items:center; } .row input[type=range] { flex:1; }
  output { width:56px; text-align:right; font-variant-numeric: tabular-nums; color:var(--muted); }
  button { margin:10px 6px 0 0; padding:6px 10px; border:1px solid var(--line); border-radius:6px; background:var(--bg); color:var(--text); cursor:pointer; }
  pre { background:var(--bg); border:1px solid var(--line); border-radius:6px; padding:8px; white-space:pre-wrap; font-size:12px; }
  .legend span { display:inline-block; width:10px; height:10px; border-radius:2px; margin:0 4px 0 10px; vertical-align:middle; }
  .note { color:var(--muted); font-size:12px; }
  #err { color:var(--accent); padding:16px; display:none; }
  @media (max-width: 760px) { main { grid-template-columns: 1fr; } aside { border-left:0; border-top:1px solid var(--line); } }
</style>
</head>
<body>
<header>
  <h1>Membrane placement: __TITLE__</h1>
  <p>Membrane normal is vertical (+z up = non-cytoplasmic side). Grey slab = hydrophobic core.</p>
</header>
<main>
  <div id="view"><div id="err">Could not load 3Dmol.js (__URL__). Check your internet connection.</div></div>
  <aside>
    <div class="legend">Residues:<span style="background:#e0a100"></span>hydrophobic<span style="background:#4a8bdf"></span>Lys/Arg<span style="background:#d6453d"></span>Asp/Glu<span style="background:#2f9e58"></span>Trp/Tyr</div>
    <label>Tilt about x (deg)</label><div class="row"><input id="rx" type="range" min="-180" max="180" step="1" value="0"><output id="rxo">0</output></div>
    <label>Tilt about y (deg)</label><div class="row"><input id="ry" type="range" min="-180" max="180" step="1" value="0"><output id="ryo">0</output></div>
    <label>Spin about z (deg)</label><div class="row"><input id="rz" type="range" min="-180" max="180" step="1" value="0"><output id="rzo">0</output></div>
    <label>Depth: move protein along z (&Aring;)</label><div class="row"><input id="dz" type="range" min="-40" max="40" step="0.5" value="0"><output id="dzo">0</output></div>
    <label>Hydrophobic thickness (&Aring;, display only)</label><div class="row"><input id="th" type="range" min="18" max="46" step="1" value="__THICK__"><output id="tho">__THICK__</output></div>
    <button id="reset">Reset to MDBuddy's guess</button><button id="flip">Flip in/out</button>
    <label>Spec lines (orientation: manual)</label>
    <pre id="yaml"></pre>
    <button id="copy">Copy</button>
    <p class="note">__NOTE__</p>
  </aside>
</main>
<script>window.MDBUDDY = __DATA__;</script>
<script src="__URL__"></script>
<script>
(function () {
  const D = window.MDBUDDY;
  const ids = ["rx", "ry", "rz", "dz", "th"];
  const el = Object.fromEntries(ids.map(k => [k, document.getElementById(k)]));
  const HYD = new Set(["ALA","VAL","LEU","ILE","MET","PHE","CYS","PRO"]);
  function matmul(A, B) { return A.map((r, i) => [0,1,2].map(j => r[0]*B[0][j] + r[1]*B[1][j] + r[2]*B[2][j])); }
  function euler(rx, ry, rz) {   // R = Rz . Ry . Rx, degrees (same as mdbuddy.orient.euler_xyz)
    const [a, b, g] = [rx, ry, rz].map(v => v * Math.PI / 180);
    const Rx = [[1,0,0],[0,Math.cos(a),-Math.sin(a)],[0,Math.sin(a),Math.cos(a)]];
    const Ry = [[Math.cos(b),0,Math.sin(b)],[0,1,0],[-Math.sin(b),0,Math.cos(b)]];
    const Rz = [[Math.cos(g),-Math.sin(g),0],[Math.sin(g),Math.cos(g),0],[0,0,1]];
    return matmul(Rz, matmul(Ry, Rx));
  }
  function state() {
    const R = matmul(euler(+el.rx.value, +el.ry.value, +el.rz.value), D.rotation);
    return { R, shift: D.shift + (+el.dz.value) };
  }
  // membrane normal and centre in the input frame: x' = R (x - origin) + (0,0,shift)
  function membrane(st) {
    const n = [st.R[2][0], st.R[2][1], st.R[2][2]];                    // R^T z
    const c = [0,1,2].map(i => D.origin[i] - st.shift * st.R[2][i]);    // origin + R^T (0,0,-shift)
    return { n, c };
  }
  const f = v => (Math.abs(v) < 5e-5 ? 0 : v).toFixed(4);
  function yaml(st) {
    const m = membrane(st);
    return "protein:\n  orientation: manual\n  manual:\n" +
      "    normal: [" + m.n.map(f).join(", ") + "]\n" +
      "    center_point: [" + m.c.map(v => v.toFixed(3)).join(", ") + "]   # Angstrom, input frame\n";
  }
  function pdbText(st) {
    const out = [];
    for (let i = 0; i < D.lines.length; i++) {
      const p = D.xyz[i], q = [p[0] - D.origin[0], p[1] - D.origin[1], p[2] - D.origin[2]];
      const x = st.R[0][0]*q[0] + st.R[0][1]*q[1] + st.R[0][2]*q[2];
      const y = st.R[1][0]*q[0] + st.R[1][1]*q[1] + st.R[1][2]*q[2];
      const z = st.R[2][0]*q[0] + st.R[2][1]*q[1] + st.R[2][2]*q[2] + st.shift;
      const s = D.lines[i];
      out.push(s.slice(0, 30) + x.toFixed(3).padStart(8) + y.toFixed(3).padStart(8) + z.toFixed(3).padStart(8) + s.slice(54));
    }
    return out.join("\n");
  }
  let viewer = null, first = true;
  function draw() {
    ids.forEach(k => document.getElementById(k + "o").textContent = el[k].value);
    const st = state();
    document.getElementById("yaml").textContent = yaml(st);
    if (!viewer) return;
    viewer.removeAllModels(); viewer.removeAllShapes();
    const m = viewer.addModel(pdbText(st), "pdb");
    m.setStyle({}, { cartoon: { color: "#c9c9c4", opacity: 0.85 } });
    m.setStyle({ resn: Array.from(HYD) }, { cartoon: { color: "#e0a100" } });
    m.setStyle({ resn: ["LYS","ARG"] }, { cartoon: { color: "#4a8bdf" }, stick: { color: "#4a8bdf", radius: 0.25 } });
    m.setStyle({ resn: ["ASP","GLU"] }, { cartoon: { color: "#d6453d" } });
    m.setStyle({ resn: ["TRP","TYR"] }, { cartoon: { color: "#2f9e58" }, stick: { color: "#2f9e58", radius: 0.25 } });
    const half = (+el.th.value) / 2, w = D.extent * 1.7;
    viewer.addBox({ center: { x: 0, y: 0, z: 0 }, dimensions: { w: w, h: w, d: 2 * half }, color: "#9a9a94", opacity: 0.12 });
    [half, -half].forEach(z => viewer.addBox({ center: { x: 0, y: 0, z: z }, dimensions: { w: w, h: w, d: 0.4 },
                                               color: z > 0 ? "#c2410c" : "#4a8bdf", opacity: 0.55 }));
    if (first) { viewer.zoomTo(); viewer.rotate(-80, "x"); first = false; }
    viewer.render();
  }
  ids.forEach(k => el[k].addEventListener("input", draw));
  document.getElementById("reset").onclick = () => { ["rx","ry","rz","dz"].forEach(k => el[k].value = 0); el.th.value = D.thickness; draw(); };
  document.getElementById("flip").onclick = () => { el.rx.value = ((+el.rx.value + 360) % 360) - 180; draw(); };
  document.getElementById("copy").onclick = () => navigator.clipboard && navigator.clipboard.writeText(document.getElementById("yaml").textContent);
  window.mdbuddyState = () => ({ ...membrane(state()), yaml: yaml(state()) });
  if (window.$3Dmol) {
    viewer = $3Dmol.createViewer(document.getElementById("view"), { backgroundColor: getComputedStyle(document.body).backgroundColor });
  } else {
    document.getElementById("err").style.display = "block";
  }
  draw();
})();
</script>
</body>
</html>
"""


def write_html(path: Path, atoms: list[Atom], record: dict, thickness_A: float = 30.0,
               title: str = "protein", note: Optional[str] = None) -> Path:
    """Write the viewer for ``atoms`` (input frame) starting from placement ``record``."""
    import math

    protein = [a for a in atoms if a.record in ("ATOM", "HETATM")]
    rotation = record.get("rotation") or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    if "origin_A" in record:
        origin = record["origin_A"]
    else:
        n = max(len(protein), 1)
        origin = [sum(a.x for a in protein) / n, sum(a.y for a in protein) / n, sum(a.z for a in protein) / n]
    extent = max((math.dist((a.x, a.y, a.z), origin) for a in protein), default=30.0)
    data = {
        "lines": [a.to_pdb_line(serial=i + 1) for i, a in enumerate(protein)],
        "xyz": [[round(a.x, 3), round(a.y, 3), round(a.z, 3)] for a in protein],
        "rotation": rotation, "origin": origin, "shift": float(record.get("z_shift_A", 0.0)),
        "thickness": round(float(thickness_A)), "extent": extent,
    }
    placement = record.get("placement") or {}
    note = note or (f"Starting point: {record.get('method', 'input frame')}"
                    + (f" - {placement['verdict']}." if placement.get("verdict") else ".")
                    + " Paste the spec lines into your spec and rebuild with --from orient.")
    html = (PAGE.replace("__DATA__", json.dumps(data, separators=(",", ":")))
            .replace("__URL__", THREEDMOL_URL).replace("__TITLE__", title)
            .replace("__THICK__", str(data["thickness"])).replace("__NOTE__", note))
    path.write_text(html)
    return path
