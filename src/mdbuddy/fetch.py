"""Get the input structure: OPM (pre-oriented in the membrane), RCSB, or a local file."""

from __future__ import annotations

import shutil
import urllib.error
import urllib.request
from pathlib import Path

from .provenance import sha256

OPM_URL = "https://opm-assets.storage.googleapis.com/pdb/{id}.pdb"
RCSB_URL = "https://files.rcsb.org/download/{ID}.pdb"
USER_AGENT = "MDBuddy (https://github.com/Mingyu1998/MDBuddy)"


class FetchError(RuntimeError):
    pass


def download(url: str, dest: Path, timeout: float = 60.0) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp, open(dest, "wb") as out:
            shutil.copyfileobj(resp, out)
    except urllib.error.HTTPError as exc:
        dest.unlink(missing_ok=True)
        raise FetchError(f"HTTP {exc.code} while downloading {url}") from exc
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        dest.unlink(missing_ok=True)
        raise FetchError(f"could not download {url}: {exc}") from exc


def source_url(source: str, pdb_id: str) -> str:
    if source == "opm":
        return OPM_URL.format(id=pdb_id.lower())
    if source == "rcsb":
        return RCSB_URL.format(ID=pdb_id.upper())
    raise ValueError(source)


def fetch_structure(spec, out_dir: Path) -> dict:
    """Place the raw input at ``out_dir/input.pdb`` and return a provenance record."""
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / "input.pdb"
    p = spec.protein
    if p["source"] == "file":
        src = spec.resolve_path(p["path"])
        if not src.exists():
            raise FetchError(f"protein.path not found: {src}")
        if src.suffix.lower() not in (".pdb", ".ent"):
            raise FetchError(f"only PDB files are supported for now (got {src.name}); "
                             "convert mmCIF with e.g. `gmx editconf` or PyMOL first")
        shutil.copyfile(src, dest)
        origin = {"source": "file", "path": str(src)}
    else:
        url = source_url(p["source"], p["id"])
        download(url, dest)
        origin = {"source": p["source"], "id": p["id"], "url": url}
        if dest.stat().st_size < 200:
            raise FetchError(f"{url} returned an almost empty file; is '{p['id']}' in {p['source'].upper()}?")
    origin["file"] = str(dest.name)
    origin["sha256"] = sha256(dest)
    return origin
