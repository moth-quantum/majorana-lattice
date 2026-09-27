"""Build explorer/explorer.html: the template with the JavaScript engine and
lattice data inlined, so the page works from a single file.

    python explorer/build.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from majorana_lattice.export import export  # noqa: E402

SIZES = ["star:4x4", "star:6x4", "star:6x6"]


def build() -> str:
    template = (ROOT / "explorer" / "template.html").read_text()
    engine = (ROOT / "js" / "engine.js").read_text()
    lattices = json.dumps(export(SIZES), separators=(",", ":"))
    return template.replace("/*ENGINE*/", engine).replace("/*LATTICES*/", lattices)


if __name__ == "__main__":
    out = ROOT / "explorer" / "explorer.html"
    out.write_text(build())
    print(f"wrote {out.relative_to(ROOT)}")
