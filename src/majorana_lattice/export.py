"""Write lattice presets as JSON for the JavaScript engine.

    python -m majorana_lattice.export out.json star:4x4 star:6x6
"""
import json
import sys

from .lattice import PRESETS


def export(specs: list[str], minify: bool = True) -> dict:
    out = {}
    for spec in specs:
        name, size = spec.split(":")
        cols, rows = (int(v) for v in size.split("x"))
        d = PRESETS[name](cols, rows).to_json()
        if minify:
            d["positions"] = [[round(a, 4), round(b, 4)] for a, b in d["positions"]]
            d["edges"] = [[a, b, l, round(x, 4), round(y, 4)] for a, b, l, x, y in d["edges"]]
            d["period"] = [round(v, 4) for v in d["period"]]
        out[f"{name.replace('_', ' ')} {cols}x{rows}"] = d
    return out


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 2:
        print(__doc__)
        return 2
    with open(argv[0], "w") as f:
        json.dump(export(argv[1:]), f, separators=(",", ":"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
