"""The browser engine (js/engine.js) must match Stim exactly.

Random sessions are run in Node, then replayed here in Stim with every random
outcome forced to match; every value is compared after every step.
"""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from majorana_lattice import MatchingCode, star
from majorana_lattice.export import export

HERE = Path(__file__).parent


@pytest.mark.skipif(shutil.which("node") is None, reason="needs Node.js")
def test_browser_engine_matches_stim(tmp_path):
    lat_file, runs_file = tmp_path / "lattices.json", tmp_path / "runs.json"
    lat_file.write_text(json.dumps(export(["star:4x4"])))
    subprocess.run(["node", str(HERE / "js" / "crosscheck_run.js"), str(lat_file), str(runs_file)], check=True)
    runs = json.loads(runs_file.read_text())
    lat = star(4, 4)
    compared = 0
    for run in runs:
        rec = run["record"]
        c = MatchingCode(lat, plaquette_outcomes=rec["plaquette_raw"])
        log, k = rec["log"], 0
        for snap in run["snaps"]:
            while k < snap["steps"]:
                st, a = log[k], log[k]["action"]
                if a == "measure_link":
                    assert c.measure_link(st["edge"], outcome=st["outcome"]) == st["outcome"], st
                elif a == "apply_link":
                    c.apply_link(st["edge"])
                elif a == "apply_pauli":
                    c.apply_pauli(st["qubit"], st["pauli"])
                elif a == "release":
                    c.release(st["qubit"])
                k += 1
            v = c.values()
            pairs = sorted([min(P.a, P.b), max(P.a, P.b), P.kind, v["pairs"][pid]] for pid, P in c.pairs.items())
            assert pairs == snap["pairs"]
            assert v["plaquettes"] == snap["plaq"]
            assert c.colouring()[1] == snap["col"]["ok"]
            compared += 1
    assert compared == 12 * 120
