"""Majoranas, fermions and anyons share one move API, in Python and in the browser engine.

For each particle: a list of the moves on offer, a method that makes one move,
a ValueError for anything not on the list, and the same lists from both engines.
"""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from majorana_lattice import MatchingCode, star
from majorana_lattice.export import export

HERE = Path(__file__).parent


def _board(seed=1):
    """Some of everything: a split dimer, a fermion, an e and an m."""
    L = star(4, 4)
    c = MatchingCode(L, seed=seed)
    z = [E for E in L.edges if E.label == "z"]
    c.apply_pauli(z[0].u, "X")                 # an e, an m and a fermion
    c.release(z[5].u)                          # two Majoranas
    return c


def test_every_particle_has_moves_and_a_move_method():
    c = _board()
    for name in ("majorana_moves", "fermion_moves", "anyon_moves", "move_majorana", "move_fermion", "move_anyon"):
        assert callable(getattr(c, name))
    assert not hasattr(c, "hop") and not hasattr(c, "hop_options")


def test_move_majorana_lands_where_the_list_says():
    c = _board()
    q = c.majoranas()[0]
    link, dest = c.majorana_moves(q)[0]
    assert c.move_majorana(q, link) == dest
    assert dest in c.majoranas()


def test_move_fermion_moves_the_fermion():
    c = _board()
    v = lambda: c.values()["pairs"]
    holder = next(pid for pid, x in v().items() if x == -1)
    q = c.pairs[holder].a
    before = sum(1 for x in v().values() if x == -1)
    link, dest = next((l, d) for l, d in c.fermion_moves(q) if v()[c.owner[d]] == 1)
    assert c.move_fermion(q, link) == dest
    assert v()[holder] == 1 and v()[c.owner[dest]] == -1
    assert sum(1 for x in v().values() if x == -1) == before


@pytest.mark.parametrize("move, listing, args", [
    ("move_majorana", "majorana_moves", lambda c: c.majoranas()[0]),
    ("move_fermion", "fermion_moves", lambda c: 0),
    ("move_anyon", "anyon_moves", lambda c: next(f for f, x in enumerate(c.values()["plaquettes"]) if x == -1)),
])
def test_a_move_not_on_the_list_is_refused(move, listing, args):
    c = _board()
    x = args(c)
    allowed = {e for e, _ in getattr(c, listing)(x)}
    bad = next(e for e in range(len(c.lat.edges)) if e not in allowed)
    log = len(c.log)
    with pytest.raises(ValueError):
        getattr(c, move)(x, bad)
    assert len(c.log) == log


@pytest.mark.skipif(shutil.which("node") is None, reason="needs Node.js")
def test_browser_engine_offers_the_same_moves(tmp_path):
    lat_file, script = tmp_path / "lattices.json", tmp_path / "moves.js"
    lat_file.write_text(json.dumps(export(["star:4x4"])))
    script.write_text("""
const { MatchingCode } = require(process.argv[2]);
const lats = JSON.parse(require("fs").readFileSync(process.argv[3], "utf8"));
const lat = { ...lats["star 4x4"], name: "star 4x4" };
const c = new MatchingCode(lat, { seed: 1 });
const z = c.edges.filter(E => E.label === "z");
c.applyPauli(z[0].u, "X");
c.release(z[5].u);
const key = o => [o.edge, o.to];
const sorted = l => l.map(key).sort((a, b) => a[0] - b[0] || a[1] - b[1]);
const out = { majorana: {}, fermion: {}, anyon: {} };
for (let q = 0; q < c.n; q++) { out.majorana[q] = sorted(c.majoranaMoves(q)); out.fermion[q] = sorted(c.fermionMoves(q)); }
lat.faces.forEach((_, f) => { out.anyon[f] = sorted(c.anyonMoves(f)); });
console.log(JSON.stringify(out));
""")
    js = json.loads(subprocess.run(["node", str(script), str(HERE.parent / "js" / "engine.js"), str(lat_file)],
                                   capture_output=True, text=True, check=True).stdout)
    c = _board()
    for q in range(c.lat.n):
        assert sorted(map(list, c.majorana_moves(q))) == js["majorana"][str(q)], ("majorana", q)
        assert sorted(map(list, c.fermion_moves(q))) == js["fermion"][str(q)], ("fermion", q)
    for f in range(len(c.lat.faces)):
        assert sorted(map(list, c.anyon_moves(f))) == js["anyon"][str(f)], ("anyon", f)


# ------------------------------------------------------------ fixing a -1 outcome
def _move(seed, fix=None, fix_rate=None):
    """Split a dimer and move its Majorana. The outcome depends only on the seed."""
    L = star(4, 4)
    c = MatchingCode(L, seed=seed)
    if fix_rate is not None:
        c.fix_rate = fix_rate
    c.release(next(E for E in L.edges if E.label == "z").u)
    q = c.majoranas()[0]
    link, _ = c.majorana_moves(q)[0]
    c.move_majorana(q, link, fix)
    return c


def _outcome(c):
    return next(s["outcome"] for s in c.log if s["action"] == "measure_link")


def _minus_one_seeds(n):
    """The first n seeds whose move measures -1."""
    out, seed = [], 0
    while len(out) < n:
        seed += 1
        if _outcome(_move(seed)) == -1:
            out.append(seed)
    return out


def test_fix_one_removes_the_fermions_and_zero_leaves_them():
    s = _minus_one_seeds(1)[0]
    assert _move(s, fix=1).summary()["fermions"] == 0
    assert _move(s, fix=0).summary()["fermions"] == 2
    assert _move(s, fix=True).summary()["fermions"] == 0     # old bool spelling
    assert _move(s, fix=False).summary()["fermions"] == 2


def test_class_default_applies_unless_the_move_says_otherwise():
    s = _minus_one_seeds(1)[0]
    assert MatchingCode(star(4, 4)).fix_rate == 1.0
    assert _move(s).summary()["fermions"] == 0                         # default is to fix
    assert _move(s, fix_rate=0).summary()["fermions"] == 2             # class default used
    assert _move(s, fix_rate=0, fix=1).summary()["fermions"] == 0      # per move wins
    assert _move(s, fix_rate=1, fix=0).summary()["fermions"] == 2


def test_fix_is_a_probability():
    seeds = _minus_one_seeds(60)
    fixed = sum(_move(s, fix=0.5).summary()["fermions"] == 0 for s in seeds)
    assert 15 < fixed < 45


def test_fixing_does_not_change_the_measurement_outcome():
    s = _minus_one_seeds(1)[0]
    for fix in (0, 0.5, 1):
        assert _outcome(_move(s, fix=fix)) == -1


def test_partial_fix_replays_exactly():
    for s in _minus_one_seeds(5):
        c = _move(s, fix=0.5)
        assert MatchingCode.replay(c.lat, c.record()).values() == c.values()


# ------------------------------------------------------------ choosing an outcome
def _link_between_dimers(c):
    return next(i for i, E in enumerate(c.lat.edges) if c.owner[E.u] != c.owner[E.v])


def test_outcome_postselects_a_random_measurement():
    for want in (1, -1):
        for seed in (1, 2, 3):
            c = MatchingCode(star(4, 4), seed=seed)
            assert c.measure_link(_link_between_dimers(c), outcome=want) == want


def test_outcome_never_changes_a_deterministic_measurement():
    c = MatchingCode(star(4, 4), seed=1)
    z = next(i for i, E in enumerate(c.lat.edges) if c.owner[E.u] == c.owner[E.v])   # a dimer's own link
    assert c.measure_link(z, outcome=-1) == 1


def test_omitting_outcome_leaves_it_random():
    seen = set()
    for s in range(1, 15):
        c = MatchingCode(star(4, 4), seed=s)
        seen.add(c.measure_link(_link_between_dimers(c)))
    assert seen == {1, -1}
