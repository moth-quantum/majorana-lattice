"""The circuit engine (Program, ops, queries, run) against the interactive one
(MatchingCode), against Stim's own flow checker, and against the physics."""
import importlib.util
import random
from pathlib import Path

import numpy as np
import pytest
import stim

from majorana_lattice import (Link, MatchingCode, Plaquette, Program, Vertex, anyon_moves, apply_link,
                              apply_pauli, create_anyons, create_majoranas, dimers, fermion_moves,
                              honeycomb, is_dimer, majorana_moves, majoranas, measure_link,
                              modified_honeycomb, move_anyon, move_fermion, move_majorana,
                              plaquette_values, prepare, readout, release, run, split_dimer, star,
                              stim_executor)
from majorana_lattice.run import Result

ROOT = Path(__file__).resolve().parent.parent


def ints(lat, xs):
    return [lat.index(x) for x in xs]


def moves_as_ints(lat, moves):
    return [(lat.index(a), lat.index(b)) for a, b in moves]


def value_on(p, P, bits):
    """A pair's value given one shot's measurement records."""
    return int(Result(p, np.array([bits], dtype=bool)).value(P)[0])


def flow_holds(circuit, P):
    """Stim's check that the pair's operator equals the parity of its records."""
    op = -P.op if P.flip else P.op
    recs = "".join(f" xor rec[{k}]" for k in sorted(P.records))
    return circuit.has_flow(stim.Flow(f"1 -> {str(op).replace('+', '')}{recs}"))


# ------------------------------------------------------------------ vacuum
@pytest.mark.parametrize("make", [star, honeycomb, modified_honeycomb])
@pytest.mark.parametrize("background", "xyz")
def test_vacuum_is_the_reference_state(make, background):
    L = make(4, 4)
    p = prepare(L, background)
    ref = MatchingCode(L, background=background, seed=0)
    for seed in range(3):
        sim = stim.TableauSimulator(seed=seed)
        sim.do(p.circuit)
        assert sim.canonical_stabilizers() == ref.sim.canonical_stabilizers()
    assert [str(a) for a in p.plaq_ops] == [str(b) for b in ref.plaq]


def test_vacuum_reads_empty():
    p = prepare(star(4, 4))
    res = run(p, shots=200, seed=1)
    assert len(dimers(p)) == 48 and p.circuit.num_qubits == 96
    assert all((v == 1).all() for v in res.pairs().values())
    assert (res.plaquettes() == 1).all()
    assert p.measurements == tuple(("prepare", f) for f in range(48))


# ------------------------------------------------- against MatchingCode
def random_step(rng, c, p):
    """The same random action on both engines; returns the new Program."""
    L = c.lat
    r = rng.random()
    if r < 0.3:
        e = rng.randrange(len(L.edges))
        c.measure_link(e)
        return measure_link(p, e)
    if r < 0.4:
        q, P = rng.randrange(L.n), rng.choice("XYZ")
        c.apply_pauli(q, P)
        return apply_pauli(p, q, P)
    if r < 0.5:
        e = rng.randrange(len(L.edges))
        c.apply_link(e)
        return apply_link(p, e)
    if r < 0.6:
        q = rng.randrange(L.n)
        c.release(q)
        return release(p, q)
    if r < 0.85 and c.majoranas():
        q = rng.choice(c.majoranas())
        mv = c.majorana_moves(q)
        if mv:
            e = rng.choice(mv)[0]
            c.move_majorana(q, e)
            return move_majorana(p, q, e)
    if r < 0.92:
        q = rng.randrange(L.n)
        e = rng.choice(c.fermion_moves(q))[0]
        c.move_fermion(q, e)
        return move_fermion(p, q, e)
    fs = [f for f in range(len(L.faces)) if c.anyon_moves(f)]
    if fs:
        f = rng.choice(fs)
        e = rng.choice(c.anyon_moves(f))[0]
        c.move_anyon(f, e)
        return move_anyon(p, f, e)
    return p


@pytest.mark.parametrize("seed", range(4))
def test_bookkeeping_matches_matching_code(seed):
    """After every step: the same pairs and operators, the same moves on offer,
    and every pair's value, evaluated on MatchingCode's own outcomes, equals
    what MatchingCode reads off its quantum state."""
    rng = random.Random(seed)
    L = star(4, 4)
    c, p = MatchingCode(L, seed=seed), prepare(L)
    for _ in range(120):
        p = random_step(rng, c, p)
        assert sorted((P.a, P.b, P.kind, P.path, str(P.op)) for P in p.pairs) == \
            sorted((P.a, P.b, P.kind, tuple(P.path), str(P.op)) for P in c.pairs.values())
        outcomes = [s["outcome"] == -1 for s in c.log if s["action"] == "measure_link"]
        bits = [False] * len(L.faces) + outcomes        # prep records never enter a pair's value
        ref = c.values()
        for pid, P in c.pairs.items():
            mine = p.pair_at(P.a)
            assert value_on(p, mine, bits) == ref["pairs"][pid]
        assert plaquette_values(p) == ref["plaquettes"]
        assert ints(L, majoranas(p)) == c.majoranas()
        q = rng.randrange(L.n)
        assert moves_as_ints(L, majorana_moves(p, q)) == c.majorana_moves(q)
        assert moves_as_ints(L, fermion_moves(p, q)) == c.fermion_moves(q)
        f = rng.randrange(len(L.faces))
        assert moves_as_ints(L, anyon_moves(p, f)) == c.anyon_moves(f)


@pytest.mark.parametrize("seed", range(3))
def test_stim_confirms_every_pair_value(seed):
    """Stim's flow checker agrees, sign included, that each tracked pair's
    operator equals the parity of its records."""
    rng = random.Random(seed)
    L = star(4, 4)
    c, p = MatchingCode(L, seed=seed), prepare(L)
    for _ in range(60):
        p = random_step(rng, c, p)
    circuit = p.circuit
    assert all(flow_holds(circuit, P) for P in p.pairs)


def test_sampled_pairs_match_their_readouts():
    """Reading every pair out with a real measurement gives, in every shot,
    the value its records predicted."""
    rng = random.Random(5)
    L = star(4, 4)
    c, p = MatchingCode(L, seed=5), prepare(L)
    for _ in range(80):
        p = random_step(rng, c, p)
    pairs = list(p.pairs)
    for k, P in enumerate(pairs):
        p = readout(p, f"pair{k}", P)
    res = run(p, shots=300, seed=2)
    for k, P in enumerate(pairs):
        assert (res[f"pair{k}"] == res.value(P)).all()
    assert len({tuple(res.value(P)) for P in pairs}) > 1          # some really are random


# ------------------------------------------------------------- physics
def load_example(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "examples" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_braid():
    b = load_example("braid")
    out = {}
    for name, stage in [("start", b.p), ("once", b.once), ("twice", b.twice)]:
        res = run(readout(readout(stage, "A", b.piA), "B", b.piB), shots=1000, seed=3)
        out[name] = ((res["A"] == 1).mean(), (res["B"] == 1).mean())
    assert out["start"] == (1.0, 1.0)
    assert all(0.4 < x < 0.6 for x in out["once"])
    assert out["twice"] == (0.0, 0.0)


def test_linking():
    m = load_example("linking")
    for linked in (False, True):
        res = run(m.experiment(linked), shots=500, seed=4)
        both_empty = (res["V"] == 1) & (res["H"] == 1)
        assert both_empty.all() if not linked else not both_empty.any()


def test_fermions_come_in_pairs():
    rng = random.Random(9)
    L = star(4, 4)
    p = prepare(L)
    stages = []
    for _ in range(60):
        p = measure_link(p, rng.randrange(len(L.edges)))
        stages.append(p)
    res = run(p, shots=300, seed=1)
    for s in stages[::10]:
        assert (res.summary(at=s)["fermions"] % 2 == 0).all()
    assert res.summary()["fermions"].max() > 0


# ----------------------------------------------------------- the API itself
def test_addresses_in_and_out():
    L = star(4, 4)
    p = split_dimer(prepare(L), Link(0, 1, "z"))
    assert majoranas(p) == [Vertex(0, 1, "D.z"), Vertex(0, 1, "U.z")]
    link, dest = majorana_moves(p, Vertex(0, 1, "D.z"))[0]
    p = move_majorana(p, Vertex(0, 1, "D.z"), link)
    assert dest in majoranas(p)
    p = create_anyons(p, Link(2, 2, "U.z"))
    assert sorted(plaquette_values(p)).count(-1) == 2
    f = next(Plaquette(*L.address(Plaquette, k).__dict__.values())
             for k, v in enumerate(plaquette_values(p)) if v == -1)
    assert anyon_moves(p, f)


def test_programs_are_values():
    L = star(4, 4)
    p0 = prepare(L)
    p1 = split_dimer(p0, Link(0, 0, "z"))
    p2 = create_anyons(p0, Link(1, 1, "z"))
    assert p0.instructions == () and len(majoranas(p0)) == 0
    assert len(majoranas(p1)) == 2 and len(majoranas(p2)) == 0
    assert p1.prep is p0.prep and p1.extends(p0) and not p1.extends(p2)


@pytest.mark.parametrize("bad", [
    lambda p: move_majorana(p, Vertex(0, 0, "D.x"), Link(0, 0, "x")),     # not a Majorana
    lambda p: split_dimer(p, Link(0, 0, "x")),                            # not a dimer
    lambda p: move_anyon(p, Plaquette(0, 0, "D"), Link(0, 0, "x")),       # no anyon there
    lambda p: create_majoranas(p, Link(0, 0, "z")),                       # already one dimer
])
def test_refused(bad):
    with pytest.raises(ValueError):
        bad(prepare(star(4, 4)))


def test_wrong_kind_of_address():
    with pytest.raises(TypeError):
        apply_pauli(prepare(star(4, 4)), Link(0, 0, "z"), "X")


def test_from_circuit_round_trip():
    rng = random.Random(3)
    L = star(4, 4)
    c, p = MatchingCode(L, seed=3), prepare(L)
    for _ in range(80):
        p = random_step(rng, c, p)
    p = readout(p, "first", p.pairs[0])
    q = Program.from_circuit(L, p.circuit)
    assert str(q.circuit) == str(p.circuit)
    assert q.measurements == p.measurements and q.readouts == p.readouts
    assert [(P.a, P.b, P.kind, P.path, str(P.op), P.flip, P.records) for P in q.pairs] == \
        [(P.a, P.b, P.kind, P.path, str(P.op), P.flip, P.records) for P in p.pairs]


def test_fix_chance_is_drawn_when_building():
    L = star(4, 4)
    def build(seed):
        rng, p = random.Random(seed), split_dimer(prepare(L), Link(1, 1, "z"))
        for _ in range(6):
            v = majoranas(p)[0]
            p = move_majorana(p, v, majorana_moves(p, v)[0][0], fix=0.5, rng=rng)
        return str(p.circuit)
    assert build(1) == build(1)
    assert len({build(s) for s in range(6)}) > 1


def test_archaeo_executor():
    pytest.importorskip("archaeo")
    from majorana_lattice import archaeo_executor
    L = star(4, 4)
    p = split_dimer(prepare(L), Link(1, 1, "z"))
    for _ in range(4):
        v = majoranas(p)[0]
        p = move_majorana(p, v, majorana_moves(p, v)[0][0])
    p = measure_link(p, next(e for e in range(len(L.edges)) if not is_dimer(p, e)))
    res = run(p, shots=200, executor=archaeo_executor())
    assert res.shots == 200
    # pairs whose value is fixed read the same through archaeo as through local Stim
    local = run(p, shots=200, executor=stim_executor(1))
    for P in p.pairs:
        a, b = res.value(P), local.value(P)
        if len(set(b)) == 1:
            assert (a == b[0]).all()
