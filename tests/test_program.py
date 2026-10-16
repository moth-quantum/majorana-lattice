"""Programs: the vacuum, the bookkeeping (checked by Stim itself), running,
replay from the circuit, and executors."""
import random

import numpy as np
import pytest
import stim

from majorana_lattice import (Link, Plaquette, Program, Vertex, anyon_moves, apply_link, apply_pauli,
                              create_anyons, create_majoranas, dimers, fermion_moves, honeycomb,
                              majorana_moves, majoranas, measure_link, modified_honeycomb,
                              move_anyon, move_fermion, move_majorana, plaquette_values, prepare,
                              readout, release, run, split_dimer, star, stim_executor)


def postselected_vacuum(L, background):
    """The vacuum by definition: the background product state, postselected
    onto every plaquette reading +1. Returns the state and the plaquettes."""
    sim = stim.TableauSimulator()
    for q in range(L.n):
        if background == "x":
            sim.h(q)
        elif background == "y":
            sim.h_yz(q)
    plaqs = []
    for f in L.faces:
        p = stim.PauliString(L.n)
        for q, P in f.pauli.items():
            p[q] = P
        if sim.peek_observable_expectation(p) == 0:
            sim.postselect_observable(p, desired_value=False)
        plaqs.append(p if sim.peek_observable_expectation(p) == 1 else -p)
    return sim, plaqs


def random_step(rng, p):
    """One random action, chosen from what the bookkeeping offers."""
    L = p.lattice
    r = rng.random()
    if r < 0.3:
        return measure_link(p, rng.randrange(len(L.edges)))
    if r < 0.4:
        return apply_pauli(p, rng.randrange(L.n), rng.choice("XYZ"))
    if r < 0.5:
        return apply_link(p, rng.randrange(len(L.edges)))
    if r < 0.6:
        return release(p, rng.randrange(L.n))
    if r < 0.85 and majoranas(p):
        v = rng.choice(majoranas(p))
        mv = majorana_moves(p, v)
        if mv:
            return move_majorana(p, v, rng.choice(mv)[0])
    if r < 0.92:
        v = rng.randrange(L.n)
        return move_fermion(p, v, rng.choice(fermion_moves(p, v))[0])
    fs = [f for f in range(len(L.faces)) if anyon_moves(p, f)]
    if fs:
        f = rng.choice(fs)
        return move_anyon(p, f, rng.choice(anyon_moves(p, f))[0])
    return p


def random_program(seed, steps, L=None):
    rng = random.Random(seed)
    p = prepare(L or star(4, 4))
    for _ in range(steps):
        p = random_step(rng, p)
    return p


def flow_holds(circuit, P):
    """Stim's check that the pair's operator equals the parity of its records."""
    op = -P.op if P.flip else P.op
    recs = "".join(f" xor rec[{k}]" for k in sorted(P.records))
    return circuit.has_flow(stim.Flow(f"1 -> {str(op).replace('+', '')}{recs}"))


# ------------------------------------------------------------------ vacuum
@pytest.mark.parametrize("make", [star, honeycomb, modified_honeycomb])
@pytest.mark.parametrize("background", "xyz")
def test_vacuum_is_the_postselected_state(make, background):
    L = make(4, 4)
    p = prepare(L, background)
    ref, plaqs = postselected_vacuum(L, background)
    for seed in range(3):
        sim = stim.TableauSimulator(seed=seed)
        sim.do(p.circuit)
        assert sim.canonical_stabilizers() == ref.canonical_stabilizers()
    assert [str(a) for a in p.plaq_ops] == [str(b) for b in plaqs]


def test_vacuum_reads_empty():
    p = prepare(star(4, 4))
    res = run(p, shots=200, seed=1)
    assert len(dimers(p)) == 48 and p.circuit.num_qubits == 96
    assert all((v == 1).all() for v in res.pairs().values())
    assert (res.plaquettes() == 1).all()
    assert p.measurements == tuple(("prepare", f) for f in range(48))
    assert all(flow_holds(p.circuit, P) for P in p.pairs)


# ------------------------------------------------------------ bookkeeping
@pytest.mark.parametrize("seed", range(4))
def test_stim_confirms_every_pair_value(seed):
    """As the program grows, Stim's flow checker agrees, sign included, that
    each tracked pair's operator equals the parity of its records."""
    rng = random.Random(seed)
    p = prepare(star(4, 4))
    for step in range(50):
        p = random_step(rng, p)
        if step % 10 == 9:
            circuit = p.circuit
            assert all(flow_holds(circuit, P) for P in p.pairs)


@pytest.mark.parametrize("seed", range(3))
def test_plaquette_values_are_what_stim_measures(seed):
    p = random_program(seed, 60)
    values = plaquette_values(p)
    for f, op in enumerate(p.plaq_ops):
        p = readout(p, f"plaq{f}", op)
    res = run(p, shots=50, seed=seed)
    for f, v in enumerate(values):
        assert (res[f"plaq{f}"] == v).all()


def test_sampled_pairs_match_their_readouts():
    """Reading every pair out with a real measurement gives, in every shot,
    the value its records predicted."""
    p = random_program(5, 80)
    pairs = list(p.pairs)
    for k, P in enumerate(pairs):
        p = readout(p, f"pair{k}", P)
    res = run(p, shots=300, seed=2)
    for k, P in enumerate(pairs):
        assert (res[f"pair{k}"] == res.value(P)).all()
    assert len({tuple(res.value(P)) for P in pairs}) > 1          # some really are random


def test_pairs_cover_every_vertex_once():
    p = random_program(7, 100)
    ends = sorted(q for P in p.pairs for q in (P.a, P.b))
    assert ends == list(range(p.lattice.n))


# ----------------------------------------------------------- the API itself
def test_addresses_in_and_out():
    L = star(4, 4)
    p = split_dimer(prepare(L), Link(0, 1, "z"))
    assert majoranas(p) == [Vertex(0, 1, "D.z"), Vertex(0, 1, "U.z")]
    link, dest = majorana_moves(p, Vertex(0, 1, "D.z"))[0]
    p = move_majorana(p, Vertex(0, 1, "D.z"), link)
    assert dest in majoranas(p)
    p = create_anyons(p, Link(2, 2, "U.z"))
    anyons = [L.address(Plaquette, k) for k, v in enumerate(plaquette_values(p)) if v == -1]
    assert len(anyons) == 2 and anyon_moves(p, anyons[0])


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
    lambda p: readout(readout(p, "a", p.pairs[0]), "a", p.pairs[1]),      # name taken
])
def test_refused(bad):
    with pytest.raises(ValueError):
        bad(prepare(star(4, 4)))


def test_wrong_kind_of_address():
    with pytest.raises(TypeError):
        apply_pauli(prepare(star(4, 4)), Link(0, 0, "z"), "X")


def test_from_circuit_round_trip():
    p = random_program(3, 80)
    p = readout(p, "first", p.pairs[0])
    q = Program.from_circuit(p.lattice, stim.Circuit(str(p.circuit)))
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


# ------------------------------------------------------------- executors
def counting(executor, calls):
    def execute(circuit, shots):
        calls.append(shots)
        return executor(circuit, shots)
    return execute


def test_a_program_carries_its_executor():
    calls, other = [], []
    p = prepare(star(4, 4), executor=counting(stim_executor(1), calls))
    p = split_dimer(p, Link(1, 1, "z"))                           # inherited by new programs
    run(p, shots=7)
    assert calls == [7]
    run(p, shots=3, executor=counting(stim_executor(2), other))   # run can override it
    assert calls == [7] and other == [3]
    run(p.with_executor(None), shots=5)                           # back to local Stim
    assert calls == [7]


def test_executor_survives_replay():
    calls = []
    p = split_dimer(prepare(star(4, 4)), Link(1, 1, "z"))
    q = Program.from_circuit(p.lattice, p.circuit, executor=counting(stim_executor(1), calls))
    run(q, shots=4)
    assert calls == [4]


def test_executor_must_return_every_record():
    with pytest.raises(ValueError):
        run(prepare(star(4, 4)), shots=2, executor=lambda c, s: np.zeros((s, 1), dtype=bool))
