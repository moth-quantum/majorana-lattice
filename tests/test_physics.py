"""The physics of arXiv:1501.07779 on the star lattice, from running programs."""
import importlib.util
import random
from pathlib import Path

import stim

from majorana_lattice import (PRESETS, Link, Vertex, apply_pauli, colouring, create_majoranas,
                              fuse, is_majorana, link_pauli, majorana_moves, measure_label, measure_link, move_majorana,
                              partner, plaquette_values, prepare, readout, route_majorana, run, star)

ROOT = Path(__file__).resolve().parent.parent


def load_example(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "examples" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_every_plaquette_commutes_with_every_link():
    for name, make in PRESETS.items():
        L = make(4, 4)
        for f in L.faces:
            W = stim.PauliString(L.n)
            for q, P in f.pauli.items():
                W[q] = P
            assert all(link_pauli(L, e).commutes(W) for e in range(len(L.edges))), name


def test_link_measurements_never_make_anyons_and_fermions_come_in_pairs():
    rng = random.Random(7)
    L = star(4, 4)
    p = prepare(L)
    stages = []
    for _ in range(300):
        p = measure_link(p, rng.randrange(len(L.edges)))
        stages.append(p)
    assert all(v == 1 for v in plaquette_values(p))
    assert all(P.records is not None for P in p.pairs)      # every pair stays determined
    res = run(p, shots=300, seed=1)
    for s in stages[::25]:
        assert (res.summary(at=s)["fermions"] % 2 == 0).all()
    assert res.summary()["fermions"].max() > 0


def test_a_single_error():
    """X on one qubit flips two plaquettes. The anyons match (e,e or m,m) with
    an even number of fermions, or differ (e,m) with an odd number."""
    L = star(4, 4)
    for q in range(0, L.n, 7):
        for P in "XYZ":
            p = apply_pauli(prepare(L), q, P)
            s = run(p, shots=20, seed=q).summary()
            assert s["e"] + s["m"] == 2
            mixed = s["e"] == s["m"]
            assert ((s["fermions"] % 2 == 1) == mixed).all()


def test_a_floquet_round_repeated_gives_the_same_outcomes():
    L = star(4, 4)
    n = sum(1 for E in L.edges if E.label == "x")
    p = measure_label(measure_label(prepare(L), "x"), "x")
    out = run(p, shots=200, seed=2).outcomes()
    first, second = out[:, -2 * n:-n], out[:, -n:]
    assert (first == second).all()
    assert len({tuple(r) for r in first}) > 1                  # and the first round is random


def test_braid():
    """One exchange puts both parities into superposition; two flip them."""
    b = load_example("braid")
    out = {}
    for name, stage in [("start", b.p), ("once", b.once), ("twice", b.twice)]:
        res = run(readout(readout(stage, "A", b.piA), "B", b.piB), shots=1000, seed=3)
        out[name] = ((res["A"] == 1).mean(), (res["B"] == 1).mean())
    assert out["start"] == (1.0, 1.0)
    assert all(0.4 < x < 0.6 for x in out["once"])
    assert out["twice"] == (0.0, 0.0)


def test_linked_and_unlinked_pairs():
    """Linked pairs both fuse to a fermion; unlinked pairs both fuse to vacuum."""
    m = load_example("linking")
    for linked in (False, True):
        res = run(m.experiment(linked), shots=500, seed=4)
        assert (res["V"] == (-1 if linked else 1)).all()
        assert (res["H"] == (-1 if linked else 1)).all()


def test_route_independence():
    """A pair from vacuum, wandered round a patch at random and fused, comes
    back empty whatever route it took."""
    L = star(10, 10)
    C = (5.0, 5.0)
    out = {L.address(Vertex, q) for q in range(L.n) if L.distance(L.address(Vertex, q), C) > 3.6}
    done = 0
    for s in range(10):
        r = random.Random(s)
        p = prepare(L)
        link = min((L.address(Link, e) for e, E in enumerate(L.edges)
                    if p.owner[E.u] != p.owner[E.v]), key=lambda l: L.distance(l, C) + 0.01 * r.random())
        E = L.edges[L.index(link)]
        a, b = partner(p, E.u), partner(p, E.v)
        p = create_majoranas(p, link)
        for _ in range(40):
            if not is_majorana(p, a):
                break               # the path collapsed to one link: the pair fused by itself
            mover = r.choice([0, 1])
            x = (a, b)[mover]
            opts = [(f, t) for f, t in majorana_moves(p, x)
                    if t not in out and t not in (a, b) and L.address(Vertex, L.edges[L.index(f)].other(L.index(x))) not in (a, b)]
            if opts:
                f, t = r.choice(opts)
                p = move_majorana(p, x, f)
                a, b = (t, b) if mover == 0 else (a, t)
        if not is_majorana(p, a):
            assert (run(p, shots=200, seed=s).value(p.pair_at(a)) == 1).all(), s
            done += 1
            continue
        for f in L.incident[L.index(b)]:
            t = L.address(Vertex, L.edges[f].other(L.index(b)))
            if t in out:
                continue
            try:
                if a != t:
                    p = route_majorana(p, a, t, avoid=out | {b})
                    a = t
                break
            except ValueError:
                continue
        e = L.edge_between(L.index(a), L.index(b))
        if e is None:
            continue                    # walled in by the patch edge
        p = fuse(p, e, "pair")
        assert (run(p, shots=200, seed=s)["pair"] == 1).all(), s
        done += 1
    assert done >= 7


def test_anyon_types_follow_the_colouring():
    """In vacuum the z links join same-colour plaquettes, so a z-link Pauli
    makes two anyons of one type."""
    L = star(4, 4)
    col, ok = colouring(prepare(L))
    assert ok
    for e, E in enumerate(L.edges):
        if E.label != "z":
            continue
        p = apply_pauli(prepare(L), E.u, "Z")
        flipped = [k for k, v in enumerate(plaquette_values(p)) if v == -1]
        assert len(flipped) == 2 and col[flipped[0]] == col[flipped[1]]
