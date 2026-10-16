"""Majoranas, fermions and anyons share one move API: a list of the moves on
offer, a function that makes one, and a ValueError for anything else."""
import random

import pytest

from majorana_lattice import (Link, Plaquette, anyon_moves, apply_pauli, colouring, fermion_moves,
                              majorana_moves, majoranas, move_anyon, move_fermion, move_majorana,
                              partner, plaquette_values, prepare, release, run, split_dimer, star)


def board():
    """Some of everything: a split dimer, a fermion, an e and an m."""
    L = star(4, 4)
    z = [E for E in L.edges if E.label == "z"]
    p = apply_pauli(prepare(L), z[0].u, "X")             # an e, an m and a fermion
    return release(p, z[5].u)                            # two Majoranas


def anyons(p):
    return [p.lattice.address(Plaquette, k) for k, v in enumerate(plaquette_values(p)) if v == -1]


# ---------------------------------------------------------------- Majoranas
def test_move_majorana_lands_where_the_list_says():
    p = board()
    v = majoranas(p)[0]
    link, dest = majorana_moves(p, v)[0]
    p = move_majorana(p, v, link)
    assert dest in majoranas(p) and v not in majoranas(p)


def test_fix_removes_the_fermions_and_no_fix_leaves_them():
    p = split_dimer(prepare(star(4, 4)), Link(1, 1, "z"))
    v = majoranas(p)[0]
    link, _ = majorana_moves(p, v)[0]
    fixed = run(move_majorana(p, v, link, fix=True), shots=300, seed=1)
    loose = run(move_majorana(p, v, link, fix=False), shots=300, seed=1)
    assert (fixed.summary()["fermions"] == 0).all()
    # without the fix, a -1 outcome leaves a fermion on each of the two new pairs
    minus = loose.outcomes()[:, -1] == -1
    assert ((loose.summary()["fermions"] == 2) == minus).all() and 0 < minus.mean() < 1


def test_fix_chance_is_a_probability():
    rng = random.Random(3)
    p = split_dimer(prepare(star(4, 4)), Link(1, 1, "z"))
    v = majoranas(p)[0]
    link, _ = majorana_moves(p, v)[0]
    fixed = sum(len(move_majorana(p, v, link, fix=0.5, rng=rng).instructions) > 2 for _ in range(200))
    assert 70 < fixed < 130


# ---------------------------------------------------------------- fermions
def test_move_fermion_moves_the_fermion():
    p = board()
    key = lambda v: tuple(sorted((v, partner(p, v)), key=p.lattice.index))
    before = run(p, shots=1, seed=0).pairs()
    holder = next(k for k, v in before.items() if v[0] == -1)
    link, dest = next((l, d) for l, d in fermion_moves(p, holder[0]) if before[key(d)][0] == 1)
    after = run(move_fermion(p, holder[0], link), shots=1, seed=0).pairs()
    assert after[holder][0] == 1 and after[key(dest)][0] == -1
    assert sum(v[0] == -1 for v in after.values()) == sum(v[0] == -1 for v in before.values())


# ------------------------------------------------------------------ anyons
def with_e():
    L = star(4, 4)
    z = [E for E in L.edges if E.label == "z"][0]
    p = apply_pauli(prepare(L), z.u, "X")
    col, _ = colouring(p)
    return p, next(f for f in anyons(p) if col[p.lattice.index(f)] == 0)


def test_moving_an_anyon_moves_only_the_anyon():
    p, f = with_e()
    fermions = run(p, shots=50, seed=1).summary()["fermions"]
    for _ in range(12):
        opts = anyon_moves(p, f)
        assert opts, f"the anyon on {f} has nowhere to go"
        before, col = anyons(p), colouring(p)[0]
        link, dest = opts[0]
        p = move_anyon(p, f, link)
        assert sorted(anyons(p), key=p.lattice.index) == \
            sorted([a for a in before if a != f] + [dest], key=p.lattice.index)
        assert colouring(p)[0] == col and col[p.lattice.index(dest)] == 0
        assert (run(p, shots=50, seed=1).summary()["fermions"] == fermions).all()
        f = dest


def test_only_allowed_plaquettes_are_offered():
    p, _ = with_e()
    col, _ = colouring(p)
    for k in range(len(p.lattice.faces)):
        f = p.lattice.address(Plaquette, k)
        opts = anyon_moves(p, f)
        if f not in anyons(p):
            assert opts == []                                 # nothing to move
        for _, g in opts:
            assert col[p.lattice.index(g)] == col[k]          # same type
            assert g not in anyons(p)                         # nothing annihilated


# ---------------------------------------------------------------- refusals
@pytest.mark.parametrize("kind", ["majorana", "fermion", "anyon"])
def test_a_move_not_on_the_list_is_refused(kind):
    p = board()
    L = p.lattice
    if kind == "majorana":
        x, listing, move = majoranas(p)[0], majorana_moves, move_majorana
    elif kind == "fermion":
        x, listing, move = 0, fermion_moves, move_fermion
    else:
        x, listing, move = anyons(p)[0], anyon_moves, move_anyon
    allowed = {L.index(e) for e, _ in listing(p, x)}
    bad = next(e for e in range(len(L.edges)) if e not in allowed)
    with pytest.raises(ValueError):
        move(p, x, bad)
