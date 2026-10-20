"""Majoranas, fermions and anyons share one move API: a list of the moves on
offer, a function that makes one, and a ValueError for anything else."""
import numpy as np
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


def one_move(fix_value):
    """Split a dimer and move one of its Majoranas with this fix_value."""
    p = split_dimer(prepare(star(4, 4)), Link(1, 1, "z"))
    v = majoranas(p)[0]
    link, _ = majorana_moves(p, v)[0]
    return move_majorana(p, v, link, fix_value=fix_value)


@pytest.mark.parametrize("fix_value, empty", [(1, 1.0), (-1, 0.0)])
def test_fix_value_one_empties_and_minus_one_fills(fix_value, empty):
    s = run(one_move(fix_value), shots=300, seed=1).summary()["fermions"]
    assert (s == (0 if empty else 2)).all()        # the new dimer and the far pair together


@pytest.mark.parametrize("fix_value", [0, 0.5, -0.8])
def test_fix_value_is_the_expected_value(fix_value):
    shots = 20000
    s = run(one_move(fix_value), shots=shots, seed=2).summary()["fermions"]
    assert set(s) <= {0, 2}
    mean = 1 - (s == 2).mean() * 2                  # the new dimer's mean value
    assert abs(mean - fix_value) < 0.03


def test_the_new_dimer_is_the_coin():
    """Whatever the measurement gave, the dimer's value is the coin's record,
    so the result is independent of the outcome."""
    p = one_move(0)
    res = run(p, shots=2000, seed=3)
    link = [k for k, m in enumerate(p.measurements) if m[0] == "link"][-1]
    coin = [k for k, m in enumerate(p.measurements) if m[0] == "coin"][-1]
    fermions = res.summary()["fermions"]
    assert ((fermions == 2) == res.bits[:, coin]).all()
    outcome = res.bits[:, link]
    assert 0.4 < outcome.mean() < 0.6 and abs(np.corrcoef(outcome, res.bits[:, coin])[0, 1]) < 0.1


def test_every_fix_uses_the_coin():
    """The same instructions for every fix_value: measure, fix, coin, flip."""
    for fix_value in (-1, 0, 0.3, 1):
        p = one_move(fix_value)
        assert [inst.tag.split(":")[0] for inst in p.instructions[-5:]] == ["link", "fix", "coin", "coin", "flip"]
        assert [inst.name for inst in p.instructions[-3:-1]] == ["R", "HERALDED_PAULI_CHANNEL_1"]
        coin = p.instructions[-2]
        assert coin.gate_args_copy() == [(1 - fix_value) / 2, 0, 0, 0]
        assert coin.targets_copy()[0].value == p.lattice.n        # the extra qubit


@pytest.mark.parametrize("bad, error", [(1.5, ValueError), (-2, ValueError), (True, TypeError), ("1", TypeError)])
def test_fix_value_must_be_a_number_from_minus_one_to_one(bad, error):
    with pytest.raises(error):
        one_move(bad)


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
