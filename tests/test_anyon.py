"""Moving e and m anyons."""
import pytest

from majorana_lattice import MatchingCode, star


def _anyons(c):
    v = c.values()["plaquettes"]
    return [f for f, x in enumerate(v) if x == -1]


def _with_e(seed=1, z_edge=0, end="u"):
    """A board with an e and an m, made by X on a qubit of a z link."""
    L = star(4, 4)
    c = MatchingCode(L, seed=seed)
    z = [E for E in L.edges if E.label == "z"][z_edge]
    c.apply_pauli(z.u if end == "u" else z.v, "X")
    col, _ = c.colouring()
    return c, next(f for f in _anyons(c) if col[f] == 0)


def test_moves_only_the_anyon():
    c, e = _with_e()
    fermions = c.summary()["fermions"]
    f = e
    for _ in range(12):
        opts = c.anyon_moves(f)
        assert opts, f"the anyon on {f} has nowhere to go"
        before = _anyons(c)
        col = c.colouring()[0]
        link, dest = opts[0]
        assert c.move_anyon(f, link) == dest
        assert sorted(_anyons(c)) == sorted([a for a in before if a != f] + [dest])
        assert c.colouring()[0] == col
        assert c.colouring()[0][dest] == 0
        assert c.summary()["fermions"] == fermions
        f = dest


def test_only_allowed_cells_are_offered():
    c, e = _with_e()
    col, _ = c.colouring()
    for f in range(len(c.lat.faces)):
        opts = c.anyon_moves(f)
        if f not in _anyons(c):
            assert opts == []          # nothing to move
        for _, g in opts:
            assert col[g] == col[f]    # same type
            assert g not in _anyons(c)  # nothing annihilated


def test_refuses_a_disallowed_step():
    c, e = _with_e()
    allowed = {link for link, _ in c.anyon_moves(e)}
    bad = next(link for _, link in c.lat.face_neighbours()[e] if link not in allowed)
    log = len(c.log)
    with pytest.raises(ValueError):
        c.move_anyon(e, bad)
    assert len(c.log) == log           # and nothing was applied


def test_replays_from_the_log():
    c, f = _with_e()
    for _ in range(5):
        link, dest = c.anyon_moves(f)[0]
        c.move_anyon(f, link)
        f = dest
    again = MatchingCode.replay(c.lat, c.record())
    assert again.values()["plaquettes"] == c.values()["plaquettes"]
