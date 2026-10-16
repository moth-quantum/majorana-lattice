"""Cell addresses: every element has one, and they follow the lattice vectors."""
import math

import pytest

from majorana_lattice import PRESETS, Link, Plaquette, Vertex, honeycomb, star

SIZES = [(4, 4), (4, 6), (6, 4), (3, 2)]


def all_addresses(L):
    return ([(Vertex, L.address(Vertex, k)) for k in range(L.n)]
            + [(Link, L.address(Link, k)) for k in range(len(L.edges))]
            + [(Plaquette, L.address(Plaquette, k)) for k in range(len(L.faces))])


@pytest.mark.parametrize("name", PRESETS)
@pytest.mark.parametrize("size", SIZES)
def test_round_trip(name, size):
    L = PRESETS[name](*size)
    for kind, a in all_addresses(L):
        assert L.address(kind, L.index(a)) == a
        assert (a.i, a.j) == L.wrap(a.i, a.j)


@pytest.mark.parametrize("name", PRESETS)
@pytest.mark.parametrize("size", SIZES)
def test_links_match_positions(name, size):
    # every link's displacement agrees with its endpoints' positions on the torus
    L = PRESETS[name](*size)
    for E in L.edges:
        (xu, yu), (xv, yv) = L.positions[E.u], L.positions[E.v]
        dx, dy = L._min_image(xv - xu - E.d[0], yv - yu - E.d[1])
        assert math.hypot(dx, dy) < 1e-9, E


@pytest.mark.parametrize("size", SIZES)
def test_wrapping(size):
    cols, rows = size
    L = star(cols, rows)
    v = Vertex(1, 1, "U.x")
    assert L.index(Vertex(1 + cols, 1, "U.x")) == L.index(v)
    assert L.index(Vertex(1 - rows // 2, 1 + rows, "U.x")) == L.index(v)
    assert L.index(Vertex(1 + rows // 2, 1 - rows, "U.x")) == L.index(v)


def test_star_cell():
    L = star(4, 4)
    ends = lambda a: {L.edges[L.index(a)].u, L.edges[L.index(a)].v}
    V = lambda i, j, n: L.index(Vertex(i, j, n))
    # links between triangles
    assert ends(Link(2, 3, "z")) == {V(2, 3, "D.z"), V(2, 3, "U.z")}
    assert ends(Link(2, 3, "x")) == {V(2, 3, "D.x"), V(3, 2, "U.x")}
    assert ends(Link(2, 3, "y")) == {V(2, 3, "D.y"), V(2, 2, "U.y")}
    # a triangle's edge joins its other two corners, and carries its label
    assert ends(Link(2, 3, "D.z")) == {V(2, 3, "D.x"), V(2, 3, "D.y")}
    assert L.edges[L.index(Link(2, 3, "D.z"))].label == "z"
    # triangle plaquettes, and the up triangle is above the down one
    for s in "DU":
        f = L.faces[L.index(Plaquette(2, 3, s))]
        assert set(f.vertices) == {V(2, 3, f"{s}.{p}") for p in "xyz"}
    assert L.coords(Plaquette(2, 3, "U"))[1] > L.coords(Plaquette(2, 3, "D"))[1]
    # one hex per cell, centred at (i + 1/3, j + 1/3)
    u, v = L.coords(Plaquette(2, 3, "hex"))
    assert abs(u - (2 + 1 / 3)) < 1e-9 and abs(v - (3 + 1 / 3)) < 1e-9
    assert len(L.faces[L.index(Plaquette(2, 3, "hex"))].vertices) == 12


def test_honeycomb_cell():
    L = honeycomb(4, 4)
    E = L.edges[L.index(Link(1, 2, "x"))]
    assert {E.u, E.v} == {L.index(Vertex(1, 2, "D")), L.index(Vertex(2, 1, "U"))}


@pytest.mark.parametrize("name", PRESETS)
def test_coords_and_nearest(name):
    L = PRESETS[name](4, 6)
    for kind, a in all_addresses(L):
        u, v = L.coords(a)
        assert abs(u - a.i) < 1.5 and abs(v - a.j) < 1.5     # near its own cell
        assert L.nearest(kind, (u, v)) == a


def test_unknown_address():
    with pytest.raises(KeyError):
        star(4, 4).index(Vertex(0, 0, "D"))          # star corners are D.x etc.
