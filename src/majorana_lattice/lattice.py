"""Trivalent lattices with x/y/z link labels, for matching codes.

A lattice is any graph where every vertex has three edges, labelled x, y or z
so that no two edges of the same label meet at a vertex (Wootton 2015,
arXiv:1501.07779). Lattices live on a torus: each edge stores the displacement
from its first vertex to its second, so faces can be traced from geometry.

Presets: honeycomb, star (every honeycomb vertex replaced by a triangle) and
modified_honeycomb (the down sites replaced by triangles, as in the paper).

Preset elements have addresses: a unit cell (i, j) in axial coordinates on
the lattice vectors a1 = (sqrt 3, 0) and a2 = (sqrt 3 / 2, 3 / 2), and a name
within the cell. Each honeycomb cell has a down site D and an up site U one
unit above it. On the star lattice, cell (i, j) holds:

  Vertex(i, j, "D.x")   corner of triangle D that the x link leaves from
                        (D.x D.y D.z U.x U.y U.z); plain honeycomb sites are
                        just "D" and "U"
  Link(i, j, "z")       link D(i, j)-U(i, j); "x" is D(i, j)-U(i+1, j-1) and
                        "y" is D(i, j)-U(i, j-1)
  Link(i, j, "D.z")     the z edge of triangle D, between its x and y corners
  Plaquette(i, j, "D")  triangle D; "U" is triangle U, "hex" the large face
                        centred at (i + 1/3, j + 1/3)

Cells wrap round the torus: (i, j + rows) is the same cell as (i + rows/2, j),
because a2 leans right.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field

LABELS = ("x", "y", "z")
_S3 = math.sqrt(3)
BASIS = ((_S3, 0.0), (_S3 / 2, 1.5))


@dataclass(frozen=True)
class Vertex:
    i: int
    j: int
    name: str


@dataclass(frozen=True)
class Link:
    i: int
    j: int
    name: str


@dataclass(frozen=True)
class Plaquette:
    i: int
    j: int
    name: str


Address = Vertex | Link | Plaquette


@dataclass
class Edge:
    id: int
    u: int
    v: int
    label: str
    d: tuple[float, float]  # displacement from u to v (unwrapped)

    def other(self, w: int) -> int:
        return self.v if w == self.u else self.u


@dataclass
class Face:
    id: int
    vertices: list[int]
    edges: list[int]
    pauli: dict[int, str]  # the plaquette operator, qubit -> Pauli


@dataclass
class Lattice:
    positions: list[tuple[float, float]]
    edges: list[Edge]
    period: tuple[float, float]
    name: str = "custom"
    cells: tuple[int, int] | None = None        # (cols, rows); presets only
    vertex_addr: list[Vertex] | None = None     # vertex id -> address
    link_addr: list[Link] | None = None         # edge id -> address
    plaquette_addr: list[Plaquette] | None = field(default=None, init=False)
    faces: list[Face] = field(default_factory=list)
    incident: list[list[int]] = field(default_factory=list)

    @property
    def n(self) -> int:
        return len(self.positions)

    def __post_init__(self):
        self.incident = [[] for _ in range(self.n)]
        for e in self.edges:
            self.incident[e.u].append(e.id)
            self.incident[e.v].append(e.id)
        self.validate()
        self.faces = self._trace_faces()
        if self.cells is not None:
            self.plaquette_addr = self._address_faces()
            self._index = {}
            for kind in (self.vertex_addr, self.link_addr, self.plaquette_addr):
                for k, a in enumerate(kind):
                    if a in self._index:
                        raise ValueError(f"address {a} is used twice")
                    self._index[a] = k

    # -------------------------------------------------------------- addresses
    def wrap(self, i: int, j: int) -> tuple[int, int]:
        """The canonical cell for (i, j): 0 <= i < cols and 0 <= j < rows."""
        cols, rows = self._need_cells()
        return _wrap(i, j, cols, rows)

    def index(self, a: Address | int) -> int:
        """Vertex, edge or face id of an address; an int is returned as is."""
        if isinstance(a, int):
            return a
        self._need_cells()
        key = type(a)(*self.wrap(a.i, a.j), a.name)
        try:
            return self._index[key]
        except KeyError:
            raise KeyError(f"{self.name} has no {a}") from None

    def address(self, kind: type, k: int) -> Address:
        """Address of vertex, edge or face id k; kind is Vertex, Link or Plaquette."""
        self._need_cells()
        return {Vertex: self.vertex_addr, Link: self.link_addr, Plaquette: self.plaquette_addr}[kind][k]

    def coords(self, a: Address) -> tuple[float, float]:
        """Centre of an element in axial lattice coordinates, near its cell."""
        k = self.index(a)
        i, j = self.wrap(a.i, a.j)
        x, y = self._centre(type(a), k)
        cx, cy = to_xy(i, j)
        dx, dy = self._min_image(x - cx, y - cy)
        return to_axial(cx + dx, cy + dy)

    def nearest(self, kind: type, uv: tuple[float, float]) -> Address:
        """The element of this kind closest to axial coordinates uv."""
        x0, y0 = to_xy(*uv)
        count = {Vertex: self.n, Link: len(self.edges), Plaquette: len(self.faces)}[kind]
        best = min(range(count), key=lambda k: math.hypot(*self._min_image(
            self._centre(kind, k)[0] - x0, self._centre(kind, k)[1] - y0)))
        return self.address(kind, best)

    def _need_cells(self) -> tuple[int, int]:
        if self.cells is None:
            raise ValueError(f"lattice {self.name!r} has no cell addresses")
        return self.cells

    def _min_image(self, dx: float, dy: float) -> tuple[float, float]:
        W, H = self.period
        return (dx - W * round(dx / W), dy - H * round(dy / H))

    def _centre(self, kind: type, k: int) -> tuple[float, float]:
        """Centre of an element in Cartesian coordinates (not wrapped)."""
        if kind is Vertex:
            return self.positions[k]
        if kind is Link:
            E = self.edges[k]
            x, y = self.positions[E.u]
            return (x + E.d[0] / 2, y + E.d[1] / 2)
        pts = self._unwrapped_face(k)
        return (sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts))

    def _unwrapped_face(self, f: int) -> list[tuple[float, float]]:
        F = self.faces[f]
        x, y = self.positions[F.vertices[0]]
        pts = []
        for w, e in zip(F.vertices, F.edges):
            pts.append((x, y))
            E = self.edges[e]
            s = 1 if E.u == w else -1
            x, y = x + s * E.d[0], y + s * E.d[1]
        return pts

    def _address_faces(self) -> list[Plaquette]:
        """A triangle takes the name of the site it replaced; any other face is
        the cell's "hex", centred at (i + 1/3, j + 1/3)."""
        out = []
        for f, F in enumerate(self.faces):
            va = [self.vertex_addr[w] for w in F.vertices]
            sites = {(a.i, a.j, a.name.split(".")[0]) for a in va}
            if len(F.vertices) == 3 and len(sites) == 1 and "." in va[0].name:
                (i, j, s), = sites
                out.append(Plaquette(i, j, s))
                continue
            u, v = to_axial(*self._centre(Plaquette, f))
            i, j = round(u - 1 / 3), round(v - 1 / 3)
            if abs(u - 1 / 3 - i) > 1e-6 or abs(v - 1 / 3 - j) > 1e-6:
                raise ValueError(f"face {f} is not centred on a hexagon")
            out.append(Plaquette(*self.wrap(i, j), "hex"))
        return out

    # ----------------------------------------------------------------- checks
    def validate(self):
        for w, inc in enumerate(self.incident):
            if len(inc) != 3:
                raise ValueError(f"vertex {w} has {len(inc)} edges, need 3")
            labels = sorted(self.edges[e].label for e in inc)
            if labels != ["x", "y", "z"]:
                raise ValueError(f"vertex {w} has labels {labels}, need x, y, z")

    def edge_between(self, a: int, b: int) -> int | None:
        for e in self.incident[a]:
            if self.edges[e].other(a) == b:
                return e
        return None

    def edge_with_label(self, w: int, label: str) -> int:
        for e in self.incident[w]:
            if self.edges[e].label == label:
                return e
        raise KeyError((w, label))

    # ------------------------------------------------------------------ faces
    def _half_edges(self, w: int):
        """Half-edges leaving w as (angle, edge id, neighbour, displacement)."""
        out = []
        for e in self.incident[w]:
            E = self.edges[e]
            d = E.d if E.u == w else (-E.d[0], -E.d[1])
            out.append((math.atan2(d[1], d[0]), e, E.other(w), d))
        out.sort()
        return out

    def _trace_faces(self) -> list[Face]:
        rot = {w: self._half_edges(w) for w in range(self.n)}
        seen: set[tuple[int, int]] = set()  # (edge, from vertex)
        faces = []
        for w0 in range(self.n):
            for _, e0, _, _ in rot[w0]:
                if (e0, w0) in seen:
                    continue
                verts, edges = [], []
                w, e = w0, e0
                while (e, w) not in seen:
                    seen.add((e, w))
                    verts.append(w)
                    edges.append(e)
                    nxt = self.edges[e].other(w)
                    # at nxt, turn to the half-edge just clockwise of the way back
                    hs = rot[nxt]
                    k = next(i for i, h in enumerate(hs) if h[1] == e and h[2] == w)
                    e = hs[(k - 1) % len(hs)][1]
                    w = nxt
                faces.append((verts, edges))
        if len(faces) * 2 != self.n:
            raise ValueError(f"found {len(faces)} faces, expected {self.n // 2} on a torus")
        result = []
        for fid, (verts, edges) in enumerate(faces):
            eset = set(edges)
            pauli = {}
            for w in verts:
                # plaquette = product of its link operators = the Pauli of the
                # one edge at each vertex that is not part of the face
                (outer,) = [e for e in self.incident[w] if e not in eset]
                pauli[w] = self.edges[outer].label.upper()
            result.append(Face(fid, verts, edges, pauli))
        return result

    def face_neighbours(self) -> dict[int, list[tuple[int, int]]]:
        """face -> [(neighbour face, shared edge)]"""
        owner: dict[int, list[int]] = {}
        for f in self.faces:
            for e in f.edges:
                owner.setdefault(e, []).append(f.id)
        nb: dict[int, list[tuple[int, int]]] = {f.id: [] for f in self.faces}
        for e, fs in owner.items():
            if len(fs) == 2 and fs[0] != fs[1]:
                nb[fs[0]].append((fs[1], e))
                nb[fs[1]].append((fs[0], e))
        return nb

    def to_json(self) -> dict:
        return {
            "name": self.name,
            "period": self.period,
            "positions": self.positions,
            "edges": [[e.u, e.v, e.label, e.d[0], e.d[1]] for e in self.edges],
            "faces": [{"vertices": f.vertices, "edges": f.edges,
                       "pauli": {str(k): v for k, v in f.pauli.items()}} for f in self.faces],
        }


# ------------------------------------------------------------------- presets
def to_xy(u: float, v: float) -> tuple[float, float]:
    """Axial lattice coordinates to Cartesian."""
    return (u * BASIS[0][0] + v * BASIS[1][0], v * BASIS[1][1])


def to_axial(x: float, y: float) -> tuple[float, float]:
    """Cartesian coordinates to axial lattice coordinates."""
    v = y / BASIS[1][1]
    return ((x - v * BASIS[1][0]) / BASIS[0][0], v)


def _wrap(i: int, j: int, cols: int, rows: int) -> tuple[int, int]:
    # wrap j first; (i, j + rows) is (i + rows/2, j) because a2 leans right
    q, j = divmod(j, rows)
    return ((i + q * (rows // 2)) % cols, j)


def _honeycomb_graph(cols: int, rows: int):
    """Honeycomb with vertical z links on a cols x rows torus.

    Returns sites [(D or U, i, j, pos)], edges [(a, b, label, disp)] and the
    period. Edge k lives in cell sites[a]'s cell and is named by its label."""
    if cols < 2 or rows < 2:
        raise ValueError("need at least 2 x 2 cells")
    period = (cols * _S3, rows * 1.5)
    if rows % 2:
        raise ValueError("rows must be even so the torus is rectangular")

    index, sites = {}, []
    for j in range(rows):
        for i in range(cols):
            for s in "DU":
                x, y = to_xy(i, j)
                y += 0.0 if s == "D" else 1.0
                x %= period[0]
                index[(s, i, j)] = len(sites)
                sites.append((s, i, j, (x, y)))

    def site(s, i, j):
        return index[(s, *_wrap(i, j, cols, rows))]

    edges = []
    for j in range(rows):
        for i in range(cols):
            D = index[("D", i, j)]
            edges.append((D, site("U", i, j), "z", (0.0, 1.0)))
            edges.append((D, site("U", i + 1, j - 1), "x", (_S3 / 2, -0.5)))
            edges.append((D, site("U", i, j - 1), "y", (-_S3 / 2, -0.5)))
    return sites, edges, period


def honeycomb(cols: int = 4, rows: int = 4) -> Lattice:
    sites, edges, period = _honeycomb_graph(cols, rows)
    E = [Edge(k, a, b, lab, d) for k, (a, b, lab, d) in enumerate(edges)]
    return Lattice([s[3] for s in sites], E, period, name=f"honeycomb {cols}x{rows}",
                   cells=(cols, rows),
                   vertex_addr=[Vertex(i, j, s) for s, i, j, _ in sites],
                   link_addr=[Link(*sites[a][1:3], lab) for a, _, lab, _ in edges])


def _triangulate(cols, rows, replace, t=0.28):
    """Replace chosen honeycomb vertices with triangles, keeping labels valid.

    A corner of a triangle takes the label of the honeycomb edge it leads to;
    the triangle edge between two corners takes the remaining third label.
    Returns the keyword arguments for Lattice, addresses included."""
    sites, edges, period = _honeycomb_graph(cols, rows)
    nbrs = {k: [] for k in range(len(sites))}
    for a, b, lab, d in edges:
        nbrs[a].append((lab, d))
        nbrs[b].append((lab, (-d[0], -d[1])))

    pos, corner, vaddr = [], {}, []
    for k, (s, i, j, (x, y)) in enumerate(sites):
        if replace(s):
            for lab, d in nbrs[k]:
                corner[(k, lab)] = len(pos)
                pos.append(((x + t * d[0]) % period[0], (y + t * d[1]) % period[1]))
                vaddr.append(Vertex(i, j, f"{s}.{lab}"))
        else:
            vid = len(pos)
            pos.append((x, y))
            vaddr.append(Vertex(i, j, s))
            for lab, _ in nbrs[k]:
                corner[(k, lab)] = vid

    def offset(k, lab):
        if not replace(sites[k][0]):
            return (0.0, 0.0)
        d = dict(nbrs[k])[lab]
        return (t * d[0], t * d[1])

    E, laddr = [], []
    for a, b, lab, d in edges:
        oa, ob = offset(a, lab), offset(b, lab)
        disp = (d[0] + ob[0] - oa[0], d[1] + ob[1] - oa[1])
        E.append(Edge(len(E), corner[(a, lab)], corner[(b, lab)], lab, disp))
        laddr.append(Link(*sites[a][1:3], lab))
    for k, (s, i, j, _) in enumerate(sites):
        if not replace(s):
            continue
        labs = [lab for lab, _ in nbrs[k]]
        for p in range(3):
            for q in range(p + 1, 3):
                la, lb = labs[p], labs[q]
                (third,) = set(LABELS) - {la, lb}
                oa, ob = offset(k, la), offset(k, lb)
                E.append(Edge(len(E), corner[(k, la)], corner[(k, lb)], third,
                              (ob[0] - oa[0], ob[1] - oa[1])))
                laddr.append(Link(i, j, f"{s}.{third}"))
    return dict(positions=pos, edges=E, period=period, cells=(cols, rows),
                vertex_addr=vaddr, link_addr=laddr)


def star(cols: int = 4, rows: int = 4) -> Lattice:
    """Star lattice (3.12.12): every honeycomb vertex becomes a triangle."""
    return Lattice(**_triangulate(cols, rows, lambda s: True), name=f"star {cols}x{rows}")


def modified_honeycomb(cols: int = 4, rows: int = 4) -> Lattice:
    """The paper's lattice: the down sites become triangles."""
    return Lattice(**_triangulate(cols, rows, lambda s: s == "D"),
                   name=f"modified honeycomb {cols}x{rows}")


PRESETS = {"honeycomb": honeycomb, "star": star, "modified_honeycomb": modified_honeycomb}


if __name__ == "__main__":
    for name, fn in PRESETS.items():
        L = fn(4, 4)
        sizes = sorted({len(f.vertices) for f in L.faces})
        print(f"{L.name}: {L.n} qubits, {len(L.edges)} links, {len(L.faces)} faces, face sizes {sizes}")
    print(json.dumps(star(2, 2).to_json())[:200], "...")
