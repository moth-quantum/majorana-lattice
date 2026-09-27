"""Trivalent lattices with x/y/z link labels, for matching codes.

A lattice is any graph where every vertex has three edges, labelled x, y or z
so that no two edges of the same label meet at a vertex (Wootton 2015,
arXiv:1501.07779). Lattices live on a torus: each edge stores the displacement
from its first vertex to its second, so faces can be traced from geometry.

Presets: honeycomb, star (every honeycomb vertex replaced by a triangle) and
modified_honeycomb (one sublattice replaced by triangles, as in the paper).
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field

LABELS = ("x", "y", "z")
_S3 = math.sqrt(3)


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
def _honeycomb_graph(cols: int, rows: int):
    """Honeycomb with vertical z links on a cols x rows torus.

    Returns sites [(sublattice, i, j, pos)] and edges [(a, b, label, disp)]."""
    if cols < 2 or rows < 2:
        raise ValueError("need at least 2 x 2 cells")
    a1, a2 = (_S3, 0.0), (_S3 / 2, 1.5)
    period = (cols * _S3, rows * 1.5)
    if rows % 2:
        raise ValueError("rows must be even so the torus is rectangular")

    def cell_pos(i, j):
        return (i * a1[0] + j * a2[0], j * a2[1])

    index, sites = {}, []
    for j in range(rows):
        for i in range(cols):
            for s in "AB":
                x, y = cell_pos(i, j)
                y += 0.0 if s == "A" else 1.0
                x %= period[0]
                index[(s, i, j)] = len(sites)
                sites.append((s, i, j, (x, y)))

    def site(s, i, j):
        # wrap j first; moving across the top shifts i because a2 leans right
        shift = 0
        while j >= rows:
            j -= rows
            shift -= rows // 2
        while j < 0:
            j += rows
            shift += rows // 2
        return index[(s, (i + shift) % cols, j)]

    edges = []
    for j in range(rows):
        for i in range(cols):
            A = index[("A", i, j)]
            edges.append((A, site("B", i, j), "z", (0.0, 1.0)))
            edges.append((A, site("B", i + 1, j - 1), "x", (_S3 / 2, -0.5)))
            edges.append((A, site("B", i, j - 1), "y", (-_S3 / 2, -0.5)))
    return sites, edges, period


def honeycomb(cols: int = 4, rows: int = 4) -> Lattice:
    sites, edges, period = _honeycomb_graph(cols, rows)
    E = [Edge(k, a, b, lab, d) for k, (a, b, lab, d) in enumerate(edges)]
    return Lattice([s[3] for s in sites], E, period, name=f"honeycomb {cols}x{rows}")


def _triangulate(cols, rows, replace, t=0.28, scale=1.0):
    """Replace chosen honeycomb vertices with triangles, keeping labels valid.

    A corner of a triangle takes the label of the honeycomb edge it leads to;
    the triangle edge between two corners takes the remaining third label."""
    sites, edges, period = _honeycomb_graph(cols, rows)
    nbrs = {k: [] for k in range(len(sites))}
    for a, b, lab, d in edges:
        nbrs[a].append((lab, d))
        nbrs[b].append((lab, (-d[0], -d[1])))

    pos, corner = [], {}
    for k, s in enumerate(sites):
        x, y = s[3]
        if replace(s[0]):
            for lab, d in nbrs[k]:
                corner[(k, lab)] = len(pos)
                pos.append(((x + t * d[0]) % period[0], (y + t * d[1]) % period[1]))
        else:
            vid = len(pos)
            pos.append((x, y))
            for lab, _ in nbrs[k]:
                corner[(k, lab)] = vid

    def offset(k, lab):
        if not replace(sites[k][0]):
            return (0.0, 0.0)
        d = dict(nbrs[k])[lab]
        return (t * d[0], t * d[1])

    E = []
    for a, b, lab, d in edges:
        oa, ob = offset(a, lab), offset(b, lab)
        disp = (d[0] + ob[0] - oa[0], d[1] + ob[1] - oa[1])
        E.append(Edge(len(E), corner[(a, lab)], corner[(b, lab)], lab, disp))
    for k, s in enumerate(sites):
        if not replace(s[0]):
            continue
        labs = [lab for lab, _ in nbrs[k]]
        for p in range(3):
            for q in range(p + 1, 3):
                la, lb = labs[p], labs[q]
                (third,) = set(LABELS) - {la, lb}
                oa, ob = offset(k, la), offset(k, lb)
                E.append(Edge(len(E), corner[(k, la)], corner[(k, lb)], third,
                              (ob[0] - oa[0], ob[1] - oa[1])))
    s = scale
    return ([(x * s, y * s) for x, y in pos],
            [Edge(e.id, e.u, e.v, e.label, (e.d[0] * s, e.d[1] * s)) for e in E],
            (period[0] * s, period[1] * s))


def star(cols: int = 4, rows: int = 4) -> Lattice:
    """Star lattice (3.12.12): every honeycomb vertex becomes a triangle."""
    pos, E, period = _triangulate(cols, rows, lambda s: True)
    return Lattice(pos, E, period, name=f"star {cols}x{rows}")


def modified_honeycomb(cols: int = 4, rows: int = 4) -> Lattice:
    """The paper's lattice: A-sublattice vertices become triangles."""
    pos, E, period = _triangulate(cols, rows, lambda s: s == "A")
    return Lattice(pos, E, period, name=f"modified honeycomb {cols}x{rows}")


PRESETS = {"honeycomb": honeycomb, "star": star, "modified_honeycomb": modified_honeycomb}


if __name__ == "__main__":
    for name, fn in PRESETS.items():
        L = fn(4, 4)
        sizes = sorted({len(f.vertices) for f in L.faces})
        print(f"{L.name}: {L.n} qubits, {len(L.edges)} links, {len(L.faces)} faces, face sizes {sizes}")
    print(json.dumps(star(2, 2).to_json())[:200], "...")
