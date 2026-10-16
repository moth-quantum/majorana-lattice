"""majorana-lattice: matching codes, anyons and Majorana modes on trivalent lattices.

Implements the matching codes of J. R. Wootton, arXiv:1501.07779, as Stim
circuits. Build a program without running it, then run it:

    from majorana_lattice import Link, prepare, run, split_dimer, star
    p = split_dimer(prepare(star(4, 4)), Link(0, 1, "z"))
    res = run(p, shots=1000)

MatchingCode is the interactive engine: it simulates as it goes, one shot at
a time, and is the reference the circuit engine is tested against.
"""
from .lattice import (Edge, Face, Lattice, Link, Plaquette, PRESETS, Vertex, honeycomb,
                      modified_honeycomb, star, to_axial, to_xy)
from .matching_code import MatchingCode, Pair, link_pauli
from .program import Program, prepare
from .ops import (apply_link, apply_pauli, create_anyons, create_majoranas, from_circuit, fuse,
                  measure_label, measure_link, move_anyon, move_fermion, move_majorana, readout,
                  release, route_majorana, split_dimer)
from .queries import (anyon_moves, colouring, dimers, fermion_moves, is_dimer, is_majorana,
                      majorana_moves, majoranas, partner, plan_route, plaquette_values)
from .run import Result, archaeo_executor, run, stim_executor

__all__ = ["Edge", "Face", "Lattice", "Link", "Plaquette", "PRESETS", "Vertex", "honeycomb",
           "modified_honeycomb", "star", "to_axial", "to_xy",
           "MatchingCode", "Pair", "link_pauli",
           "Program", "prepare",
           "apply_link", "apply_pauli", "create_anyons", "create_majoranas", "from_circuit", "fuse",
           "measure_label", "measure_link", "move_anyon", "move_fermion", "move_majorana", "readout",
           "release", "route_majorana", "split_dimer",
           "anyon_moves", "colouring", "dimers", "fermion_moves", "is_dimer", "is_majorana",
           "majorana_moves", "majoranas", "partner", "plan_route", "plaquette_values",
           "Result", "archaeo_executor", "run", "stim_executor"]
__version__ = "0.1.0"
