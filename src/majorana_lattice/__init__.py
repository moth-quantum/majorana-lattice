"""majorana-lattice: matching codes, anyons and Majorana modes on trivalent lattices.

Implements the matching codes of J. R. Wootton, arXiv:1501.07779, simulated
with Stim. Start with:

    from majorana_lattice import MatchingCode, star
    code = MatchingCode(star(4, 4), seed=1)
"""
from .lattice import Edge, Face, Lattice, PRESETS, honeycomb, modified_honeycomb, star
from .matching_code import MatchingCode, Pair, link_pauli

__all__ = ["Edge", "Face", "Lattice", "PRESETS", "honeycomb", "modified_honeycomb", "star",
           "MatchingCode", "Pair", "link_pauli"]
__version__ = "0.1.0"
