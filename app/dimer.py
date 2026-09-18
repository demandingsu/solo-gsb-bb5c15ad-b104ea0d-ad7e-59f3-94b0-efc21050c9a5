"""Core primer-dimer alignment logic.

Conventions
-----------
* Primer sequences are stored 5' -> 3'.
* A dimer forms antiparallel, so primer B is reversed before alignment:
  ``rev_b = seq_b[::-1]`` (index 0 of ``rev_b`` is therefore B's 3' end).
* ``offset`` is the index in ``rev_b`` that aligns with ``seq_a[0]`` minus
  zero, i.e. ``rev_b`` index ``j`` pairs with ``seq_a`` index ``i`` where
  ``j = i - offset``.  Valid offsets range from ``-(len(rev_b) - 1)`` to
  ``len(seq_a) - 1``.
* For self-dimers only offsets >= 0 are enumerated: swapping two identical
  molecules maps offset ``o`` to ``-o``, so negative offsets duplicate
  physical dimers already seen.
"""
from __future__ import annotations

from dataclasses import dataclass

# IUPAC degenerate base symbols -> the concrete bases they represent.
IUPAC_BASES: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "U": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

_WATSON_CRICK = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Risk reason codes reported by the API.
REASON_MAX_CONSECUTIVE = "MAX_CONSECUTIVE_COMPLEMENTARITY"
REASON_THREE_PRIME = "THREE_PRIME_END_COMPLEMENTARITY"


def are_complementary(base_a: str, base_b: str) -> bool:
    """Return True when two IUPAC bases can form a Watson-Crick pair.

    Degenerate bases are compared as sets: the pair is complementary if at
    least one concrete base from each set can pair (A-T or C-G).
    """
    set_a = IUPAC_BASES[base_a]
    set_b = IUPAC_BASES[base_b]
    return any(_WATSON_CRICK[a] in set_b for a in set_a)


@dataclass(frozen=True)
class Alignment:
    """Statistics for one antiparallel alignment of two primers."""

    offset: int
    max_consecutive: int  # longest run of consecutive complementary bases
    three_prime_consecutive: int  # longest complementary run touching an aligned 3' end


def evaluate_alignment(seq_a: str, rev_b: str, offset: int) -> Alignment:
    """Evaluate a single alignment of ``seq_a`` against reversed ``seq_b``."""
    n, m = len(seq_a), len(rev_b)
    lo = max(0, offset)
    hi = min(n, offset + m)  # exclusive
    max_run = 0
    three_prime_run = 0
    run = 0
    run_touches_3prime = False
    for i in range(lo, hi):
        j = i - offset
        if are_complementary(seq_a[i], rev_b[j]):
            run += 1
            # i == n - 1: 3' end of primer A; j == 0: 3' end of primer B.
            if i == n - 1 or j == 0:
                run_touches_3prime = True
        else:
            max_run = max(max_run, run)
            if run_touches_3prime:
                three_prime_run = max(three_prime_run, run)
            run = 0
            run_touches_3prime = False
    max_run = max(max_run, run)
    if run_touches_3prime:
        three_prime_run = max(three_prime_run, run)
    return Alignment(
        offset=offset,
        max_consecutive=max_run,
        three_prime_consecutive=three_prime_run,
    )


def _rank(alignment: Alignment) -> tuple[int, int, int]:
    """Danger ranking; higher wins.

    Ties in (max_consecutive, three_prime_consecutive) are broken by the
    fixed rule "smallest offset wins" (offset enters negated), which makes
    the reported alignment unique and deterministic.
    """
    return (
        alignment.max_consecutive,
        alignment.three_prime_consecutive,
        -alignment.offset,
    )


def best_alignment(seq_a: str, seq_b: str, self_dimer: bool = False) -> Alignment:
    """Return the most dangerous alignment across all offsets."""
    rev_b = seq_b[::-1]
    n, m = len(seq_a), len(rev_b)
    lo = 0 if self_dimer else -(m - 1)
    best: Alignment | None = None
    for offset in range(lo, n):
        candidate = evaluate_alignment(seq_a, rev_b, offset)
        if best is None or _rank(candidate) > _rank(best):
            best = candidate
    assert best is not None  # ranges are never empty for non-empty primers
    return best


def assess(
    alignment: Alignment,
    max_consecutive_threshold: int,
    three_prime_threshold: int,
) -> list[str]:
    """Return the list of triggered risk reason codes (empty means safe)."""
    reasons: list[str] = []
    if alignment.max_consecutive >= max_consecutive_threshold:
        reasons.append(REASON_MAX_CONSECUTIVE)
    if alignment.three_prime_consecutive >= three_prime_threshold:
        reasons.append(REASON_THREE_PRIME)
    return reasons


def screen_primers(
    primers: list[tuple[str, str]],
    max_consecutive_threshold: int,
    three_prime_threshold: int,
) -> list[dict]:
    """Screen every self-dimer and pairwise heterodimer.

    ``primers`` is a list of ``(name, sequence)`` tuples in submission
    order.  Returns one result dict per pair, deterministically ordered:
    for each primer its self-dimer first, then heterodimers with the
    primers that follow it.
    """
    results: list[dict] = []
    for i, (name_a, seq_a) in enumerate(primers):
        results.append(
            _screen_pair(
                name_a,
                seq_a,
                name_a,
                seq_a,
                "self",
                max_consecutive_threshold,
                three_prime_threshold,
            )
        )
        for j in range(i + 1, len(primers)):
            name_b, seq_b = primers[j]
            results.append(
                _screen_pair(
                    name_a,
                    seq_a,
                    name_b,
                    seq_b,
                    "hetero",
                    max_consecutive_threshold,
                    three_prime_threshold,
                )
            )
    return results


def _screen_pair(
    name_a: str,
    seq_a: str,
    name_b: str,
    seq_b: str,
    kind: str,
    max_consecutive_threshold: int,
    three_prime_threshold: int,
) -> dict:
    alignment = best_alignment(seq_a, seq_b, self_dimer=(kind == "self"))
    reasons = assess(alignment, max_consecutive_threshold, three_prime_threshold)
    return {
        "primer_a": name_a,
        "primer_b": name_b,
        "kind": kind,
        "offset": alignment.offset,
        "max_consecutive": alignment.max_consecutive,
        "three_prime_consecutive": alignment.three_prime_consecutive,
        "flagged": bool(reasons),
        "reasons": reasons,
    }
