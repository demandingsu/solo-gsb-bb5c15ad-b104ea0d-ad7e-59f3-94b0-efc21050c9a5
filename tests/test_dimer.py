"""Unit tests for the core dimer-alignment logic."""
from app.dimer import (
    REASON_MAX_CONSECUTIVE,
    REASON_THREE_PRIME,
    are_complementary,
    assess,
    best_alignment,
    evaluate_alignment,
    screen_primers,
)


class TestComplementarity:
    def test_watson_crick_pairs(self):
        assert are_complementary("A", "T")
        assert are_complementary("T", "A")
        assert are_complementary("C", "G")
        assert are_complementary("G", "C")

    def test_non_pairs(self):
        assert not are_complementary("A", "A")
        assert not are_complementary("A", "C")
        assert not are_complementary("G", "T")

    def test_degenerate_bases_match_by_set_overlap(self):
        # R = A/G can pair with T (via A) and with C (via G).
        assert are_complementary("R", "T")
        assert are_complementary("R", "C")
        assert not are_complementary("R", "A")
        # N matches anything; W (A/T) pairs with W (A/T).
        assert are_complementary("N", "G")
        assert are_complementary("W", "W")
        assert not are_complementary("S", "W")  # CG vs AT share no pair

    def test_u_treated_as_t(self):
        assert are_complementary("A", "U")


class TestEvaluateAlignment:
    def test_perfect_reverse_complement(self):
        # rev("GCAT") == "TACG"; offset 0 aligns all four bases.
        aln = evaluate_alignment("ATGC", "TACG", 0)
        assert aln.max_consecutive == 4
        # The run touches both aligned 3' ends (A's last base, B's rev index 0).
        assert aln.three_prime_consecutive == 4

    def test_offset_shifts_overlap(self):
        # rev_b = "TACG"; offset 1 aligns seq_a[1:] ("ATGC") with all of rev_b.
        aln = evaluate_alignment("AATGC", "TACG", 1)
        assert aln.max_consecutive == 4
        assert aln.three_prime_consecutive == 4  # run reaches seq_a's 3' end
        # Without the shift only the first base pairs.
        aln0 = evaluate_alignment("AATGC", "TACG", 0)
        assert aln0.max_consecutive == 1

    def test_no_complementarity(self):
        aln = evaluate_alignment("AAAA", "AAAA", 0)
        assert aln.max_consecutive == 0
        assert aln.three_prime_consecutive == 0

    def test_three_prime_run_requires_end_involvement(self):
        # Internal run only: complementary bases away from any 3' end.
        aln = evaluate_alignment("TAT", "ATA", 0)
        # positions: T/A ok, A/T ok, T/A ok -> full run touches both ends.
        assert aln.three_prime_consecutive == 3
        aln2 = evaluate_alignment("GATG", "CTAC", 0)
        # G/C ok, A/T ok, T/A ok, G/C ok -> all four; ends included.
        assert aln2.max_consecutive == 4
        aln3 = evaluate_alignment("CATG", "GTAC", 0)
        # C/G ok, A/T ok, T/A ok, G/C ok
        assert aln3.max_consecutive == 4

    def test_internal_run_not_counted_as_three_prime(self):
        # "AATT" vs rev "AATT": offset 0 -> A/A x, A/A x, T/T x, T/T x: 0.
        # Build a case with an internal complementary stretch only:
        # seq_a = "CGAT", rev_b = "GCXX" is invalid IUPAC; use "TAGC":
        # offset 1: i=1 G vs rev[0]=T no; i=2 A vs rev[1]=A no; i=3 T vs rev[2]=G no
        aln = evaluate_alignment("CGAT", "TAGC", 1)
        assert aln.max_consecutive == 0
        # offset 0: C/T no, G/A no, A/G no, T/C no -> 0
        assert evaluate_alignment("CGAT", "TAGC", 0).max_consecutive == 0
        # offset -1: i=0 C vs rev[1]=A no; i=1 G vs rev[2]=G no; i=2 A vs rev[3]=C no
        assert evaluate_alignment("CGAT", "TAGC", -1).max_consecutive == 0
        # A genuinely internal run: seq_a="ACGT", rev_b="TGCA" offset 0:
        # A/T ok, C/G ok, G/C ok, T/A ok -> touches ends; instead shift:
        aln4 = evaluate_alignment("ACGT", "TGCA", 1)
        # i=1 C vs rev[0]=T no; i=2 G vs rev[1]=G no; i=3 T vs rev[2]=C no
        assert aln4.max_consecutive == 0


class TestBestAlignment:
    def test_finds_most_dangerous_offset(self):
        aln = best_alignment("ATGC", "GCAT")
        assert aln.offset == 0
        assert aln.max_consecutive == 4
        assert aln.three_prime_consecutive == 4

    def test_self_dimer_uses_non_negative_offsets_only(self):
        # For identical molecules offsets o and -o are the same physical
        # dimer; the search must stay in o >= 0.
        aln = best_alignment("AT", "AT", self_dimer=True)
        assert aln.offset >= 0
        assert aln.max_consecutive == 2

    def test_tie_broken_by_smallest_offset(self):
        # "AAAA" vs "AATT" (rev "TTAA"): offsets 0, 1 and 2 all give
        # (max_consecutive=2, three_prime_consecutive=2); the fixed rule
        # "smallest offset wins" must pick offset 0.
        aln = best_alignment("AAAA", "AATT")
        assert aln.offset == 0
        assert aln.max_consecutive == 2
        assert aln.three_prime_consecutive == 2

    def test_three_prime_breaks_max_consecutive_tie(self):
        # "TAT" vs "TAT" (rev "TAT"): offset -1 gives a 2-base internal run
        # (three_prime=0); offset 1 gives a 2-base run reaching primer A's
        # 3' end (three_prime=2). Equal max run -> the 3'-end run must win.
        aln = best_alignment("TAT", "TAT")
        assert aln.offset == 1
        assert aln.max_consecutive == 2
        assert aln.three_prime_consecutive == 2

    def test_single_base_primers(self):
        aln = best_alignment("A", "T")
        assert aln.max_consecutive == 1
        assert aln.three_prime_consecutive == 1


class TestAssess:
    def test_reasons_triggered_at_threshold(self):
        aln = best_alignment("ATGC", "GCAT")
        assert assess(aln, 4, 3) == [REASON_MAX_CONSECUTIVE, REASON_THREE_PRIME]
        assert assess(aln, 5, 5) == []
        assert assess(aln, 4, 5) == [REASON_MAX_CONSECUTIVE]
        assert assess(aln, 5, 4) == [REASON_THREE_PRIME]


class TestScreenPrimers:
    def test_pair_enumeration_order_and_count(self):
        results = screen_primers(
            [("p1", "ATGC"), ("p2", "GCAT"), ("p3", "AAAA")], 4, 3
        )
        # 3 self-dimers + 3 heterodimers, deterministic order.
        kinds = [(r["primer_a"], r["primer_b"], r["kind"]) for r in results]
        assert kinds == [
            ("p1", "p1", "self"),
            ("p1", "p2", "hetero"),
            ("p1", "p3", "hetero"),
            ("p2", "p2", "self"),
            ("p2", "p3", "hetero"),
            ("p3", "p3", "self"),
        ]

    def test_flagging(self):
        results = screen_primers([("fwd", "ATGC"), ("rev", "GCAT")], 4, 3)
        hetero = next(r for r in results if r["kind"] == "hetero")
        assert hetero["flagged"] is True
        assert hetero["reasons"] == [
            REASON_MAX_CONSECUTIVE,
            REASON_THREE_PRIME,
        ]
