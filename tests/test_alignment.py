"""Unit tests for needleman_wunsch() in src/align.py"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from align import needleman_wunsch

def identity(a1, a2):
    matches = sum(1 for x, y in zip(a1, a2) if x == y and x != '-')
    return 100 * matches / len(a1)

def test_identical_sequences_are_100_percent():
    a1, a2, score = needleman_wunsch("ACGT", "ACGT")
    assert a1 == "ACGT"
    assert a2 == "ACGT"
    assert identity(a1, a2) == 100.0

def test_single_mismatch_detected():
    a1, a2, score = needleman_wunsch("ACGT", "ACCT")
    assert identity(a1, a2) == 75.0

def test_completely_different_sequences_score_low():
    a1, a2, score = needleman_wunsch("AAAA", "TTTT")
    assert identity(a1, a2) == 0.0

def test_insertion_produces_gap():
    # seq2 has an extra base inserted in the middle
    a1, a2, score = needleman_wunsch("ACGT", "ACXGT")
    assert '-' in a1 or '-' in a2
    assert len(a1) == len(a2)

def test_empty_vs_empty():
    a1, a2, score = needleman_wunsch("", "")
    assert a1 == "" and a2 == ""

def test_alignment_output_always_equal_length():
    a1, a2, score = needleman_wunsch("AAGCT", "AGCTTT")
    assert len(a1) == len(a2)

def test_real_data_shape_sanity():
    # Mirrors the real scale check: two sequences differing by one
    # substitution should be scored higher than two sharing nothing in common
    a1, a2, score_similar = needleman_wunsch("ACGTACGTACGT", "ACGTACGTACGA")
    a3, a4, score_different = needleman_wunsch("ACGTACGTACGT", "TTTTTTTTTTTT")
    assert score_similar > score_different
