"""Unit tests for translate() in src/build_proteins.py"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from build_proteins import translate, find_orf, read_fasta

def test_single_codon():
    assert translate("ATG") == "M"

def test_known_short_peptide():
    # ATG GTG TTC -> Met-Val-Phe
    assert translate("ATGGTGTTC") == "MVF"

def test_stops_at_stop_codon():
    # ATG then TAA (stop) then more sequence that should be ignored
    assert translate("ATGTAAGGGGGG") == "M"

def test_handles_all_three_stop_codons():
    for stop in ["TAA", "TAG", "TGA"]:
        assert translate("ATG" + stop + "AAA") == "M"

def test_empty_sequence():
    assert translate("") == ""

def test_incomplete_trailing_codon_ignored():
    # 4 nt - the first full codon (ATG) still translates; the trailing
    # single "A" is not a full codon and is silently dropped, not an error
    assert translate("ATGA") == "M"
    assert translate("ATGAAAT") == "MK"  # trailing "T" alone is dropped

def test_find_orf_locates_first_atg():
    seq = "GGGGCCCATGAAATTT"
    assert find_orf(seq) == "ATGAAATTT"

def test_find_orf_no_atg_returns_minus_one_slice():
    # find() returns -1 if absent; seq[-1:] is the last character -
    # documenting this known edge-case behavior rather than hiding it
    seq = "GGGCCC"
    result = find_orf(seq)
    assert result == seq[-1:]

def test_read_fasta_parses_multi_record_file(tmp_path):
    fasta_content = ">seq1\nATGAAA\n>seq2\nTTTGGG\nCCC\n"
    f = tmp_path / "test.fasta"
    f.write_text(fasta_content)
    result = read_fasta(str(f))
    assert result == {"seq1": "ATGAAA", "seq2": "TTTGGGCCC"}
