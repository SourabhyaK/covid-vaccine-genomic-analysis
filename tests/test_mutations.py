"""Unit tests for apply_novavax_mutations() in src/build_proteins.py"""
import sys, os
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from build_proteins import apply_novavax_mutations

def make_fake_wt(length=1273):
    """Build a minimal fake WT sequence with RRAR at 682-685 and K,V at 986,987
    (1-indexed), padded with 'A' elsewhere - enough to exercise the mutation
    logic without needing the full 1273-residue real sequence in every test."""
    seq = list("A" * length)
    seq[681:685] = list("RRAR")   # positions 682-685
    seq[985] = 'K'                # position 986
    seq[986] = 'V'                # position 987
    return "".join(seq)

def test_furin_site_knocked_out():
    wt = make_fake_wt()
    result = apply_novavax_mutations(wt)
    assert result[681:685] == "QQAQ"

def test_2p_mutation_applied():
    wt = make_fake_wt()
    result = apply_novavax_mutations(wt)
    assert result[985] == 'P'  # K986P
    assert result[986] == 'P'  # V987P

def test_sequence_length_unchanged():
    wt = make_fake_wt()
    result = apply_novavax_mutations(wt)
    assert len(result) == len(wt)

def test_only_expected_positions_change():
    wt = make_fake_wt()
    result = apply_novavax_mutations(wt)
    diffs = [i for i in range(len(wt)) if wt[i] != result[i]]
    # RRAR -> QQAQ only actually changes R->Q at positions 682,683,685 (0-indexed
    # 681,682,684) since the 'A' at position 684 (0-indexed 683) is unchanged;
    # plus K986P/V987P (0-indexed 985,986) = 5 positions total
    assert sorted(diffs) == [681, 682, 684, 985, 986]

def test_raises_on_wrong_furin_site():
    wt = list(make_fake_wt())
    wt[681] = 'X'  # corrupt the expected RRAR motif
    with pytest.raises(AssertionError):
        apply_novavax_mutations("".join(wt))

def test_raises_on_wrong_2p_site():
    wt = list(make_fake_wt())
    wt[985] = 'X'  # corrupt the expected K at position 986
    with pytest.raises(AssertionError):
        apply_novavax_mutations("".join(wt))
