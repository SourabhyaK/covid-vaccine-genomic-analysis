"""
Build all 4 real vaccine antigen (spike) protein sequences from verified sources.
See README.md for full data provenance and citations.

Refactored into pytest-testable functions: translate(), read_fasta(), find_orf(),
and apply_novavax_mutations() are each independently unit tested in tests/.
"""

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def translate(cds):
    """Translate a nucleotide CDS to protein using the standard genetic code.
    Stops at the first stop codon (not included in output)."""
    protein = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        protein.append(aa)
    return "".join(protein)

def read_fasta(path):
    seqs = {}; name = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]; seqs[name] = []
            elif line:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}

def find_orf(seq):
    """Find the first ATG and return the sequence from there to the end."""
    idx = seq.find("ATG")
    return seq[idx:]

def apply_novavax_mutations(wt_protein):
    """Construct the Novavax NVX-CoV2373 antigen from the WT reference spike
    by applying the two documented, cited mutations:
      - furin cleavage site knockout: positions 682-685, RRAR -> QQAQ
      - 2P stabilizing substitution: K986P, V987P

    Raises AssertionError if the WT reference does not match the expected
    residues at those positions (a safety check against using a stale or
    wrong reference sequence).
    """
    seq = list(wt_protein)
    if "".join(seq[681:685]) != "RRAR":
        raise AssertionError(
            f"Expected RRAR at positions 682-685, found {''.join(seq[681:685])} "
            "- refusing to apply furin site mutation to an unexpected sequence"
        )
    seq[681:685] = list("QQAQ")

    if seq[985] != 'K' or seq[986] != 'V':
        raise AssertionError(
            f"Expected K,V at positions 986-987, found {seq[985]},{seq[986]} "
            "- refusing to apply 2P mutation to an unexpected sequence"
        )
    seq[985] = 'P'  # K986P
    seq[986] = 'P'  # V987P

    return "".join(seq)


if __name__ == "__main__":
    mrna = read_fasta("../data/pfizer_moderna_real.fasta")
    wt = list(read_fasta("../data/spike_wt.fasta").values())[0]

    names = list(mrna.keys())
    pfizer_cds = find_orf(mrna[names[0]])
    moderna_cds = find_orf(mrna[names[1]])

    pfizer_protein = translate(pfizer_cds)
    moderna_protein = translate(moderna_cds)

    az_protein = wt  # AstraZeneca: unmodified WT spike
    novavax_protein = apply_novavax_mutations(wt)

    proteins = {
        "Pfizer_BNT162b2": pfizer_protein,
        "Moderna_mRNA1273": moderna_protein,
        "AstraZeneca_ChAdOx1": az_protein,
        "Novavax_NVXCoV2373": novavax_protein,
    }

    with open("../data/all_4_vaccine_spike_proteins.fasta", "w") as f:
        for name, seq in proteins.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

    for name, seq in proteins.items():
        print(f"{name}: {len(seq)} aa")

    print("\nPfizer  aa 986-987:", pfizer_protein[985:987], "(expect PP)")
    print("Moderna aa 986-987:", moderna_protein[985:987], "(expect PP)")
    print("Pfizer  aa 682-685:", pfizer_protein[681:685], "(expect RRAR)")
    print("Moderna aa 682-685:", moderna_protein[681:685], "(expect RRAR)")
