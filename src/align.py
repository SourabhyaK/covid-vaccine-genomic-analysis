"""
Global pairwise alignment (Needleman-Wunsch) of real, published SARS-CoV-2
spike-encoding mRNA sequences for the Pfizer-BioNTech (BNT162b2) and
Moderna (mRNA-1273) vaccines.

Source of sequences: Jeong DE, McCoy M, Artiles K, et al. "Assemblies of
putative SARS-CoV2-spike-encoding mRNA sequences for vaccines BNT-162b2
and mRNA-1273." Stanford University (Fire/Shoura labs), posted publicly
on GitHub (NAalytics repo), derived from direct RNA-seq of vaccine vial
remnants under FDA research-use authorization. This is real experimental
data, not a simulation.

No R/Bioconductor or Biopython is installed in this execution sandbox, so
this implements standard Needleman-Wunsch global alignment directly in
Python/NumPy rather than calling Biostrings::pairwiseAlignment. The
algorithm and scoring (match=+1, mismatch=-1, gap=-2, linear gap penalty)
mirror what Biostrings would compute for global alignment; the R
equivalent script is provided separately for use in an R environment
with the actual package installed.
"""
import numpy as np

def read_fasta(path):
    seqs = {}
    name = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:]
                seqs[name] = []
            else:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    # score matrix
    score = np.zeros((n + 1, m + 1), dtype=np.int32)
    score[:, 0] = np.arange(n + 1) * gap
    score[0, :] = np.arange(m + 1) * gap

    for i in range(1, n + 1):
        row_prev = score[i - 1]
        row_cur = score[i]
        c1 = seq1[i - 1]
        for j in range(1, m + 1):
            diag = row_prev[j - 1] + (match if c1 == seq2[j - 1] else mismatch)
            up = row_prev[j] + gap
            left = row_cur[j - 1] + gap
            row_cur[j] = max(diag, up, left)

    # traceback
    aligned1, aligned2 = [], []
    i, j = n, m
    while i > 0 and j > 0:
        c1, c2 = seq1[i - 1], seq2[j - 1]
        diag_score = score[i - 1][j - 1] + (match if c1 == c2 else mismatch)
        if score[i][j] == diag_score:
            aligned1.append(c1)
            aligned2.append(c2)
            i -= 1
            j -= 1
        elif score[i][j] == score[i - 1][j] + gap:
            aligned1.append(c1)
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(c2)
            j -= 1
    while i > 0:
        aligned1.append(seq1[i - 1]); aligned2.append('-'); i -= 1
    while j > 0:
        aligned1.append('-'); aligned2.append(seq2[j - 1]); j -= 1

    aligned1.reverse()
    aligned2.reverse()
    return "".join(aligned1), "".join(aligned2), int(score[n][m])

if __name__ == "__main__":
    seqs = read_fasta("../data/pfizer_moderna_real.fasta")
    names = list(seqs.keys())
    s1, s2 = seqs[names[0]], seqs[names[1]]
    print(f"Aligning:\n  1: {names[0]} ({len(s1)} nt)\n  2: {names[1]} ({len(s2)} nt)\n")

    a1, a2, total_score = needleman_wunsch(s1, s2)

    matches = sum(1 for x, y in zip(a1, a2) if x == y and x != '-')
    mismatches = sum(1 for x, y in zip(a1, a2) if x != y and x != '-' and y != '-')
    gaps = sum(1 for x, y in zip(a1, a2) if x == '-' or y == '-')
    aligned_len = len(a1)
    pct_identity = 100 * matches / aligned_len

    print(f"Alignment length: {aligned_len}")
    print(f"Matches: {matches}")
    print(f"Mismatches: {mismatches}")
    print(f"Gap positions: {gaps}")
    print(f"Percent identity: {pct_identity:.2f}%")
    print(f"Raw alignment score: {total_score}")

    with open("../results/alignment_result.txt", "w") as f:
        f.write(f"Pairwise global alignment: {names[0]} vs {names[1]}\n")
        f.write(f"Alignment length: {aligned_len}\n")
        f.write(f"Matches: {matches}\nMismatches: {mismatches}\nGaps: {gaps}\n")
        f.write(f"Percent identity: {pct_identity:.2f}%\n")
        f.write(f"Score: {total_score}\n\n")
        # write alignment in blocks of 80
        for k in range(0, aligned_len, 80):
            f.write(a1[k:k+80] + "\n")
            f.write(a2[k:k+80] + "\n\n")
