# Comparative Genomic & Clinical Analysis of COVID-19 Vaccine Platforms

A sequence-level and clinical comparison of four COVID-19 vaccines spanning
three distinct platforms: mRNA-lipid nanoparticle (Pfizer, Moderna), viral
vector (AstraZeneca), and recombinant protein subunit (Novavax).

## Background

All four vaccines target the same antigen — the SARS-CoV-2 spike protein —
but reach it through completely different mechanisms, and each manufacturer
made distinct engineering choices to the antigen itself. This project asks:
how similar are these vaccines actually at the sequence level, and do the
differences that exist line up with real, documented design decisions?

## Data sources (all real, no synthetic or placeholder sequences)

| Vaccine | Sequence source | Provenance |
|---|---|---|
| Pfizer BNT162b2 | Direct RNA-seq of vaccine vial remnants | Jeong, McCoy, Artiles et al. (Stanford, Fire/Shoura labs), published via NAalytics GitHub repo |
| Moderna mRNA-1273 | Direct RNA-seq of vaccine vial remnants | Same study as above |
| AstraZeneca ChAdOx1-S | Reference spike (UniProt P0DTC2), unmodified | Confirmed unmodified WT design per Watanabe et al. 2021, *ACS Central Science* |
| Novavax NVX-CoV2373 | Reference spike (UniProt P0DTC2) + 2 documented mutations | Furin site knockout and 2P substitution per Bangaru et al. 2020 and Tian et al. 2021, *Science* |

AstraZeneca and Novavax sequences are reconstructed from the verified WT
reference plus the exact, multiply-cited published mutations — not guessed.
The mutation positions were confirmed against the actual reference sequence
before being applied (see `build_proteins.py`).

## Methods

1. **Translation & sequence assembly** (`build_proteins.py`) — translated
   the real Pfizer/Moderna mRNA to protein using the standard genetic code;
   constructed AstraZeneca (WT) and Novavax (WT + furin knockout + 2P) from
   the verified reference.
2. **Pairwise alignment** (`align.py`, `analysis.R`) — Needleman-Wunsch
   global alignment (Python prototype) and Biostrings `pairwiseAlignment`
   with BLOSUM62 (R, for the actual submission) across all 6 vaccine pairs.
3. **Nucleotide-level analysis** — GC content and codon usage comparison
   between the Pfizer and Moderna mRNA sequences, which encode a near-
   identical protein but were independently codon-optimized.
4. **Clinical context** — published Phase 3 trial efficacy for each vaccine,
   pulled from the original pivotal trial papers.

## Key findings

- **Protein level:** 99.6–100% identity across all 4 vaccines. Pfizer and
  Moderna are **100% identical** at the protein level.
- **Nucleotide level:** Pfizer and Moderna mRNA are only **89.9% identical**
  despite encoding the identical protein — the gap is entirely due to each
  company's independent codon-optimization choices, confirmed by direct GC
  content (56.9% vs 62.1%) and codon-usage comparison. Both companies favor
  **CTG** as their single most-used codon for leucine (105 vs 103 uses), but
  diverge sharply after that — Moderna's second-most-common codon, **AGC**
  (serine, 96 uses), doesn't appear in Pfizer's top 5 at all. This is direct
  evidence of two independently-run codon-optimization algorithms converging
  on the dominant choice but diverging on secondary ones, rather than any
  biological difference in the antigen itself.
- **Design differences map to real, confirmed mutations:** translating the
  actual sequenced mRNA independently reproduced the published 2P mutation
  (positions 986–987 = PP) and intact furin site (682–685 = RRAR) for both
  mRNA vaccines, without those facts being assumed going in — a built-in
  validation that the pipeline is correct.
- AstraZeneca is the only platform using the **unmodified** wild-type spike
  (no 2P), while Novavax is the only one with **both** a mutated furin site
  and the 2P substitution.

## Clinical efficacy (for context — see limitations)

| Vaccine | Efficacy | Source |
|---|---|---|
| Pfizer BNT162b2 | 95.0% | Polack et al. 2020, *NEJM* |
| Moderna mRNA-1273 | 94.1% | Baden et al. 2021, *NEJM* |
| AstraZeneca ChAdOx1-S | 74.0% | Falsey et al. 2021, *NEJM* (US trial) |
| Novavax NVX-CoV2373 | 89.7% | Heath et al. 2021, *NEJM* (UK trial) |

## Limitations (read before citing the efficacy numbers above)

- **These efficacy trials are not directly comparable.** They were run in
  different countries, different time periods, and against different
  dominant circulating variants — the AstraZeneca and Novavax trial
  populations in particular were not matched to the Pfizer/Moderna U.S.
  trial population. A higher or lower number does not mean one vaccine is
  simply "better" than another in a head-to-head sense.
  - AstraZeneca: a separate UK/Brazil trial reported ~70.4% with a
    different dosing regimen; the 74.0% here is specifically the later
    US Phase 3 trial.
- **Sequence identity is not the same as immunogenicity.** Two spike
  sequences being 99.6% identical doesn't mean they produce equally
  strong immune responses — presentation context (mRNA in a cell vs.
  purified protein with an adjuvant) matters as much as the sequence.
- **AstraZeneca and Novavax sequences here are reconstructed**, not
  directly sequenced the way Pfizer and Moderna were. They are built from
  a verified reference plus specific, well-documented published mutations,
  but this is a different evidentiary standard than direct sequencing of
  the physical vaccine product.
- **No variant/strain analysis.** This project compares the original
  Wuhan-based vaccine designs only; it doesn't address how well any of
  these spike sequences match later variants (Delta, Omicron, etc.).

## How to reproduce

```bash
# Python prototype (sequence construction + alignment, no R required)
cd scripts
python3 build_proteins.py
python3 align.py

# Full analysis (requires R + Bioconductor)
Rscript analysis.R
```

Requires `Biostrings` (and optionally `msa` for the multi-sequence
alignment figure): see comments at the top of `analysis.R`.

**Verified reproducibility:** every number in this README was independently
cross-checked two ways — a from-scratch Needleman-Wunsch implementation in
Python, and `Biostrings::pairwiseAlignment` with BLOSUM62 in R — on two
separate machines. Both produced identical pairwise identity values
(100.00% / 99.84% / 99.76% / 99.61%) and identical GC content (56.91% /
62.10%), confirming the pipeline is correct rather than an artifact of one
implementation's quirks.

Environment used for the R analysis: R 4.3.3, Bioconductor `Biostrings`
(installed via `BiocManager`). Run `sessionInfo()` after loading the
packages to capture exact versions for your own environment.

## Why this project

Beyond the sequence comparison itself, this project is a first step toward
a longer-term interest in applying computational methods to antigen and
therapeutic design — understanding how small, deliberate sequence changes
(like the 2P stabilizing mutation or furin site knockout found here) affect
a molecule's behavior is the same underlying question that shows up in
ML-driven drug discovery and protein engineering more broadly.

## Repository structure

```
covid-vaccine-genomic-analysis/
├── src/                 # source code
│   ├── build_proteins.py  # sequence construction + mutation verification
│   ├── align.py            # Python Needleman-Wunsch prototype
│   └── analysis.R          # full Biostrings-based analysis (the actual submission)
├── tests/               # pytest unit tests
│   ├── test_translation.py # translation + FASTA parsing logic
│   ├── test_mutations.py   # Novavax mutation construction + safety checks
│   └── test_alignment.py   # Needleman-Wunsch correctness on known cases
├── data/                # real sequence data
│   ├── pfizer_moderna_real.fasta
│   ├── spike_wt.fasta
│   └── all_4_vaccine_spike_proteins.fasta
├── figures/             # output visualizations
├── results/             # generated analysis outputs (CSV, alignment text)
├── requirements.txt     # Python dependencies
├── environment.yml      # full conda environment (Python + R/Bioconductor)
└── README.md
```

## Running the tests

```bash
pip install -r requirements.txt
cd tests
pytest -v
```

22 unit tests cover: codon translation correctness (including stop-codon
handling and incomplete trailing codons), the Novavax mutation construction
(verifying the furin-site and 2P substitutions land on the exact expected
positions, and that the function safely refuses to run if the input sequence
doesn't match what's expected), and the Needleman-Wunsch alignment
implementation (identity scoring, gap handling, and relative scoring
sanity checks).

## Reproducing the full environment

```bash
conda env create -f environment.yml
conda activate covid-vaccine-genomic-analysis
```

This installs both the Python stack and R/Bioconductor (Biostrings, msa) in
one step, rather than installing each manually.
