# Comparative Genomic Analysis of COVID-19 Vaccine Platforms
# Sourabhya Kasireddy
#
# Run this in R (RStudio, Posit Cloud, or your course environment) where
# Bioconductor/Biostrings is installed. Install once with:
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install("Biostrings")

library(Biostrings)

# ---- 1. Load the 4 real vaccine spike antigen sequences -------------------
# ../data/all_4_vaccine_spike_proteins.fasta contains:
#   - Pfizer_BNT162b2   : translated from real sequenced mRNA
#                         (Jeong, McCoy, Artiles et al., Stanford/Fire lab,
#                          NAalytics GitHub repo)
#   - Moderna_mRNA1273  : translated from real sequenced mRNA (same source)
#   - AstraZeneca_ChAdOx1 : unmodified WT spike, UniProt P0DTC2
#                           (per Watanabe et al. 2021, ACS Cent Sci)
#   - Novavax_NVXCoV2373  : WT spike + furin site knockout (RRAR->QQAQ,
#                           positions 682-685) + 2P mutation (K986P/V987P)
#                           (per Bangaru et al. 2020; Tian et al. 2021, Science)

proteins <- readAAStringSet("../data/all_4_vaccine_spike_proteins.fasta")
print(proteins)

# ---- 2. Pairwise global alignment (all 6 combinations) ---------------------
names_vec <- names(proteins)
pairs <- combn(names_vec, 2, simplify = FALSE)

results <- data.frame(seq1 = character(), seq2 = character(),
                       pct_identity = numeric(), score = numeric())

for (p in pairs) {
  aln <- pairwiseAlignment(proteins[[p[1]]], proteins[[p[2]]],
                            type = "global",
                            substitutionMatrix = "BLOSUM62",
                            gapOpening = 10, gapExtension = 0.5)
  pid <- pid(aln, type = "PID1")  # percent identity
  results <- rbind(results, data.frame(
    seq1 = p[1], seq2 = p[2],
    pct_identity = round(pid, 2),
    score = round(score(aln), 1)
  ))
}

print(results)
write.csv(results, "../results/r_pairwise_identity_results.csv", row.names = FALSE)

# ---- 3. GC content / codon usage (nucleotide level, Pfizer vs Moderna) -----
mrna <- readDNAStringSet("../data/pfizer_moderna_real.fasta")
freqs <- letterFrequency(mrna, letters = c("G", "C"))
gc_content <- 100 * rowSums(freqs) / width(mrna)
print(data.frame(vaccine = names(mrna), gc_pct = round(gc_content, 2)))
# Codon usage table for each mRNA (from first ATG)
get_cds <- function(seq) {
  char_seq <- as.character(seq)
  start <- regexpr("ATG", char_seq)[1]
  DNAString(substr(char_seq, start, nchar(char_seq)))
}
for (i in seq_along(mrna)) {
  cds <- as.character(get_cds(mrna[[i]]))
  n <- nchar(cds)
  n <- n - (n %% 3)  # trim to a whole number of codons
  codon_list <- substring(cds, seq(1, n, by = 3), seq(3, n, by = 3))
  tbl <- table(codon_list)
  cat("\nTop 5 codons -", names(mrna)[i], ":\n")
  print(head(sort(tbl, decreasing = TRUE), 5))
}

# ---- 4. Multiple sequence alignment across all 4 (for visualization) ------
# Requires the msa package: BiocManager::install("msa")
if (requireNamespace("msa", quietly = TRUE)) {
  library(msa)
  alignment <- msa(proteins, method = "ClustalOmega")
  print(alignment, show = "complete")
  msaPrettyPrint(alignment, output = "pdf", file = "../results/msa_4_vaccines.pdf",
                  askForOverwrite = FALSE, verbose = FALSE)
} else {
  message("Install the 'msa' package for the full 4-way visual alignment: ",
          "BiocManager::install('msa')")
}
