library(GenomicFeatures)
library(AnnotationDbi)
library(GenomicRanges)
library(dplyr)

source("Non3nExonFinder_v1.0.R")

tx_ids_tia1 <- c(
  "NM_001351508.2", "NM_001351509.2", "NM_001351510.2", "NM_001351511.1",
  "NM_001351512.1", "NM_001351513.1", "NM_001351514.2", "NM_001351515.2"
)

gtf_file <- "GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
dir.create("results", showWarnings = FALSE)

# Run Non3nExonFinder
tia1_res <- analyze_common_non3n_exons(tx_ids_tia1, txdb_human)
stopifnot(length(tia1_res$dropped_tx_ids) == 0)

# Extract the selected transcripts directly from the RefSeq GTF
escaped_ids <- gsub("\\.", "\\\\.", tx_ids_tia1)
pattern <- paste0('transcript_id "(', paste(escaped_ids, collapse = "|"), ')"')

con <- gzfile(gtf_file, open = "rt")
gtf_lines <- character()

repeat {
  x <- readLines(con, n = 100000)
  if (!length(x)) break
  gtf_lines <- c(gtf_lines, x[grepl(pattern, x)])
}
close(con)

writeLines(gtf_lines, "results/TIA1_raw_GTF_subset.gtf")

tia1_raw <- read.delim(
  "results/TIA1_raw_GTF_subset.gtf",
  header = FALSE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE
)

colnames(tia1_raw) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
tia1_raw$tx_id <- sub('.*transcript_id "([^"]+)".*', "\\1", tia1_raw$attribute)

# Reconstruct exon coordinates and exon rank from the raw GTF
raw_exons <- tia1_raw %>%
  filter(feature == "exon") %>%
  transmute(
    tx_id, chr, start, end, strand,
    exon_rank = as.integer(sub('.*exon_number "([0-9]+)".*', "\\1", attribute))
  ) %>%
  distinct() %>%
  arrange(tx_id, exon_rank)

# Reconstruct coding intervals from CDS and stop_codon features
raw_coding <- tia1_raw %>%
  filter(feature %in% c("CDS", "stop_codon")) %>%
  select(tx_id, chr, start, end, strand) %>%
  distinct()

coding_gr <- GRanges(
  seqnames = raw_coding$chr,
  ranges = IRanges(start = raw_coding$start, end = raw_coding$end),
  strand = raw_coding$strand
)
coding_gr$tx_id <- raw_coding$tx_id

# Reduce coding intervals within each transcript before calculating overlap
coding_list <- split(coding_gr, coding_gr$tx_id)
coding_list <- lapply(coding_list, reduce)
coding_reduced <- unlist(GRangesList(coding_list), use.names = FALSE)
coding_reduced$tx_id <- rep(names(coding_list), lengths(coding_list))

exon_gr <- GRanges(
  seqnames = raw_exons$chr,
  ranges = IRanges(start = raw_exons$start, end = raw_exons$end),
  strand = raw_exons$strand
)
exon_gr$tx_id <- raw_exons$tx_id

hits <- findOverlaps(exon_gr, coding_reduced, ignore.strand = FALSE)
same_tx <- exon_gr$tx_id[queryHits(hits)] == coding_reduced$tx_id[subjectHits(hits)]
hits <- hits[same_tx]

raw_exons$cds_bp <- 0L
if (length(hits)) {
  ov <- pintersect(exon_gr[queryHits(hits)], coding_reduced[subjectHits(hits)])
  bp_by_exon <- tapply(width(ov), queryHits(hits), sum)
  raw_exons$cds_bp[as.integer(names(bp_by_exon))] <- as.integer(bp_by_exon)
}

# Derive transcript-specific exon annotations independently
raw_exons <- raw_exons %>%
  group_by(tx_id) %>%
  mutate(
    max_exon_rank = max(exon_rank),
    is_internal_exon = exon_rank > 1L & exon_rank < max_exon_rank
  ) %>%
  ungroup() %>%
  mutate(
    exon_len = end - start + 1L,
    region_type = case_when(
      cds_bp == 0L ~ "UTR_only",
      cds_bp == exon_len ~ "CDS_only",
      cds_bp > 0L & cds_bp < exon_len ~ "mixed",
      TRUE ~ "unknown"
    ),
    cds_len_mod3 = if_else(cds_bp > 0L, cds_bp %% 3L, NA_integer_),
    is_non3n_cds = cds_bp > 0L & cds_bp %% 3L != 0L,
    exon_key = paste(chr, start, end, strand, sep = ":")
  )

# Compare raw-GTF annotations with Non3nExonFinder output
app_exons <- tia1_res$exon_table %>%
  select(tx_id, exon_rank, chr, start, end, strand, cds_bp, region_type, cds_len_mod3, is_non3n_cds) %>%
  rename_with(~ paste0("app_", .x), -c(tx_id, exon_rank))

raw_compare <- raw_exons %>%
  select(tx_id, exon_rank, chr, start, end, strand, cds_bp, region_type, cds_len_mod3, is_non3n_cds) %>%
  rename_with(~ paste0("raw_", .x), -c(tx_id, exon_rank))

comparison <- full_join(raw_compare, app_exons, by = c("tx_id", "exon_rank"))

same_value <- function(x, y) {
  (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
}

comparison <- comparison %>%
  mutate(
    coordinates_match = same_value(raw_chr, app_chr) &
      same_value(raw_start, app_start) &
      same_value(raw_end, app_end) &
      same_value(raw_strand, app_strand),
    coding_bp_match = same_value(raw_cds_bp, app_cds_bp),
    region_type_match = same_value(raw_region_type, app_region_type),
    mod3_match = same_value(raw_cds_len_mod3, app_cds_len_mod3),
    non3n_match = same_value(raw_is_non3n_cds, app_is_non3n_cds)
  )

exon_validation <- data.frame(
  metric = c(
    "Exon coordinates / rank", "Coding contribution", "CDS/UTR classification",
    "Modulo-three classification", "Non-3n classification"
  ),
  matched = c(
    sum(comparison$coordinates_match), sum(comparison$coding_bp_match),
    sum(comparison$region_type_match), sum(comparison$mod3_match),
    sum(comparison$non3n_match)
  ),
  total = nrow(comparison)
)

# Reconstruct candidate sets from the raw GTF
n_selected <- length(tx_ids_tia1)

raw_common <- raw_exons %>%
  group_by(exon_key, chr, start, end, strand) %>%
  summarise(
    n_transcripts = n_distinct(tx_id),
    all_cds_only = all(region_type == "CDS_only"),
    all_non3n = all(is_non3n_cds),
    all_internal = all(is_internal_exon),
    .groups = "drop"
  ) %>%
  filter(n_transcripts == n_selected)

raw_sets <- list(
  "All common exons" = raw_common$exon_key,
  "Common CDS-only exons" = raw_common$exon_key[raw_common$all_cds_only],
  "Common non-3n exons" = raw_common$exon_key[raw_common$all_cds_only & raw_common$all_non3n],
  "High-confidence candidates" = raw_common$exon_key[
    raw_common$all_cds_only & raw_common$all_non3n & raw_common$all_internal
  ]
)

app_sets <- list(
  "All common exons" = tia1_res$common_all$exon_key,
  "Common CDS-only exons" = tia1_res$common_cds_only$exon_key,
  "Common non-3n exons" = tia1_res$common_non3n$exon_key,
  "High-confidence candidates" = tia1_res$high_confidence$exon_key
)

candidate_validation <- data.frame(
  stage = names(raw_sets),
  raw_gtf_n = lengths(raw_sets),
  app_n = lengths(app_sets),
  exact_match = mapply(setequal, raw_sets, app_sets),
  row.names = NULL
)

# Save validation outputs
write.csv(data.frame(tx_id = tx_ids_tia1), "results/TIA1_transcripts.csv", row.names = FALSE)
write.csv(comparison, "results/TIA1_exon_validation.csv", row.names = FALSE)
write.csv(exon_validation, "results/TIA1_exon_validation_summary.csv", row.names = FALSE)
write.csv(candidate_validation, "results/TIA1_candidate_validation.csv", row.names = FALSE)

print(exon_validation)
print(candidate_validation)

stopifnot(
  all(exon_validation$matched == exon_validation$total),
  all(candidate_validation$exact_match)
)

