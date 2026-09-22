library(GenomicFeatures)
library(AnnotationDbi)
library(IRanges)
library(dplyr)

source("Non3nExonFinder_v1.0.R")

HUMAN_TXDB <- "txdb_human.sqlite"
MOUSE_TXDB <- "txdb_mouse.sqlite"
HUMAN_GTF <- "GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
MOUSE_GTF <- "GCF_000001635.27_GRCm39_genomic.gtf.gz"

N_GENES <- 40L
HUMAN_SEED <- 260921L
MOUSE_SEED <- 260922L
BENCHMARK_SIZES <- c(1L, 5L, 10L, 20L, 40L)
BENCHMARK_REPS <- 10L
VALID_NM_PATTERN <- "^NM_[0-9]+(\\.[0-9]+)?$"
OUTPUT_DIR <- "results_LARGE"

dir.create(OUTPUT_DIR, showWarnings = FALSE)

txdb_human <- AnnotationDbi::loadDb(HUMAN_TXDB)
txdb_mouse <- AnnotationDbi::loadDb(MOUSE_TXDB)

# Build the eligible NM_ transcript catalog
build_nm_catalog <- function(txdb) {
  genes <- AnnotationDbi::keys(txdb, keytype = "GENEID")
  suppressMessages(
    AnnotationDbi::select(txdb, keys = genes, keytype = "GENEID", columns = c("GENEID", "TXNAME"))
  ) %>%
    filter(!is.na(GENEID), !is.na(TXNAME), grepl(VALID_NM_PATTERN, TXNAME)) %>%
    distinct(GENEID, TXNAME) %>%
    arrange(GENEID, TXNAME)
}

human_nm <- build_nm_catalog(txdb_human)
mouse_nm <- build_nm_catalog(txdb_mouse)

human_multi <- human_nm %>% count(GENEID, name = "n_tx") %>% filter(n_tx >= 2L) %>% arrange(GENEID)
mouse_multi <- mouse_nm %>% count(GENEID, name = "n_tx") %>% filter(n_tx >= 2L) %>% arrange(GENEID)

# Exclude TIA1 from the human sampling pool
TIA1_TX <- c(
  "NM_001351508.2", "NM_001351509.2", "NM_001351510.2", "NM_001351511.1",
  "NM_001351512.1", "NM_001351513.1", "NM_001351514.2", "NM_001351515.2"
)

tia1_gene_ids <- human_nm %>% filter(TXNAME %in% TIA1_TX) %>% distinct(GENEID) %>% pull(GENEID)
if (!length(tia1_gene_ids)) stop("TIA1 could not be mapped in the human TxDb.")

set.seed(HUMAN_SEED)
human_genes <- human_multi %>%
  filter(!GENEID %in% tia1_gene_ids) %>%
  slice_sample(n = N_GENES) %>%
  arrange(GENEID)

set.seed(MOUSE_SEED)
mouse_genes <- mouse_multi %>%
  slice_sample(n = N_GENES) %>%
  arrange(GENEID)

human_tx <- human_nm %>% filter(GENEID %in% human_genes$GENEID) %>% arrange(GENEID, TXNAME)
mouse_tx <- mouse_nm %>% filter(GENEID %in% mouse_genes$GENEID) %>% arrange(GENEID, TXNAME)

validation_cohort <- bind_rows(
  human_genes %>% mutate(species = "Human"),
  mouse_genes %>% mutate(species = "Mouse")
) %>%
  transmute(species, gene = GENEID, n_transcripts = n_tx) %>%
  arrange(species, gene)

validation_transcripts <- bind_rows(
  human_tx %>% mutate(species = "Human"),
  mouse_tx %>% mutate(species = "Mouse")
) %>%
  transmute(species, gene = GENEID, transcript = TXNAME) %>%
  arrange(species, gene, transcript)

write.csv(validation_cohort, file.path(OUTPUT_DIR, "validation_cohort.csv"), row.names = FALSE)
write.csv(validation_transcripts, file.path(OUTPUT_DIR, "validation_transcripts.csv"), row.names = FALSE)

# Run Non3nExonFinder for each validation gene
run_app_batch <- function(selected_genes, tx_map, txdb, species) {
  exon_out <- list()
  candidate_out <- list()
  summary_out <- list()

  for (i in seq_len(nrow(selected_genes))) {
    gene_i <- selected_genes$GENEID[i]
    tx_i <- tx_map %>% filter(GENEID == gene_i) %>% pull(TXNAME) %>% unique() %>% sort()

    res <- analyze_common_non3n_exons(tx_i, txdb)
    if (length(res$dropped_tx_ids)) stop(species, " ", gene_i, ": one or more transcripts were dropped.")

    exon_out[[i]] <- res$exon_table %>% mutate(species = species, gene = gene_i)

    candidate_out[[i]] <- bind_rows(
      res$common_all %>% transmute(exon_key, stage = "All common exons"),
      res$common_cds_only %>% transmute(exon_key, stage = "Common CDS-only"),
      res$common_non3n %>% transmute(exon_key, stage = "Common non-3n"),
      res$high_confidence %>% transmute(exon_key, stage = "High-confidence")
    ) %>%
      mutate(species = species, gene = gene_i)

    summary_out[[i]] <- data.frame(
      species = species, gene = gene_i, n_transcripts = length(tx_i),
      n_exon_instances = nrow(res$exon_table),
      n_common = nrow(res$common_all),
      n_cds_only = nrow(res$common_cds_only),
      n_non3n = nrow(res$common_non3n),
      n_high_confidence = nrow(res$high_confidence)
    )
  }

  list(
    exons = bind_rows(exon_out),
    candidates = bind_rows(candidate_out),
    summary = bind_rows(summary_out)
  )
}

human_app <- run_app_batch(human_genes, human_tx, txdb_human, "Human")
mouse_app <- run_app_batch(mouse_genes, mouse_tx, txdb_mouse, "Mouse")

app_exons <- bind_rows(human_app$exons, mouse_app$exons)
app_candidates <- bind_rows(human_app$candidates, mouse_app$candidates)
app_summary <- bind_rows(human_app$summary, mouse_app$summary)

# Extract exon, CDS, and stop_codon records directly from the original GTF
extract_gtf_records <- function(gtf_file, tx_ids, chunk_size = 200000L) {
  con <- gzfile(gtf_file, open = "rt")
  on.exit(close(con), add = TRUE)
  out <- list()
  k <- 0L

  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (!length(lines)) break

    lines <- lines[nzchar(lines) & substr(lines, 1L, 1L) != "#"]
    if (!length(lines)) next

    fields <- strsplit(lines, "\t", fixed = TRUE)
    fields <- fields[lengths(fields) >= 9L]
    if (!length(fields)) next

    feature <- vapply(fields, `[[`, character(1), 3L)
    keep <- feature %in% c("exon", "CDS", "stop_codon")
    fields <- fields[keep]
    if (!length(fields)) next

    attributes <- vapply(fields, `[[`, character(1), 9L)
    tx <- sub('.*transcript_id "([^"]+)".*', "\\1", attributes)
    keep <- tx %in% tx_ids
    if (!any(keep)) next

    fields <- fields[keep]
    attributes <- attributes[keep]
    tx <- tx[keep]
    k <- k + 1L

    out[[k]] <- data.frame(
      chr = vapply(fields, `[[`, character(1), 1L),
      feature = vapply(fields, `[[`, character(1), 3L),
      start = as.integer(vapply(fields, `[[`, character(1), 4L)),
      end = as.integer(vapply(fields, `[[`, character(1), 5L)),
      strand = vapply(fields, `[[`, character(1), 7L),
      attributes = attributes,
      tx_id = tx,
      stringsAsFactors = FALSE
    )
  }

  bind_rows(out)
}

human_raw <- extract_gtf_records(HUMAN_GTF, unique(human_tx$TXNAME))
mouse_raw <- extract_gtf_records(MOUSE_GTF, unique(mouse_tx$TXNAME))

if (length(setdiff(unique(human_tx$TXNAME), unique(human_raw$tx_id)))) {
  stop("One or more human validation transcripts were not recovered from the raw GTF.")
}
if (length(setdiff(unique(mouse_tx$TXNAME), unique(mouse_raw$tx_id)))) {
  stop("One or more mouse validation transcripts were not recovered from the raw GTF.")
}

# Reconstruct transcript-level exon annotations independently
get_exon_rank <- function(x) {
  out <- rep(NA_integer_, length(x))
  has_rank <- grepl('exon_number "', x, fixed = TRUE)
  out[has_rank] <- as.integer(sub('.*exon_number "([0-9]+)".*', "\\1", x[has_rank]))
  out
}

build_raw_exons <- function(raw_gtf) {
  exons <- raw_gtf %>%
    filter(feature == "exon") %>%
    mutate(exon_rank = get_exon_rank(attributes)) %>%
    select(tx_id, chr, start, end, strand, exon_rank) %>%
    distinct()

  if (anyNA(exons$exon_rank)) stop("A selected exon lacks a numeric exon_number in the raw GTF.")

  coding <- raw_gtf %>%
    filter(feature %in% c("CDS", "stop_codon")) %>%
    select(tx_id, chr, start, end, strand) %>%
    distinct()

  cds_bp <- integer(nrow(exons))

  for (i in seq_len(nrow(exons))) {
    e <- exons[i, ]
    x <- coding %>%
      filter(tx_id == e$tx_id, chr == e$chr, strand == e$strand, start <= e$end, end >= e$start)

    if (nrow(x)) {
      reduced <- IRanges::reduce(IRanges::IRanges(start = x$start, end = x$end))
      overlap_start <- pmax(e$start, IRanges::start(reduced))
      overlap_end <- pmin(e$end, IRanges::end(reduced))
      cds_bp[i] <- sum(pmax(0L, overlap_end - overlap_start + 1L))
    }
  }

  exons %>%
    group_by(tx_id) %>%
    mutate(max_exon_rank = max(exon_rank)) %>%
    ungroup() %>%
    mutate(
      exon_len = end - start + 1L,
      cds_bp = cds_bp,
      region_type = case_when(
        cds_bp == 0L ~ "UTR_only",
        cds_bp == exon_len ~ "CDS_only",
        cds_bp > 0L & cds_bp < exon_len ~ "mixed",
        TRUE ~ "unknown"
      ),
      cds_len_mod3 = if_else(cds_bp > 0L, cds_bp %% 3L, NA_integer_),
      is_non3n_cds = cds_bp > 0L & cds_bp %% 3L != 0L,
      is_internal_exon = exon_rank > 1L & exon_rank < max_exon_rank,
      exon_key = paste(chr, start, end, strand, sep = ":")
    )
}

human_raw_exons <- build_raw_exons(human_raw) %>%
  left_join(human_tx %>% transmute(tx_id = TXNAME, gene = GENEID), by = "tx_id") %>%
  mutate(species = "Human")

mouse_raw_exons <- build_raw_exons(mouse_raw) %>%
  left_join(mouse_tx %>% transmute(tx_id = TXNAME, gene = GENEID), by = "tx_id") %>%
  mutate(species = "Mouse")

raw_exons <- bind_rows(human_raw_exons, mouse_raw_exons)
if (anyNA(raw_exons$gene)) stop("A raw-GTF exon could not be assigned to its selected gene.")

# Reconstruct the four candidate sets from the raw GTF
build_raw_candidates <- function(raw_exons, selected_genes, species) {
  out <- list()

  for (i in seq_len(nrow(selected_genes))) {
    gene_i <- selected_genes$GENEID[i]
    expected_n <- selected_genes$n_tx[i]
    x <- raw_exons %>% filter(gene == gene_i)

    if (n_distinct(x$tx_id) != expected_n) {
      stop(species, " ", gene_i, ": transcript count differs between the selected cohort and raw GTF.")
    }

    common <- x %>%
      group_by(exon_key, chr, start, end, strand) %>%
      summarise(
        n_tx = n_distinct(tx_id),
        all_cds_only = all(region_type == "CDS_only"),
        all_non3n = all(is_non3n_cds),
        all_internal = all(is_internal_exon),
        .groups = "drop"
      ) %>%
      filter(n_tx == expected_n)

    out[[i]] <- bind_rows(
      common %>% transmute(exon_key, stage = "All common exons"),
      common %>% filter(all_cds_only) %>% transmute(exon_key, stage = "Common CDS-only"),
      common %>% filter(all_cds_only, all_non3n) %>% transmute(exon_key, stage = "Common non-3n"),
      common %>% filter(all_cds_only, all_non3n, all_internal) %>% transmute(exon_key, stage = "High-confidence")
    ) %>%
      mutate(species = species, gene = gene_i)
  }

  bind_rows(out)
}

raw_candidates <- bind_rows(
  build_raw_candidates(human_raw_exons, human_genes, "Human"),
  build_raw_candidates(mouse_raw_exons, mouse_genes, "Mouse")
)

# Compare transcript-level annotations
equal_na <- function(x, y) {
  z <- (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
  z[is.na(z)] <- FALSE
  z
}

raw_compare <- raw_exons %>%
  transmute(
    species, gene, tx_id, exon_rank, raw_present = TRUE,
    raw_chr = chr, raw_start = start, raw_end = end, raw_strand = strand,
    raw_cds_bp = cds_bp, raw_region_type = region_type, raw_mod3 = cds_len_mod3,
    raw_non3n = is_non3n_cds, raw_internal = is_internal_exon
  )

app_compare <- app_exons %>%
  transmute(
    species, gene, tx_id, exon_rank, app_present = TRUE,
    app_chr = chr, app_start = start, app_end = end, app_strand = strand,
    app_cds_bp = cds_bp, app_region_type = region_type, app_mod3 = cds_len_mod3,
    app_non3n = is_non3n_cds, app_internal = is_internal_exon
  )

exon_validation <- full_join(raw_compare, app_compare, by = c("species", "gene", "tx_id", "exon_rank")) %>%
  mutate(
    row_match = !is.na(raw_present) & !is.na(app_present),
    coordinates_match = row_match & equal_na(raw_chr, app_chr) & equal_na(raw_start, app_start) &
      equal_na(raw_end, app_end) & equal_na(raw_strand, app_strand),
    coding_bp_match = row_match & equal_na(raw_cds_bp, app_cds_bp),
    region_type_match = row_match & equal_na(raw_region_type, app_region_type),
    mod3_match = row_match & equal_na(raw_mod3, app_mod3),
    non3n_match = row_match & equal_na(raw_non3n, app_non3n),
    internal_match = row_match & equal_na(raw_internal, app_internal)
  ) %>%
  arrange(species, gene, tx_id, exon_rank)

# Compare candidate sets gene by gene
compare_candidate_sets <- function(selected_genes, raw_candidates, app_candidates, species) {
  stages <- c("All common exons", "Common CDS-only", "Common non-3n", "High-confidence")
  out <- list()
  k <- 0L

  for (gene_i in selected_genes$GENEID) {
    for (stage_i in stages) {
      raw_set <- raw_candidates %>%
        filter(species == !!species, gene == gene_i, stage == stage_i) %>%
        pull(exon_key) %>% unique()

      app_set <- app_candidates %>%
        filter(species == !!species, gene == gene_i, stage == stage_i) %>%
        pull(exon_key) %>% unique()

      k <- k + 1L
      out[[k]] <- data.frame(
        species = species, gene = gene_i, stage = stage_i,
        n_raw = length(raw_set), n_app = length(app_set),
        exact_match = setequal(raw_set, app_set),
        missing_from_app = paste(setdiff(raw_set, app_set), collapse = ";"),
        extra_in_app = paste(setdiff(app_set, raw_set), collapse = ";")
      )
    }
  }

  bind_rows(out)
}

candidate_validation <- bind_rows(
  compare_candidate_sets(human_genes, raw_candidates, app_candidates, "Human"),
  compare_candidate_sets(mouse_genes, raw_candidates, app_candidates, "Mouse")
)

# Validation summaries
metric_summary <- bind_rows(
  exon_validation %>% group_by(species) %>% summarise(metric = "Exon coordinates / rank", matched = sum(coordinates_match), total = n(), .groups = "drop"),
  exon_validation %>% group_by(species) %>% summarise(metric = "Coding contribution", matched = sum(coding_bp_match), total = n(), .groups = "drop"),
  exon_validation %>% group_by(species) %>% summarise(metric = "CDS/UTR classification", matched = sum(region_type_match), total = n(), .groups = "drop"),
  exon_validation %>% group_by(species) %>% summarise(metric = "Modulo-three classification", matched = sum(mod3_match), total = n(), .groups = "drop"),
  exon_validation %>% group_by(species) %>% summarise(metric = "Non-3n classification", matched = sum(non3n_match), total = n(), .groups = "drop"),
  exon_validation %>% group_by(species) %>% summarise(metric = "Internal-exon status", matched = sum(internal_match), total = n(), .groups = "drop")
) %>%
  mutate(concordance_pct = 100 * matched / total)

candidate_summary <- candidate_validation %>%
  group_by(species, stage) %>%
  summarise(matched = sum(exact_match), total = n(), .groups = "drop") %>%
  mutate(concordance_pct = 100 * matched / total)

exon_mismatches <- exon_validation %>%
  filter(!coordinates_match | !coding_bp_match | !region_type_match | !mod3_match | !non3n_match | !internal_match)

candidate_mismatches <- candidate_validation %>% filter(!exact_match)

gene_summary <- app_summary %>%
  left_join(exon_mismatches %>% count(species, gene, name = "n_exon_mismatches"), by = c("species", "gene")) %>%
  left_join(candidate_mismatches %>% count(species, gene, name = "n_candidate_mismatches"), by = c("species", "gene")) %>%
  mutate(
    n_exon_mismatches = coalesce(n_exon_mismatches, 0L),
    n_candidate_mismatches = coalesce(n_candidate_mismatches, 0L),
    validation_passed = n_exon_mismatches == 0L & n_candidate_mismatches == 0L
  ) %>%
  arrange(species, gene)

# Computational performance benchmark
benchmark_gene_row <- human_multi %>%
  filter(n_tx >= max(BENCHMARK_SIZES)) %>%
  mutate(distance = abs(n_tx - max(BENCHMARK_SIZES))) %>%
  arrange(distance, GENEID) %>%
  slice(1L)

if (!nrow(benchmark_gene_row)) stop("No human gene contains enough NM_ transcripts for the benchmark.")

benchmark_gene <- benchmark_gene_row$GENEID[1L]
benchmark_tx <- human_nm %>% filter(GENEID == benchmark_gene) %>% arrange(TXNAME) %>% pull(TXNAME)
benchmark_runs <- list()
k <- 0L

for (n_i in BENCHMARK_SIZES) {
  tx_i <- benchmark_tx[seq_len(n_i)]
  invisible(analyze_common_non3n_exons(tx_i, txdb_human))

  for (rep_i in seq_len(BENCHMARK_REPS)) {
    gc(verbose = FALSE)
    elapsed <- system.time(invisible(analyze_common_non3n_exons(tx_i, txdb_human)))[["elapsed"]]
    k <- k + 1L
    benchmark_runs[[k]] <- data.frame(
      gene = benchmark_gene, n_transcripts = n_i, replicate = rep_i, elapsed_seconds = elapsed
    )
  }
}

benchmark_runs <- bind_rows(benchmark_runs)
benchmark_summary <- benchmark_runs %>%
  group_by(gene, n_transcripts) %>%
  summarise(
    median_seconds = median(elapsed_seconds),
    q1_seconds = unname(quantile(elapsed_seconds, 0.25)),
    q3_seconds = unname(quantile(elapsed_seconds, 0.75)),
    .groups = "drop"
  )

# Save outputs
write.csv(exon_validation, file.path(OUTPUT_DIR, "exon_level_validation.csv"), row.names = FALSE)
write.csv(candidate_validation, file.path(OUTPUT_DIR, "candidate_set_validation.csv"), row.names = FALSE)
write.csv(gene_summary, file.path(OUTPUT_DIR, "gene_level_validation_summary.csv"), row.names = FALSE)
write.csv(metric_summary, file.path(OUTPUT_DIR, "exon_validation_summary.csv"), row.names = FALSE)
write.csv(candidate_summary, file.path(OUTPUT_DIR, "candidate_validation_summary.csv"), row.names = FALSE)
write.csv(benchmark_runs, file.path(OUTPUT_DIR, "performance_benchmark.csv"), row.names = FALSE)
write.csv(benchmark_summary, file.path(OUTPUT_DIR, "performance_benchmark_summary.csv"), row.names = FALSE)
write.csv(exon_mismatches, file.path(OUTPUT_DIR, "QC_exon_mismatches.csv"), row.names = FALSE)
write.csv(candidate_mismatches, file.path(OUTPUT_DIR, "QC_candidate_mismatches.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))

print(metric_summary)
print(candidate_summary)
print(benchmark_summary)

stopifnot(
  nrow(exon_mismatches) == 0L,
  nrow(candidate_mismatches) == 0L,
  all(gene_summary$validation_passed)
)
