# Non3nExonFinder v1.0
# Required files: txdb_human.sqlite, txdb_mouse.sqlite
# Required packages: GenomicFeatures, AnnotationDbi, dplyr, shiny, DT

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(dplyr)
  library(shiny)
  library(DT)
})

VALID_NM_PATTERN <- "^NM_[0-9]+(\\.[0-9]+)?$"

txdb_human <- AnnotationDbi::loadDb("txdb_human.sqlite")
txdb_mouse <- AnnotationDbi::loadDb("txdb_mouse.sqlite")

annotation_label <- list(
  human = "Human RefSeq annotation (GRCh38.p14; local TxDb)",
  mouse = "Mouse RefSeq annotation (GRCm39; local TxDb)"
)

# Input and annotation utilities

parse_tx_ids <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(character(0))

  x <- gsub("[\u2028\u2029]", "\n", x)
  ids <- unlist(strsplit(x, "[,[:space:]]+"))
  ids <- gsub("[\"']", "", ids)
  unique(trimws(ids[nzchar(ids)]))
}

require_columns <- function(txdb, cols) {
  missing <- setdiff(cols, columns(txdb))
  if (length(missing)) {
    stop(
      "This TxDb does not contain required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

validate_nm_ids <- function(tx_ids) {
  bad <- tx_ids[!grepl(VALID_NM_PATTERN, tx_ids)]
  if (length(bad)) {
    stop(
      "Only NCBI RefSeq protein-coding transcript accessions beginning with NM_ are accepted. ",
      "Invalid input: ", paste(bad, collapse = ", ")
    )
  }
  invisible(TRUE)
}

get_same_gene_nm <- function(tx_id, txdb) {
  validate_nm_ids(tx_id)

  if (!all(c("TXNAME", "GENEID") %in% columns(txdb))) {
    stop("GENEID mapping is not available in this TxDb.")
  }

  tx_gene <- suppressMessages(
    AnnotationDbi::select(
      txdb, keys = tx_id, keytype = "TXNAME",
      columns = c("TXNAME", "GENEID")
    )
  )

  gene_ids <- unique(na.omit(tx_gene$GENEID))
  if (!length(gene_ids)) stop("No GeneID was found for ", tx_id, ".")
  if (length(gene_ids) > 1L) {
    stop(
      "The transcript maps to more than one GeneID in this annotation: ",
      paste(gene_ids, collapse = ", ")
    )
  }

  all_tx <- suppressMessages(
    AnnotationDbi::select(
      txdb, keys = gene_ids, keytype = "GENEID",
      columns = c("GENEID", "TXNAME")
    )
  )

  sort(unique(all_tx$TXNAME[grepl(VALID_NM_PATTERN, all_tx$TXNAME)]))
}

# Transcript-level exon and CDS annotation

query_selected_exons <- function(tx_ids, txdb) {
  needed <- c(
    "TXNAME", "EXONCHROM", "EXONSTART", "EXONEND",
    "EXONSTRAND", "EXONRANK", "EXONID"
  )
  require_columns(txdb, needed)

  suppressMessages(
    AnnotationDbi::select(
      txdb, keys = tx_ids, keytype = "TXNAME", columns = needed
    )
  ) %>%
    transmute(
      tx_id = TXNAME,
      exon_id = EXONID,
      chr = as.character(EXONCHROM),
      start = as.integer(EXONSTART),
      end = as.integer(EXONEND),
      strand = as.character(EXONSTRAND),
      exon_rank = as.integer(EXONRANK)
    ) %>%
    filter(
      !is.na(tx_id), !is.na(chr), !is.na(start), !is.na(end),
      !is.na(strand), !is.na(exon_rank)
    ) %>%
    distinct()
}

query_selected_cds <- function(tx_ids, txdb) {
  needed <- c(
    "TXNAME", "CDSCHROM", "CDSSTART", "CDSEND",
    "CDSSTRAND", "CDSID", "CDSPHASE"
  )
  require_columns(txdb, needed)

  suppressMessages(
    AnnotationDbi::select(
      txdb, keys = tx_ids, keytype = "TXNAME", columns = needed
    )
  ) %>%
    transmute(
      tx_id = TXNAME,
      cds_id = CDSID,
      chr = as.character(CDSCHROM),
      start = as.integer(CDSSTART),
      end = as.integer(CDSEND),
      strand = as.character(CDSSTRAND),
      cds_phase = suppressWarnings(as.integer(CDSPHASE))
    ) %>%
    filter(
      !is.na(tx_id), !is.na(chr), !is.na(start), !is.na(end),
      !is.na(strand)
    ) %>%
    distinct()
}

calculate_exon_cds_table <- function(tx_ids, txdb) {
  tx_ids <- unique(trimws(tx_ids))
  validate_nm_ids(tx_ids)

  tx_check <- suppressMessages(
    AnnotationDbi::select(
      txdb, keys = tx_ids, keytype = "TXNAME", columns = "TXNAME"
    )
  )

  valid <- sort(unique(na.omit(tx_check$TXNAME)))
  dropped <- setdiff(tx_ids, valid)

  if (!length(valid)) {
    stop("None of the supplied NM_ transcript accessions were found in the selected annotation.")
  }

  exons <- query_selected_exons(valid, txdb)
  cds <- query_selected_cds(valid, txdb)

  tx_max_rank <- exons %>%
    group_by(tx_id) %>%
    summarise(max_exon_rank = max(exon_rank, na.rm = TRUE), .groups = "drop")

  exon_cds_bp <- integer(nrow(exons))

  for (i in seq_len(nrow(exons))) {
    e <- exons[i, ]

    csub <- cds %>%
      filter(
        tx_id == e$tx_id,
        chr == e$chr,
        strand == e$strand,
        start <= e$end,
        end >= e$start
      )

    if (!nrow(csub)) {
      exon_cds_bp[i] <- 0L
      next
    }

    ov_start <- pmax(e$start, csub$start)
    ov_end <- pmin(e$end, csub$end)

    intervals <- data.frame(
      start = ov_start,
      end = ov_end,
      width = pmax(0L, ov_end - ov_start + 1L)
    ) %>%
      filter(width > 0L) %>%
      distinct(start, end, .keep_all = TRUE)

    exon_cds_bp[i] <- sum(intervals$width)
  }

  exons <- exons %>%
    mutate(
      exon_len = end - start + 1L,
      cds_bp = exon_cds_bp,
      cds_len_mod3 = if_else(cds_bp > 0L, cds_bp %% 3L, NA_integer_),
      region_type = case_when(
        cds_bp == 0L ~ "UTR_only",
        cds_bp == exon_len ~ "CDS_only",
        cds_bp > 0L & cds_bp < exon_len ~ "mixed",
        TRUE ~ "unknown"
      ),
      is_non3n_cds = cds_bp > 0L & cds_bp %% 3L != 0L
    ) %>%
    left_join(tx_max_rank, by = "tx_id") %>%
    mutate(
      is_internal_exon = exon_rank > 1L & exon_rank < max_exon_rank,
      exon_key = paste(chr, start, end, strand, sep = ":")
    )

  list(
    valid_tx_ids = valid,
    dropped_tx_ids = dropped,
    exon_table = exons
  )
}

# Common exon filtering

summarize_common_exons <- function(exon_df, valid_tx_ids) {
  n_all <- length(valid_tx_ids)

  common <- exon_df %>%
    group_by(exon_key, chr, start, end, strand) %>%
    filter(n_distinct(tx_id) == n_all) %>%
    summarise(
      n_tx = n_distinct(tx_id),
      exon_len = first(exon_len),
      exon_ranks = paste(paste0(tx_id, ":", exon_rank), collapse = "; "),
      cds_bp_by_tx = paste(paste0(tx_id, ":", cds_bp), collapse = "; "),
      region_type_by_tx = paste(paste0(tx_id, ":", region_type), collapse = "; "),
      all_have_cds = all(cds_bp > 0L),
      all_cds_only = all(region_type == "CDS_only"),
      all_internal = all(is_internal_exon),
      all_non3n = all(cds_bp > 0L & cds_bp %% 3L != 0L),
      any_non3n = any(cds_bp > 0L & cds_bp %% 3L != 0L),
      same_cds_bp = n_distinct(cds_bp) == 1L,
      cds_len_mod3 = if (
        n_distinct(cds_bp) == 1L && !is.na(first(cds_bp)) && first(cds_bp) > 0L
      ) first(cds_bp) %% 3L else NA_integer_,
      cds_bp = if (
        n_distinct(cds_bp) == 1L && !is.na(first(cds_bp))
      ) first(cds_bp) else NA_integer_,
      .groups = "drop"
    ) %>%
    arrange(chr, start, end)

  common_cds_only <- common %>% filter(all_cds_only)
  common_non3n <- common %>% filter(all_cds_only, all_non3n)
  high_confidence <- common_non3n %>% filter(all_internal)

  list(
    common_all = common,
    common_cds_only = common_cds_only,
    common_non3n = common_non3n,
    high_confidence = high_confidence
  )
}

analyze_common_non3n_exons <- function(tx_ids_raw, txdb) {
  base <- calculate_exon_cds_table(tx_ids_raw, txdb)
  common <- summarize_common_exons(base$exon_table, base$valid_tx_ids)
  c(base, common)
}

# Shiny interface

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { background-color: #F7F9FB; }
      .custom-title {
        font-size: 32px; font-weight: 600; text-align: center;
        color: #6A1B9A; margin-top: 22px; margin-bottom: 6px;
      }
      .subtitle {
        text-align: center; font-size: 16px; color: #6C757D;
        margin-bottom: 24px;
      }
      .info-box {
        background: #FFFFFF; border: 1px solid #E2E6EA;
        border-left: 4px solid #6A1B9A; border-radius: 4px;
        padding: 12px 14px; margin-top: 14px; margin-bottom: 14px;
      }
      .guide-section { max-width: 950px; padding: 12px 22px 28px 22px; }
      .guide-section h3 { color: #4A4A4A; margin-top: 12px; }
      .guide-section h4 { margin-top: 24px; }
      .annotation-text { color: #6C757D; font-size: 12px; }
      .btn-primary, .btn-default { margin-top: 4px; }
    "))
  ),

  tags$div(
    class = "custom-title",
    icon("dna", lib = "font-awesome"),
    " Non3nExonFinder"
  ),
  tags$div(
    class = "subtitle",
    "Transcript-aware identification of shared non-triplet coding exons"
  ),

  tabsetPanel(
    id = "main_tabs",

    tabPanel(
      "Analysis",
      br(),
      sidebarLayout(
        sidebarPanel(
          radioButtons(
            "species", "Species",
            choices = c("Human" = "human", "Mouse" = "mouse"),
            selected = "human", inline = TRUE
          ),

          textAreaInput(
            "tx_ids",
            "NCBI RefSeq NM_ transcript IDs",
            placeholder = paste(
              "Enter one or more NM_ accessions",
              "Example:",
              "NM_001351508.2, NM_001351509.2",
              sep = "\n"
            ),
            rows = 9
          ),

          fluidRow(
            column(6, actionButton("same_gene", "Get same-gene NM_", width = "100%")),
            column(6, actionButton("run", "Run analysis", width = "100%"))
          ),

          tags$div(
            class = "info-box",
            tags$b("Candidate definition"),
            tags$p(
              style = "margin-top: 8px; margin-bottom: 6px;",
              "A common non-3n candidate is an exon with identical genomic boundaries and strand ",
              "across all selected transcripts, classified as CDS-only in every transcript, and with ",
              "a coding length not divisible by three in every transcript."
            ),
            tags$p(
              style = "margin-bottom: 0;",
              "High-confidence candidates additionally require the exon to be internal in every selected transcript."
            )
          ),

          tags$b("Annotation in use"),
          tags$div(class = "annotation-text", textOutput("annotation_info"))
        ),

        mainPanel(
          tags$h4("Analysis summary"),
          verbatimTextOutput("summary"),
          tags$hr(),

          tabsetPanel(
            tabPanel(
              "High-confidence candidates", br(),
              downloadButton("download_high", "Download CSV"),
              br(), br(), DTOutput("table_high")
            ),
            tabPanel(
              "Common non-3n CDS-only", br(),
              downloadButton("download_non3n", "Download CSV"),
              br(), br(), DTOutput("table_non3n")
            ),
            tabPanel(
              "Common CDS-only exons", br(),
              downloadButton("download_cds_only", "Download CSV"),
              br(), br(), DTOutput("table_cds_only")
            ),
            tabPanel(
              "All common exons", br(),
              downloadButton("download_common", "Download CSV"),
              br(), br(), DTOutput("table_common")
            ),
            tabPanel(
              "Per-transcript exon details", br(),
              downloadButton("download_details", "Download CSV"),
              br(), br(), DTOutput("table_details")
            )
          )
        )
      )
    ),

    tabPanel(
      "User Guide",
      tags$div(
        class = "guide-section",
        tags$h3("How to use Non3nExonFinder"),
        tags$ol(
          tags$li(
            tags$b("Select a species: "),
            "Choose Human or Mouse. The current implementation uses fixed local RefSeq-derived ",
            "TxDb annotations based on GRCh38.p14 for human and GRCm39 for mouse."
          ),
          tags$li(
            tags$b("Enter RefSeq transcripts: "),
            "Paste one or more NCBI RefSeq protein-coding transcript accessions beginning with NM_. ",
            "Accessions may be separated by spaces, commas, or line breaks."
          ),
          tags$li(
            tags$b("Optional same-gene retrieval: "),
            "Enter a transcript and click 'Get same-gene NM_' to retrieve NM_ transcripts mapped to ",
            "the same Gene ID in the selected local TxDb annotation."
          ),
          tags$li(
            tags$b("Run the analysis: "),
            "Click 'Run analysis'. Exon/CDS overlap is calculated separately for each transcript, ",
            "followed by CDS/UTR classification and coding-length divisibility assessment."
          ),
          tags$li(
            tags$b("Review candidate tables: "),
            "Results are organized from the high-confidence candidate set to all common exons. ",
            "Transcript-specific exon annotations are available in the final tab."
          ),
          tags$li(
            tags$b("Export results: "),
            "Each result table can be downloaded as a CSV file."
          )
        ),

        tags$h4("How candidates are defined"),
        tags$ul(
          tags$li(
            tags$b("Common exon: "),
            "chromosome, genomic start, genomic end, and strand are identical across all selected transcripts."
          ),
          tags$li(
            tags$b("CDS-only: "),
            "the entire exon is covered by coding sequence in every selected transcript."
          ),
          tags$li(
            tags$b("Non-triplet: "),
            "the transcript-specific coding contribution is not divisible by three."
          ),
          tags$li(
            tags$b("High-confidence candidate: "),
            "a common CDS-only non-triplet exon that is internal in every selected transcript."
          )
        ),

        tags$h4("Important notes"),
        tags$ul(
          tags$li(
            "Non3nExonFinder is an upstream annotation utility and does not design guide RNAs or ",
            "guarantee that deletion of a candidate exon will produce a functional knockout."
          ),
          tags$li(
            "Downstream consequences, including premature termination and nonsense-mediated decay, ",
            "should be evaluated separately."
          ),
          tags$li(
            "Results depend on the selected transcript set and annotation version."
          ),
          tags$li(
            "The current implementation is restricted to human and mouse RefSeq protein-coding ",
            "transcripts represented by NM_ accessions."
          )
        )
      )
    ),

    tabPanel(
      "About",
      tags$div(
        class = "guide-section",
        tags$h3("About Non3nExonFinder"),
        tags$p(
          "Non3nExonFinder is an interactive R/Shiny-based tool for prioritizing shared non-triplet ",
          "coding exons across user-defined human or mouse NCBI RefSeq protein-coding transcripts."
        ),
        tags$p(
          "The program evaluates exon and coding-sequence annotations independently for each transcript, ",
          "identifies exons with exact shared genomic boundaries, and sequentially filters candidates by ",
          "CDS-only status, coding-length divisibility, and internal-exon status."
        ),
        tags$p(
          "Independent validation showed complete concordance with reconstructed RefSeq GTF annotations ",
          "across 80 human and mouse genes and 6,848 transcript-specific exon instances, with exact ",
          "recovery of the corresponding candidate sets."
        ),
        tags$h4("Scope"),
        tags$p(
          "The tool is intended to support transcript-aware exon prioritization before downstream ",
          "CRISPR guide design and experimental validation."
        )
      )
    ),

    tabPanel(
      "Contact",
      tags$div(
        class = "guide-section",
        tags$h3("Contact Information"),
        tags$p("For questions, bug reports, or issues related to Non3nExonFinder, please contact:"),
        tags$p(
          tags$strong("Sung-Yeon Lee"),
          tags$br(),
          tags$a(
            href = "mailto:sungyeonlee0711@gmail.com",
            "sungyeonlee0711@gmail.com"
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  txdb_selected <- reactive({
    if (identical(input$species, "mouse")) txdb_mouse else txdb_human
  })

  output$annotation_info <- renderText({
    annotation_label[[input$species]]
  })

  parsed_ids <- reactive(parse_tx_ids(input$tx_ids))

  observeEvent(input$same_gene, {
    ids <- parsed_ids()
    validate(need(length(ids) >= 1L, "Enter at least one NM_ transcript ID first."))

    nm <- tryCatch(
      get_same_gene_nm(ids[1], txdb_selected()),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 8)
        NULL
      }
    )

    if (!is.null(nm)) {
      updateTextAreaInput(session, "tx_ids", value = paste(nm, collapse = "\n"))
      showNotification(
        paste(length(nm), "NM_ transcripts loaded from the same GeneID."),
        type = "message"
      )
    }
  })

  res <- eventReactive(input$run, {
    ids <- parsed_ids()
    validate(need(length(ids) > 0L, "Enter at least one NM_ transcript ID."))

    tryCatch(
      analyze_common_non3n_exons(ids, txdb_selected()),
      error = function(e) validate(need(FALSE, conditionMessage(e)))
    )
  })

  output$summary <- renderPrint({
    r <- res()

    list(
      species = input$species,
      annotation = annotation_label[[input$species]],
      n_input_transcripts = length(parsed_ids()),
      n_valid_transcripts = length(r$valid_tx_ids),
      valid_tx_ids = r$valid_tx_ids,
      dropped_tx_ids = r$dropped_tx_ids,
      n_common_exons = nrow(r$common_all),
      n_common_cds_only = nrow(r$common_cds_only),
      n_common_non3n = nrow(r$common_non3n),
      n_high_confidence_internal_non3n = nrow(r$high_confidence)
    )
  })

  output$table_high <- renderDT({
    datatable(
      res()$high_confidence,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$table_non3n <- renderDT({
    datatable(
      res()$common_non3n,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$table_cds_only <- renderDT({
    datatable(
      res()$common_cds_only,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$table_common <- renderDT({
    datatable(
      res()$common_all,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$table_details <- renderDT({
    details <- res()$exon_table %>%
      select(
        tx_id, exon_rank, chr, start, end, strand,
        exon_len, cds_bp, cds_len_mod3, region_type,
        is_non3n_cds, is_internal_exon
      ) %>%
      arrange(tx_id, exon_rank)

    datatable(
      details,
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$download_high <- downloadHandler(
    filename = function() paste0("Non3nExonFinder_high_confidence_", Sys.Date(), ".csv"),
    content = function(file) write.csv(res()$high_confidence, file, row.names = FALSE)
  )

  output$download_non3n <- downloadHandler(
    filename = function() paste0("Non3nExonFinder_common_non3n_", Sys.Date(), ".csv"),
    content = function(file) write.csv(res()$common_non3n, file, row.names = FALSE)
  )

  output$download_cds_only <- downloadHandler(
    filename = function() paste0("Non3nExonFinder_common_CDS_only_", Sys.Date(), ".csv"),
    content = function(file) write.csv(res()$common_cds_only, file, row.names = FALSE)
  )

  output$download_common <- downloadHandler(
    filename = function() paste0("Non3nExonFinder_all_common_exons_", Sys.Date(), ".csv"),
    content = function(file) write.csv(res()$common_all, file, row.names = FALSE)
  )

  output$download_details <- downloadHandler(
    filename = function() paste0("Non3nExonFinder_exon_details_", Sys.Date(), ".csv"),
    content = function(file) write.csv(res()$exon_table, file, row.names = FALSE)
  )
}

shinyApp(ui, server)
