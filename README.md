# Non3nExonFinder

Non3nExonFinder is an R/Shiny tool for identifying shared non-triplet coding exons across user-selected human or mouse NCBI RefSeq protein-coding transcripts.

The tool compares exon and CDS annotations across selected transcripts and identifies shared CDS-only exons with coding lengths not divisible by three.

## Features

- Human and mouse RefSeq NM_ transcripts
- Transcript-specific exon/CDS annotation
- Exact shared exon-boundary detection
- CDS/UTR classification
- Coding-length modulo-three analysis
- Internal-exon filtering
- Same-gene NM_ transcript retrieval
- Interactive Shiny interface
- CSV export

## Requirements

R 4.3.3

Required packages:

- shiny
- GenomicFeatures
- AnnotationDbi
- dplyr
- DT

Additional packages used for validation:

- GenomicRanges
- IRanges

## Annotation

Human:
- GRCh38.p14
- NCBI accession: GCF_000001405.40

Mouse:
- GRCm39
- NCBI accession: GCF_000001635.27

The program uses local RefSeq-derived TxDb databases:

- txdb_human.sqlite
- txdb_mouse.sqlite

The TxDb files can be generated using `build_TxDb.R`.

The original RefSeq GTF files are not included in this repository.

## Usage

Place the TxDb files in the same directory as `Non3nExonFinder_v1.0.R` and run:

source("Non3nExonFinder_v1.0.R")

Enter one or more RefSeq NM_ transcript accessions in the Shiny interface and run the analysis.

## Candidate definition

High-confidence candidates are shared exons that:

1. have identical genomic boundaries across all selected transcripts;
2. are CDS-only in all selected transcripts;
3. have coding lengths not divisible by three in all selected transcripts; and
4. are internal exons in all selected transcripts.

Non3nExonFinder is intended for upstream exon prioritization. Functional consequences of exon deletion should be evaluated separately.

## Validation

Independent validation was performed using the original NCBI RefSeq GTF annotations.

Validation scripts are provided in the `validation` directory:

- `validation/TIA1_independent_validation.R`
- `validation/large_scale_validation_and_benchmark.R`

The large-scale validation included 40 human and 40 mouse genes, comprising 456 transcripts and 6,848 transcript-specific exon instances.

Random seeds:

- Human: 260921
- Mouse: 260922

No discrepancies were detected between Non3nExonFinder and the independently reconstructed annotations or candidate sets.

Validation output files are provided in the `results` and `results_LARGE` directories.
Random seeds:

- Human: 260921
- Mouse: 260922

No discrepancies were detected between Non3nExonFinder and the independently reconstructed annotations or candidate sets.

Validation output files are provided in the `results` directory.

## Repository contents

```text
Non3nExonFinder/
├── Non3nExonFinder_v1.0.R
├── build_TxDb.R
├── validation/
├── results/
└── results_LARGE/

## Citation

Lee S-Y. Non3nExonFinder: transcript-aware identification of shared non-triplet coding exons.

Citation information will be updated after publication.
