library(GenomicFeatures)
library(AnnotationDbi)
library(txdbmaker)

txdb_human <- txdbmaker::makeTxDbFromGFF(
  "GCF_000001405.40_GRCh38.p14_genomic.gtf.gz",
  format = "gtf"
)

txdb_mouse <- txdbmaker::makeTxDbFromGFF(
  "GCF_000001635.27_GRCm39_genomic.gtf.gz",
  format = "gtf"
)

AnnotationDbi::saveDb(txdb_human, "txdb_human.sqlite")
AnnotationDbi::saveDb(txdb_mouse, "txdb_mouse.sqlite")