# 1. Configurar biblioteca personal
dir.create("~/R/library", recursive = TRUE, showWarnings = FALSE)
.libPaths(c("~/R/library", .libPaths()))

# 2. Instalar paquetes (solo si no están instalados)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib = "~/R/library")

BiocManager::install(c("Rsamtools", "GenomicFeatures", "GenomicAlignments", "org.Ce.eg.db"),
                     lib = "~/R/library")
pacman::p_load(Rsamtools,GenomicFeatures,GenomicAlignments,org.Ce.eg.db)
gtfFile = "celegans.gff3"
txdb = makeTxDbFromGFF(gtfFile, format="gff3")
genes = exonsBy(txdb, by="gene")
indir = getwd()
files = list.files(indir, pattern = '*sorted.bam')
bamLst = BamFileList(files, index=character(),obeyQname=TRUE)
PRJNA701955 = summarizeOverlaps(features = genes,
read=bamLst,
mode="Union",
singleEnd=FALSE,
ignore.strand=TRUE,
fragments=FALSE)
save(PRJNA701955,file="PRJNA701955.rda")
