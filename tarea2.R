# Este script permite construir un SummarizedExperiment a paertir de un rda con
# una matriz de conteos.

library(Biobase)
library(SummarizedExperiment)
pacman::p_load("org.Ce.eg.db")
library(dplyr)
setwd("C:/Users/adria/Desktop/Master/Bioinformática estadística")

load("PRJNA701955.rda")

###############################
# Creamos la variable fenotípica a partir de el % de oxigeno
###############################

csv_rowdata<-read.csv("SraRunTable.csv", sep=",")


csv_rowdata <- csv_rowdata %>%
  mutate(Group = ifelse(oxygen_exposure == "21%", "control", "hypoxia"))
write.csv(csv_rowdata, 'SraRunTable.csv', row.names = FALSE)

#----------------------------------------------------------------------
conteo<-assay(PRJNA701955)
colnames(conteo) <- gsub(".sorted.bam", "", colnames(conteo), fixed = TRUE)

row_data <- DataFrame(
  Symbol = rownames(conteo),
  ENSEMBL = mapIds(org.Ce.eg.db, 
                   keys = rownames(conteo),
                   keytype = "SYMBOL",
                   column = "ENSEMBL"),
  ENTREZID = mapIds(org.Ce.eg.db, 
                  keys = rownames(conteo),
                  keytype = "SYMBOL",
                  column = "ENTREZID"),
  GO = mapIds(org.Ce.eg.db, 
              keys = rownames(conteo),
              keytype = "SYMBOL",
              column = "GO"),
  WORMBASE = mapIds(org.Ce.eg.db, 
                    keys = rownames(conteo),
                    keytype = "SYMBOL",
                    column = "WORMBASE"),
  GENENAME = mapIds(org.Ce.eg.db, 
                    keys = rownames(conteo),
                    keytype = "SYMBOL",
                    column = "GENENAME"),
  stringsAsFactors = FALSE
)

rownames(conteo) <- row_data$ENTREZID

metadatos<-read.csv("SraRunTable.csv", sep=",")
col_data <- metadatos %>%
  dplyr::select(
    Run, 
    LibraryLayout,
    genotype,
    oxygen_exposure,
    Group,
    source_name,
    Instrument,
    LibrarySelection,
    LibrarySource
  ) %>%
  dplyr::rename(
    SampleID = Run,
    Condition = oxygen_exposure,
    Tissue = source_name
  ) %>%
  as.data.frame()

rownames(col_data) <- col_data$SampleID

PRJNA701955_SE <- SummarizedExperiment(
  assays   = list(counts = conteo),
  rowData  = row_data,
  colData  = col_data
)

library(edgeR)

colData(PRJNA701955_SE)$Group<-factor(colData(PRJNA701955_SE)$Group, levels = 
                                        c("control", "hypoxia"))
colData(PRJNA701955_SE)$mutant <- ifelse(colData(PRJNA701955_SE)$genotype == "N2", "WT", "MUT")
colData(PRJNA701955_SE)$mutant <- factor(colData(PRJNA701955_SE)$mutant, levels = c("WT", "MUT"))
colData(PRJNA701955_SE)

annot <- as.data.frame(rowData(PRJNA701955_SE))

table(colData(PRJNA701955_SE)[,"Group"])
elim<-which(is.na(colData(PRJNA701955_SE)$"Group")|is.na(colData(PRJNA701955_SE)
                                                         $"mutant"))
if (length(elim) > 0) {
  PRJNA701955_SE <- PRJNA701955_SE[, -elim]
}

count_matrix <- assay(PRJNA701955_SE, "counts")
dge<-DGEList(counts = count_matrix)
to_keep = rowSums(cpm(dge) > 0.5) >= 9
dge = dge[to_keep,keep.lib.sizes=FALSE]

#Clasico

# Condición

exact<-DGEList(counts=dge$counts, group=PRJNA701955_SE$Group)

dge.c<-estimateCommonDisp(exact)
et.c<-exactTest(dge.c)
tt_exact<-topTags(et.c, n=Inf, adjust.method = "BH", sort.by = "none")$table

# Mutantes

exact_mut<-DGEList(counts=dge$counts, group=PRJNA701955_SE$mutant)

dge.c_mut<-estimateCommonDisp(exact_mut)
et.c_mut<-exactTest(dge.c_mut)
tt_exact_mut<-topTags(et.c_mut)


#GLM

condicion<-colData(PRJNA701955_SE)$Group
mutante<-colData(PRJNA701955_SE)$mutant

design<-model.matrix(~condicion*mutante, data = colData(PRJNA701955_SE))

dge.glm<-estimateGLMCommonDisp(dge, design)

colnames(design)<-make.names(colnames(design))
contrasts<-makeContrasts(contrast1=condicionhypoxia,
                         contrast2=mutanteMUT,
                         contrast3=condicionhypoxia.mutanteMUT,
                         levels=design)

fit=glmFit(dge.glm,design=design)

lrt1<-glmLRT(fit, contrast = contrasts[,"contrast1"])
tt_glm_hyp<-topTags(lrt1, n=Inf, sort.by = "none", adjust.method = "BH")$table

tt_annotated_glm_hyp <- merge(tt_glm_hyp, annot, by = "row.names", 
                              all.x = TRUE)
# Mutantes

lrt2<-glmLRT(fit, contrast = contrasts[,"contrast2"])
tt_glm_mut<-topTags(lrt2, n=Inf, sort.by = "none", adjust.method = "BH")$table

tt_annotated_glm_mut <- merge(tt_glm_mut, annot, by = "row.names", 
                              all.x = TRUE)

# Interacción

lrt3<-glmLRT(fit, contrast = contrasts[,"contrast3"])
tt_glm_int<-topTags(lrt3, n=Inf, sort.by = "none", adjust.method = "BH")$table

tt_annotated_glm_int <- merge(tt_glm_int, annot, by = "row.names", 
                              all.x = TRUE)

#Colecciones de genes

library(celegans.db)
library(GO.db)
library(GSEABase)

#GO
frame <- toTable(org.Ce.egGO)
goframeData <- data.frame(
        go_id = frame$go_id,
        Evidence = frame$Evidence,
        gene_id = frame$gene_id,
        stringsAsFactors = FALSE
      )
goFrame <- GOFrame(goframeData, organism = "Caenorhabditis elegans")
goAllFrame <- GOAllFrame(goFrame)
gscCe1 <- GeneSetCollection(goAllFrame, setType = GOCollection())
gscCe<-geneIds(gscCe1)
names(gscCe) <- paste0(
      names(gscCe), "_",
      Term(GOTERM[names(gscCe)])
    )

save(gscCe,file = "gscCe.rda")

#KEGG
kegg_annotations <- as.list(org.Ce.egPATH2EG)
kegg_genesets <- lapply(names(kegg_annotations), function(pathway_id) {
  genes <- kegg_annotations[[pathway_id]]
  if (length(genes) > 0) {
    GeneSet(
      geneIds = genes,
      setName = pathway_id,  
      collectionType = KEGGCollection(),
      organism = "Caenorhabditis elegans"
    )
  }
})

gscKEGG0 <- GeneSetCollection(kegg_genesets)
gscKEGG <- geneIds(gscKEGG0)
save(gscKEGG, file = "gscKEGG_Celegans.rda")

# Análisis de sobrerepresentación

# GO
elim<-which(sapply(gscCe, length)<5 | sapply(gscCe, length)>100)
gscCe<-gscCe[-elim]

library(tami)

selected_genes1<-which(tt_annotated_glm_hyp$FDR<0.05)
entrez_selected<-as.character(tt_annotated_glm_hyp[selected_genes1,"ENTREZID"])

ora_hyp_glm_go<-ora(entrez_selected, gsc=gscCe)

# Mutantes

tt_upregulated_glm_mut<-subset(tt_annotated_glm_mut, logFC>0)

selected_genes2<-which(tt_upregulated_glm_mut$FDR<0.05)
entrez_selected2<-as.character(tt_upregulated_glm_mut[selected_genes2,"ENTREZID"])

ora_mut_glm_go<-ora(entrez_selected2, gsc=gscCe)
sum(ora_mut_glm_go$rawp<0.05)

max(ora_mut_glm_go$OR)

# Interaccion

selected_genes3<-which(tt_annotated_glm_int$FDR<0.05)
entrez_selected3<-as.character(tt_annotated_glm_int[selected_genes3,"ENTREZID"])
ora_int_glm_go<-ora(entrez_selected3, gsc=gscCe)


# KEGG

# Hypoxia
selected_genes1<-which(tt_annotated_glm_hyp$FDR<0.05)
entrez_selected<-as.character(tt_annotated_glm_hyp[selected_genes1,"ENTREZID"])
ora_hyp_glm_kegg<-ora(entrez_selected, gsc=gscKEGG)

sum(ora_hyp_glm_kegg$rawp<0.05)

# Mutantes

tt_upregulated_glm_mut<-subset(tt_annotated_glm_mut, logFC>0)

selected_genes2<-which(tt_upregulated_glm_mut$FDR<0.05)
entrez_selected2<-as.character(tt_upregulated_glm_mut[selected_genes2,"ENTREZID"])

ora_mut_glm_kegg<-ora(entrez_selected2, gsc=gscKEGG)
sum(ora_mut_glm_kegg$rawp<0.05)

max(ora_mut_glm_kegg$OR)
# Interaccion

selected_genes3<-which(tt_annotated_glm_int$FDR<0.05)
entrez_selected3<-as.character(tt_annotated_glm_int[selected_genes3,"ENTREZID"])
ora_int_glm_kegg<-ora(entrez_selected3, gsc=gscKEGG)



# Análisis de grupos de genes

# Hipótesis autocontenida (para ver si hay algo significativo)

library(tami)

# Hipoxia/normoxia
eset<-PRJNA701955_SE
PRJNA701955_self_group <- GeneSetTest(x =eset ,y="Group",gsc=gscCe,
                    test = edgercommon,association="pvalue",correction="BH",
                               GeneNullDistr = "randomization",
                               GeneSetNullDistr ="self-contained",
                               alternative="less",nmax = 1000,
                               id = "ENTREZID",descriptive=mean)
save(PRJNA701955_self_group, file = "PRJNA701955_self_group.rda")
SE_self_group_go<-PRJNA701955_self_group@GeneSetStat

self_group_go<-data.frame(PRJNA701955_self_group@GeneSetData,PRJNA701955_self_group@GeneSetStat)

# Análisis de mutantes

PRJNA701955_self_mut <- GeneSetTest(x =eset ,y="mutant",gsc=gscCe,
                                      test = edgercommon,association="pvalue",correction="BH",
                                      GeneNullDistr = "randomization",
                                      GeneSetNullDistr ="self-contained",
                                      alternative="less",nmax = 1000,
                                      id = "ENTREZID",descriptive=mean)
save(PRJNA701955_self_mut, file = "PRJNA701955_self_mut.rda")


# Hipótesis competitiva ()
# Hipoxia/normoxia
PRJNA701955_comp_group <- GeneSetTest(x =eset ,y="Group",gsc=gscCe,
                                      test = edgercommon,association="pvalue",
                                      correction="BH",
                                      GeneNullDistr = "randomization",
                                      GeneSetNullDistr ="competitive",
                                      alternative="less",nmax = 1000,
                                      id = "ENTREZID",descriptive=mean)

# Análisis de mutantes
PRJNA701955_comp_mutant <- GeneSetTest(x =eset ,y="mutant",gsc=gscCe,
                                      test = edgercommon,association="pvalue",
                                      correction="BH",
                                      GeneNullDistr = "randomization",
                                      GeneSetNullDistr ="competitive",
                                      alternative="less",nmax = 1000,
                                      id = "ENTREZID",descriptive=mean)


# KEGG

# Hipoxia/normoxia
eset<-PRJNA701955_SE
PRJNA701955_self_group_k <- GeneSetTest(x =eset ,y="Group",gsc=gscKEGG,
                                      test = edgercommon,association="pvalue",
                                      correction="BH",
                                      GeneNullDistr = "randomization",
                                      GeneSetNullDistr ="self-contained",
                                      alternative="less",nmax = 1000,
                                      id = "ENTREZID",descriptive=mean)
save(PRJNA701955_self_group_k, file = "PRJNA701955_self_group_k.rda")


# Análisis de mutantes

PRJNA701955_self_mut_k <- GeneSetTest(x =eset ,y="mutant",gsc=gscKEGG,
                                    test = edgercommon,association="pvalue",
                                    correction="BH",
                                    GeneNullDistr = "randomization",
                                    GeneSetNullDistr ="self-contained",
                                    alternative="less",nmax = 1000,
                                    id = "ENTREZID",descriptive=mean)
save(PRJNA701955_self_mut_k, file = "PRJNA701955_self_mut_k.rda")


# Hipótesis competitiva
# Hipoxia/normoxia
PRJNA701955_comp_group_k <- GeneSetTest(x =eset ,y="Group",gsc=gscKEGG,
                                      test = edgercommon,association="pvalue",
                                      correction="BH",
                                      GeneNullDistr = "randomization",
                                      GeneSetNullDistr ="competitive",
                                      alternative="less",nmax = 1000,
                                      id = "ENTREZID",descriptive=mean)

# Análisis de mutantes
PRJNA701955_comp_mutant_k<- GeneSetTest(x =eset ,y="mutant",gsc=gscKEGG,
                                       test = edgercommon,association="pvalue",
                                       correction="BH",
                                       GeneNullDistr = "randomization",
                                       GeneSetNullDistr ="competitive",
                                       alternative="less",nmax = 1000,
                                       id = "ENTREZID",descriptive=mean)

