### Finalized Analysis of MULTI-ATAC Pilot Redo Experiment ###

### Set Up Environment ------------------------------------------------------

library(ArchR)
library(deMULTIplex2)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggpubr)
'%ni%' <- Negate('%in%')

source("/Volumes/DannySSD/MULTI_ATAC/CustomFunctions.R")

addArchRThreads(threads = 1) 
addArchRGenome("hg38")

wd <- "/Volumes/DannySSD/MULTI_ATAC/20211102-PilotRedo/ArchR4"
dir.create(wd)
setwd(wd)

dir.create("Plots")


### Create Arrow Files ------------------------------------------------------

inputFiles <- c(
  Ill_labeled = "../cellranger-atac/gDNA3/fragments.tsv.gz"
)

ArrowFiles <- vector()

for (i in 1:length(inputFiles)) {
  arrow <- createArrowFiles(
    inputFiles = inputFiles[i],
    sampleNames = names(inputFiles)[i],
    minTSS = 2,
    minFrags = 100,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE)
  ArrowFiles[i] <- arrow
  rm(arrow)
}

ArrowFiles <- c("Ill_labeled.arrow")



### Create ArchR Project ----------------------------------------------------

proj_1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "proj_1",
  copyArrows = TRUE
)

thresh.TSS <- 5
thresh.nFrag <- 800

ggplot(as.data.frame(proj_1@cellColData[proj_1$nFrags > 100 & proj_1$TSSEnrichment > 3,]), aes(x = log10(nFrags), y = TSSEnrichment)) + 
  geom_density_2d_filled(bins = 10) + 
  geom_vline(xintercept = log10(thresh.nFrag)) + 
  geom_hline(yintercept = thresh.TSS) + 
  facet_grid(~Sample) + 
  theme_minimal() + 
  theme(legend.position = "none")
ggplot(as.data.frame(proj_1@cellColData[proj_1$nFrags > 100 & proj_1$TSSEnrichment > 3,]), aes(x = log10(nFrags), y = TSSEnrichment)) + 
  geom_point(size = 0.1, alpha = 0.1) + 
  geom_vline(xintercept = log10(thresh.nFrag)) + 
  geom_hline(yintercept = thresh.TSS) + 
  facet_grid(~Sample) + 
  theme_minimal() + 
  theme(legend.position = "none")
ggplot(as.data.frame(proj_1@cellColData), aes(x = log10(nFrags), y = TSSEnrichment, color = TSSEnrichment > (-13*log10(nFrags) + 45))) + 
  geom_point(size = 0.5, alpha = 0.1) + 
  facet_grid(~Sample) + 
  theme_minimal() + 
  theme(legend.position = "none")


proj_1 <- proj_1[proj_1$nFrags > thresh.nFrag & proj_1$TSSEnrichment > thresh.TSS,]

proj_1 <- addIterativeLSI(proj_1)
proj_1 <- addClusters(proj_1, dimsToUse = 1:30)
proj_1 <- addUMAP(proj_1, dimsToUse = 1:30)

g <- lapply(c("Sample", "Clusters", "log10(nFrags)", "TSSEnrichment"), function(x) {
  plotEmbedding(proj_1, name = x, plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle(x) 
}) %>% plot_grid(plotlist = ., ncol = 2)

pdfPlot(g, file = "Plots/proj1_Metadata.pdf", width = 12, height = 8)


## Add Alternate cell barcode naming schemes

proj_1$CB_CellRanger <- gsub(".*#", "", proj_1$cellNames)
proj_1$CB <- gsub("-.*", "", proj_1$CB_CellRanger)
proj_1$CB_ReverseComp <- revComplement(proj_1$CB)



### deMULTIplex2 Sample Classification --------------------------------------

dir.create(paste(wd, "SampleClassification",sep="/"))
data("multiseq_oligos")

bar.ref.rev <- multiseq_oligos[1:3]
names(bar.ref.rev) <- paste("Bar", 1:3, sep="")

cellIDs <- proj_1$CB_ReverseComp

## Create read & barcode tables

readTable_3 <- readTags(dir = "/Volumes/DannySSD/MULTI_ATAC/20211102-PilotRedo/MULTI/trimmed",
                        barcode = "MULTI-ATAC", 
                        assay = "ATAC",
                        name = "MULTI3")

saveRDS(readTable_3, "SampleClassification/readTable_3.rds")

# readTable_3 <- readRDS("SampleClassification/readTable_3.rds")

barTable_3 <- alignTags(readTable_3, 
                        bar.ref.rev,
                        filter.cells = cellIDs)

saveRDS(barTable_3, "SampleClassification/barTable_3.rds")

rm(readTable_3)
gc()

# barTable_3 <- readRDS("SampleClassification/barTable_3.rds")

tagHist(barTable_3) + geom_vline(xintercept = 100)
tagHist(barTable_3, minUMI = 100)

barTable_3 <- barTable_3[Matrix::rowSums(barTable_3) >= 100,]

## Run classification

res_3 <- demultiplexTags(barTable_3,
                         plot.path = "SampleClassification/",
                         plot.name = "Ill_Rev",
                         plot.diagnostics = T)

saveRDS(res_3, "SampleClassification/res_3.rds")




### Incorporate Sample Classifications --------------------------------------

barDF <- lapply(1, function(x) { 
  res <- list(res_3,res_5)[[x]]
  bar <- list(barTable_3,barTable_5)[[x]] %>% as.matrix() %>% as.data.frame()
  bar$nUMI <- rowSums(bar)
  cos <- do.call(cbind, lapply(1:length(res$df_list), function(x) {
    df <- data.frame(res$df_list[[x]]$cos.umi)
  }))
  colnames(cos) <- paste("CosUMI_Bar",1:length(res$df_list), sep="")
  
  df <- data.frame(res$assign_table,
                   res$umap,
                   bar,
                   cos)
  rownames(df) <- paste(names(inputFiles)[x], "#", revComplement(rownames(df)), "-1", sep = "")      
  df
})
barDF <- do.call(rbind, barDF)

barDF$Library <- gsub("#.*","",rownames(barDF))
colnames(barDF)[1:4] <- c("BarcodeAssign","BarcodeCount","DropletType","Call")

barDF$MULTI <- factor(barDF$Call,
                      levels = c("Bar1", "Bar2", "Bar3", "multiplet", "negative"),
                      labels = c("Bar1", "Bar2", "Bar3", "Doublet", "Negative"))

barDF$Donor <- factor(barDF$Call,
                      levels = c("Bar1", "Bar2", "Bar3", "multiplet", "negative"),
                      labels = c("DonorA", "DonorB", "DonorC", "Doublet", "Negative"))

barDF$Source <- factor(barDF$Call,
                       levels = c("Bar1", "Bar2", "Bar3", "multiplet", "negative"),
                       labels = c("Hemacare", "Cytologics", "Promab", "Doublet", "Negative"))

cells <- rownames(barDF)

proj_2 <- subsetArchRProject(proj_1,
                             outputDirectory = "proj_2",
                             cells = cells,
                             dropCells = T, force = T)

proj_2@cellColData <- cbind(proj_2@cellColData, barDF)

proj_2 <- addIterativeLSI(proj_2, force = T)
proj_2 <- addClusters(proj_2, dimsToUse = 1:30, force = T)
proj_2 <- addUMAP(proj_2, dimsToUse = 1:30, force = T)
proj_2 <- addImputeWeights(proj_2, dimsToUse = 1:30)


### Vireo Classifications ---------------------------------------------------

dir.create("vireo")
write.table(proj_2$CB_CellRanger[proj_2$Sample == 'Ill_labeled'], "vireo/cells_3.tsv", quote = F, row.names = F, col.names = F, sep = "\t")

## Run Vireo To Classify Cells by Genotype

vireoFiles <- list(
  Ill_labeled = "vireo/3/vireo_LargeVCF/donor_ids.tsv"
)

vireoDF <- lapply(1, function(x) {
  tmp <- read.table(vireoFiles[[x]], header = T, row.names = 1)
  # rownames(tmp) <- paste(names(vireoFiles)[x], "#", rownames(tmp), "-1", sep = "")
  rownames(tmp) <- paste(names(vireoFiles)[x], "#", rownames(tmp), sep = "")
  tmp
})
vireoDF <- do.call(rbind, vireoDF)

proj_2$Vireo <- vireoDF[proj_2$cellNames, "donor_id"]

table(proj_2$Vireo)

plotEmbedding(proj_2, name = "Vireo")



### Accessibility Analysis --------------------------------------------------

pathToMacs2 <- "/Users/dannyconrad/opt/miniconda3/bin/macs2"

proj_2 <- addGroupCoverages(proj_2, groupBy = "Clusters", force = T)
proj_2 <- addReproduciblePeakSet(
  ArchRProj = proj_2, 
  groupBy = "Clusters", 
  pathToMacs2 = pathToMacs2, 
  force = T)
proj_2 <- addPeakMatrix(proj_2, force = T)

## Plot Metadata
g <- plot_grid(plotlist = 
                 lapply(c("Sample", "Clusters", "MULTI", "Vireo", "log10(nUMI)", "log10(nFrags)", "TSSEnrichment", "FRIP"), function(x) {
                   plotEmbedding(proj_2, name = x, plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle(x) 
                 }),
               ncol = 4)


### Remove Low-Quality Cells ------------------------------------------------

## Identify low quality clusters
metadata <- c("Clusters", "Donor", "nUMI", "nFrags", "TSSEnrichment", "NucleosomeRatio", "PromoterRatio", "FRIP")

lapply(metadata, function(x) {
  if (x %in% c('nUMI','nFrags')) { x <- paste("log10(",x,")",sep="") }
  plotEmbedding(proj_2, name = x, plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle(x) 
}) %>% plot_grid(plotlist = ., nrow = 2)

lapply(c("Donor","Vireo"), function(x) {
  plotEmbedding(proj_2, name = x, plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle(x) 
}) %>% plot_grid(plotlist = ., nrow = 1)

lapply(metadata[3:8], function(x) {
  g <- proj_2@cellColData %>% as.data.frame %>% ggplot(aes(x = Donor, y = .data[[x]], fill = Donor)) + geom_boxplot(outliers = F)
  if (x %in% c('nUMI','nFrags')) { g <- g + scale_y_log10() }
  g
}) %>% plot_grid(plotlist = ., nrow = 2)

lapply(metadata[3:8], function(x) {
  g <- proj_2@cellColData %>% as.data.frame %>% ggplot(aes(x = Clusters, y = .data[[x]], fill = Clusters)) + geom_boxplot(outliers = F)
  if (x %in% c('nUMI','nFrags')) { g <- g + scale_y_log10() }
  g
}) %>% plot_grid(plotlist = ., nrow = 2)


plotEmbedding(proj_2, name = "FRIP >= 0.4", plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void()

cells <- proj_2$cellNames[proj_2$FRIP >= 0.4 & proj_2$Clusters %ni% c("C1","C4","C5","C17","C18")]


## Filter identified cells from project
plotEmbedding(proj_2, name = 'proj_2$cellNames %in% cells', plotAs = 'points', size = 1) + theme_void() + ggtitle("Cells Passing QC")

proj_3 <- subsetArchRProject(proj_2,
                             outputDirectory = "./proj_3",
                             cells = cells,
                             dropCells = T, 
                             force = T)

proj_3 <- addIterativeLSI(proj_3, force = T)
proj_3 <- addHarmony(proj_3, name = "Harmony", groupBy = "Donor", force = T)
proj_3 <- addClusters(proj_3, force = T, reducedDims = "Harmony")
proj_3 <- addUMAP(proj_3, force = T, reducedDims = "Harmony")

## Plot Metadata
g <- plot_grid(plotlist = 
                 lapply(c("Sample", "Clusters", "MULTI", "Vireo", "log10(nUMI)", "log10(nFrags)", "TSSEnrichment", "FRIP"), function(x) {
                   plotEmbedding(proj_3, name = x, plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle(x) 
                 }),
               ncol = 4)


### Compare Vireo & MULTI-ATAC Classifications ------------------------------

colors_bar <- c("#FAC056", "#7A41A8", "#45C845","grey","black")
names(colors_bar) <- c("Bar1","Bar2","Bar3","multiplet","negative")
colors_donor <- c("#FAC056", "#7A41A8", "#45C845","grey","black")
names(colors_donor) <- c("DonorA","DonorB","DonorC","Doublet","Unclassified")

tmp <- cbind(UMAP1 = proj_3@embeddings$UMAP$df$`Harmony#UMAP_Dimension_1`,
             UMAP2 = proj_3@embeddings$UMAP$df$`Harmony#UMAP_Dimension_2`,
             proj_3@cellColData %>% as.data.frame) %>% 
  mutate(MULTIATAC = factor(Call, levels = names(colors_bar), labels = names(colors_donor)),
         Vireo = factor(Vireo, levels = sort(unique(Vireo))[c(3,1,2,4,5)], labels = names(colors_donor)))

classDF <- tmp %>% dplyr::select(MULTIATAC, Vireo) %>% 
  group_by(MULTIATAC, Vireo, .drop = F) %>% summarize(Freq = n()) %>% 
  # mutate(Fill = Freq)
  # mutate(Fill = log1p(Freq))
  # mutate(Fill = Freq / sum(Freq))
  # group_by(Vireo) %>% mutate(Fill = Freq / sum(Freq))
  group_by(Vireo) %>% mutate(ColSum = sum(Freq)) %>% group_by(MULTIATAC) %>% mutate(RowSum = sum(Freq)) %>% ungroup() %>% mutate(Fill = Freq / (RowSum + ColSum - Freq))


g <- ggplot(classDF, aes(x = Vireo, y = MULTIATAC, fill = Fill, label = Freq)) + 
  geom_tile() + geom_text() + 
  theme_DC + theme(line = element_blank()) + labs(y = "MULTI-ATAC", fill = "Fraction") +
  scale_fill_viridis_c(option = "E", alpha = 0.7)
save_plot(g, "Plots/proj3_ClassificationMatrix.pdf", w = 4, h = 3, show.legend = "right")
g <- ggplot(classDF, aes(x = Vireo, y = MULTIATAC, fill = Fill, label = Freq)) + 
  geom_tile() + 
  # geom_text(size = 2) + 
  shadowtext::geom_shadowtext(size = 2, color = 'black', bg.color = 'white') +
  theme_DC + theme(line = element_blank()) + labs(y = "MULTI-ATAC", fill = "Fraction") +
  scale_fill_viridis_c(option = "E", alpha = 0.7)
save_plot(g, "Plots/proj3_ClassificationMatrix_v2.pdf", w = 2.75, h = 2, show.legend = "right")
g <- ggplot(classDF, aes(x = MULTIATAC, y = Vireo, fill = Fill, label = Freq)) + 
  geom_tile() + 
  # geom_text(size = 2) + 
  # shadowtext::geom_shadowtext(size = 2.75, color = 'black', bg.color = 'white', fontface = 'bold') +
  shadowtext::geom_shadowtext(size = 2.75, bg.color = 'black', color = 'white', fontface = 'bold', br.r = 0) +
  theme_DC + theme(line = element_blank(), axis.text.x = element_text(angle=-30,vjust=0,hjust=0.25)) + 
  labs(x = "MULTI-ATAC", fill = "Fraction") + scale_fill_viridis_c(option = "E", alpha = 0.7)
save_plot(g, "Plots/proj3_ClassificationMatrix_v2.pdf", w = 2.75, h = 2, show.legend = "right")



g <- ggplot(tmp, aes(x = UMAP1, y = UMAP2, color = MULTIATAC)) + 
  ggrastr::rasterise(geom_point(stroke = 0, size = 0.75), dpi = 600, scale = 1) +
  labs(title = "MULTI-ATAC", color = NULL) + scale_color_manual(values = colors_donor) + 
  guides(color = guide_legend(override.aes = list(size=3))) + 
  theme_DC + theme_UMAP + theme(legend.background = element_blank())
save_plot(g, "Plots/proj3_DonorClassifications_MULTIATAC.pdf", w = 2, h = 2, show.legend = c(0.85,0.225))
g <- ggplot(tmp, aes(x = UMAP1, y = UMAP2, color = Vireo)) + 
  ggrastr::rasterise(geom_point(stroke = 0, size = 0.75), dpi = 600, scale = 1) +
  labs(title = "Vireo", color = NULL) + scale_color_manual(values = colors_donor) + 
  guides(color = guide_legend(override.aes = list(size=3))) + 
  theme_DC + theme_UMAP + theme(legend.background = element_blank())
save_plot(g, "Plots/proj3_DonorClassifications_Vireo.pdf", w = 2, h = 2, show.legend = c(0.85,0.225))



### Independent Doublet Identification --------------------------------------

## ArchR Doublet Scoring
proj_3 <- addDoubletScores(proj_3, knnMethod = "LSI")
plot_grid(plotEmbedding(proj_3, name = "DoubletEnrichment", size = 1, plotAs = 'points') + theme_void(),
          plotEmbedding(proj_3, name = "Call", size = 1, plotAs = 'points') + theme_void(),
          plotEmbedding(proj_3, name = "Vireo", size = 1, plotAs = 'points') + theme_void(),
          nrow = 1)

## AMULET
dir.create("../AMULET/")

csv <- read.csv("/Volumes/DannySSD/MULTI_ATAC/20211102-PilotRedo/cellranger-atac/gDNA3/singlecell.csv")
csv <- csv[match(proj_3$CB_CellRanger, csv$barcode),] 

proj_3@cellColData <- proj_3@cellColData[,1:39]
proj_3@cellColData <- cbind(proj_3@cellColData, csv)

csv$is__cell_barcode <- 1
write.table(csv, "../AMULET/singlecell_3.csv", row.names = F, col.names = T, quote = F, sep = ",")

## in Terminal:
#AMULET.sh /path/to/fragments.tsv.gz /path/to/singlecell.csv /path/to/human_autosomes.txt /path/to/repeatfilter.bed /path/to/output/ /path/to/shellscript/
AMULET="/Volumes/DannySSD/Software/AMULET/"
AUTO="/Volumes/DannySSD/Software/AMULET/human_autosomes.txt"
REP="/Volumes/DannySSD/Software/AMULET/blacklist_repeats_segdups_rmsk_hg38.bed"

OUT="/Volumes/DannySSD/MULTI_ATAC/20211102-PilotRedo/AMULET"
outs1="/Volumes/DannySSD/MULTI_ATAC/20211102-PilotRedo/cellranger-atac/gDNA3/"

mkdir $OUT/3
$AMULET/AMULET.sh $outs1/fragments.tsv.gz $OUT/singlecell_3.csv $AUTO $REP $OUT/3 $AMULET

## Import AMULET multiplet calls
amulet <- paste("Ill_labeled#", read.table("/Volumes/DannySSD/MULTI_ATAC/20211102-PilotRedo/AMULET/3/MultipletBarcodes_01.txt")$V1, sep = "")
proj_3$AMULET <- "Singlet"
proj_3@cellColData[amulet,"AMULET"] <- "Multiplet"
proj_3$AMULET <- factor(proj_3$AMULET, levels = c("Singlet","Multiplet"))

plotEmbedding(proj_3, name = "AMULET", size = 1, labelAsFactors = F, labelMeans = F) + theme_void()
plotEmbedding(proj_3, name = "Donor == 'Doublet'", size = 1, labelAsFactors = F, labelMeans = F) + theme_void()
plotEmbedding(proj_3, name = "is__cell_barcode == 1", size = 1, labelAsFactors = F, labelMeans = F) + theme_void()

table(proj_3@cellColData[,c("Sample","AMULET")])

## Compare Doublets Classified by MULTI-ATAC to those by Vireo & AMULET

# Just Vireo
proj_3$DoubletCall <- proj_3@cellColData %>% as.data.frame %>%
  mutate(DoubletCall = case_when(
    MULTI == "Doublet" & Vireo == "doublet" ~ "Both",
    MULTI == "Doublet" ~ "MULTI-ATAC Only",
    Vireo == "doublet" ~ "Vireo Only",
    TRUE ~ "Singlet"
  )) %>% .$DoubletCall

# plotEmbedding(proj_3, name = "DoubletCall", size = 1) + theme_void()

g <- proj_3@cellColData %>% as.data.frame %>% mutate(DoubletCall = factor(DoubletCall, levels = c("Both","MULTI-ATAC Only","Vireo Only", "Singlet"))) %>%
  ggplot(aes(x = as.numeric(DoubletCall), y = log10(nFrags), fill = DoubletCall)) + 
  # geom_boxplot(outliers = F) +
  geom_violin() +
  # geom_text(data = proj_3@cellColData %>% as.data.frame %>% group_by(DoubletCall) %>% summarize(Count = n()), mapping = aes(x = as.numeric(DoubletCall), y = 3, label = Count)) +
  theme_DC + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + labs(fill = NULL, x = "Doublet Call") + #scale_y_log10() + 
  stat_compare_means(method = 't.test', size = 2,
                     comparisons = list(c(1,2),c(1,3),c(2,4),c(3,4)) %>% rev) 
save_plot(g, "./Plots/proj3_DoubletClassificationComparison_nFrags.pdf", w = 3, h = 2, show.legend = 'right')

g <- proj_3@cellColData %>% as.data.frame %>% mutate(DoubletCall = factor(DoubletCall, levels = c("Both","MULTI-ATAC Only","Vireo Only", "Singlet"))) %>%
  ggplot(aes(x = as.numeric(DoubletCall), y = log10(DoubletEnrichment + 1), fill = DoubletCall)) + 
  # geom_boxplot(outliers = F) +
  geom_violin() +
  # geom_text(data = proj_3@cellColData %>% as.data.frame %>% group_by(DoubletCall) %>% summarize(Count = n()), mapping = aes(x = as.numeric(DoubletCall), y = 3, label = Count)) +
  theme_DC + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + labs(fill = NULL, x = "Doublet Call") + scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  stat_compare_means(method = 't.test', size = 2,
                     comparisons = list(c(1,2),c(1,3),c(2,4),c(3,4)) %>% rev) 
save_plot(g, "./Plots/proj3_DoubletClassificationComparison_DoubletEnrichment.pdf", w = 3, h = 2, show.legend = 'right')


## Vireo & AMULET
proj_3$DoubletCall <- proj_3@cellColData %>% as.data.frame %>%
  mutate(DoubletCall = case_when(
    MULTI == "Doublet" & Vireo == "doublet" & AMULET == "Multiplet" ~ "Consensus Doublet",
    MULTI == "Doublet" & Vireo == "doublet" ~ "MULTI-ATAC + Vireo",
    MULTI == "Doublet" & AMULET == "Multiplet" ~ "MULTI-ATAC + AMULET",
    Vireo == "doublet" & AMULET == "Multiplet" ~ "Vireo + AMULET",
    MULTI == "Doublet" ~ "MULTI-ATAC Only",
    Vireo == "doublet" ~ "Vireo Only",
    AMULET == "Multiplet" ~ "AMULET Only",
    TRUE ~ "Singlet"
  )) %>% .$DoubletCall


proj_3$SingletCall <- proj_3@cellColData %>% as.data.frame %>% 
  mutate(SingletCall = case_when(
    MULTI %ni% c("Doublet","Negative") & Vireo %ni% c("doublet","unassigned") & AMULET == "Singlet" ~ "Consensus Singlet",
    MULTI %ni% c("Doublet","Negative") & Vireo %ni% c("doublet","unassigned") ~ "MULTI-ATAC + Vireo",
    MULTI %ni% c("Doublet","Negative") & AMULET == "Singlet" ~ "MULTI-ATAC + AMULET",
    Vireo %ni% c("doublet","unassigned") & AMULET == "Singlet" ~ "Vireo + AMULET",
    MULTI %ni% c("Doublet","Negative") ~ "MULTI-ATAC Only",
    Vireo %ni% c("doublet","unassigned") ~ "Vireo Only",
    AMULET == "Singlet" ~ "AMULET Only",
    TRUE ~ "Doublet"
  )) %>% .$SingletCall

## Plot Venn Diagrams
tmp <- proj_3@cellColData %>% as.data.frame() %>% filter(MULTI != "Negative", Vireo != "unassigned")
tmp <- proj_3@cellColData %>% as.data.frame() #%>% filter(MULTI != "Negative", Vireo != "unassigned")


g <- ggVennDiagram::ggVennDiagram(list(`MULTI-ATAC` = rownames(tmp)[tmp$MULTI == "Doublet"],
                                       Vireo = rownames(tmp)[tmp$Vireo == "doublet"],
                                       AMULET = rownames(tmp)[tmp$AMULET == "Multiplet"]),
                                  label_alpha = 0, label = "count", label_size = 3, set_size = 2.5, edge_size = 0.5) + paletteer::scale_fill_paletteer_c("ggthemes::Red-Gold") + 
  labs(x = NULL, y = NULL, fill = "Count", title = "Doublet Classification") +
  theme_DC + coord_cartesian(clip = 'off', xlim = c(-6,10)) + theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())#, plot.title.position = 0.6)
save_plot(g, file = "Plots/proj3_DoubletClassificationVenn_v2.pdf", w = 2.25, h = 2, show.legend = c(0.95,0.5))

g <- ggVennDiagram::ggVennDiagram(list(`MULTI-ATAC` = rownames(tmp)[tmp$MULTI != "Doublet"],
                                       Vireo = rownames(tmp)[tmp$Vireo != "doublet"],
                                       AMULET = rownames(tmp)[tmp$AMULET == "Singlet"]),
                                  label_alpha = 0, label = "count", label_size = 3, set_size = 2.5, edge_size = 0.5) + paletteer::scale_fill_paletteer_c("ggthemes::Red-Gold") + 
  labs(x = NULL, y = NULL, fill = "Count", title = "Singlet Classification") +
  theme_DC + coord_cartesian(clip = 'off', xlim = c(-6,10)) + theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())#, plot.title.position = 0.6)
save_plot(g, file = "Plots/proj3_SingletClassificationVenn_v2.pdf", w = 2.25, h = 2, show.legend = c(0.95,0.5))





### Export R Objects --------------------------------------------------------

dir.create("Export")

sessionOut("Export/MultidonorPBMCPilot_Packages.tsv")

## Single-cell metadata
saveRDS(proj_3@cellColData, file = "Export/MultidonorPBMCPilot_SingleCellMetadata.rds")










