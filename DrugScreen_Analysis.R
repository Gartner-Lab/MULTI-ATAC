### Finalized Analysis of Multiome MULTI-ATAC Dataset ###

### Set Up Environment ---------------------------------------------------------

library(ArchR)
library(Seurat)
library(Signac)
library(deMULTIplex2)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
# library(ComplexHeatmap)
library(pheatmap)


# library(GenomicRanges)
library(rtracklayer)

library(msigdbr)
library(fgsea)

library(tibble)
library(data.table)
library(tidyr)
library(dplyr)

source("/Volumes/DannySSD/MULTI_ATAC/CustomFunctions.R")

addArchRThreads(threads = 1) 
addArchRGenome("hg38")

wd <- "/Volumes/DannySSD/MULTI_ATAC/20231114_PROTACExperiment2/ArchR"
dir.create(wd)
setwd(wd)

dir.create("Plots")


### Define Colors Palettes --------------------------------------------------

library(RColorBrewer)
library(paletteer)

colors <- c(
  paletteer_d("RColorBrewer::PuOr",n = 11)[c(8:10)],
  paletteer_d("RColorBrewer::BrBG",n = 11)[c(8:10)],
  paletteer_d("RColorBrewer::BrBG",n = 11)[c(4:2)]
)
# colors %>% scales::show_col()

drug <- c("DMSO(-)","DMSO(+)","MS177","AU-15330","dCBP-1","EPZ-6438","BRM014","GNE-781")
drugdose <- c("DMSO(-)","DMSO(+)", paste(rep(c("MS177","AU-15330","dCBP-1","EPZ-6438","BRM014","GNE-781"), each = 3), c("10nM","100nM","1µM")))
celltype <- c("T","Myeloid","B")

# Classification Type
col_multi <- c("#208A42","#D51F26","#272E6A")
names(col_multi) <- c("Singlet","Multiplet","Negative")

# Drug
col_drug <- c("#D6604DFF","grey80","grey65",colors[c(2,3,5,6,8,9)]) 
names(col_drug) <- c("Activation","DMSO(-)","DMSO(+)","MS177","EPZ-6438","AU-15330","BRM014","dCBP-1","GNE-781")
col_drug2 <- c("#D6604DFF","grey80","grey65",colors[c(2,3,5,6,8,9)]) 
names(col_drug2) <- c("Activation","NegCtrl","PosCtrl","MS177","EPZ-6438","AU-15330","BRM014","dCBP-1","GNE-781")

# Dose
col_dose <- c("grey80","grey65",viridis::magma(n = 5)[4:2])
names(col_dose) <- c("DMSO(-)","DMSO(+)","10nM","100nM","1µM")
col_dose2 <- c("grey80","grey65",viridis::magma(n = 5)[4:2])
names(col_dose2) <- c("NegCtrl","PosCtrl","10nM","100nM","1µM")

# Drug Target
col_target <- c("grey80","grey65",colors[c(2,5,8)])
names(col_target) <- c("DMSO(-)", "DMSO(+)", "PRC2", "SWI/SNF", "p300/CBP")

# Plate Side
col_side <- c("#984ea3","#ff7f00")
names(col_side) <- c("L","R")

# Major Cell Types
col_cellmaj <- ArchRPalettes$grove[as.character(1:3)]
names(col_cellmaj) <- names(table(proj_4$CellType_Major))

# Minor Cell Types
col_cellmin <- ArchRPalettes$circus[as.character(1:(length(unique(proj_4$CellType_Minor))+1))[-3]]
names(col_cellmin) <- names(table(proj_4$CellType_Minor))

# Hierarchical Clusters
col_clust <- paletteer::paletteer_d("ggsci::category10_d3", n = 10) %>% as.character()
names(col_clust) <- 1:10


### Create Arrow Files ---------------------------------------------------------

# User cellranger-arc output files (lane 3 was lost due to microfluidic clog)
inputFiles <- c(
  L1 = "/Volumes/DannySSD/MULTI_ATAC/20231114_PROTACExperiment2/outs1/atac_fragments.tsv.gz",
  L2 = "/Volumes/DannySSD/MULTI_ATAC/20231114_PROTACExperiment2/outs2/atac_fragments.tsv.gz",
  L4 = "/Volumes/DannySSD/MULTI_ATAC/20231114_PROTACExperiment2/outs4/atac_fragments.tsv.gz"
)

ArrowFiles <- vector()

for (i in 1:length(inputFiles)) {
  ArrowFiles[i] <- createArrowFiles(
    inputFiles = inputFiles[i],
    sampleNames = names(inputFiles)[i],
    minTSS = 2, 
    minFrags = 100, 
    addTileMat = TRUE,
    addGeneScoreMat = TRUE)
}

ArrowFiles <- c("L1.arrow", "L2.arrow", "L4.arrow")

### Create ArchR Project -------------------------------------------------------

proj_1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "proj_1",
  copyArrows = TRUE
)

g <- ggplot(as.data.frame(proj_1@cellColData), aes(x = log10(nFrags), y = TSSEnrichment, color = TSSEnrichment > (-13*log10(nFrags) + 45))) + 
  ggrastr::rasterise(geom_point(size = 1, alpha = 0.3, stroke = 0, color = "grey"), dpi = 600, scale = 1) +
  geom_abline(slope = -13, intercept = 45) + 
  geom_density_2d(bins = 100) + 
  labs(color = "PassFilter") + 
  theme_DC +
  scale_color_manual(values = c("black","dodgerblue")) +
  facet_grid(.~Sample) + 
  ylim(c(0,quantile(as.data.frame(proj_1@cellColData)$TSSEnrichment, 0.999)))

save_plot(g, "./Plots/proj1_PassFilter_Density.pdf", w = 5, h = 2.5, show.legend = "none")

# import scRNA
GEX <- import10xFeatureMatrix(
  input = c("../outs1/raw_feature_bc_matrix.h5",
            "../outs2/raw_feature_bc_matrix.h5",
            "../outs4/raw_feature_bc_matrix.h5"),
  names = names(inputFiles) 
)

proj_1 <- addGeneExpressionMatrix(input = proj_1, seRNA = GEX, force = TRUE)

### Filter By ATAC & GEX QC ----------------------------------------------------

## nFrags, TSS enrichment, and nUMI
## Use permissive/generous filtering because of epigenetic drugs used in the experiment

g <- ggplot(as.data.frame(proj_1@cellColData), aes(x = log10(nFrags), y = TSSEnrichment, color = TSSEnrichment > (-13*log10(nFrags) + 45) & Gex_nUMI > 10)) + 
  # geom_point(alpha = 0.3, size = 0.5) + 
  ggrastr::rasterise(geom_point(size = 1, alpha = 0.3, stroke = 0), dpi = 600, scale = 1) +
  geom_abline(slope = -13, intercept = 45) + geom_vline(xintercept = 2) + geom_hline(yintercept = 2) + 
  labs(color = "PassFilter", title = "TSSEnrichment > (-13*log10(nFrags) + 45) & Gex_nUMI > 10") + 
  theme_DC +
  scale_color_manual(values = c("black","dodgerblue")) +
  facet_grid(~Sample) +
  ylim(c(0,quantile(as.data.frame(proj_1@cellColData)$TSSEnrichment, 0.999)))

save_plot(g, "./Plots/proj1_PassFilter.pdf", w = 5, h = 2.5, show.legend = "none")

cells <- proj_1$cellNames[proj_1$TSSEnrichment > (-13*log10(proj_1$nFrags) + 45) &
                            proj_1$Gex_nUMI > 10]

proj_2 <- subsetArchRProject(proj_1,
                             cells = cells,
                             outputDirectory = "proj_2",
                             dropCells = F, 
                             force = T)
# proj_2 <- readRDS("proj_2/Save-ArchR-Project.rds")

## Compute Dim Reduction and UMAP embedding with proper scaling

proj_2 <- dimReduce(proj_2, scale = "Col", default = "NoScale", dims = 2:30)


### Add CellRanger Metadata ----------------------------------------------------

csv1 <- read.csv("/Volumes/DannySSD/MULTI_ATAC/20231114_PROTACExperiment2/outs1/per_barcode_metrics.csv")
csv2 <- read.csv("/Volumes/DannySSD/MULTI_ATAC/20231114_PROTACExperiment2/outs2/per_barcode_metrics.csv")
csv4 <- read.csv("/Volumes/DannySSD/MULTI_ATAC/20231114_PROTACExperiment2/outs4/per_barcode_metrics.csv")
rownames(csv1) <- paste("L1#", csv1$barcode, sep = "")
rownames(csv2) <- paste("L2#", csv2$barcode, sep = "")
rownames(csv4) <- paste("L4#", csv4$barcode, sep = "")
csv <- rbind(csv1,csv2,csv4)
csv <- csv[proj_2$cellNames,]
csv$Library <- gsub("#.*", "", rownames(csv))

rm(csv1,csv2,csv4)

## Alternate cell barcodes
# For Multiome, barcode in fragment file is "GEX barcode", not "ATAC barcode"

proj_2$CB_CellRangerGEX <- gsub(".*#", "", rownames(proj_2@cellColData))
proj_2$CB_CellRangerATAC <- csv$atac_barcode
proj_2$CB <- gsub("-.*", "", proj_2$CB_CellRangerATAC)
proj_2$CB_ReverseComp <- revComplement(proj_2$CB)

## ATAC data
proj_2$ATAC_Reads <- csv$atac_raw_reads
plotEmbedding(proj_2, name = "ATAC_Reads > 25000", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("ATAC_Reads > 25000")

proj_2$ATAC_PercentMito <- csv$atac_mitochondrial_reads / csv$atac_raw_reads * 100
plotEmbedding(proj_2, name = "ATAC_PercentMito", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("ATAC_PercentMito")

## Plot Embedding & Metadata
g <- plot_grid(plotlist = lapply(c("Clusters", "log10(nFrags)", "TSSEnrichment", "log10(Gex_nUMI)", "Gex_MitoRatio", "ATAC_PercentMito"), function(x) {
  g <- plotEmbedding(proj_2, name = x, embedding = "UMAP", plotAs = 'points', size = 0.75, labelAsFactors=F, labelMeans=T) + labs(x = "UMAP1", y = "UMAP2") +
    theme_DC + theme_UMAP + ggtitle(x)
  if (x == "Clusters") { g <- g + theme(legend.position = 'none') }
  g
}), ncol = 3) 
save_plot(g, "./Plots/proj2_QC.pdf", w = 6, h = 4, show.legend = "none")


### MULTI-ATAC Classification --------------------------------------------------

dir.create(paste(wd, "SampleClassification",sep="/"))

bar.ref <- read.csv("/Users/dannyconrad/Library/CloudStorage/Box-Box/Data/BarcodeSetRedesign/Final/MULTI-ATAC_Barcodes_96.csv")
bar.ref.rev <- revComplement(bar.ref$Barcode_Sequence,cbLength = 8)
names(bar.ref.rev) <- bar.ref$Well_Position

# Lane 1
readTable_1 <- readTags(dir = "/Volumes/GDRIVEPRO/CAT_Server_Data/RawData/MULTIATAC_PROTAC_Multiome2/",
                        name = "Bar1_PreSPRI", 
                        barcode.type = "MULTI-ATAC", 
                        assay = "Multiome")
barTable_1 <- alignTags(read_table = readTable_1,
                        tag.ref = bar.ref.rev,
                        filter.cells = proj_2$CB_ReverseComp[proj_2$Sample == "L1"])
res_1 <- demultiplexTags(barTable_1, 
                         plot.diagnostics = T, 
                         plot.path = "SampleClassification/",
                         plot.name = "L1pre")

# Lane 2
readTable_2 <- readTags(dir = "/Volumes/GDRIVEPRO/CAT_Server_Data/RawData/MULTIATAC_PROTAC_Multiome2/",
                        name = "Bar2_PreSPRI", 
                        barcode.type = "MULTI-ATAC", 
                        assay = "Multiome")
barTable_2 <- alignTags(read_table = readTable_2,
                        tag.ref = bar.ref.rev,
                        filter.cells = proj_2$CB_ReverseComp[proj_2$Sample == "L2"])
res_2 <- demultiplexTags(barTable_2, 
                         plot.diagnostics = T, 
                         plot.path = "SampleClassification/",
                         plot.name = "L2pre")

# Lane 4
readTable_4 <- readTags(dir = "/Volumes/GDRIVEPRO/CAT_Server_Data/RawData/MULTIATAC_PROTAC_Multiome2/",
                        name = "Bar4_PreSPRI", 
                        barcode.type = "MULTI-ATAC", 
                        assay = "Multiome")
barTable_4 <- alignTags(read_table = readTable_4,
                        tag.ref = bar.ref.rev,
                        filter.cells = proj_2$CB_ReverseComp[proj_2$Sample == "L4"])
res_4 <- demultiplexTags(barTable_4, 
                         plot.diagnostics = T, 
                         plot.path = "SampleClassification/",
                         plot.name = "L4pre")

# saveRDS(readTable_1, "SampleClassification/readTable_1pre.rds")
# saveRDS(readTable_2, "SampleClassification/readTable_2pre.rds")
# saveRDS(readTable_4, "SampleClassification/readTable_4pre.rds")
# saveRDS(barTable_1, "SampleClassification/barTable_1pre.rds")
# saveRDS(barTable_2, "SampleClassification/barTable_2pre.rds")
# saveRDS(barTable_4, "SampleClassification/barTable_4pre.rds")
# saveRDS(res_1, "SampleClassification/res_1pre.rds")
# saveRDS(res_2, "SampleClassification/res_2pre.rds")
# saveRDS(res_4, "SampleClassification/res_4pre.rds")

# barTable_1 <- readRDS("SampleClassification/barTable_1pre.rds")
# barTable_2 <- readRDS("SampleClassification/barTable_2pre.rds")
# barTable_4 <- readRDS("SampleClassification/barTable_4pre.rds")
# res_1 <- readRDS("SampleClassification/res_1pre.rds")
# res_2 <- readRDS("SampleClassification/res_2pre.rds")
# res_4 <- readRDS("SampleClassification/res_4pre.rds")

calls <- list(res_1$final_assign, res_2$final_assign, res_4$final_assign)
saveRDS(calls, "SampleClassification/calls.rds")

tagCallHeatmap(barTable_1, res_1$final_assign)
tagCallHeatmap(barTable_2, res_2$final_assign)
tagCallHeatmap(barTable_4, res_4$final_assign)

g <- tagCallHeatmap(barTable_1, res_1$final_assign) + scale_fill_viridis_c(option="H")
g[[2]][[2]] <- NULL
pdfPlot(g, width = 15, height = 15, "./Plots/res1_HeatmapV2.pdf")

g <- cbind(res_1$umap,res_1$assign_table) %>% 
  {
    ggplot(., aes(x = UMAP_1, y = UMAP_2, color = droplet_type, size = droplet_type)) +
      scale_size_manual(values = c(0.5,0.5,1)) + 
      ggrastr::rasterise(geom_point(alpha = 0.6, stroke = 0), dpi = 600, scale = 1) +
      scale_color_manual(values = c("#AD72D6FF","#5AAE61FF","grey10")) +
      theme_DC + theme_UMAP +
      labs(color = "Droplet Type", x = "UMAP1", y = "UMAP2", title = "deMULTIplex2 Classifications") + #, subtitle = "Lane 1 of 3") + 
      ggrepel::geom_label_repel(
        data = . %>% group_by(droplet_type, final_assign) %>% summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) %>% filter(droplet_type == "singlet"),
        aes(label = final_assign), 
        label.padding = 0.08, max.overlaps = 20, color = 'black', alpha = 0.8, size = 1)
  }
save_plot(g, "./Plots/res1_UMAP.pdf", w = 2, h = 2, show.legend = 'none')

g <- cbind(res_1$umap,res_1$assign_table) %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = final_assign, size = droplet_type, alpha = droplet_type)) +
  scale_size_manual(values = c(0.5,0.5,1)) + scale_alpha_manual(values = c(0.2,0.2,0.5)) + 
  ggrastr::rasterise(geom_point(stroke = 0), dpi = 600, scale = 1) +
  scale_color_manual(values = c(colorRampPalette(c("navyblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1", "maroon"))(96),"black","grey")) + 
  theme_DC + theme_UMAP +
  labs(x = "UMAP1", y = "UMAP2", title = "deMULTIplex2 Classifications") #, subtitle = "Lane 1 of 3") 
save_plot(g, "./Plots/res1_UMAP_Rainbow.pdf", w = 2, h = 2, show.legend = 'none')

g <- cbind(res_1$umap,res_1$assign_table) %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = final_assign, size = droplet_type, alpha = droplet_type)) +
  scale_size_manual(values = c(0,0,1)) + scale_alpha_manual(values = c(0,0,0.5)) + 
  ggrastr::rasterise(geom_point(stroke = 0), dpi = 600, scale = 1) +
  scale_color_manual(values = c(colorRampPalette(c("navyblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1", "maroon"))(96),"black","grey")) + 
  theme_DC + theme_UMAP +
  labs(x = "UMAP1", y = "UMAP2", title = "deMULTIplex2 Classifications", subtitle = "Lane 1 of 3")
save_plot(g, "./Plots/res1_UMAP_Rainbow_SingletsOnly.pdf", w = 2, h = 2, show.legend = 'none')


### Merge Classifications ------------------------------------------------------

proj_2$MULTI[proj_2$Sample == "L1"] <- calls[[1]][as.vector(proj_2@cellColData[proj_2$Sample == "L1",]$CB_ReverseComp)]
proj_2$MULTI[proj_2$Sample == "L2"] <- calls[[2]][as.vector(proj_2@cellColData[proj_2$Sample == "L2",]$CB_ReverseComp)]
proj_2$MULTI[proj_2$Sample == "L4"] <- calls[[3]][as.vector(proj_2@cellColData[proj_2$Sample == "L4",]$CB_ReverseComp)]
proj_2$MULTI[is.na(proj_2$MULTI)] <- "negative"

proj_2$tagUMI[proj_2$Sample == "L1"] <- barTable_1[proj_2$CB_ReverseComp[proj_2$Sample == "L1"],] %>% Matrix::rowSums()
proj_2$tagUMI[proj_2$Sample == "L2"] <- barTable_2[proj_2$CB_ReverseComp[proj_2$Sample == "L2"],] %>% Matrix::rowSums()
proj_2$tagUMI[proj_2$Sample == "L4"] <- barTable_4[proj_2$CB_ReverseComp[proj_2$Sample == "L4"],] %>% Matrix::rowSums()

plotEmbedding(proj_2, name = "log10(tagUMI)", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("log10(tagUMI)")

rm(list = ls(pattern = "barTable"))
rm(list = ls(pattern = "res_"))

## Breakdown Classifications into Singlets, Multiplets, Negatives
cbind(Singlets = lapply(calls, function(x) { sum(x %ni% c("multiplet","negative")) }),
      Multiplets = lapply(calls, function(x) { sum(x == "multiplet") }),
      Negatives = lapply(calls, function(x) { sum(x == "negative") })) %>% `rownames<-`(paste("L",c(1:2,4),sep=""))

proj_2$MULTItype <- factor(proj_2$MULTI, levels = c(names(bar.ref.rev),"multiplet","negative"), labels = c(rep("Singlet",96),"Multiplet","Negative"))

plotEmbedding(proj_2, name = "MULTItype", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("MULTItype")
plotEmbedding(proj_2, name = "MULTItype == 'Singlet'", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("Singlets") + theme(legend.position = "none") + scale_color_manual(values = c("grey","blue"))
plotEmbedding(proj_2, name = "MULTItype == 'Negative'", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("Negatives") + theme(legend.position = "none") + scale_color_manual(values = c("grey","blue"))
plotEmbedding(proj_2, name = "MULTItype == 'Multiplet'", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("Multiplets") + theme(legend.position = "none") + scale_color_manual(values = c("grey","blue"))


### Translate Classifications to Sample Metadata -------------------------------

metadata <- data.frame(Barcode = names(bar.ref.rev),
                       Drug = rep(c("DMSO(-)","DMSO(+)","MS177","AU-15330","dCBP-1","EPZ6438","BRM014","GNE-781"), each = 12),
                       Dose = c(rep(c("DMSO(-)","DMSO(+)"), each = 12), rep(c("10nM","100nM","1µM"), 24)),
                       Type = c(rep("Control", 24), rep("PROTAC", 36), rep("Inhibitor", 36)),
                       Target = c(rep(c("DMSO(-)","DMSO(+)"), each = 12), rep(rep(c("PRC2","SWI/SNF","p300/CBP"), each = 12), 2)),
                       Replicate = c(rep(1:12, 2), rep(1:4, 6, each = 3)),
                       Plate = rep(c("A","B"), 8, each = 6),
                       PlateColumn = rep(c(1:3,10:12), 16),
                       PlateSide = rep(c("L","R"), each = 3, 16))

matrix(metadata$Barcode, nrow = 8, ncol = 12, byrow = T)
matrix(metadata$Dose, nrow = 8, ncol = 12, byrow = T)
matrix(metadata$Drug, nrow = 8, ncol = 12, byrow = T)
matrix(metadata$Type, nrow = 8, ncol = 12, byrow = T)
matrix(metadata$Target, nrow = 8, ncol = 12, byrow = T)
matrix(metadata$Replicate, nrow = 8, ncol = 12, byrow = T)
matrix(metadata$Plate, nrow = 8, ncol = 12, byrow = T)
matrix(metadata$PlateColumn, nrow = 8, ncol = 12, byrow = T)
matrix(metadata$PlateSide, nrow = 8, ncol = 12, byrow = T)

metadata <- rbind(metadata, "multiplet", "negative")

proj_2$Drug <- metadata$Drug[match(proj_2$MULTI, metadata$Barcode)]
proj_2$Dose <- metadata$Dose[match(proj_2$MULTI, metadata$Barcode)]
proj_2$Type <- metadata$Type[match(proj_2$MULTI, metadata$Barcode)]
proj_2$Target <- metadata$Target[match(proj_2$MULTI, metadata$Barcode)]
proj_2$Replicate <- metadata$Replicate[match(proj_2$MULTI, metadata$Barcode)]
proj_2$Plate <- metadata$Plate[match(proj_2$MULTI, metadata$Barcode)]
proj_2$PlateColumn <- metadata$PlateColumn[match(proj_2$MULTI, metadata$Barcode)]
proj_2$PlateSide <- metadata$PlateSide[match(proj_2$MULTI, metadata$Barcode)]

proj_2@cellColData$Drug <- factor(proj_2$Drug,
                                  levels = c("DMSO(-)","DMSO(+)","MS177","AU-15330","dCBP-1","EPZ6438","BRM014","GNE-781"),
                                  labels = c("DMSO(-)","DMSO(+)","MS177","AU-15330","dCBP-1","EPZ6438","BRM014","GNE-781"))
proj_2@cellColData$Dose <- factor(proj_2$Dose,
                                  levels = c("DMSO(-)","DMSO(+)","10nM","100nM","1µM"),
                                  labels = c("DMSO(-)","DMSO(+)","10nM","100nM","1µM"))

## Sample Metadata Embeddings
plotEmbedding(proj_2, name = "Plate", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("Plate") + theme(legend.position = "bottom")
plotEmbedding(proj_2, name = "Target", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("Target") + theme(legend.position = "bottom")
plotEmbedding(proj_2, name = "Type", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("Type") + theme(legend.position = "bottom")
plotEmbedding(proj_2, name = "Drug", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("Drug") + theme(legend.position = "bottom")
plotEmbedding(proj_2, name = "Dose", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("Dose") + theme(legend.position = "bottom")


### Remove Doublets & Negatives ------------------------------------------------

## Cluster Composition
plotEmbedding(proj_2, embedding = "UMAP", name = 'MULTI %in% c("multiplet","negative",NA)') + theme_void()
plotEmbedding(proj_2, embedding = "UMAP", name = 'Clusters', size = 0.1) + theme_void() + ggtitle("Clusters") + theme(legend.position = "none")

ggplot(as.data.frame(proj_2@cellColData), aes(x = factor(Clusters, levels = paste("C",1:length(unique(Clusters)),sep="")), fill = MULTI == 'multiplet')) + 
  geom_bar(position = 'fill') + theme_minimal() + xlab("Cluster") + ggtitle("Fraction Doublets") # C17 is all doublets
ggplot(as.data.frame(proj_2@cellColData), aes(x = factor(Clusters, levels = paste("C",1:length(unique(Clusters)),sep="")), fill = MULTI == 'negative')) + 
  geom_bar(position = 'fill') + theme_minimal() + xlab("Cluster") + ggtitle("Fraction Negatives") # C13 about 50% negatives
plotGroups(proj_2, groupBy = 'Clusters', name = "log10(nFrags)", plotAs = "violin") + theme_classic() + ggtitle("nFrags") # C13 also has lowest nFrags


cells <- proj_2$cellNames[proj_2$MULTItype == 'Singlet' & 
                            proj_2$Clusters %ni% c("C13","C17")]
plotEmbedding(proj_2, name = 'proj_2$cellNames %in% cells', embedding = "UMAP", size = 1) + theme_void() + ggtitle("Keep Cells") + scale_color_manual(values = c("black","dodgerblue"))     

g <- ggplot(as.data.frame(proj_2@cellColData), aes(x = log10(nFrags), y = log10(Gex_nUMI), color = proj_2$cellNames %in% cells)) + 
  ggrastr::rasterise(geom_point(size = 1, alpha = 0.3, stroke = 0), dpi = 600, scale = 1) +
  labs(color = "Keep Cells", title = "Singlets") + theme_DC + 
  scale_color_manual(values = c("black","dodgerblue")) + facet_grid(~Sample)
save_plot(g, "./Plots/proj2_KeepCells.pdf", w = 5, h = 2.5, show.legend = 'none')

proj_3 <- subsetArchRProject(proj_2,
                             cells = cells,
                             outputDirectory = "proj_3",
                             dropCells = F, force = T)
# proj_3 <- readRDS("proj_3/Save-ArchR-Project.rds")

proj_3 <- dimReduce(proj_3, scale = "Col", default = "NoScale", dims = 2:30)
proj_3 <- addImputeWeights(proj_3, reducedDims = "IterativeLSI", dimsToUse = 2:30, corCutOff = 1, scaleDims = F)

## Plot Embedding
p1 <- plotEmbedding(proj_3, name = "Clusters", embedding = "UMAP", size = 1, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + ggtitle("Clusters")
p2 <- plotEmbedding(proj_3, name = "log10(nFrags)", embedding = "UMAP", size = 1, labelAsFactors=F, labelMeans=T, plotAs="points") + theme_void() + ggtitle("ATAC nFrags")
p3 <- plotEmbedding(proj_3, name = "log10(Gex_nUMI)", embedding = "UMAP", size = 1, labelAsFactors=F, labelMeans=T, plotAs="points") + theme_void() + ggtitle("GEX nUMI")
plot_grid(p1,p2,p3, nrow = 1)

# Plot Sample Metadata
plotEmbedding(proj_3, name = "Plate", plotAs = 'points', size = 0.5) + theme_void() + ggtitle("Plate") + theme(legend.position = "bottom")
plotEmbedding(proj_3, name = "Type", plotAs = 'points', size = 0.5) + theme_void() + ggtitle("Type") + theme(legend.position = "bottom")
plot_grid(plotEmbedding(proj_3, name = "Target", plotAs = 'points', size = 0.5, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle("Target") + theme(legend.position = "bottom") + scale_color_manual(values = c("lightgrey","darkgrey",colors[c(9,3,6)])),
          plotEmbedding(proj_3, name = "Drug", plotAs = 'points', size = 0.5, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle("Drug") + theme(legend.position = "bottom") + scale_color_manual(values = colors_ordered),
          plotEmbedding(proj_3, name = "Dose", plotAs = 'points', size = 0.5, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle("Dose") + theme(legend.position = "bottom") + scale_color_manual(values = c("#802A07","#FFBB80","#E96200","lightgrey","darkgrey")),
          nrow = 1)


### ATAC Peak Calling ----------------------------------------------------------

pathToMacs2 <- "/Users/dannyconrad/opt/miniconda3/bin/macs2"

proj_3 <- addGroupCoverages(proj_3, groupBy = "Clusters", force = T)
proj_3 <- addReproduciblePeakSet(
  ArchRProj = proj_3, 
  groupBy = "Clusters", 
  pathToMacs2 = pathToMacs2, 
  force = T)
proj_3 <- addPeakMatrix(proj_3, force = T)

## Plot Metadata
g <- plot_grid(plotlist = lapply(c("log10(nFrags)", "TSSEnrichment", "FRIP", "log10(Gex_nUMI)"), function(x) {
  plotEmbedding(proj_3, name = x, embedding = "UMAP", plotAs = 'points', size = 0.75) + labs(x = "UMAP1", y = "UMAP2") + theme_DC + theme_UMAP + ggtitle(x) 
}), ncol = 2) 

save_plot(g, "./Plots/proj3_QC.pdf", w = 4, h = 4, show.legend = "none")

g <- plot_grid(plotlist = lapply(c("Clusters", "log10(nFrags)", "TSSEnrichment", "FRIP", "log10(Gex_nUMI)", "log10(tagUMI)"), function(x) {
  g <- plotEmbedding(proj_3, name = x, embedding = "UMAP", plotAs = 'points', size = 0.75, labelAsFactors=F, labelMeans=T) + labs(x = "UMAP1", y = "UMAP2") +
    theme_DC + theme_UMAP + ggtitle(x)
  if (x == "Clusters") { g <- g + theme(legend.position = 'none') }
  g
}), ncol = 3) 
save_plot(g, "./Plots/proj3_QC.pdf", w = 6, h = 4, show.legend = "none")

g <- plotEmbedding(proj_3, name = "Drug", embedding = "UMAP", plotAs = 'points', size = 0.75, labelAsFactors=F, labelMeans=F) + labs(x = "UMAP1", y = "UMAP2") +
  theme_DC + theme_UMAP + ggtitle("Drug") + scale_color_manual(values = col_drug)
save_plot(g, "./Plots/proj3_Drug.pdf", w = 2, h = 2, show.legend = "none")
g <- plotEmbedding(proj_3, name = "Plate", embedding = "UMAP", plotAs = 'points', size = 0.75, labelAsFactors=F, labelMeans=F) + labs(x = "UMAP1", y = "UMAP2") +
  theme_DC + theme_UMAP + ggtitle("Drug")
save_plot(g, "./Plots/proj3_Plate.pdf", w = 2, h = 2, show.legend = "none")
g <- plotEmbedding(proj_3, name = "PlateSide", embedding = "UMAP", plotAs = 'points', size = 0.75, labelAsFactors=F, labelMeans=F) + labs(x = "UMAP1", y = "UMAP2") +
  theme_DC + theme_UMAP + ggtitle("Drug") + scale_color_manual(values = col_side)
save_plot(g, "./Plots/proj3_PlateSide.pdf", w = 2, h = 2, show.legend = "none")

## FRIP
plotGroups(proj_3, groupBy = "Clusters", name = "FRIP", plotAs = 'vioin')
ggplot(as.data.frame(proj_3@cellColData), aes(x = Clusters, fill = FRIP <= 0.4)) + geom_bar(position = 'fill') + theme_minimal() + scale_fill_manual(values = c('grey','black'))
ggplot(as.data.frame(proj_3@cellColData), aes(x = Dose, fill = FRIP <= 0.4)) + geom_bar(position = 'fill') + theme_minimal() + scale_fill_manual(values = c('grey','black')) + facet_wrap(~Drug, scales = "free")
plotEmbedding(proj_3, name = 'FRIP <= 0.5', embedding = "UMAP", size = 1) + theme_void() + ggtitle("FRIP <= 0.5") + theme(legend.position = "none")               



### Cell Type Annotation ----------------------------------------------------

## Major Cell Type Annotation
plotGene(proj_3, "CD3G") # T
plotGene(proj_3, "BLK") # B
plotGene(proj_3, "FTL") # Myeloid

g <- plotGene(proj_3, c("CD3G", "BLK", "FTL"), "Major Cell Type Markers")
save_plot(g, "./Plots/proj3_MajorMarkers.pdf", w = 6, h = 4, show.legend = "none")

proj_3$CellType_Major <- "LowQuality"
proj_3$CellType_Major[proj_3$Clusters %in% c("C20")] <- "B"
proj_3$CellType_Major[proj_3$Clusters %in% paste("C",c(1,2,4,5,6,7),sep="")] <- "Myeloid"
proj_3$CellType_Major[proj_3$Clusters %in% paste("C",c(11:19,24),sep="")] <- "T"

table(proj_3$CellType_Major)

g <- plotEmbedding(proj_3, name = "CellType_Major", size = 0.75, labelAsFactors=F, labelMeans=T) + theme_DC + theme_UMAP + labs(x = "UMAP1", y = "UMAP2", title = "Major CellTypes")
save_plot(g, "./Plots/proj3_CellTypes.pdf", w = 2, h = 2, show.legend = "none")

## Compare Data Quality by Cell Type + Drug
ggplot(as.data.frame(proj_3@cellColData), aes(x = Drug, y = FRIP, fill = Drug)) + 
  geom_boxplot(outliers = F) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  facet_wrap(~CellType_Major, nrow = 1) + ggtitle("FRIP") + scale_fill_manual(values = col_drug)
ggplot(as.data.frame(proj_3@cellColData), aes(x = Drug, y = log10(nFrags), fill = Drug)) + 
  geom_boxplot(outliers = F) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  facet_wrap(~CellType_Major, nrow = 1) + ggtitle("log10(nFrags)") + scale_fill_manual(values = col_drug)
ggplot(as.data.frame(proj_3@cellColData), aes(x = Drug, y = log10(Gex_nUMI), fill = Drug)) + 
  geom_boxplot(outliers = F) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  facet_wrap(~CellType_Major, nrow = 1) + ggtitle("log10(nUMI)") + scale_fill_manual(values = col_drug)



### Remove Low Quality Clusters ---------------------------------------------

cells <- proj_3$cellNames[proj_3$CellType_Major != 'LowQuality']
plotEmbedding(proj_3, name = 'proj_3$cellNames %in% cells', size = 1) + theme_void() + ggtitle("Keep Cells") + theme(legend.position = "none") + scale_color_manual(values = c("black","dodgerblue"))        

g <- ggplot(as.data.frame(proj_3@cellColData), aes(x = log10(nFrags), y = log10(Gex_nUMI), color = proj_3$cellNames %in% cells)) + 
  ggrastr::rasterise(geom_point(size = 1, alpha = 0.3, stroke = 0), dpi = 600, scale = 1) +
  labs(color = "Keep Cells", title = "Good Quality Cells") + theme_DC + 
  scale_color_manual(values = c("black","dodgerblue")) + facet_grid(~Sample)
save_plot(g, "./Plots/proj3_KeepCells.pdf", w = 5, h = 2.5, show.legend = 'none')

proj_4 <- subsetArchRProject(proj_3,
                             cells = cells,
                             outputDirectory = "proj_4",
                             dropCells = T, force = T)
# proj_4 <- readRDS("proj_4/Save-ArchR-Project.rds")

proj_4 <- dimReduce(proj_4, scale = "Col", default = "NoScale", dims = c(2:30))


# Quality metrics different between cells from replicates on each plate side, and these plate side differences are captured in LSI4
g <- df %>% select(CellType_Major, PlateSide, Plate, nFrags, FRIP, Gex_nUMI, TSSEnrichment) %>%
  pivot_longer(cols = 4:7, names_to = "Metric", values_to = "Value") %>%
  group_by(CellType_Major, PlateSide, Plate, Metric) %>% summarize(Value = mean(Value)) %>%
  ggplot(aes(x = PlateSide, y = Value)) +
  geom_point(mapping = aes(color = PlateSide, shape = CellType_Major)) + 
  geom_line(mapping = aes(group=paste(CellType_Major, Plate), linetype = Plate), color = 'black') +
  stat_compare_means(method = 't.test', paired = T, size = 2) +
  facet_wrap(~Metric, scales = 'free', nrow = 1) + 
  scale_color_manual(values = col_side) + labs(shape = "Cell Type", color = "Plate Side", title = "QC Metrics by Plate Side", y = "Mean") +
  theme_DC
save_plot(g, "./Plots/PlateSide_QC.pdf", w = 5, h = 2, show.legend = 'right')

g <- data.frame(df, proj_4@reducedDims$IterativeLSI$matSVD) %>%
  ggplot(aes(x = log10(nFrags), y = LSI4, color = PlateSide)) + 
  geom_point() + 
  facet_grid(~CellType_Major) +
  scale_color_manual(values = col_side) +
  theme_DC 
g <- data.frame(df, proj_4@reducedDims$IterativeLSI$matSVD) %>%
  ggplot(aes(x = log10(nFrags), y = LSI4, color = PlateSide)) + 
  ggrastr::rasterise(geom_point(size = 0.75, color = 'grey80', alpha = 0.75, stroke = 0), dpi = 600, scale = 1) +
  geom_density2d(size = 0.5, bins = 8) +
  facet_grid(~CellType_Major) +
  scale_color_manual(values = col_side) + labs(color = "Plate Side", title = "LSI4 Captures Experimental Plate Side") +
  theme_DC 
save_plot(g, "./Plots/PlateSide_LSI4.pdf", w = 5, h = 2, show.legend = 'right')

## thus, exclude LSI4 from downstream embeddings (only affects visualization and minor cell type annotation)

proj_4 <- dimReduce(proj_4, scale = "Col", default = "NoScale", dims = c(2:3,5:30))
proj_4 <- addImputeWeights(proj_4, reducedDims = "IterativeLSI", dimsToUse = c(2:3,5:30), corCutOff = 1, scaleDims = F)


lapply(1:length(unique(proj_4$Clusters)), function(x) {
  cbind(proj_4@embeddings$UMAP$df, proj_4@cellColData) %>% as.data.frame() %>% 
    ggplot(aes(x = .[,1], y = .[,2], color = Clusters == paste("C",x,sep=""), alpha = Clusters == paste("C",x,sep=""))) + 
    geom_point(size = 0.5, stroke = 0) + 
    theme_DC + theme_UMAP + theme(legend.position = 'none') + labs(x = "UMAP1", y = "UMAP2", title = x) +
    scale_alpha_manual(values = c(0.2,0.8)) + scale_color_manual(values = c("grey","purple"))
}) %>% plot_grid(plotlist = .)

plotEmbedding(proj_4, name = "nFrags >= 1000", plotAs = 'points', size = 0.75, labelAsFactors=F, labelMeans=T) + labs(x = "UMAP1", y = "UMAP2") + theme_DC + theme_UMAP
plotEmbedding(proj_4, name = "Clusters != 'C3'", plotAs = 'points', size = 0.75, labelAsFactors=F, labelMeans=T) + labs(x = "UMAP1", y = "UMAP2") + theme_DC + theme_UMAP

cells <- proj_4$cellNames[proj_4$nFrags >= 1000 & proj_4$Clusters != "C3"]

proj_4 <- proj_4[cells]

proj_4 <- dimReduce(proj_4, scale = "Col", default = "NoScale", dims = c(2:3,5:30))
proj_4 <- addImputeWeights(proj_4, reducedDims = "IterativeLSI", dimsToUse = c(2:3,5:30), corCutOff = 1, scaleDims = F)


g <- plot_grid(plotlist = lapply(c("Clusters", "log10(nFrags)", "TSSEnrichment", "FRIP", "log10(Gex_nUMI)", "log10(tagUMI)"), function(x) {
  g <- plotEmbedding(proj_4, name = x, embedding = "UMAP", plotAs = 'points', size = 0.75, labelAsFactors=F, labelMeans=T) + labs(x = "UMAP1", y = "UMAP2") +
    theme_DC + theme_UMAP + ggtitle(x)
  if (x == "Clusters") { g <- g + theme(legend.position = 'none') }
  g
}), ncol = 3) 
save_plot(g, "./Plots/proj4_QC.pdf", w = 6, h = 4, show.legend = "none")

g <- plot_grid(plotEmbedding(proj_4, name = "Clusters", labelAsFactors = F, labelMeans = T, size = 0.5) + labs(x = "UMAP1", y = "UMAP2", title = "Clusters") + theme_DC + theme_UMAP + theme(legend.position = 'none'),
               plotEmbedding(proj_4, name = "CellType_Major", labelAsFactors = F, labelMeans = T, size = 0.5) + labs(x = "UMAP1", y = "UMAP2", title = "Major Cell Types") + theme_DC + theme_UMAP + theme(legend.position = 'none') + scale_color_manual(values = col_cellmaj),
               plotEmbedding(proj_4, name = "Drug", labelAsFactors = F, labelMeans = F, size = 0.5) + labs(x = "UMAP1", y = "UMAP2", title = "Drug") + theme_DC + theme_UMAP + theme(legend.position = "right") + scale_color_manual(values = col_drug),
               plotEmbedding(proj_4, name = "PlateSide", labelAsFactors = F, labelMeans = F, size = 0.5) + labs(x = "UMAP1", y = "UMAP2", title = "Plate Side") + theme_DC + theme_UMAP + theme(legend.position = "right") + scale_color_manual(values = col_side),
               nrow = 2)
save_plot(g, "./Plots/proj4_Meta.pdf", w = 4, h = 4, show.legend = "none")


### Cell Subtype Annotation -------------------------------------------------

## Major Cell Types
g <- plotGene(proj_4, c("CD3G", "BLK", "FTL"), "Major Cell Type Markers")
save_plot(g, "./Plots/proj4_MajorMarkers.pdf", w = 6, h = 4, show.legend = "none")

g <- plotEmbedding(proj_4, name = "CellType_Major", size = 1, labelAsFactors=F, labelMeans=T) + labs(x = "UMAP1", y = "UMAP2", title = "Major Cell Types") + theme_DC + theme_UMAP + scale_color_manual(values = col_cellmaj)
save_plot(g, "./Plots/proj4_CellTypes.pdf", w = 2, h = 2, show.legend = "none")

## Assign Minor Cell Types
plotGene(proj_4, c("FOXP3", "IL2RA"), "T-Reg Cells") %>% save_plot("./Plots/proj4_MinorMarkers_TReg.pdf", w = 4, h = 4, show.legend = 'none')
plotGene(proj_4, c("NCR1", "GNLY"), "NK Cells") %>% save_plot("./Plots/proj4_MinorMarkers_NK.pdf", w = 4, h = 4, show.legend = 'none')
plotGene(proj_4, c("CD8A", "IFNG"), "CD8 T Cells") %>% save_plot("./Plots/proj4_MinorMarkers_CD8.pdf", w = 4, h = 4, show.legend = 'none')
plotGene(proj_4, c("CD4", "LEF1"), "CD4 T Cells") %>% save_plot("./Plots/proj4_MinorMarkers_CD4.pdf", w = 4, h = 4, show.legend = 'none')
plotGene(proj_4, c("CPQ","FTL"), "Monocytes") %>% save_plot("./Plots/proj4_MinorMarkers_Mono.pdf", w = 4, h = 4, show.legend = 'none')
plotGene(proj_4, c("FLT3","CD86"), "Dendritic Cells") %>% save_plot("./Plots/proj4_MinorMarkers_DC.pdf", w = 4, h = 4, show.legend = 'none')

proj_4 <- addClusters(proj_4, corCutOff = 1, scaleDims = F, dimsToUse = c(2:3,5:30), resolution = 1, force = T)

g <- plotEmbedding(proj_4, name = "Clusters", size = 0.5, labelAsFactors=F, labelMeans=T) + labs(x = "UMAP1", y = "UMAP2", title = "Clusters") + theme_DC + theme_UMAP
save_plot(g, "./Plots/proj4_Clusters.pdf", w = 3, h = 3, show.legend = "none")

proj_4$CellType_Minor <- proj_4$CellType_Major

proj_4$CellType_Minor[proj_4$Clusters %in% c("C18","C14")] <- "Myeloid SWI/SNF"
proj_4$CellType_Minor[proj_4$Clusters %in% c("C16")] <- "Myeloid p300/CBP"
proj_4$CellType_Minor[proj_4$Clusters %in% c("C8")] <- "T SWI/SNF"
proj_4$CellType_Minor[proj_4$Clusters %in% c("C1")] <- "T PRC2"

proj_4$CellType_Minor[proj_4$Clusters %in% c("C6")] <- "B"
proj_4$CellType_Minor[proj_4$Clusters %in% c("C9")] <- "NK"
proj_4$CellType_Minor[proj_4$Clusters %in% c("C10","C11","C7")] <- "Naive T"
proj_4$CellType_Minor[proj_4$Clusters %in% c("C5")] <- "CD8+ Tmem"
proj_4$CellType_Minor[proj_4$Clusters %in% c("C3","C4")] <- "CD4+ Tmem"
proj_4$CellType_Minor[proj_4$Clusters %in% c("C2")] <- "Treg"
proj_4$CellType_Minor[proj_4$Clusters %in% c("C12")] <- "Monocyte Stim"
proj_4$CellType_Minor[proj_4$Clusters %in% c("C13","C15")] <- "Monocyte Unstim"
proj_4$CellType_Minor[proj_4$Clusters %in% c("C17")] <- "DC"


plotEmbedding(proj_4, name = 'CellType_Major', size = 0.1, labelAsFactors = F, labelMeans = T) + 
  theme_void() + ggtitle("Major Cell Types") + 
  scale_color_manual(values = col_cellmaj)
plotEmbedding(proj_4, name = 'CellType_Minor', size = 0.1, labelAsFactors = F, labelMeans = T) + 
  theme_void() + ggtitle("Minor Cell Types") + 
  scale_color_manual(values = col_cellmin)


### Create Group BigWig Files -----------------------------------------------

proj_4$DrugDose <- proj_4$Drug
proj_4$DrugDose[proj_4$Type != "Control"] <- paste(proj_4$DrugDose, proj_4$Dose)[proj_4$Type != "Control"]

proj_4$CellDrugDose <- paste(proj_4$CellType_Major, proj_4$DrugDose)

proj_4 <- addGroupCoverages(proj_4, groupBy = "CellDrugDose")
# getGroupBW(proj_4, groupBy = "CellDrugDose", normMethod = "ReadsInTSS")
# getGroupBW(proj_4, groupBy = "CellDrugDose", normMethod = "nFrags")
getGroupBW(proj_4, groupBy = "CellDrugDose", normMethod = "nCells")
getGroupBW(proj_4, groupBy = "CellDrugDose", normMethod = "None")



### chromVAR Motif Deviations -----------------------------------------------

# Add Motif Peak Annotations
proj_4 <- addMotifAnnotations(ArchRProj = proj_4, motifSet = "cisbp", annoName = "Motif", force = T) ## cis-bp v1

proj_4 <- addBgdPeaks(proj_4)
proj_4 <- addDeviationsMatrix(
  ArchRProj = proj_4,
  peakAnnotation = "Motif", force = T
)
getVarDeviations(proj_4, plot = TRUE, name = "MotifMatrix")


### UMAP Embeddings ---------------------------------------------------------

df <- data.frame(UMAP1 = proj_4@embeddings$UMAP$df[,1], 
                 UMAP2 = proj_4@embeddings$UMAP$df[,2],
                 Cell = proj_4$cellNames,
                 proj_4@cellColData)

df$Clusters <- factor(df$Clusters, levels = paste("C",1:length(unique(df$Clusters)),sep=""))
df$Drug <- factor(df$Drug, levels = c("DMSO(-)","DMSO(+)","MS177","EPZ-6438","AU-15330","BRM014","dCBP-1","GNE-781"))
df$Target <- factor(df$Target, levels = c("DMSO(-)", "DMSO(+)", "PRC2", "SWI/SNF", "p300/CBP"))
df$PlateColumn <- factor(df$PlateColumn, levels = 1:12)
df$CellType_Major <- factor(df$CellType_Major, levels = c("T","Myeloid","B","Other"))

## Clusters
g <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Clusters)) + 
  ggrastr::rasterise(geom_point(size = 0.5, alpha = 0.75, stroke = 0), dpi = 600, scale = 1) +
  ggrepel::geom_label_repel(
    data = df %>% group_by(Clusters) %>% summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), 
    aes(label = Clusters),
    max.overlaps = 20, size = 2, alpha = 0.8, label.padding = 0.1) +
  theme_DC + theme_UMAP +
  scale_color_manual(values = paletteer::paletteer_c("grDevices::Dark 3", 24)) +
  ggtitle("Clusters") 
save_plot(g,  "./Plots/proj4_ClusterUMAP.pdf", w = 2, h = 2, show.legend = "none")

## Drugs + Doses
g <- plot_grid(plotlist = lapply(levels(df$Drug), function(x) {
  ggplot(df, aes(x = UMAP1, y = UMAP2)) + 
    ggrastr::rasterise(geom_point(color = 'grey90', size = 0.5, alpha = 0.3, stroke = 0), dpi = 600, scale = 1) +
    ggrastr::rasterise(geom_point(data = df[df$Drug == x,], mapping = aes(x = UMAP1, y = UMAP2, color = Dose), size = 0.75, alpha = 0.75, stroke = 0), dpi = 600, scale = 1) +
    theme_DC + theme_UMAP + theme(legend.position = 'none') + ggtitle(x) +
    scale_color_manual(values = col_dose)
}), nrow = 2, byrow = F)
save_plot(g, "./Plots/proj4_DrugDoseUMAPs.pdf", w = 8, h = 4, show.legend = "none")

## Drugs
g <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Drug)) + 
  ggrastr::rasterise(geom_point(size = 0.5, alpha = 0.75, stroke = 0), dpi = 600, scale = 1) +
  scale_color_manual(values = col_drug) +
  theme_DC + theme_UMAP + ggtitle("Drug") +
  guides(colour = guide_legend(override.aes = list(size=3, shape = 16)), linetype = 'none')
save_plot(g, "./Plots/proj4_DrugUMAP.pdf", w = 2.5, h = 2, show.legend = 'right')

## Major Cell Type
g <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = CellType_Major)) + 
  ggrastr::rasterise(geom_point(size = 0.5, alpha = 0.75, stroke = 0), dpi = 600, scale = 1) +
  geom_label(
    data = df %>% group_by(CellType_Major) %>% summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), 
    aes(label = CellType_Major),
    size = 3, alpha = 0.9) +
  theme_DC + theme_UMAP + ggtitle("Major Cell Types") +
  scale_color_manual(values = col_cellmaj)
save_plot(g, "./Plots/proj4_MajorCellTypeUMAP.pdf", w = 2, h = 2, show.legend = 'none')

## Minor Cell Type
g <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = CellType_Minor)) + 
  ggrastr::rasterise(geom_point(size = 0.5, alpha = 0.75, stroke = 0), dpi = 600, scale = 1) +
  ggrepel::geom_label_repel(
    data = df %>% group_by(CellType_Minor) %>% summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), 
    aes(label = CellType_Minor),
    max.overlaps = 20, size = 2, alpha = 0.8, label.padding = 0.1) +
  theme_DC + theme_UMAP + ggtitle("Minor Cell Types") +
  scale_color_manual(values = col_cellmin) 
save_plot(g, "./Plots/proj4_MinorCellTypeUMAP.pdf", w = 2, h = 2, show.legend = 'none')

### Basic Characterization --------------------------------------------------

## Plot Quality Metrics Across Conditions & Cell Types
df %>% dupControl %>% filter(Dose != "DMSO(-)") %>% group_by(CellType_Major, Drug, Dose, Replicate) %>% 
  summarize(TSSEnrichment = mean(TSSEnrichment)) %>% 
  ggplot(aes(x = Drug, y = TSSEnrichment, fill = Dose)) + 
  geom_boxplot(outliers = F) + facet_grid(CellType_Major~.) + scale_fill_manual(values = col_dose) + theme_DC +
  stat_compare_means(method = "anova", label = "p.signif", size = 2, vjust = 0)
df %>% dupControl %>% filter(Dose != "DMSO(-)") %>% group_by(CellType_Major, Drug, Dose, Replicate) %>% 
  summarize(nFrags = mean(nFrags)) %>% 
  ggplot(aes(x = Drug, y = log10(nFrags), fill = Dose)) + 
  geom_boxplot(outliers = F) + facet_grid(CellType_Major~.) + scale_fill_manual(values = col_dose) + theme_DC +
  stat_compare_means(method = "anova", label = "p.signif", size = 2, vjust = 0)
df %>% dupControl %>% filter(Dose != "DMSO(-)") %>% group_by(CellType_Major, Drug, Dose, Replicate) %>% 
  summarize(Gex_nUMI = mean(Gex_nUMI)) %>% 
  ggplot(aes(x = Drug, y = log10(Gex_nUMI), fill = Dose)) + 
  geom_boxplot(outliers = F) + facet_grid(CellType_Major~.) + scale_fill_manual(values = col_dose) + theme_DC +
  stat_compare_means(method = "anova", label = "p.signif", size = 2, vjust = 0)
df %>% dupControl %>% filter(Dose != "DMSO(-)") %>% group_by(CellType_Major, Drug, Dose, Replicate) %>% 
  summarize(FRIP = mean(FRIP)) %>% 
  ggplot(aes(x = Drug, y = FRIP, fill = Dose)) + 
  geom_boxplot(outliers = F) + facet_grid(CellType_Major~.) + scale_fill_manual(values = col_dose) + theme_DC +
  stat_compare_means(method = "anova", label = "p.signif", size = 2, vjust = 0)
df %>% dupControl %>% filter(Dose != "DMSO(-)") %>% group_by(CellType_Major, Drug, Dose, Replicate) %>% 
  summarize(Length = mean(Length)) %>% 
  ggplot(aes(x = Drug, y = Length, fill = Dose)) + 
  geom_boxplot(outliers = F) + facet_grid(CellType_Major~.) + scale_fill_manual(values = col_dose) + theme_DC +
  stat_compare_means(method = "anova", label = "p.signif", size = 2, vjust = 0)


## Cluster Composition

plotEmbedding(proj_4, name = "Clusters", plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = T) + theme_void() + ggtitle("Clusters") + theme(legend.position = "none")
# ggplot(df, aes(x = Clusters, fill = Target)) + geom_bar(position = 'fill') + theme_minimal() + scale_fill_manual(values = col_target)
g <- ggplot(df, aes(x = Clusters, fill = Drug)) + geom_bar(position = 'fill') + 
  theme_minimal() + scale_fill_manual(values = col_drug) + theme_DC

save_plot(g, "./Plots/proj4_ClusterDrugComposition.pdf", w = 5, h = 2.5, show.legend = "right")


## Overall Recovery Stats

table(proj_4$MULTI) %>% mean
table(proj_4$MULTI) %>% sd
table(proj_4$Drug) %>% mean
table(proj_4$Drug) %>% sd
table(proj_4$DrugDose[proj_4$Type != "Control"]) %>% mean
table(proj_4$DrugDose[proj_4$Type != "Control"]) %>% sd
table(proj_4$CellType_Major) / length(proj_4$CellType_Major) * 100

df$MULTI %>% table
df$DrugDose %>% table

g <- df %>% group_by(MULTI) %>% summarize(Count = n()) %>% 
  ggplot(aes(x = Count)) + geom_histogram(color = 'white', fill = 'purple', bins = 15) + 
  labs(x = "Nuclei Recovered per Replicate Well", y = "Count") + theme_DC
save_plot(g, "./Plots/RecoveryPerWell.pdf", w = 2, h = 2, show.legend = 'none')
g <- df %>% group_by(MULTI, Plate) %>% summarize(Count = n()) %>% 
  ggplot(aes(x = Count, fill = Plate)) + geom_histogram(color = 'white', bins = 15, position = 'identity', alpha = 0.5) + 
  labs(x = "Nuclei Recovered per Replicate Well", y = "Count") + theme_DC + scale_fill_manual(values = c("red","blue"))
save_plot(g, "./Plots/RecoveryPerWell_ByPlate.pdf", w = 2, h = 2, show.legend = c(0.85,0.85))
g <- df %>% group_by(Drug, MULTI) %>% summarize(Count = n()) %>% 
  ggplot(aes(x = Drug, y = Count, fill = Drug)) + geom_boxplot(alpha = 1) + scale_fill_manual(values = col_drug) + 
  labs(x = "Drug", y = "Nuclei Recovered per Replicate Well", fill = NULL) + 
  theme_DC + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_plot(g, "./Plots/RecoveryPerWell_ByDrug.pdf", w = 3, h = 2, show.legend = "right")

## Compare Recovery Between Drugs/Doses

df %>% group_by(Drug, Dose) %>% summarize(Count = n(), CountPerWell = n()/length(unique(Replicate))) %>% tidyr::complete(Dose, fill = list(Count = 0, CountPerWell = 0)) %>%
  ggplot(aes(x = Dose, y = Drug, fill = Count, label = Count)) + 
  geom_tile() + geom_text() + 
  ggtitle("Total Cells Recovered per Condition") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  paletteer::scale_fill_paletteer_c("ggthemes::Red-Gold")
df %>% group_by(Drug, Dose) %>% summarize(Count = n(), CountPerWell = (n()/length(unique(Replicate))) %>% round(digits =)) %>% tidyr::complete(Dose, fill = list(Count = 0, CountPerWell = 0)) %>%
  ggplot(aes(x = Dose, y = Drug, fill = CountPerWell, label = CountPerWell)) + 
  geom_tile() + geom_text() + 
  ggtitle("Mean Cells Recovered per Replicate") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  paletteer::scale_fill_paletteer_c("ggthemes::Red-Gold")


g <- plot_grid(df %>% mutate(PlateRow = str_sub(MULTI, 1, 1)) %>% group_by(Plate, PlateRow, PlateColumn, Drug) %>% summarize(Count = n()) %>%
                 ggplot(aes(x = PlateColumn, y = PlateRow, fill = Count, label = Count)) +
                 geom_tile() + geom_text() +
                 paletteer::scale_fill_paletteer_c("grDevices::Red-Purple", direction = -1) +
                 ggnewscale::new_scale_fill() +
                 geom_label(data = df %>% mutate(PlateRow = str_sub(MULTI, 1, 1)) %>% group_by(Plate, PlateRow, Drug) %>% summarize(n()) ,
                            mapping = aes(y = PlateRow, label = Drug, fill = Drug), alpha = 0.5, x = 6.5, inherit.aes = F) +
                 scale_fill_manual(values = col_drug) +
                 facet_wrap(~Plate, nrow = 1) + scale_y_discrete(limits = rev) + theme_DC + theme(legend.position = 'none') + scale_x_discrete(drop = FALSE),
               ggplot(df %>% group_by(Plate, PlateColumn) %>% summarize(Count = n(), DMSO = sum(Type == "Control")), aes(x = PlateColumn, y = Count)) +
                 geom_col(fill = 'white', color = 'black') + geom_col(mapping=aes(y=DMSO), fill = 'grey60', width = 0.7, color = 'black') +
                 geom_point(size = 2, mapping = aes(y = Count / DMSO * 100)) +
                 facet_wrap(~Plate, nrow = 1) + theme_DC + scale_x_discrete(drop = FALSE) + labs(y = "Total Singlets"),
               nrow = 2)
save_plot(g, "./Plots/proj4_WellColumnRecovery.pdf", w = 8, h = 6, )


# Break down by minor cell type
ggplot(df, aes(x = Dose, fill = CellType_Minor)) +
  geom_bar(position = 'fill') +
  facet_wrap(~Drug, nrow = 2, dir = "v", scales = 'free_x') + theme_bw() +
  theme(axis.text.x = element_text(angle = -25, hjust = 0)) +
  scale_fill_manual(values = col_cellmin)

ggplot(df[df$CellType_Major %ni% c("B","Other"),], aes(x = Dose, fill = CellType_Minor)) +
  geom_bar(position = 'fill') +
  facet_grid(CellType_Major~Drug, scales = 'free_x', space = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(angle = -25, hjust = 0)) +
  scale_fill_manual(values = col_cellmin)

df %>% 
  dupControl %>% filter(Drug %ni% c("DMSO(-)","DMSO(+)")) %>%
  filter(CellType_Minor %ni% grep("Drug|PRC2|p300|SWI", CellType_Minor, value = T)) %>%
  group_by(CellType_Minor, Drug, Dose, ColScaleFactor) %>% summarize(Count = n()) %>% mutate(Scaled = Count * ColScaleFactor) %>% summarize(Count = mean(Count), Scaled = mean(Scaled)) %>%
  # ggplot(aes(x = CellType_Minor, y = Count, fill = Dose)) +
  ggplot(aes(x = CellType_Minor, y = Scaled, fill = Dose)) +
  geom_col(position = 'dodge') +
  facet_wrap(~Drug, nrow = 2, dir = "v") +
  theme_DC + theme(axis.text.x = element_text(angle=-25, hjust=0)) +
  scale_fill_manual(values = col_dose)


## Relative Frequencies

g <- df %>% dupControl %>%
  group_by(Drug, Dose, Replicate, CellType_Major) %>% summarize(Count = n()) %>% mutate(Freq = Count/sum(Count)) %>% mutate(Drug = factor(Drug, levels = names(col_drug))) %>%
  {
    ggplot(., aes(x = Dose, y = Freq, color = CellType_Major)) +
      geom_boxplot(position = "identity", outliers = F, linewidth = 0.3, alpha = 0) +
      geom_line(data = . %>% group_by(Drug, Dose, CellType_Major) %>% summarize(Freq = mean(Freq), Count = mean(Count)), mapping = aes(group = CellType_Major), linewidth = 0.5, linetype = "11") +
      facet_wrap(~Drug, scales = 'free_x', dir = "v", nrow = 2) + 
      scale_color_manual(values = col_cellmaj) + ylim(c(0,0.85)) +
      scale_x_discrete(labels = c("Neg","Pos","10nM","100nM","1µM")) +
      theme_DC + labs(y = "Frequency", color = "Cell Type") + theme(axis.text.x = element_text(colour = col_dose), legend.position = 'none')
  }
save_plot(g, "./Plots/proj4_CellTypeFrequencies.pdf", w = 5, h = 4, show.legend = "none")

g <- lapply(drug[2:8], function(d) { 
  df %>% filter(Drug %in% c("DMSO(-)","DMSO(+)", d)) %>% mutate(Drug = case_when(d == "DMSO(+)" ~ "Activation", TRUE ~ d)) 
}) %>% do.call(rbind,.) %>%
  filter(!(Dose == "DMSO(-)" & Drug != "Activation")) %>%
  group_by(Drug, Dose, Replicate, CellType_Major) %>% summarize(Count = n()) %>% mutate(Freq = Count/sum(Count)) %>%
  group_by(Drug, Dose, CellType_Major) %>% summarize(Freq = mean(Freq), Count = mean(Count)) %>% mutate(Drug = factor(Drug, levels = names(col_drug))) %>%
  ggplot(aes(x = Dose, y = Freq, group = paste(CellType_Major, Drug), linetype = CellType_Major, color = Drug)) +
  geom_line(linewidth = 0.75) + #facet_wrap(~Drug) +
  scale_color_manual(values = col_drug) + scale_linetype_manual(values = c("solid","22","11")) + 
  ylim(c(0,0.8)) + scale_x_discrete(labels = c("Neg","Pos","10nM","100nM","1µM")) +
  theme_DC + labs(y = "Frequency", linetype = "Cell Type", color = "Treatment") + theme(axis.text.x = element_text(colour = col_dose))
save_plot(g, "./Plots/proj4_CellTypeFrequencies_v2.pdf", w = 3, h = 2.5, show.legend = "right")


### Signac/Seurat -----------------------------------------------------------

dir.create("../Seurat/")

# Load Cell Ranger ARC Metadata
md.1 <- read.table(
  file = "/Volumes/DannySSD/MULTI_ATAC/20231114_PROTACExperiment2/outs1/per_barcode_metrics.csv",
  stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1
)
md.2 <- read.table(
  file = "/Volumes/DannySSD/MULTI_ATAC/20231114_PROTACExperiment2/outs2/per_barcode_metrics.csv",
  stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1
)
md.4 <- read.table(
  file = "/Volumes/DannySSD/MULTI_ATAC/20231114_PROTACExperiment2/outs4/per_barcode_metrics.csv",
  stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1
)

md.1 <- md.1[proj_4$CB_CellRangerGEX[proj_4$Sample == "L1"],]
md.2 <- md.2[proj_4$CB_CellRangerGEX[proj_4$Sample == "L2"],]
md.4 <- md.4[proj_4$CB_CellRangerGEX[proj_4$Sample == "L4"],]

# Create Fragment Objects
frags.1 <- CreateFragmentObject(
  path = "/Volumes/DannySSD/MULTI_ATAC/20231114_PROTACExperiment2/outs1/atac_fragments.tsv.gz",
  cells = rownames(md.1)
)
frags.2 <- CreateFragmentObject(
  path = "/Volumes/DannySSD/MULTI_ATAC/20231114_PROTACExperiment2/outs2/atac_fragments.tsv.gz",
  cells = rownames(md.2)
)
frags.4 <- CreateFragmentObject(
  path = "/Volumes/DannySSD/MULTI_ATAC/20231114_PROTACExperiment2/outs4/atac_fragments.tsv.gz",
  cells = rownames(md.4)
)

# Quantify peaks in each dataset
L1_atac <- FeatureMatrix(
  fragments = frags.1,
  features = proj_4@peakSet,
  cells = rownames(md.1)
)
L2_atac <- FeatureMatrix(
  fragments = frags.2,
  features = proj_4@peakSet,
  cells = rownames(md.2)
)
L4_atac <- FeatureMatrix(
  fragments = frags.4,
  features = proj_4@peakSet,
  cells = rownames(md.4)
)

# Create Assays & Seurat Objects
L1_atac <- CreateChromatinAssay(L1_atac, fragments = frags.1)
L1 <- CreateSeuratObject(L1_atac, assay = "ATAC", meta.data=md.1)
L1_rna <- Read10X_h5("../outs1/raw_feature_bc_matrix.h5")$`Gene Expression`
L1_rna <- L1_rna[, match(df$CB_CellRangerGEX[df$Sample == "L1"], colnames(L1_rna))] %>% CreateAssayObject()
L1[["RNA"]] <- L1_rna

L2_atac <- CreateChromatinAssay(L2_atac, fragments = frags.2)
L2 <- CreateSeuratObject(L2_atac, assay = "ATAC", meta.data=md.2)
L2_rna <- Read10X_h5("../outs2/raw_feature_bc_matrix.h5")$`Gene Expression`
L2_rna <- L2_rna[, match(df$CB_CellRangerGEX[df$Sample == "L2"], colnames(L2_rna))] %>% CreateAssayObject()
L2[["RNA"]] <- L2_rna

L4_atac <- CreateChromatinAssay(L4_atac, fragments = frags.4)
L4 <- CreateSeuratObject(L4_atac, assay = "ATAC", meta.data=md.4)
L4_rna <- Read10X_h5("../outs4/raw_feature_bc_matrix.h5")$`Gene Expression`
L4_rna <- L4_rna[, match(df$CB_CellRangerGEX[df$Sample == "L4"], colnames(L4_rna))] %>% CreateAssayObject()
L4[["RNA"]] <- L4_rna

L1$dataset <- 'L1'
L2$dataset <- 'L2'
L4$dataset <- 'L4'

saveRDS(L1, "../Seurat/L1.rds")
saveRDS(L2, "../Seurat/L2.rds")
saveRDS(L4, "../Seurat/L4.rds")

# Merge All Datasets
seu <- merge(
  x = L1,
  y = list(L2,L4),
  add.cell.ids = c("L1", "L2", "L4")
)
saveRDS(seu, "../Seurat/seu.rds")

# Dimensionality Reduction & UMAP Embedding on ATAC Assay
seu <- RunTFIDF(seu)
seu <- FindTopFeatures(seu, min.cutoff = 20)
seu <- RunSVD(seu)
seu <- RunUMAP(seu, dims = 2:30, reduction = 'lsi', reduction.name = "umap.atac")
DimPlot(seu, group.by = 'dataset', reduction = "umap.atac", pt.size = 0.1)

# Add Metadata from ArchR Project
seu@meta.data <- cbind(seu@meta.data, proj_4@cellColData[gsub("_","#",colnames(seu)),colnames(proj_4@cellColData) %ni% colnames(seu@meta.data)])

data.frame(seu@reductions$umap.atac@cell.embeddings, seu@meta.data) %>% ggplot(aes(.[,1],.[,2],color = log10(nFrags))) + geom_point(size = 0.5, alpha = 0.5) + scale_color_viridis_c(option = "B") + theme_classic()
data.frame(seu@reductions$umap.atac@cell.embeddings, seu@meta.data) %>% ggplot(aes(.[,1],.[,2],color = FRIP)) + geom_point(size = 0.5, alpha = 0.5) + scale_color_viridis_c(option = "B") + theme_classic()
data.frame(seu@reductions$umap.atac@cell.embeddings, seu@meta.data) %>% ggplot(aes(.[,1],.[,2],color = log10(Gex_nUMI))) + geom_point(size = 0.5, alpha = 0.5) + scale_color_viridis_c(option = "B") + theme_classic()
data.frame(seu@reductions$umap.atac@cell.embeddings, seu@meta.data) %>% ggplot(aes(.[,1],.[,2],color = Drug)) + geom_point(size = 1) + scale_color_manual(values = c("lightgrey","darkgrey",colors[c(2,5,8,3,6,9)])) + theme_classic()

# Get EnsDb Gene Annotations & Change to UCSC Style
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(seu) <- annotations
seu[['GeneScore']] <- CreateAssayObject(counts = GeneActivity(seu))

DefaultAssay(seu) <- 'GeneScore'
seu <- NormalizeData(
  object = seu,
  scale.factor = median(seu$nCount_GeneScore)
)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, reduction.name = "pca.score")
seu <- FindNeighbors(seu, reduction = "pca.score", dims = 1:30)
seu <- FindClusters(seu, reduction = "pca.score", resolution = 0.8, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca.score", dims = 1:30, reduction.name = "umap.score")
DimPlot(seu, group.by = 'dataset', reduction = "umap.score", pt.size = 0.1)

DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, reduction.name = "pca.rna")
seu <- FindNeighbors(seu, reduction = "pca.rna", dims = 1:30)
seu <- FindClusters(seu, reduction = "pca.rna", resolution = 0.8, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca.rna", dims = 1:30, reduction.name = "umap.rna")
DimPlot(seu, group.by = 'dataset', reduction = "umap.rna", pt.size = 0.1)

seu <- SCTransform(seu)
DefaultAssay(seu) <- 'ATAC'

saveRDS(seu, "../Seurat/seu.rds")

rm(L1,L2,L4,
   L1_atac, L2_atac, L4_atac,
   L1_rna, L2_rna, L4_rna,
   frags.1, frags.2, frags.4,
   md.1, md.2, md.4,
   annotations)



### ArchR Markers - CellDrugDose --------------------------------------------

proj_4$CellDrugDose2 <- proj_4@cellColData %>% as.data.frame %>% 
  mutate(CellDrugDose2 = case_when(
    Drug == "DMSO(+)" ~ paste(CellType_Major, "PosCtrl", sep=" "),
    Drug == "DMSO(-)" ~ paste(CellType_Major, "NegCtrl", sep=" "),
    TRUE ~ paste(CellType_Major, Drug, Dose, sep=" "))) %>% 
  .$CellDrugDose2

## Marker Computation by Drug Dose & Cell Type
dir.create("./MarkerSets/CellDrugDose2")

# bias <- c("TSSEnrichment", "log10(nFrags)")
bias <- c("FRIP", "log10(Gex_nUMI)")

## Peaks
markerPeaks_DrugDoseT <- getMarkerFeatures(
  proj_4,
  groupBy = "CellDrugDose2",
  useMatrix = "PeakMatrix",
  testMethod = "wilcoxon",
  maxCells = 1000,
  bias = bias,
  useGroups = grep("^B|^M", unique(proj_4$CellDrugDose2), value = T, invert = T),
  bgdGroups = grep("T PosCtrl", unique(proj_4$CellDrugDose2), value = T)
)
saveRDS(markerPeaks_DrugDoseT, "./MarkerSets/CellDrugDose2/markerPeaksT_vPos.rds")
markerPeaks_DrugDoseM <- getMarkerFeatures(
  proj_4,
  groupBy = "CellDrugDose2",
  useMatrix = "PeakMatrix",
  testMethod = "wilcoxon",
  maxCells = 1000,
  bias = bias,
  useGroups = grep("^B|^T", unique(proj_4$CellDrugDose2), value = T, invert = T),
  bgdGroups = grep("Myeloid PosCtrl", unique(proj_4$CellDrugDose2), value = T)
)
saveRDS(markerPeaks_DrugDoseM, "./MarkerSets/CellDrugDose2/markerPeaksM_vPos.rds")

rm(list = ls(pattern = "marker.*DrugDose."))

# Use comparisons against the Positive Control
markerPeaks_DrugDoseT <- readRDS("./MarkerSets/CellDrugDose2/markerPeaksT_vPos.rds")
markerPeaks_DrugDoseM <- readRDS("./MarkerSets/CellDrugDose2/markerPeaksM_vPos.rds")


### Dose-Dependent Markers -------------------------------------------------

## Perform linear regression between mean values per replicate and the drug dose

dir.create("./RepTest")

lm_atac <- repTest("ATAC", bins = 200, repsumm = "median", test = 'lm', log = F)
lm_score <- repTest("Score", bins = 50, repsumm = "median", test = 'lm', log = F)
lm_rna <- repTest("RNA", bins = 50, repsumm = "median", test = 'lm', log = F)
lm_rna2 <- repTest("RNA", bins = 50, repsumm = "median", test = 'lm', log = T)
lm_sct <- repTest("SCT", bins = 50, repsumm = "median", test = 'lm', log = F)
lm_sct2 <- repTest("SCT", bins = 50, repsumm = "median", test = 'lm', log = T)
lm_dev <- repTest("Dev1", bins = 3, repsumm = "median", test = 'lm', log = F)
lm_dev2 <- repTest("Dev2", bins = 3, repsumm = "median", test = 'lm', log = F)

saveRDS(lm_atac, "./RepTest/lm_atac.rds")
saveRDS(lm_score, "./RepTest/lm_score.rds")
saveRDS(lm_rna, "./RepTest/lm_rna.rds")
saveRDS(lm_rna2, "./RepTest/lm_rna2.rds")
saveRDS(lm_sct, "./RepTest/lm_sct.rds")
saveRDS(lm_sct2, "./RepTest/lm_sct2.rds")
saveRDS(lm_dev, "./RepTest/lm_dev.rds")
saveRDS(lm_dev2, "./RepTest/lm_dev2.rds")

markers_lm <- 
  lapply(c("Myeloid","T"), function(c) {
    df <- lapply(levels(df$Drug)[c(1,3:8)], function(d) { 
      peaks <- data.frame(Gene = colnames(lm_atac$p[,-1:-2]),
                          P = lm_atac$p %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                          Coef = lm_atac$coef %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                          Log2FC = lm_atac$log2fc.Max %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                          Mean = lm_atac$summ.Max %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                          MarkerType = "Peaks")
      genes <- data.frame(Gene = colnames(lm_rna$p[,-1:-2]) %>% fixGene(new),
                          P = lm_rna$p %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                          Coef = lm_rna$coef %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                          Log2FC = lm_rna$log2fc.Max %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                          Mean = lm_rna$summ.Max %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                          MarkerType = "Genes")
      scores <- data.frame(Gene = colnames(lm_score$p[,-1:-2]) %>% fixGene(new),
                           P = lm_score$p %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                           Coef = lm_score$coef %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                           Log2FC = lm_score$log2fc.Max %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                           Mean = lm_score$summ.Max %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                           MarkerType = "Scores")
      tmp <- rbind(peaks,genes,scores)
      tmp$Drug <- d
      tmp
    }) %>% do.call(rbind,.)
    df$CellType <- c
    df
  }) %>% do.call(rbind, .)


markers_lm_sig <-
  markers_lm %>% 
  filter(P < 0.01, Log2FC > 1) %>%
  mutate(FDR = P,
         Log2FC = Log2FC * sign(Coef),
         Treatment = factor(Drug, levels = names(col_drug)[c(2,4:9)], labels = names(col_drug)[c(1,4:9)]),
         Target = factor(Treatment, levels = c("Activation","MS177","EPZ6438","AU-15330","BRM014","dCBP-1","GNE-781"), labels = c("Activation","PRC2","PRC2","SWI/SNF","SWI/SNF","p300/CBP","p300/CBP")),
         Type = factor(sign(Coef), levels = c(-1,1), labels = c("Down","Up")),
         name = Gene) %>%
  select(-Gene, -Drug)

g <- markers_lm_sig %>%
  mutate(CellType = factor(CellType,levels = c("T","Myeloid")),
         MarkerType = factor(MarkerType,levels = c("Peaks","Scores","Genes")),
         Type = factor(Type, levels = c("Up","Down"))) %>%
  ggplot(aes(y = Treatment, fill = Type)) + geom_bar(position = position_dodge2(reverse = TRUE, padding = 0), width = 0.7) + 
  facet_grid(CellType~MarkerType, scales = 'free_x', axes = 'all_y') + scale_y_discrete(drop=FALSE, limits = rev) + 
  paletteer::scale_fill_paletteer_d("ggsci::alternating_igv", direction = -1) + theme_DC + theme(axis.text.y = element_text(color = rev(col_drug[-2:-3]))) + 
  labs(x = "Marker Count (p < 0.01 & Log2FC > 1)")
save_plot(g, "./Plots/Markers_All_Count.pdf", w = 5.5, h = 2, show.legend = 'right')


## Examples
# data.frame(Peak = resT$hclust$labels, 
#            Group = resT$cut,
#            Var = colVars(resT$heatmap)) %>% group_by(Group) %>% slice_max(order_by = Var, n = 1) %>% .$Peak

x <- "chr3-4978588-4979088" 

g <- lm_atac$rep[,c("CellType", "Drug", "Dose", "Replicate", x)] %>% {colnames(.)[5] <- "Feature";.} %>% 
  dupControl(controls = "DMSO(+)", drugs = drug[c(1,3:8)]) %>% filter(CellType == "T", Drug %in% c("DMSO(-)","MS177","AU-15330","dCBP-1")) %>% 
  mutate(Group = factor(Drug, levels = c("DMSO(-)","MS177","AU-15330","dCBP-1"), labels = c("Activation","Drug","Drug","Drug"))) %>%
  mutate(Drug = factor(Drug, levels = c("DMSO(-)","MS177","AU-15330","dCBP-1"), labels = c("Activation","MS177","AU-15330","dCBP-1"))) %>%
  {
    ggplot(., aes(x = as.numeric(Dose), y = Feature)) +
      geom_jitter(data = . %>% filter(Group == "Activation"), mapping = aes(color = Dose), width = 0.2, size = 1, show.legend = F) + 
      scale_color_manual(values = c(col_dose[1:2])) + guides(color = 'none') +
      ggnewscale::new_scale_color() +
      geom_jitter(data = . %>% filter(Group == "Drug" & Dose != "DMSO(+)"), mapping = aes(color = Drug), width = 0.2, size = 1, show.legend = F) + 
      geom_smooth(mapping = aes(color = Drug, fill = Drug), method = "lm", alpha = 0.2, show.legend = T) + 
      scale_color_manual(values = col_drug) + scale_fill_manual(values = col_drug) + 
      theme_DC + scale_x_continuous(breaks = 1:5, labels = levels(df$Dose)) + theme(axis.text.x = element_text(color = col_dose)) + guides(fill = 'none') +
      labs(color = "Treatment", x = "Dose", y = "Mean Accessibility", title = "Example Peak", subtitle = x) + #title = paste("Peak:",x)) + 
      theme(axis.text.x = element_text(angle = -25, hjust = 0)) +
      guides(color = guide_legend(override.aes = list(size = 2, alpha = 0)))
  }
save_plot(g, "./Plots/RepTest_ExamplePeak.pdf", w = 2.5, h = 2, show.legend = 'right')


x <- "MX1"

g <- lm_score$rep[,c("CellType", "Drug", "Dose", "Replicate", x)] %>% {colnames(.)[5] <- "Feature";.} %>% 
  dupControl(controls = "DMSO(+)", drugs = drug[c(1,3:8)]) %>% filter(CellType == "T", Drug %in% c("DMSO(-)","MS177","AU-15330","dCBP-1")) %>% 
  mutate(Group = factor(Drug, levels = c("DMSO(-)","MS177","AU-15330","dCBP-1"), labels = c("Activation","Drug","Drug","Drug"))) %>%
  mutate(Drug = factor(Drug, levels = c("DMSO(-)","MS177","AU-15330","dCBP-1"), labels = c("Activation","MS177","AU-15330","dCBP-1"))) %>%
  {
    ggplot(., aes(x = as.numeric(Dose), y = Feature)) +
      geom_jitter(data = . %>% filter(Group == "Activation"), mapping = aes(color = Dose), width = 0.2, size = 1, show.legend = F) + 
      scale_color_manual(values = c(col_dose[1:2])) + guides(color = 'none') +
      ggnewscale::new_scale_color() +
      geom_jitter(data = . %>% filter(Group == "Drug" & Dose != "DMSO(+)"), mapping = aes(color = Drug), width = 0.2, size = 1, show.legend = F) + 
      geom_smooth(mapping = aes(color = Drug, fill = Drug), method = "lm", alpha = 0.2, show.legend = T) +
      scale_color_manual(values = col_drug) + scale_fill_manual(values = col_drug) + 
      theme_DC + scale_x_continuous(breaks = 1:5, labels = levels(df$Dose)) + theme(axis.text.x = element_text(color = col_dose)) + guides(fill = 'none') +
      labs(color = "Treatment", x = "Dose", y = "Mean Accessibility", title = "Example Gene", subtitle = x) + #title = paste("Gene Score:",x)) + 
      theme(axis.text.x = element_text(angle = -25, hjust = 0)) +
      guides(color = guide_legend(override.aes = list(size = 2, alpha = 0)))
  }
save_plot(g, "./Plots/RepTest_ExampleGene.pdf", w = 2.5, h = 2, show.legend = 'right')

x <- "BACH2"

lm_rna$rep[,c("CellType", "Drug", "Dose", "Replicate", x)] %>% {colnames(.)[5] <- "Feature";.} %>% 
  dupControl(controls = "DMSO(+)", drugs = drug[c(1,3:8)]) %>% filter(CellType == "T", Drug %in% c("DMSO(-)","MS177","AU-15330","dCBP-1")) %>% 
  mutate(Group = factor(Drug, levels = c("DMSO(-)","MS177","AU-15330","dCBP-1"), labels = c("Activation","Drug","Drug","Drug"))) %>%
  mutate(Drug = factor(Drug, levels = c("DMSO(-)","MS177","AU-15330","dCBP-1"), labels = c("Activation","MS177","AU-15330","dCBP-1"))) %>%
  {
    ggplot(., aes(x = as.numeric(Dose), y = Feature)) +
      geom_jitter(data = . %>% filter(Group == "Activation"), mapping = aes(color = Dose), width = 0.2, size = 1, show.legend = F) + 
      scale_color_manual(values = c(col_dose[1:2])) + guides(color = 'none') +
      ggnewscale::new_scale_color() +
      geom_jitter(data = . %>% filter(Group == "Drug" & Dose != "DMSO(+)"), mapping = aes(color = Drug), width = 0.2, size = 1, show.legend = F) + 
      geom_smooth(mapping = aes(color = Drug, fill = Drug), method = "lm", alpha = 0.2, show.legend = T) +
      scale_color_manual(values = col_drug) + scale_fill_manual(values = col_drug) + 
      theme_DC + scale_x_continuous(breaks = 1:5, labels = levels(df$Dose)) + theme(axis.text.x = element_text(color = col_dose)) + guides(fill = 'none') +
      labs(color = "Treatment", x = "Dose", y = "Mean Expression", title = paste("Gene:",x))
  }


### GSEA Analysis - Correlated Marker Genes ---------------------------------

## Determine if gsea is appropriate
data.frame(P = lm_rna$p[10,-1:-2] %>% t, C = lm_rna$coef[10,-1:-2] %>% t) %>% mutate(P = -log10(P) * sign(C), Rank = rank(P, ties.method = "first")) %>% rownames_to_column() %>% # Myeloid AU-15330
  { ggplot(., aes(x = Rank, y = P, label = rowname)) + geom_point() + ggrepel::geom_label_repel(data = . %>% group_by(P > 0) %>% slice_max(order_by = abs(P), n = 10), max.overlaps = 20) }
data.frame(P = lm_rna$p[16,-1:-2] %>% t, C = lm_rna$coef[16,-1:-2] %>% t) %>% mutate(P = -log10(P) * sign(C), Rank = rank(P, ties.method = "first")) %>% rownames_to_column() %>% # T MS177
  { ggplot(., aes(x = Rank, y = P, label = rowname)) + geom_point() + ggrepel::geom_label_repel(data = . %>% group_by(P > 0) %>% slice_max(order_by = abs(P), n = 10), max.overlaps = 20) }

data.frame(P = lm_score$p[10,-1:-2] %>% t, C = lm_score$coef[10,-1:-2] %>% t) %>% mutate(P = -log10(P) * sign(C), Rank = rank(P, ties.method = "first")) %>% rownames_to_column() %>% # Myeloid AU-15330
  { ggplot(., aes(x = Rank, y = P, label = rowname)) + geom_point() + ggrepel::geom_label_repel(data = . %>% group_by(P > 0) %>% slice_max(order_by = abs(P), n = 10), max.overlaps = 20) }
data.frame(P = lm_score$p[16,-1:-2] %>% t, C = lm_score$coef[16,-1:-2] %>% t) %>% mutate(P = -log10(P) * sign(C), Rank = rank(P, ties.method = "first")) %>% rownames_to_column() %>% # T MS177
  { ggplot(., aes(x = Rank, y = P, label = rowname)) + geom_point() + ggrepel::geom_label_repel(data = . %>% group_by(P > 0) %>% slice_max(order_by = abs(P), n = 10), max.overlaps = 20) }

## Rank Genes by P value of linear regression fit (most significant positive --> most significant negative fit)

sets <- list(Hallmarks = msigdbr(species = "Homo sapiens", category = "H") %>% split(x = .$gene_symbol, f = .$gs_name),
             Immune = msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB") %>% split(x = .$gene_symbol, f = .$gs_name),
             CellType = msigdbr(species = "Homo sapiens", category = "C8") %>% split(x = .$gene_symbol, f = .$gs_name),
             Reactome = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% split(x = .$gene_symbol, f = .$gs_name),
             KEGG = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% split(x = .$gene_symbol, f = .$gs_name),
             TF_Targets = msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>% split(x = .$gene_symbol, f = .$gs_name),
             GO_BP =  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% split(x = .$gene_symbol, f = .$gs_name),
             GO_MF =  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>% split(x = .$gene_symbol, f = .$gs_name),
             GO_CC =  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC") %>% split(x = .$gene_symbol, f = .$gs_name))
sets <- unlist(sets, recursive = F)

gsea_lm_rna <- lapply(c("Myeloid","T"), function(c) {
  df <- lapply(levels(df$Drug)[c(1,3:8)], function(d) { 
    genes <- data.frame(P = lm_rna$p %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                        C = lm_rna$coef %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                        Gene = colnames(lm_rna$p[,-1:-2]))
    genes$Gene <- fixGene(genes$Gene, new)
    genes <- genes %>%
      mutate(CorrP = -log10(P) * sign(C)) %>%
      arrange(-CorrP) %>% filter(CorrP != 0) %>% select(Gene, CorrP) %>% tibble::deframe()
    res <- fgseaMultilevel(sets, stats = genes, scoreType = "std", maxSize = 500, nproc = 0)
    res$Collection <- gsub("\\..*", "", res$pathway)
    res$pathway <- gsub(".*\\.", "", res$pathway)
    res$Drug <- d
    res <- res %>% filter(!is.na(log2err), !is.na(padj))
    res
  }) %>% do.call(rbind,.)
  df$CellType <- c
  df
}) %>% do.call(rbind, .)
saveRDS(gsea_lm_rna, "./RepTest/gsea_lm_rna.rds")

gsea_lm_score <- lapply(c("Myeloid","T"), function(c) {
  df <- lapply(levels(df$Drug)[c(1,3:8)], function(d) { 
    genes <- data.frame(P = lm_score$p %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                        C = lm_score$coef %>% filter(CellType == c, Drug == d) %>% ungroup %>% select(-CellType, -Drug) %>% t,
                        Gene = colnames(lm_score$p[,-1:-2])) %>%
      mutate(Gene = fixGene(Gene, new),
             CorrP = -log10(P) * sign(C)) %>%
      arrange(-CorrP) %>% filter(CorrP != 0) %>% select(Gene, CorrP) %>% tibble::deframe()
    res <- fgseaMultilevel(sets, stats = genes, scoreType = "std", maxSize = 500, nproc = 0)
    res$Collection <- gsub("\\..*", "", res$pathway)
    res$pathway <- gsub(".*\\.", "", res$pathway)
    res$Drug <- d
    res <- res %>% filter(!is.na(log2err), !is.na(padj))
    res
  }) %>% do.call(rbind,.)
  df$CellType <- c
  df
}) %>% do.call(rbind, .)
saveRDS(gsea_lm_score, "./RepTest/gsea_lm_score.rds")


## Plot Heatmaps of Correlated Terms
terms <- 
  gsea_lm_rna %>%
  filter(Collection %in% c("GO_BP","Hallmarks")) %>%
  filter(padj < 0.001) %>%
  filter(Drug != 'DMSO(-)') %>%
  arrange(padj) %>% {.[!duplicated(.$pathway),]} %>%
  group_by(CellType, Drug, NES>0) %>% slice_max(order_by = -padj, n = 3, with_ties = F) %>% .$pathway %>% unique

gsea_lm_rna %>%
  filter(pathway %in% terms, padj < 0.01) %>%
  filter(Drug != 'DMSO(-)') %>%
  # mutate(NES = -log10(padj) * sign(NES)) %>%
  # mutate(pathway = gsub("(?<=.{40}).+","...",pathway,perl=T)) %>% .$pathway # keep first 40 characters, replace rest with ellipses
  mutate(pathway = gsub("^(.{35})(.*)(.{10})$", "\\1...\\3",pathway)) %>% # keep first 35 characters and final 10, replace middle with ellipses
  # mutate(pathway = gsub("^(.{30})(.*)_.*_([^_]*$)", "\\1...\\3",pathway)) %>%  # keep first 30 characters and last word, replace middle with ellipses
  ungroup() %>% dplyr::select(pathway, CellType, Drug, NES) %>%
  as.data.frame() %>% pivot_wider(id_cols = 1, names_from = c(CellType, Drug), values_from = NES, values_fill = 0) %>% 
  tibble::column_to_rownames(var = "pathway") %>% 
  pheatmap(scale = 'none', display_numbers = F, border_color = NA, fontsize = 7, fontsize_row = 7,
           show_colnames = F, treeheight_col = 10, treeheight_row = 0,
           color = paletteer::paletteer_c("grDevices::Blue-Red 3", 100, direction = 1), breaks = seq(-max(abs(.)),max(abs(.)),length = 100),
           annotation_col = markers_pos %>% dplyr::select(CellType, Treatment) %>% distinct %>% mutate(rowname = paste(CellType, Treatment, sep="_")) %>% tibble::column_to_rownames(),
           annotation_colors = list(Treatment = col_drug[-1:-3], CellType = col_cellmaj[-1]), annotation_legend = T,
           filename = if (plot) {"./Plots/RepTest_RNA_GSEA_Heatmap.pdf"} else { NA }, width = 5, height = 5, main = "GSEA\nRNA Markers")

terms <- 
  gsea_lm_score %>%
  filter(Collection %in% c("GO_BP","Hallmarks")) %>%
  filter(padj < 0.001) %>%
  filter(Drug != 'DMSO(-)') %>%
  arrange(padj) %>% {.[!duplicated(.$pathway),]} %>% 
  group_by(CellType, Drug, NES>0) %>% slice_max(order_by = -padj, n = 3, with_ties = F) %>% .$pathway %>% unique

gsea_lm_score %>%
  filter(pathway %in% terms, padj < 0.01) %>%
  filter(Drug != 'DMSO(-)') %>%
  # mutate(NES = -log10(padj) * sign(NES)) %>%
  # mutate(pathway = gsub("(?<=.{45}).+","...",pathway,perl=T)) %>%
  mutate(pathway = gsub("^(.{35})(.*)(.{10})$", "\\1...\\3",pathway)) %>% # keep first 35 characters and final 10, replace middle with ellipses
  # mutate(pathway = gsub("^(.{30})(.*)_.*_([^_]*$)", "\\1...\\3",pathway)) %>%  # keep first 30 characters and last word, replace middle with ellipses
  ungroup() %>% dplyr::select(pathway, CellType, Drug, NES) %>%
  as.data.frame() %>% pivot_wider(id_cols = 1, names_from = c(CellType, Drug), values_from = NES, values_fill = 0) %>% 
  tibble::column_to_rownames(var = "pathway") %>% 
  pheatmap(scale = 'none', display_numbers = F, border_color = NA, fontsize = 7, fontsize_row = 7,
           show_colnames = F, treeheight_col = 10, treeheight_row = 0, 
           color = paletteer::paletteer_c("grDevices::Blue-Red 3", 100, direction = 1), breaks = seq(-max(abs(.)),max(abs(.)),length = 100),
           annotation_col = markers_pos %>% dplyr::select(CellType, Treatment) %>% distinct %>% mutate(rowname = paste(CellType, Treatment, sep="_")) %>% tibble::column_to_rownames(),
           annotation_colors = list(Treatment = col_drug[-1:-3], CellType = col_cellmaj[-1]), annotation_legend = T,
           filename = if (plot) {"./Plots/RepTest_Score_GSEA_Heatmap.pdf"} else { NA }, width = 5, height = 5, main = "GSEA\nGene Scores")


### Plot Gene Sets ----------------------------------------------------------

## only TNFA and IFNA
set <- sets[grep("Hallmark.*(TNFA|Alpha)", names(sets), ignore.case = T, value = T)][c(2,1)]
g <- plotGeneSet2(set, log2fc = T, facets = c("Type","CellType"), hline = 0, nrow = 1)
save_plot(g, "./Plots/RepTest_Combined_GeneSets_Log2FC_v2.pdf", w = 4, h = 2, show.legend = "right")
g <- plotGeneSet2(set, log2fc = F, facets = c("Type","CellType"), hline = NULL, nrow = 1)
save_plot(g, "./Plots/RepTest_Combined_GeneSets_Scaled_v2.pdf", w = 4, h = 2, show.legend = "right")


### Dose-Responsive Marker Heatmaps -----------------------------------------

nTerm <- 8

## Differential Peaks/Motifs - T Cells

# Use only the peaks with strong linear fit
# res <- plotHeatmap(lm_atac, assay = "Peaks", cellType = "T", p = 0.01, log2fc = 1, cutCol = 6, returnRes = T)
res <- plotHeatmap(lm_atac, assay = "Peaks", cellType = "T", p = 0.01, log2fc = 1, cutCol = 0, returnRes = T, filename = "./Plots/RepTest_Peaks_Marker_Heatmap_T.pdf", width = 6, height = 3)

## Differential Peaks/Motifs - Myeloid Cells

# res <- plotHeatmap(lm_atac, assay = "Peaks", cellType = "Myeloid", p = 0.01, log2fc = 1, cutCol = 8, returnRes = T)
res <- plotHeatmap(lm_atac, assay = "Peaks", cellType = "Myeloid", p = 0.01, log2fc = 1, cutCol = 0, returnRes = T, filename = "./Plots/RepTest_Peaks_Marker_Heatmap_M.pdf", width = 6, height = 3)

## Gene Score Enrichment - T Cells

# res <- plotHeatmap(lm_score, assay = "Scores", cellType = "T", p = 0.01, log2fc = 1, cutCol = 8, returnRes = T)
res <- plotHeatmap(lm_score, assay = "Scores", cellType = "T", p = 0.01, log2fc = 1, cutCol = 0, returnRes = T, filename = "./Plots/RepTest_Scores_Marker_Heatmap_T.pdf", width = 6, height = 3)

## Gene Score Enrichment - Myeloid Cells

# res <- plotHeatmap(lm_score, assay = "Scores", cellType = "Myeloid", p = 0.01, log2fc = 1, cutCol = 5, returnRes = T)
res <- plotHeatmap(lm_score, assay = "Scores", cellType = "Myeloid", p = 0.01, log2fc = 1, cutCol = 0, returnRes = T, filename = "./Plots/RepTest_Scores_Marker_Heatmap_M.pdf", width = 6, height = 3)

## Gene Expression Enrichment - T Cells

# res <- plotHeatmap(lm_rna, assay = "RNA", cellType = "T", p = 0.01, log2fc = 1, cutCol = 6, returnRes = T)
res <- plotHeatmap(lm_rna, assay = "RNA", cellType = "T", p = 0.01, log2fc = 1, cutCol = 0, returnRes = T, filename = "./Plots/RepTest_RNA_Marker_Heatmap_T.pdf", width = 6, height = 3)

## Gene Expression Enrichment - Myeloid Cells

# res <- plotHeatmap(lm_rna, assay = "RNA", cellType = "Myeloid", p = 0.01, log2fc = 1,  cutCol = 5, returnRes = T)
res <- plotHeatmap(lm_rna, assay = "RNA", cellType = "Myeloid", p = 0.01, log2fc = 1, cutCol = 0, returnRes = T, filename = "./Plots/RepTest_RNA_Marker_Heatmap_M.pdf", width = 6, height = 3)


### Motif vs TF Expression Correlation --------------------------------------

devrnacor <- lapply(c("T","Myeloid"), function(c) {
  lapply(drug[3:8], function(d) {
    tmp_dev <- lm_dev %>%
      .[c("p","coef")] %>%
      { .$p[,-1:-2] <- -log10(.$p[,-1:-2]) * sign(.$coef[,-1:-2]); .$p } %>%
      pivot_longer(cols = 3:ncol(.), names_to = "Gene", values_to = "Dev_P") %>%
      filter(CellType == c, Drug == d) %>%
      mutate(Gene = gsub("_.*","",Gene))
    tmp_rna <- lm_rna %>%
      .[c("p","coef")] %>%
      lapply(function(x) {
        x %>% ungroup %>% {.[,c("CellType","Drug",grep(paste(unique(tmp_dev$Gene), collapse = "$|^"),colnames(x),value=T))]}
      }) %>%
      { .$p[,-1:-2] <- -log10(.$p[,-1:-2]) * sign(.$coef[,-1:-2]); .$p } %>%
      pivot_longer(cols = 3:ncol(.), names_to = "Gene", values_to = "RNA_P") %>%
      filter(CellType == c, Drug == d) 
    tmp <- join(tmp_dev, tmp_rna) %>% filter(!is.na(RNA_P))
  }) %>% do.call(rbind,.)
}) %>% do.call(rbind, .)


g <- devrnacor %>%
  filter(Drug == "MS177" & CellType == "T") %>%
  # filter(Drug == "MS177" & CellType == "Myeloid") %>%
  filter(Gene %ni% grep("ZNF",Gene,value=T)) %>%
  mutate(Drug = factor(Drug, levels(df$Drug)),
         CellType = factor(CellType, c("T","Myeloid"))) %>%
  ggplot(aes(x = Dev_P, y = RNA_P, color = Drug, label = Gene)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') + geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey40') + 
  geom_point(mapping = aes(alpha = abs(RNA_P) > 2 & abs(Dev_P) > 2), size = 0.75, show.legend = F) + 
  ggrepel::geom_label_repel(data = . %>% filter(RNA_P > 2, abs(Dev_P) > 2) %>% group_by(Drug, CellType, RNA_P > 0, Dev_P > 0) %>% slice_max(order_by = abs(RNA_P) * abs(Dev_P), n = 5), color = 'black', max.overlaps = 20, size = 1.5, nudge_x = 0, label.padding = 0.09, alpha = 0.75) +
  ggrepel::geom_label_repel(data = . %>% filter(RNA_P < -2, abs(Dev_P) > 2) %>% group_by(Drug, CellType, RNA_P > 0, Dev_P > 0) %>% slice_max(order_by = abs(RNA_P) * abs(Dev_P), n = 5), color = 'black', max.overlaps = 20, size = 1.5, nudge_x = 0, label.padding = 0.09, alpha = 0.75) +
  scale_x_continuous(breaks = c(-7,7), labels = c("Down","Up"),limits=c(-14,14)) +
  scale_y_continuous(breaks = c(-3,3), labels = c("Down","Up"),limits=c(-6,6)) +
  labs(x = "Motif Accessibility", y = "TF Expression") +
  facet_wrap(Drug~., axes = 'all', nrow = 1) +
  scale_color_manual(values = col_drug) + 
  theme_DC + theme(axis.ticks = element_blank())
save_plot(g, "./Plots/DevRNACorrelation_Motif_MS177_T.pdf", w = 1.5, h = 1.5, show.legend = 'none')
# save_plot(g, "./Plots/DevRNACorrelation_Motif_MS177_M.pdf", w = 1.5, h = 1.5, show.legend = 'none')

g <- devrnacor %>%
  # filter(Drug == "MS177" & CellType == "T") %>%
  filter(Drug %in% c("MS177","EPZ6438") & CellType == "T") %>%
  filter(Gene %ni% grep("ZNF",Gene,value=T)) %>%
  mutate(Drug = factor(Drug, levels(df$Drug)),
         CellType = factor(CellType, c("T","Myeloid"))) %>%
  ggplot(aes(x = Dev_P, y = RNA_P, color = Drug, label = Gene)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') + geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey40') + 
  geom_point(mapping = aes(alpha = abs(RNA_P) > 2 & abs(Dev_P) > 2, shape = Drug), size = 1.5, show.legend = T) + 
  ggrepel::geom_label_repel(data = . %>% filter(abs(RNA_P) > 2, abs(Dev_P) > 2) %>% group_by(Drug, CellType, RNA_P > 0, Dev_P > 0) %>% slice_max(order_by = abs(RNA_P) * abs(Dev_P), n = 5), color = 'black', max.overlaps = 20, min.segment.length = 0, size = 2, nudge_x = 0, label.padding = 0.09, alpha = 1) +
  scale_x_continuous(breaks = c(-7,7), labels = c("Down","Up"),limits=c(-14,14)) +
  scale_y_continuous(breaks = c(-3,3), labels = c("Down","Up"),limits=c(-6,6)) +
  scale_shape_manual(values = c(16,15)) + scale_alpha_manual(values = c(0.1,1)) +
  labs(x = "Motif Accessibility", y = "TF Expression", title = "Gene-Motif Correlation", subtitle = "T Cells", color = NULL, shape = NULL) + guides(alpha = 'none', color = guide_legend(override.aes = list(size = 2.5))) +
  scale_color_manual(values = col_drug) + 
  theme_DC + theme(axis.ticks = element_blank(), legend.background = element_blank())
save_plot(g, "./Plots/DevRNACorrelation_Motif_MS177_T_v2.pdf", w = 2, h = 2, show.legend = c(0.2,0.95))

g <- devrnacor %>%
  # devrnacor2 %>%
  filter(Drug == "BRM014" & CellType == "T") %>%
  # filter(Drug == "BRM014" & CellType == "Myeloid") %>%
  filter(Gene %ni% grep("ZNF",Gene,value=T)) %>%
  mutate(Drug = factor(Drug, levels(df$Drug)),
         CellType = factor(CellType, c("T","Myeloid"))) %>%
  ggplot(aes(x = Dev_P, y = RNA_P, color = Drug, label = Gene)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') + geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey40') + 
  geom_point(mapping = aes(alpha = abs(RNA_P) > 2 & abs(Dev_P) > 2), size = 0.75, show.legend = F) + 
  ggrepel::geom_label_repel(data = . %>% filter(RNA_P > 2, abs(Dev_P) > 2) %>% group_by(Drug, CellType, RNA_P > 0, Dev_P > 0) %>% slice_max(order_by = abs(RNA_P) * abs(Dev_P), n = 5), color = 'black', max.overlaps = 20, size = 1.5, nudge_x = 0, label.padding = 0.09, alpha = 0.75) +
  ggrepel::geom_label_repel(data = . %>% filter(RNA_P < -2, abs(Dev_P) > 2) %>% group_by(Drug, CellType, RNA_P > 0, Dev_P > 0) %>% slice_max(order_by = abs(RNA_P) * abs(Dev_P), n = 5), color = 'black', max.overlaps = 20, size = 1.5, nudge_x = 0, label.padding = 0.09, alpha = 0.75) +
  scale_x_continuous(breaks = c(-7,7), labels = c("Down","Up"),limits=c(-14,14)) +
  scale_y_continuous(breaks = c(-3,3), labels = c("Down","Up"),limits=c(-6,6)) +
  labs(x = "Motif Accessibility", y = "TF Expression") +
  facet_wrap(Drug~., axes = 'all', nrow = 1) +
  scale_color_manual(values = col_drug) + 
  theme_DC + theme(axis.ticks = element_blank())
# save_plot(g, "./Plots/DevRNACorrelation_Motif_BRM014_T.pdf", w = 1.5, h = 1.5, show.legend = 'none')
save_plot(g, "./Plots/DevRNACorrelation_Motif_BRM014_M.pdf", w = 1.5, h = 1.5, show.legend = 'none')

g <- devrnacor %>%
  # devrnacor2 %>%
  filter(Drug %in% c("BRM014","AU-15330") & CellType == "T") %>%
  # filter(Drug %in% c("BRM014","AU-15330") & CellType == "Myeloid") %>%
  filter(Gene %ni% grep("ZNF",Gene,value=T)) %>%
  mutate(Drug = factor(Drug, levels(df$Drug)),
         CellType = factor(CellType, c("T","Myeloid"))) %>%
  ggplot(aes(x = Dev_P, y = RNA_P, color = Drug, label = Gene)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') + geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey40') + 
  geom_point(mapping = aes(alpha = abs(RNA_P) > 2 & abs(Dev_P) > 2, shape = Drug), size = 0.75, show.legend = F) + 
  ggrepel::geom_label_repel(data = . %>% filter(RNA_P > 2, abs(Dev_P) > 2) %>% group_by(CellType, RNA_P > 0, Dev_P > 0) %>% slice_max(order_by = abs(RNA_P) * abs(Dev_P), n = 5), color = 'black', max.overlaps = 20, size = 1.5, nudge_x = 0, label.padding = 0.09, alpha = 0.75) +
  ggrepel::geom_label_repel(data = . %>% filter(RNA_P < -2, abs(Dev_P) > 2) %>% group_by(CellType, RNA_P > 0, Dev_P > 0) %>% slice_max(order_by = abs(RNA_P) * abs(Dev_P), n = 5), color = 'black', max.overlaps = 20, size = 1.5, nudge_x = 0, label.padding = 0.09, alpha = 0.75) +
  scale_x_continuous(breaks = c(-7,7), labels = c("Down","Up"),limits=c(-14,14)) +
  scale_y_continuous(breaks = c(-3,3), labels = c("Down","Up"),limits=c(-6,6)) +
  scale_shape_manual(values = c(15,16)) +
  labs(x = "Motif Accessibility", y = "TF Expression") +
  # facet_wrap(Drug~., axes = 'all', nrow = 1) +
  scale_color_manual(values = col_drug) + 
  theme_DC +  theme(axis.ticks = element_blank())
save_plot(g, "./Plots/DevRNACorrelation_Motif_SWISNF_T.pdf", w = 1.5, h = 1.5, show.legend = 'none')
save_plot(g, "./Plots/DevRNACorrelation_Motif_SWISNF_M.pdf", w = 1.5, h = 1.5, show.legend = 'none')

g <- devrnacor %>%
  filter(Drug %in% c("BRM014","AU-15330") & CellType == "T") %>%
  filter(Gene %ni% grep("ZNF",Gene,value=T)) %>%
  mutate(Drug = factor(Drug, levels(df$Drug)),
         CellType = factor(CellType, c("T","Myeloid"))) %>%
  ggplot(aes(x = Dev_P, y = RNA_P, color = Drug, label = Gene)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') + geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey40') + 
  geom_point(mapping = aes(alpha = abs(RNA_P) > 2 & abs(Dev_P) > 2, shape = Drug), size = 1.5, show.legend = T) + 
  ggrepel::geom_label_repel(data = . %>% filter(abs(RNA_P) > 2, abs(Dev_P) > 2) %>% group_by(Drug, CellType, RNA_P > 0, Dev_P > 0) %>% slice_max(order_by = abs(RNA_P) * abs(Dev_P), n = 5), color = 'black', max.overlaps = 20, min.segment.length = 0, size = 2, nudge_x = 0, label.padding = 0.09, alpha = 1) +
  scale_x_continuous(breaks = c(-7,7), labels = c("Down","Up"),limits=c(-14,14)) +
  scale_y_continuous(breaks = c(-3,3), labels = c("Down","Up"),limits=c(-6,6)) +
  scale_shape_manual(values = c(16,15)) + scale_alpha_manual(values = c(0.1,1)) +
  labs(x = "Motif Accessibility", y = "TF Expression", title = "Gene-Motif Correlation", subtitle = "T Cells", color = NULL, shape = NULL) + guides(alpha = 'none', color = guide_legend(override.aes = list(size = 2.5))) +
  scale_color_manual(values = col_drug) + 
  theme_DC + theme(axis.ticks = element_blank(), legend.background = element_blank())
save_plot(g, "./Plots/DevRNACorrelation_Motif_SWISNF_T_v2.pdf", w = 2, h = 2, show.legend = c(0.2,0.95))


### Investigate Marker Peaks ------------------------------------------------

## Overlap with Enhancer Sites
enh <- read.table("../FeatureSets/FANTOM5 Enhancers/F5.hg38.enhancers.bed", header = FALSE, fill = TRUE) %>% mutate(width = abs(V3-V2))
enh <- GRanges(seqnames = enh$V1, IRanges(start = enh$V2, end = enh$V3), strand = "*")

g <- markers_lm_sig %>% filter(MarkerType == "Peaks") %>%
  separate(name, into = c("seqnames","start","end"), sep = "-") %>% mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  group_by(CellType, Treatment, Type) %>% 
  mutate(nDAR = n(),
         Overlap = GRanges(seqnames, IRanges(start,end)) %>% overlapsAny(enh)) %>% 
  group_by(CellType, Treatment, Type, nDAR) %>% summarize(nOverlaps = sum(Overlap)) %>%
  mutate(CellType = factor(CellType,levels = c("T","Myeloid")),
         Type = factor(Type, levels = c("Up","Down"))) %>%
  ggplot(aes(y = Treatment, x = nOverlaps/nDAR, fill = Type)) + geom_col(position = position_dodge2(reverse = TRUE, padding = 0), width = 0.7) + 
  facet_grid(~CellType, scales = 'free_x', axes = 'all_y') + scale_y_discrete(drop=FALSE, limits = rev) + 
  paletteer::scale_fill_paletteer_d("ggsci::alternating_igv", direction = -1) +
  theme_DC + theme(axis.text.y = element_text(color = rev(col_drug[-2:-3]))) +
  labs(x = "Fraction DARs Overlapping with FANTOM5 Enhancer Set")
save_plot(g, "./Plots/Markers_Peaks_EnhancerOverlap.pdf", w = 3, h = 1.5, show.legend = 'right')

g <- markers_lm_sig %>% filter(MarkerType == "Peaks") %>%
  separate(name, into = c("seqnames","start","end"), sep = "-") %>% mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  group_by(CellType, Treatment, Type) %>% 
  mutate(nDAR = n(),
         Overlap = GRanges(seqnames, IRanges(start,end)) %>% overlapsAny(enh)) %>% 
  group_by(CellType, Treatment, Type, nDAR) %>% summarize(nOverlaps = sum(Overlap)) %>%
  mutate(CellType = factor(CellType,levels = c("T","Myeloid")),
         Type = factor(Type, levels = c("Up","Down")),
         Fraction = nOverlaps/nDAR) %>%
  {
    res <- lapply(seq(2, nrow(.), by = 2), function(x) {
      test <- prop.test(as.matrix(.[c((x-1):x), c(5,4)])) %>% broom::tidy()
      data.frame(CellType = .[x,1],
                 Treatment = .[x,2],
                 Type = .[x,3],
                 start = "Up", 
                 end = "Down",
                 Fraction = max(.[(x-1):x,6]),
                 p = test$p.value)
    }) %>% do.call(rbind,.) %>%
      mutate(p = symnum(p,
                        symbols   = c("***","**","*",".","n.s."),
                        cutpoints = c(0,  .001,.01,.05, .1, 1),
                        corr      = FALSE
      ))
    ggplot(., aes(x = Type, y = Fraction, fill = Type)) +
      geom_col(width = 1) +
      geom_signif(data = res,
                  mapping = aes(x = Type, xmin = start, xmax = end,
                                annotations = p),
                  y_position = max(res$Fraction) + 0.02,
                  textsize = 1.75, vjust = -0.25, size = 0.15,
                  manual = T) +
      facet_grid(Treatment~CellType, scales = 'free_x', axes = 'all_y', switch = 'y') +
      scale_x_discrete(expand = expand_scale(add=0.75), limits = rev) + coord_flip(clip = 'off') +
      paletteer::scale_fill_paletteer_d("ggsci::alternating_igv", direction = -1) +
      theme_DC + theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       panel.spacing.y = unit(0, "null"),
                       strip.text.y.left = element_text(angle = 0)) +
      labs(y = "Fraction DARs Overlapping with FANTOM5 Enhancer Set", x = NULL)
  } 
save_plot(g, "./Plots/Markers_Peaks_EnhancerOverlap_v2.pdf", w = 3, h = 1.5, show.legend = 'right')


## Overlap With Insulator Sites
# CTCFBSDB predicted CTCF binding sites (lift to hg38 from hg18 using UCSC Genome Browser)
ctcf <- read.table("../FeatureSets/CTCFBSDB Predicted/allcomp_hg38_lifted.bed", sep = '\t')
ctcf <- ctcf[ctcf$V1 %in% paste("chr",c(1:22,"X","Y"), sep = ""),]
ctcf <- GRanges(seqnames = ctcf$V1, IRanges(start = ctcf$V2, end = ctcf$V3), strand = "*")

g <- markers_lm_sig %>% filter(MarkerType == "Peaks") %>%
  separate(name, into = c("seqnames","start","end"), sep = "-") %>% mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  group_by(CellType, Treatment, Type) %>% 
  mutate(nDAR = n(),
         Overlap = GRanges(seqnames, IRanges(start,end)) %>% overlapsAny(ctcf)) %>% 
  group_by(CellType, Treatment, Type, nDAR) %>% summarize(nOverlaps = sum(Overlap)) %>%
  mutate(CellType = factor(CellType,levels = c("T","Myeloid")),
         Type = factor(Type, levels = c("Up","Down"))) %>%
  ggplot(aes(y = Treatment, x = nOverlaps/nDAR, fill = Type)) + geom_col(position = position_dodge2(reverse = TRUE, padding = 0), width = 0.7) + 
  facet_grid(~CellType, scales = 'free_x', axes = 'all_y') + scale_y_discrete(drop=FALSE, limits = rev) + 
  paletteer::scale_fill_paletteer_d("ggsci::alternating_igv", direction = -1) +
  theme_DC + theme(axis.text.y = element_text(color = rev(col_drug[-2:-3]))) +
  labs(x = "Fraction DARs Overlapping with CTCFBSDB CTCF Sites")
save_plot(g, "./Plots/Markers_Peaks_CTCFOverlap.pdf", w = 3, h = 1.5, show.legend = 'right')


g <- markers_lm_sig %>% filter(MarkerType == "Peaks") %>%
  separate(name, into = c("seqnames","start","end"), sep = "-") %>% mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  group_by(CellType, Treatment, Type) %>% 
  mutate(nDAR = n(),
         Overlap = GRanges(seqnames, IRanges(start,end)) %>% overlapsAny(ctcf)) %>% 
  group_by(CellType, Treatment, Type, nDAR) %>% summarize(nOverlaps = sum(Overlap)) %>%
  mutate(CellType = factor(CellType,levels = c("T","Myeloid")),
         Type = factor(Type, levels = c("Up","Down")),
         Fraction = nOverlaps/nDAR) %>%
  {
    res <- lapply(seq(2, nrow(.), by = 2), function(x) {
      test <- prop.test(as.matrix(.[c((x-1):x), c(5,4)])) %>% broom::tidy()
      data.frame(CellType = .[x,1],
                 Treatment = .[x,2],
                 Type = .[x,3],
                 start = "Up", 
                 end = "Down",
                 Fraction = max(.[(x-1):x,6]),
                 p = test$p.value)
    }) %>% do.call(rbind,.) %>%
      mutate(p = symnum(p,
                        symbols   = c("***","**","*",".","n.s."),
                        cutpoints = c(0,  .001,.01,.05, .1, 1),
                        corr      = FALSE
      ))
    ggplot(., aes(x = Type, y = Fraction, fill = Type)) +
      geom_col(width = 1) +
      geom_signif(data = res,
                  mapping = aes(x = Type, xmin = start, xmax = end,
                                annotations = p),
                  y_position = max(res$Fraction) + 0.005,
                  textsize = 1.75, vjust = -0.25, size = 0.15,
                  manual = T) +
      facet_grid(Treatment~CellType, scales = 'free_x', axes = 'all_y', switch = 'y') +
      scale_x_discrete(expand = expand_scale(add=0.75), limits = rev) + coord_flip(clip = 'off') +
      paletteer::scale_fill_paletteer_d("ggsci::alternating_igv", direction = -1) +
      theme_DC + theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       panel.spacing.y = unit(0, "null"),
                       strip.text.y.left = element_text(angle = 0)) +
      labs(y = "Fraction DARs Overlapping with CTCFBSDB CTCF Sites", x = NULL)
  } 
save_plot(g, "./Plots/Markers_Peaks_CTCFOverlap_v2.pdf", w = 3, h = 1.5, show.legend = 'right')


### Compare Drug Peaks to Activation Peaks ------------------------------------

# Activation = DMSO+ vs DMSO-
# Drug = Drug vs DMSO+
# Drug Reduces Activation -> Positive Activation Peaks Close, Negative Activation Peaks Open
# Drug Amplifies Activation -> Positive Activation Peaks More Open, Negative Activation Peaks More Closed
# Drug-Unique Effects -> Unique From Activation Peaks
# Plot per drug:
# • Inverted Markers
# • Amplified Markers
# • Unique Markers

## T Cells
activ <- markers_lm_sig %>%  
  filter(CellType == "T",
         MarkerType == "Peaks",
         Treatment == "Activation",
         abs(Log2FC) > 1,
  ) %>% dplyr::select(name,Type)

g <- markers_lm_sig %>%
  filter(CellType == "T",
         MarkerType == "Peaks", 
         Treatment != "Activation",
         abs(Log2FC) > 1,
  ) %>%
  mutate(Type = case_when(name %in% activ$name & Type == activ$Type[match(name, activ$name)] ~ "Amplified",
                          name %in% activ$name & Type != activ$Type[match(name, activ$name)] ~ "Inverted",
                          name %ni% activ$name ~ "Unique")) %>%
  group_by(Treatment, MarkerType, CellType) %>% 
  summarize(Drug = n(), Activation = nrow(activ), Amplified = sum(Type == "Amplified"), Inverted = sum(Type == "Inverted"), Numerator = Amplified+Inverted) %>%
  mutate(Amplified = Amplified/min(Drug,Activation),
         Inverted = Inverted/min(Drug,Activation)) %>%
  dplyr::select(-Drug,-Activation) %>% melt(id.vars = c(1:3,6), variable.name = "Type", value.name = "Freq") %>%
  mutate(Treatment = factor(Treatment, levels = levels(df$Drug)[3:8])) %>%
  {
    ggplot(., aes(x = Treatment, y = Freq, fill = Type)) + 
      geom_col(position = "stack") + 
      geom_text(data = . %>% group_by(MarkerType,Treatment) %>% summarize(Total = mean(Numerator), Freq = sum(Freq)) %>% filter(Total > 0), 
                mapping = aes(label = Total, fill = NULL, y = Freq + 0.04), size = 3, hjust = 0) + coord_flip() + scale_x_discrete(limits = rev, drop = F) +
      scale_fill_manual(values = paletteer::paletteer_d("RColorBrewer::PRGn", n = 10)[c(2,9)]) +
      ylim(c(0,1)) + theme_DC + ggtitle("Activation Peaks (T)") + ylab("Overlap")
  }
save_plot(g, "./Plots/ActivationPeakComparison_OverlapCoef_T.pdf", w = 2, h = 1.5, show.legend = c(0.8,0.7))

## Myeloid Cells
activ <- markers_lm_sig %>%
  filter(CellType == "Myeloid", 
         MarkerType == "Peaks",
         Treatment == "Activation",
         abs(Log2FC) > 1,
  ) %>% dplyr::select(name,Type)

g <- markers_lm_sig %>%
  filter(CellType == "Myeloid",
         MarkerType == "Peaks", 
         Treatment != "Activation",
         abs(Log2FC) > 1,
  ) %>%
  mutate(Type = case_when(name %in% activ$name & Type == activ$Type[match(name, activ$name)] ~ "Amplified",
                          name %in% activ$name & Type != activ$Type[match(name, activ$name)] ~ "Inverted",
                          name %ni% activ$name ~ "Unique")) %>%
  group_by(Treatment, MarkerType, CellType) %>% 
  summarize(Drug = n(), Activation = nrow(activ), Amplified = sum(Type == "Amplified"), Inverted = sum(Type == "Inverted"), Numerator = Amplified+Inverted) %>%
  mutate(Amplified = Amplified/min(Drug,Activation),
         Inverted = Inverted/min(Drug,Activation)) %>%
  dplyr::select(-Drug,-Activation) %>% melt(id.vars = c(1:3,6), variable.name = "Type", value.name = "Freq") %>%
  mutate(Treatment = factor(Treatment, levels = levels(df$Drug)[3:8])) %>%
  {
    ggplot(., aes(x = Treatment, y = Freq, fill = Type)) + 
      geom_col(position = "stack") + 
      geom_text(data = . %>% group_by(MarkerType,Treatment) %>% summarize(Total = mean(Numerator), Freq = sum(Freq)) %>% filter(Total > 0), 
                mapping = aes(label = Total, fill = NULL, y = Freq + 0.04), size = 3, hjust = 0) + coord_flip() + scale_x_discrete(limits = rev, drop = F) +
      scale_fill_manual(values = paletteer::paletteer_d("RColorBrewer::PRGn", n = 10)[c(2,9)]) +
      ylim(c(0,1)) + theme_DC + ggtitle("Activation Peaks (Myeloid)") + ylab("Overlap")
  }
save_plot(g, "./Plots/ActivationPeakComparison_OverlapCoef_Myeloid.pdf", w = 2, h = 1.5, show.legend = c(0.8,0.7))


## example coverage plots of the different "classes" of markers compared to activation

activ <- markers_lm_sig %>%
  filter(CellType == "T", 
         MarkerType == "Peaks",
         Treatment == "Activation",
         abs(Log2FC) > 1,
         Mean > 0.1) %>% dplyr::select(name,Type)

g <- markers_lm_sig %>% 
  filter(Treatment == "MS177",
         CellType == "T",
         MarkerType == "Peaks") %>%
  mutate(Comp = case_when(name %in% activ$name & Type == activ$Type[match(name, activ$name)] ~ "Amplified",
                          name %in% activ$name & Type != activ$Type[match(name, activ$name)] ~ "Inverted",
                          name %ni% activ$name ~ "Unique")) %>%
  group_by(Comp) %>% slice_max(order_by = -log10(FDR) * abs(Log2FC), n = 1) %>%
  dplyr::select(name, Type, Comp, Treatment) %>%
  apply(.,1,function(x) {
    g <- CovPlot(seu[,seu$CellType_Major == "T" & seu$Drug %in% c("DMSO(-)","DMSO(+)","MS177")], x[1], group.by = "Dose", 2000, 2000, col_dose, ratio = c(7,1,1))
    if (x[3] != "Unique") { g[[1]] <- g[[1]] + ggtitle(paste(x[4], "-", x[3], " Activation Peak (", x[2], ")", sep="")) } 
    else { g[[1]] <- g[[1]] + ggtitle(paste(x[4], "-", x[3], " Peak (", x[2], ")",sep="")) }
    g
  }) %>% plot_grid(plotlist = ., byrow = F, nrow = 1)
save_plot(g, "./Plots/ActivationPeakComparison_ExamplePeaks_MS177.pdf", w = 15, h = 5, show.legend = 'none')

g <- markers_lm_sig %>% 
  filter(Treatment == "BRM014",
         CellType == "T",
         MarkerType == "Peaks") %>%
  mutate(Comp = case_when(name %in% activ$name & Type == activ$Type[match(name, activ$name)] ~ "Amplified",
                          name %in% activ$name & Type != activ$Type[match(name, activ$name)] ~ "Inverted",
                          name %ni% activ$name ~ "Unique")) %>%
  group_by(Comp) %>% slice_max(order_by = -log10(FDR) * abs(Log2FC), n = 1) %>%
  dplyr::select(name, Type, Comp, Treatment) %>%
  apply(.,1,function(x) {
    g <- CovPlot(seu[,seu$CellType_Major == "T" & seu$Drug %in% c("DMSO(-)","DMSO(+)","BRM014")], x[1], group.by = "Dose", 2000, 2000, col_dose, ratio = c(7,1,1))
    if (x[3] != "Unique") { g[[1]] <- g[[1]] + ggtitle(paste(x[4], "-", x[3], " Activation Peak (", x[2], ")", sep="")) } 
    else { g[[1]] <- g[[1]] + ggtitle(paste(x[4], "-", x[3], " Peak (", x[2], ")",sep="")) }
    g
  }) %>% plot_grid(plotlist = ., byrow = F, nrow = 1)
save_plot(g, "./Plots/ActivationPeakComparison_ExamplePeaks_BRM014.pdf", w = 15, h = 5, show.legend = 'none')


### Activation vs Drug Scoring ----------------------------------------------

## PseudoBulk Counts
sePeakRep <- se(proj_4, useMatrix = "PeakMatrix", groupBy = "paste(CellDrugDose,Replicate,sep='_')")
seScoreRep <- se(proj_4, useMatrix = "GeneScoreMatrix", groupBy = "paste(CellDrugDose,Replicate,sep='_')")
seGeneRep <- se(proj_4, useMatrix = "GeneExpressionMatrix", groupBy = "paste(CellDrugDose,Replicate,sep='_')")

g <- drugScore(seGeneRep, markers_lm_sig, "Genes")
save_plot(g, "./Plots/ActivationDrugScores_GeneReps_ByTarget.pdf", w = 4, h = 2, show.legend = 'right')

g <- drugScore(seScoreRep, markers_lm_sig, "Scores")
save_plot(g, "./Plots/ActivationDrugScores_ScoreReps_ByTarget.pdf", w = 4, h = 2, show.legend = 'right')

g <- drugScore(sePeakRep, markers_lm_sig, "Peaks")
save_plot(g, "./Plots/ActivationDrugScores_PeakReps_ByTarget.pdf", w = 4, h = 2, show.legend = 'right')


### Per-Cell Feature Counts -------------------------------------------------

ccre <- BigBedFile("/Volumes/DannySSD/MULTI_ATAC/20231114_PROTACExperiment2/FeatureSets/encodeCcreCombined.bb")
ccre <- import.bb(ccre)
ccre <- split(ccre, ccre$ucscLabel)

for (elem in names(ccre)) {
  proj_4 <- addFeatureCounts(proj_4, features = ccre[[elem]], name = paste("ENCODE_cCRE",elem,sep="_"))
}
df <- join(df,proj_4@cellColData[,c("Sample","CB", paste("ENCODE_cCRE_",names(ccre),"Ratio", sep=""))] %>% as.data.frame)

g <- plot_grid(plotlist = lapply(names(ccre) %>% rev(), function(elem) {
  plotEmbedding(proj_4, name = paste("ENCODE_cCRE_",elem,"Ratio", sep=""), embedding = "UMAP", size = 1, plotAs = 'points', labelAsFactors = F, labelMeans = F) + 
    theme_DC + theme_UMAP + labs(title = elem, x = "UMAP1", y = "UMAP2")
}), nrow = 1)
save_plot(g, "./Plots/proj4_CCREUMAPs.pdf", w = 13, h = 2.5)

# just distal enhancers
g <- plotEmbedding(proj_4, name = "ENCODE_cCRE_enhDRatio", embedding = "UMAP", size = 1, plotAs = 'points', labelAsFactors = F, labelMeans = F) + 
  theme_DC + theme_UMAP + labs(title = "Fraction of Counts in Distal Enhancers", x = "UMAP1", y = "UMAP2")
save_plot(g, "./Plots/proj4_CCREUMAP_enhD.pdf", w = 2.25, h = 2, show.legend = 'right')
g <- df %>% filter(Type == "Control") %>% 
  ggplot(aes(x = CellType_Major, y = ENCODE_cCRE_enhDRatio, color = CellType_Major, fill = CellType_Major)) + 
  geom_violin(size = 1.5, alpha = 0.3, show.legend = F) +
  geom_boxplot(width = 0.3, fill = NA, outliers = F, size = 1) + 
  scale_color_manual(values = col_cellmaj) + scale_fill_manual(values = col_cellmaj) + 
  theme_DC + labs(x = "Cell Type", y = "Fraction Counts in Distal Enhancers", title = "Control Cells Only") +
  stat_compare_means(method = "t.test", comparisons = list(c("Myeloid","T"),c("Myeloid","B")), label = "p.signif", label.y = c(0.345,0.345), vjust = 0.5, show.legend = F) +
  coord_cartesian(clip = 'off')
save_plot(g, "./Plots/proj4_enhD_violin.pdf", w = 1.5, h = 2, show.legend = 'none')



### Subset Major Cell Types -------------------------------------------------

proj_T <- subsetArchRProject(proj_4, 
                             cells = proj_4$cellNames[proj_4$CellType_Major == "T"],
                             outputDirectory = "proj_T", force = T)
proj_B <- subsetArchRProject(proj_4, 
                             cells = proj_4$cellNames[proj_4$CellType_Major == "B"],
                             outputDirectory = "proj_B", force = T)
proj_M <- subsetArchRProject(proj_4, 
                             cells = proj_4$cellNames[proj_4$CellType_Major == "Myeloid"],
                             outputDirectory = "proj_M", force = T)

proj_T <- dimReduce(proj_T, scale = "Col", default = "NoScale", dims = 2:30)
proj_B <- dimReduce(proj_B, scale = "Col", default = "NoScale", dims = 2:30)
proj_M <- dimReduce(proj_M, scale = "Col", default = "NoScale", dims = 2:30)

plot_grid(plotEmbedding(proj_T, name = "CellType_Minor", embedding = "UMAP", size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + scale_color_manual(values = col_cellmin) + ggtitle("CellType_Minor"),
          plotEmbedding(proj_T, name = "Clusters", embedding = "UMAP", size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + ggtitle("Clusters"),
          plotEmbedding(proj_T, name = "FRIP", embedding = "UMAP", plotAs = 'points', size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + ggtitle("FRIP"),
          plotEmbedding(proj_T, name = "log10(nFrags)", embedding = "UMAP", plotAs = 'points', size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + ggtitle("log10(nFrags)"),
          plotEmbedding(proj_T, name = "log10(Gex_nUMI)", embedding = "UMAP", plotAs = 'points', size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + ggtitle("log10(Gex_nUMI)"),
          # plotEmbedding(proj_T, name = "Phase", embedding = "UMAP", size = 1, labelAsFactors=F, labelMeans=F) + theme_void() + scale_color_manual(values = col_cc) + ggtitle("Phase"),
          plotEmbedding(proj_T, name = "Drug", embedding = "UMAP", size = 1, labelAsFactors=F, labelMeans=F) + theme_void() + scale_color_manual(values = col_drug) + ggtitle("Drug"),
          nrow = 2)

plot_grid(plotEmbedding(proj_M, name = "CellType_Minor", embedding = "UMAP", size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + scale_color_manual(values = col_cellmin) + ggtitle("CellType_Minor"),
          plotEmbedding(proj_M, name = "Clusters", embedding = "UMAP", size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + ggtitle("Clusters"),
          plotEmbedding(proj_M, name = "FRIP", embedding = "UMAP", plotAs = 'points', size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + ggtitle("FRIP"),
          plotEmbedding(proj_M, name = "log10(nFrags)", embedding = "UMAP", plotAs = 'points', size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + ggtitle("log10(nFrags)"),
          plotEmbedding(proj_M, name = "log10(Gex_nUMI)", embedding = "UMAP", plotAs = 'points', size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + ggtitle("log10(Gex_nUMI)"),
          # plotEmbedding(proj_M, name = "Phase", embedding = "UMAP", size = 1, labelAsFactors=F, labelMeans=F) + theme_void() + scale_color_manual(values = col_cc) + ggtitle("Phase"),
          plotEmbedding(proj_M, name = "Drug", embedding = "UMAP", size = 1, labelAsFactors=F, labelMeans=F) + theme_void() + scale_color_manual(values = col_drug) + ggtitle("Drug"),
          nrow = 2)

plot_grid(plotEmbedding(proj_B, name = "CellType_Minor", embedding = "UMAP", size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + scale_color_manual(values = col_cellmin) + ggtitle("CellType_Minor"),
          plotEmbedding(proj_B, name = "Clusters", embedding = "UMAP", size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + ggtitle("Clusters"),
          plotEmbedding(proj_B, name = "FRIP", embedding = "UMAP", plotAs = 'points', size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + ggtitle("FRIP"),
          plotEmbedding(proj_B, name = "log10(nFrags)", embedding = "UMAP", plotAs = 'points', size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + ggtitle("log10(nFrags)"),
          plotEmbedding(proj_B, name = "log10(Gex_nUMI)", embedding = "UMAP", plotAs = 'points', size = 0.5, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + ggtitle("log10(Gex_nUMI)"),
          # plotEmbedding(proj_B, name = "Phase", embedding = "UMAP", size = 1, labelAsFactors=F, labelMeans=F) + theme_void() + scale_color_manual(values = col_cc) + ggtitle("Phase"),
          plotEmbedding(proj_B, name = "Drug", embedding = "UMAP", size = 1, labelAsFactors=F, labelMeans=F) + theme_void() + scale_color_manual(values = col_drug) + ggtitle("Drug"),
          nrow = 2)


## LSI Component Correlations Between Treatments

plotLSIHeatmap("individual", replicates = F, celltypes = c("T","Myeloid"))
plotLSIHeatmap("individual", replicates = T, celltypes = c("T","Myeloid"))
plotLSIHeatmap("full", replicates = F, celltypes = c("T","Myeloid"), dims = c(2:3,5:30))

lapply(c("T","Myeloid"), function(c) {
  plotLSIHeatmap("individual", replicates = F, celltypes = c, title = paste("Treatments - ", c, " Cells",sep="")) %>% .[[1]] %>%
    pdfPlot(., width = 5, height = 3, file = paste("./Plots/LSI_Heatmaps_DrugDose_",c, ".pdf",sep=""))
  # pdfPlot(., width = 4, height = 2.5, file = paste("./Plots/LSI_Heatmaps_DrugDose_",c, ".pdf",sep=""))
  # save_plot(., paste("./Plots/LSI_Heatmaps_DrugDose_",c, ".pdf",sep=""), w = 4, h = 4, show.legend = 'none')
  plotLSIHeatmap("individual", replicates = T, celltypes = c, title = paste("Control Replicates - ", c, " Cells",sep="")) %>% .[[1]] %>%
    pdfPlot(., width = 5, height = 3, file = paste("./Plots/LSI_Heatmaps_DMSOReplicate_",c, ".pdf",sep=""))
  # pdfPlot(., width = 4, height = 2.5, file = paste("./Plots/LSI_Heatmaps_DMSOReplicate_",c, ".pdf",sep=""))
  # save_plot(., paste("./Plots/LSI_Heatmaps_DMSOReplicate_",c, ".pdf",sep=""), w = 4, h = 4, show.legend = 'none')
})


plotLSIHeatmap("individual", replicates = T, celltypes = "T", title = "Control Replicates - T Cells", fontsize = 5, filename = "./Plots/LSI_Heatmaps_DMSOReplicate_T_v2.pdf", width = 4, height = 2.5)
plotLSIHeatmap("individual", replicates = T, celltypes = "Myeloid", title = "Control Replicates - Myeloid Cells", fontsize = 5, filename = "./Plots/LSI_Heatmaps_DMSOReplicate_Myeloid_v2.pdf", width = 4, height = 2.5)

plotLSIHeatmap("individual", replicates = F, celltypes = "T", title = "Treatments - T Cells", fontsize = 5,  filename = "./Plots/LSI_Heatmaps_Treatments_T.pdf", width = 4, height = 2.5)
plotLSIHeatmap("individual", replicates = F, celltypes = "Myeloid", title = "Treatments - Myeloid Cells", fontsize = 5,  filename = "./Plots/LSI_Heatmaps_Treatments_Myeloid.pdf", width = 4, height = 2.5)

plotLSIHeatmap("individual", replicates = F, celltypes = "T", doses = names(col_dose)[-3:-4], title = "Highest Dose - T Cells", fontsize = 5) # , filename = "./Plots/LSI_Heatmaps_Treatments_T.pdf", width = 4, height = 2.5)
plotLSIHeatmap("individual", replicates = F, celltypes = "Myeloid", doses = names(col_dose)[-3:-4], title = "Highest Dose - Myeloid Cells", fontsize = 5) # , filename = "./Plots/LSI_Heatmaps_Treatments_T.pdf", width = 4, height = 2.5)

plotLSIHeatmap("individual", replicates = F, celltypes = "T", doses = names(col_dose)[-4:-5], title = "Lowest Dose - T Cells", fontsize = 5) # , filename = "./Plots/LSI_Heatmaps_Treatments_T.pdf", width = 4, height = 2.5)
plotLSIHeatmap("individual", replicates = F, celltypes = "Myeloid", doses = names(col_dose)[-4:-5], title = "Lowest Dose - Myeloid Cells", fontsize = 5) # , filename = "./Plots/LSI_Heatmaps_Treatments_T.pdf", width = 4, height = 2.5)


## Control Replicates

# set.seed(1)
# proj_4$Shuffled <- df[match(df$Cell,proj_4$cellNames),] %>% group_by(Drug,Dose,CellType_Major) %>% mutate(Shuffled = sample(Replicate, replace = F)) %>% .$Shuffled
# proj_4$Shuffled <- df[match(df$Cell,proj_4$cellNames),] %>% group_by(Drug,Dose,CellType_Major) %>% mutate(Shuffled = sample(Replicate, replace = T)) %>% .$Shuffled
# 
# df$Shuffled <- proj_4$Shuffled


set.seed(1)
proj_4$Shuffled <- df[match(df$Cell,proj_4$cellNames),] %>% group_by(Drug,Dose,CellType_Major) %>% mutate(Shuffled = sample(Replicate, replace = F)) %>% .$Shuffled
proj_4$Shuffled <- df[match(df$Cell,proj_4$cellNames),] %>% group_by(Drug,Dose,CellType_Major) %>% mutate(Shuffled = sample(Replicate, replace = T)) %>% .$Shuffled

set.seed(1)
proj_4$Shuffled <- proj_4@cellColData %>% as.data.frame %>% group_by(CellType_Major, Drug, Dose) %>% mutate(Shuffled = sample(Replicate, replace = F)) %>% .$Shuffled

# df$Shuffled <- proj_4$Shuffled

# Gene Score, average SD per Dose
g <- cbind(seu@meta.data %>% select(CellType_Major, Drug, Dose, Replicate), t(seu@assays$GeneScore$counts[seu@assays$GeneScore@var.features,])) %>%
  dupControl %>% filter(Dose != "DMSO(-)", CellType_Major != "B") %>%
  group_by(CellType_Major, Drug, Dose, Replicate) %>% summarize_all(mean) %>% select(-Replicate) %>% 
  pivot_longer(c(-CellType_Major,-Drug,-Dose), names_to = "Gene", values_to = "Mean") %>% 
  group_by(CellType_Major, Drug, Dose, Gene) %>% summarize(SD = sd(Mean)) %>%
  group_by(CellType_Major, Drug, Gene) %>% mutate(SD = SD / SD[Dose == "DMSO(+)"], Log2 = log2(SD)) %>%
  filter(is.finite(Log2)) %>%
  group_by(Dose) %>% summarize(mean = mean(SD), error = sd(SD)) %>%
  ggplot(aes(x = Dose, y = mean, fill = Dose)) +
  geom_col() +
  scale_color_manual(values = col_dose) + scale_fill_manual(values = col_dose) + 
  theme_DC + labs(x = "Dose", y = "Inter-Replicate Variability", title = "2000 Most-Variable Gene Scores")
save_plot(g, "./Plots/DoseVariability_GeneScores_byDose.pdf", w = 2, h = 2, show.legend = 'none')

# Gene Score, SD per Drug+Dose
g <- cbind(seu@meta.data %>% select(CellType_Major, Drug, Dose, Replicate), t(seu@assays$GeneScore$counts[seu@assays$GeneScore@var.features,])) %>%
  dupControl %>% filter(Dose != "DMSO(-)", CellType_Major != "B") %>%
  group_by(CellType_Major, Drug, Dose, Replicate) %>% summarize_all(mean) %>% select(-Replicate) %>% 
  pivot_longer(c(-CellType_Major,-Drug,-Dose), names_to = "Gene", values_to = "Mean") %>% 
  group_by(CellType_Major, Drug, Dose, Gene) %>% summarize(SD = sd(Mean)) %>%
  group_by(CellType_Major, Drug, Gene) %>% mutate(SD = SD / SD[Dose == "DMSO(+)"], Log2 = log2(SD)) %>%
  filter(is.finite(Log2)) %>%
  # group_by(Dose) %>% summarize(mean = mean(SD), error = sd(SD)) %>%
  ggplot(aes(x = as.numeric(Dose), y = SD, color = Drug, fill = Drug)) +
  # geom_smooth(method = 'lm', alpha = 0.2) +
  geom_smooth(method = 'loess', span = 1.5) + 
  scale_color_manual(values = col_drug) + scale_fill_manual(values = col_drug) + 
  theme_DC + labs(x = "Dose", y = "Inter-Replicate Variability", title = "2000 Most-Variable Gene Scores") + theme(axis.text.x = element_blank())
save_plot(g, "./Plots/DoseVariability_GeneScores_byDrugDose.pdf", w = 2, h = 2, show.legend = 'right')
save_plot(g, "./Plots/DoseVariability_GeneScores_byDrugDose_loess.pdf", w = 2, h = 2, show.legend = 'right')

# RNA, average SD per Dose
g <- cbind(seu@meta.data %>% select(CellType_Major, Drug, Dose, Replicate), t(seu@assays$RNA$counts[seu@assays$RNA@var.features,])) %>%
  dupControl %>% filter(Dose != "DMSO(-)", CellType_Major != "B") %>%
  group_by(CellType_Major, Drug, Dose, Replicate) %>% summarize_all(mean) %>% select(-Replicate) %>% 
  pivot_longer(c(-CellType_Major,-Drug,-Dose), names_to = "Gene", values_to = "Mean") %>% 
  group_by(CellType_Major, Drug, Dose, Gene) %>% summarize(SD = sd(Mean)) %>%
  group_by(CellType_Major, Drug, Gene) %>% mutate(SD = SD / SD[Dose == "DMSO(+)"], Log2 = log2(SD)) %>%
  filter(is.finite(Log2)) %>%
  group_by(Dose) %>% summarize(mean = mean(SD), error = sd(SD)) %>%
  ggplot(aes(x = Dose, y = mean, fill = Dose)) +
  geom_col() +
  scale_color_manual(values = col_dose) + scale_fill_manual(values = col_dose) + 
  theme_DC + labs(x = "Dose", y = "Inter-Replicate Variability", title = "2000 Most-Variable Genes")
save_plot(g, "./Plots/DoseVariability_RNA_byDose.pdf", w = 2, h = 2, show.legend = 'none')

# RNA, SD per Drug+Dose
g <- cbind(seu@meta.data %>% select(CellType_Major, Drug, Dose, Replicate), t(seu@assays$RNA$counts[seu@assays$RNA@var.features,])) %>%
  dupControl %>% filter(Dose != "DMSO(-)", CellType_Major != "B") %>%
  group_by(CellType_Major, Drug, Dose, Replicate) %>% summarize_all(mean) %>% select(-Replicate) %>% 
  pivot_longer(c(-CellType_Major,-Drug,-Dose), names_to = "Gene", values_to = "Mean") %>% 
  group_by(CellType_Major, Drug, Dose, Gene) %>% summarize(SD = sd(Mean)) %>%
  group_by(CellType_Major, Drug, Gene) %>% mutate(SD = SD / SD[Dose == "DMSO(+)"], Log2 = log2(SD)) %>%
  filter(is.finite(Log2)) %>%
  # group_by(Dose) %>% summarize(mean = mean(SD), error = sd(SD)) %>%
  ggplot(aes(x = as.numeric(Dose), y = SD, color = Drug, fill = Drug)) +
  # geom_smooth(method = 'lm', alpha = 0.2) +
  geom_smooth(method = 'loess', span = 1.5) + 
  scale_color_manual(values = col_drug) + scale_fill_manual(values = col_drug) + 
  theme_DC + labs(x = "Dose", y = "Inter-Replicate Variability", title = "2000 Most-Variable Genes") + theme(axis.text.x = element_blank())
save_plot(g, "./Plots/DoseVariability_RNA_byDrugDose.pdf", w = 2, h = 2, show.legend = 'right')
save_plot(g, "./Plots/DoseVariability_RNA_byDrugDose_loess.pdf", w = 2, h = 2, show.legend = 'right')



### Vignettes ---------------------------------------------------------------


# Type I Interferon Response ----------------------------------------------

lapply(c(RIGI="DDX58",MDA5="IFIH1",LGP2="DHX58", # RIG-I-Like Receptors (RLRs; dsRNA-binding)
         "MAVS","TBK1","IKBKE", # downstream of RLRs
         "IRF3","IRF7","NFKB1", # transcribe IFNa/b and antiviral genes'
         "IFNB1","IFNAR1","IFNAR2", # IFNb and type I IFN receptors
         "STAT1","STAT2","IRF9", # transcribe ISGs (ISGF3)
         "OAS1","OAS2","OAS3", # oligoadenylate synthetases
         "MX1","MX2",PKR="EIF2AK2",
         "ISG15","IFIT1","IFITM1",
         "IFI27","IFI44","IFI6",
         "USP18","ADAR","METTL3"
# )[1:15],
)[16:30],         
function(x) {
  plotEmbedding(proj_4, colorBy = "GeneExpressionMatrix", name = x, plotAs = 'points', size = 1) + theme_void() + ggtitle(x) + labs(color="GeneExp")
}) %>% 
  plot_grid(plotlist=.,ncol=3)

genes <- c(RIGI="DDX58",MDA5="IFIH1",LGP2="DHX58", # dsRNA-sensing
           "CGAS","ZBP1","AIM2", # DNA-sensing
           "IFNB1","STAT1","STAT2", # signaling pathway
           "ISG15","ADAR","IFI6", # ISGs
           "OAS3","MX1",PKR="EIF2AK2") # RNA-responding

tmp <- lapply(genes, function(x) {
  plotEmbedding(proj_4, colorBy = "GeneExpressionMatrix", name = x, plotAs = 'points', size = 1) + theme_DC + theme_UMAP + 
    labs(color="GeneExp", x = "UMAP1", y = "UMAP2", title = fixGene(x, new))
})

g <- plot_grid(plotlist=tmp, ncol = 3)
save_plot(g, "./Plots/InterferonGenes.pdf", w = 6, h = 8, show.legend = "none")

g <- plot_grid(plotlist=lapply(tmp[12:13], function(x) x + theme(legend.position = "bottom")), ncol = 2)
save_plot(g5, "./Plots/InterferonGenes_select_v2.pdf", w = 3, h = 2, show.legend = "none")

markers_lm_sig %>% filter(name %in% genes, Type == "Up", Target == "SWI/SNF")
markers_lm_sig %>% filter(Type == "Up", Target == "SWI/SNF") %>% {genes %in% .$name}
markers_lm_sig %>% filter(Type == "Up", Target == "SWI/SNF") %>% {fixGene(genes, new) %in% .$name}
markers_lm %>% filter(P < 0.05, Coef > 0, Drug %in% c("AU-15330","BRM014"), MarkerType != "Peaks") %>% {genes %in% .$Gene}
markers_lm %>% filter(P < 0.05, Coef > 0, Drug %in% c("AU-15330","BRM014"), MarkerType == "Genes") %>% {fixGene(genes, new) %in% .$Gene}

markers_lm %>% filter(P < 0.01, Coef > 0.05, Drug %in% c("AU-15330","BRM014"), MarkerType == "Genes") %>% group_by(Gene) %>% summarize(Count = n()) %>%
  filter(Count >= 4)


rbind(data.frame(gsea_lm_rna, Type = "RNA"),
      data.frame(gsea_lm_score, Type = "Score")) %>%
  # g <- rbind(data.frame(gsea_cor_rna2, Type = "RNA"),
  #            data.frame(gsea_cor_score2, Type = "Score")) %>%
  filter(Drug %in% c("BRM014","AU-15330")) %>% 
  filter(padj < 0.01, NES > 1) %>%
  filter(Collection %in% c("Hallmark","GO_BP")) %>%
  filter(pathway %in% grep("interferon|viral|innate", pathway, value = T, ignore.case = T)) %>%
  .$leadingEdge %>% unlist %>% table %>% sort(decreasing = T) %>% head(n = 30)

gsea_lm_score %>%
  filter(padj < 0.01, NES > 0) %>%
  filter(Drug %in% c("BRM014","AU-15330")) %>%
  filter(Collection == "GO_BP") %>%
  # filter(Collection == "Reactome") %>%
  arrange(padj) %>%
  group_by(pathway) %>% summarize(padj = prod(padj), NES = mean(NES)) %>% 
  # filter(NES > 2) %>%
  arrange(padj) %>%
  .$pathway -> terms1

gsea_lm_rna %>%
  filter(padj < 0.01, NES > 0) %>%
  filter(Drug %in% c("BRM014","AU-15330")) %>%
  filter(Collection == "GO_BP") %>%
  # filter(Collection == "Reactome") %>%
  arrange(padj) %>%
  group_by(pathway) %>% summarize(padj = prod(padj), NES = mean(NES)) %>% 
  # filter(NES > 2) %>%
  arrange(padj) %>%
  .$pathway -> terms2

# ggVennDiagram::ggVennDiagram(list(terms1, terms2))

g <- rbind(data.frame(gsea_lm_rna, Type = "RNA"),
           data.frame(gsea_lm_score, Type = "Score")) %>%
  filter(Drug %in% c("BRM014","AU-15330")) %>% 
  filter(padj < 0.01, NES > 0) %>% 
  filter(pathway %in% c(terms1, terms2)) %>%
  mutate(Group = case_when(
    pathway %in% intersect(terms1,terms2) ~ "Both",
    pathway %in% setdiff(terms1,terms2) ~ "Score Only",
    pathway %in% setdiff(terms2,terms1) ~ "RNA Only")) %>%
  filter(Group != "RNA Only") %>%
  group_by(pathway, Group, Type) %>% summarize(NES = mean(NES), padj = mean(padj)) %>% mutate(sum = prod(-log10(padj))) %>%
  group_by(Group, Type) %>% slice_max(order_by = sum, n = 10) %>%
  ggplot(aes(y = reorder(pathway, sum), x = -log10(padj), fill = NES)) + geom_col() +  facet_grid(Group~Type, scales = "free_y", space = 'free_y') +
  theme_DC + labs(y = "Gene Set", title = "Up-Regulated by SWI/SNF Perturbation") + scale_fill_viridis_c(option = "D")
save_plot(g, "./Plots/SWISNF_Responsive_GeneSets.pdf", w = 6, h = 3, show.legend = 'right')






sets$Hallmarks.HALLMARK_INTERFERON_ALPHA_RESPONSE

markers_lm_sig %>% filter(MarkerType != "Peaks") %>% group_by(CellType, MarkerType, Treatment, Type) %>%
  summarize(Fraction = sum(name %in% sets$Hallmarks.HALLMARK_INTERFERON_ALPHA_RESPONSE)/length(sets$Hallmarks.HALLMARK_INTERFERON_ALPHA_RESPONSE)) %>%
  mutate(Fraction = case_when(Type == "Down" ~ -Fraction, TRUE ~ Fraction)) %>%
  ggplot(aes(x = Treatment, y = Fraction, fill = Type)) + geom_col(position = 'identity') + 
  facet_grid(MarkerType~CellType) + theme_DC

g <- markers_lm_sig %>% filter(MarkerType == "Scores",  Treatment != "Activation") %>% group_by(CellType, MarkerType, Treatment, Type) %>%
  summarize(Fraction = sum(name %in% sets$Hallmarks.HALLMARK_INTERFERON_ALPHA_RESPONSE)/length(sets$Hallmarks.HALLMARK_INTERFERON_ALPHA_RESPONSE)) %>%
  mutate(Fraction = case_when(Type == "Down" ~ -Fraction, TRUE ~ Fraction)) %>%
  ggplot(aes(x = Treatment, y = Fraction, alpha = Type, fill = Treatment)) + geom_col(position = 'identity') + geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') +
  facet_grid(MarkerType~CellType) + theme_DC + scale_fill_manual(values = col_drug) + scale_alpha_manual(values = c(0.5,1)) + 
  labs(title = "Interferon Alpha Response")
save_plot(g, "./Plots/InterferonGenes_FractionDE_Scores.pdf")


g <- markers_lm_sig %>% filter(MarkerType == "Scores", Treatment != "Activation", Type == "Up") %>% 
  group_by(CellType, MarkerType, Treatment, Type) %>%
  summarize(Count = sum(name %in% sets$Hallmarks.HALLMARK_INTERFERON_ALPHA_RESPONSE), 
            Fraction = Count/length(unique(sets$Hallmarks.HALLMARK_INTERFERON_ALPHA_RESPONSE))) %>%
  ggplot(aes(x = Treatment, y = Count, fill = Treatment)) + geom_col(position = 'identity') +
  facet_grid(MarkerType~CellType) + theme_DC + scale_fill_manual(values = col_drug) +
  labs(y = "Interferon Alpha Response Genes Upregulated", x = "Drug") + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0))

g <- markers_lm_sig %>% 
  # filter(MarkerType == "Scores", Treatment != "Activation", Type == "Up") %>% group_by(CellType, MarkerType, Treatment, Type) %>%
  filter(MarkerType != "Peaks", Treatment != "Activation", Type == "Up") %>% group_by(Treatment) %>%
  # filter(MarkerType == "Scores", Treatment != "Activation", Type == "Up") %>% group_by(Treatment) %>%
  filter(!duplicated(name)) %>%
  summarize(Count = sum(name %in% sets$Hallmarks.HALLMARK_INTERFERON_ALPHA_RESPONSE),
            Fraction = Count/length(unique(sets$Hallmarks.HALLMARK_INTERFERON_ALPHA_RESPONSE))) %>%
  ggplot(aes(x = Treatment, y = Fraction, fill = Treatment)) + geom_col(position = 'identity') +
  theme_DC + scale_fill_manual(values = col_drug) +
  # labs(y = "Interferon Alpha Response Genes Upregulated", x = "Drug") +
  labs(y = "Fraction Up-Regulated", x = "Drug", title = "Type I IFN Genes") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0))
save_plot(g, "./Plots/InterferonGenes_FractionUp.pdf", w = 1.5, h = 2.5, show.legend = 'none')

g <- markers_lm_sig %>% 
  filter(MarkerType != "Peaks", Treatment != "Activation", Type == "Up") %>% group_by(Treatment) %>%
  filter(!duplicated(name)) %>%
  summarize(Count = sum(name %in% sets$Hallmarks.HALLMARK_INTERFERON_ALPHA_RESPONSE),
            Fraction = Count/length(unique(sets$Hallmarks.HALLMARK_INTERFERON_ALPHA_RESPONSE))) %>%
  ggplot(aes(x = Fraction, y = Treatment, fill = Treatment)) + geom_col(position = 'identity') +
  theme_DC + scale_fill_manual(values = col_drug) +
  scale_y_discrete(limits = rev) +
  labs(x = "Fraction Up-Regulated", y = "Drug", title = "Type I IFN Genes")
save_plot(g, "./Plots/InterferonGenes_FractionUp_v2.pdf", w = 3, h = 2, show.legend = 'none')



# NF-kB Signalling --------------------------------------------------------

# look at MS177 peaks specifically
resT <- plotHeatmap(lm_atac, assay = "Peaks", cellType = "T", drugsSig = drug[c(1,3)], drugsPlot = drug[c(1:3)], p = 0.01, log2fc = 1, cutCol = 7, returnRes = T, filename = "./Plots/RepTest_Peaks_Marker_Heatmap_T_select_v3.pdf", width = 6, height = 3)


markeridx <- which((proj_4@peakSet %>% as.data.frame(row.names = NULL) %>% mutate(name = paste(seqnames, start, end, sep='-')) %>% .$name) %in%
                     (names(resT$cut) %>% gsub(":","-",.)))

assays(markerPeaks_DrugDoseT)$Pval[,] <- 0
assays(markerPeaks_DrugDoseT)$Pval[markeridx,] <- resT$cut

diffmotifsT <- lapply(1:max(assays(markerPeaks_DrugDoseT)$Pval), function(c) {
  tf <- peakAnnoEnrichment(
    seMarker = markerPeaks_DrugDoseT[,1],
    ArchRProj = proj_4,
    peakAnnotation = "Motif",
    cutOff = paste("Pval ==", c)
  )
  tf <- data.frame(TF = rownames(tf),
                   mlog10Padj = assay(tf)[,1])
  tf <- tf[order(tf$mlog10Padj, decreasing = T),]
  tf$Rank <- seq_len(nrow(tf))
  tf$Cluster <- c
  tf
}) %>% do.call(rbind, .)

diffmotifsT <- diffmotifsT %>% 
  mutate(RNA = rowSums(seu@assays$RNA@counts[,seu$CellType_Major == "T"])[gsub("_.*","",TF)], # count RNA UMIs by Cell Type
         Cluster = factor(Cluster, levels = resT$cut[resT$hc$order] %>% unique, labels = resT$cut[resT$hc$order] %>% unique %>% as.character),
         TF = gsub("_.*","",TF))

g <- diffmotifsT %>% 
  filter(mlog10Padj >= -log10(0.05)) %>% 
  filter(RNA >= 1) %>%
  filter(TF %ni% grep("ZNF",TF,value=T)) %>%
  group_by(Cluster) %>% slice_max(order_by = -Rank, n = nTerm, with_ties = F) %>%
  mutate(TF = tidytext::reorder_within(TF, mlog10Padj, within = Cluster)) %>%
  ggplot(aes(x = TF, y = mlog10Padj, size = log10(RNA), color = Cluster)) + 
  geom_col(mapping=aes(fill = Cluster, alpha = log10(RNA)), width = 0.5, color = NA) +
  scale_color_manual(values = col_clust) + scale_fill_manual(values = col_clust) + scale_size_continuous(range = c(1,3)) + scale_alpha_continuous(range = c(0.3,1)) +
  facet_grid(.~Cluster, space = 'free_x', scales = 'free_x', drop = T) + labs(y = "-log10(pAdj)", x = NULL) +
  tidytext::scale_x_reordered() + coord_cartesian(clip = 'off') + expand_limits(y = 0) +
  theme_DC + theme(axis.text.x = element_text(angle = -90, hjust=0), strip.text = element_blank())
save_plot(g, "./Plots/RepTest_Peaks_Marker_Heatmap_T_DiffMotifs.pdf", w = 6, h = 2.5, show.legend = "right")

diffmotifsT %>%
  filter(Cluster %ni% c(2,3,4,7)) %>% 
  filter(mlog10Padj > -log10(0.01)) %>%
  pivot_wider(id_cols = "Cluster", names_from = "TF", values_from = mlog10Padj, values_fill = 0) %>% column_to_rownames(var = "Cluster") %>% #{.[,c(2,1)]} %>% #{.[,c(2,1)]} %>%
  pheatmap(scale = 'none', display_numbers = F, border_color = NA, treeheight_col = 0, treeheight_row = 0, cluster_cols = F, fontsize = 10, show_rownames = F, #main = "MS177-Responsive Peaks", 
           color = paletteer::paletteer_c("pals::ocean.dense", 100, direction = -1),
           annotation_row = data.frame(Cluster = levels(diffmotifsT$Cluster), rowname = levels(diffmotifsT$Cluster)) %>% column_to_rownames(), #annotation_names_row = F,
           annotation_colors = list(Cluster = col_clust), annotation_legend = F, 
           breaks = seq(0,max(abs(.)),length = 100),
           filename = if (plot) { "./Plots/MS177_Responsive_Motifs_v3.pdf" } else { NA }, width = 5, h = 3)



## What NF-kB Related Gene Sets are enriched in up-/down-regulated genes?

g <- rbind(data.frame(gsea_lm_rna %>% filter(pathway %in% grep("NF.kB|NFkB", pathway, value= T, ignore.case=T)) %>% filter(padj < 0.05), Type = "RNA"),
           data.frame(gsea_lm_score %>% filter(pathway %in% grep("NF.kB|NFkB", pathway, value= T, ignore.case=T)) %>% filter(padj < 0.05), Type = "Score")) %>%
  filter(Drug == "MS177") %>% filter(Collection != "Immune") %>% filter(CellType == "T") %>%
  filter(pathway != "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>%
  mutate(pathway = gsub("REACTOME_|HALLMARK_","",pathway)) %>%
  group_by(pathway) %>% mutate(sum = sum(-log10(padj))) %>%
  ggplot(aes(y = reorder(pathway, sum), x = -log10(padj), fill = Type)) + geom_col(position = 'dodge', show.legend = F) + facet_wrap(~Type, axes = 'all', axis.labels = 'margins') + 
  theme_DC + labs(y = "Reactome Set", title = "MS177-Treated T Cells")
save_plot(g, "./Plots/MS177_Responsive_GeneSets.pdf", w = 6, h = 1.5, show.legend = 'right')

g <- rbind(data.frame(gsea_lm_rna %>% filter(pathway %in% grep("NF.kB|NFkB", pathway, value= T, ignore.case=T)) %>% filter(padj < 0.05), Type = "RNA"),
           data.frame(gsea_lm_score %>% filter(pathway %in% grep("NF.kB|NFkB", pathway, value= T, ignore.case=T)) %>% filter(padj < 0.05), Type = "Score")) %>%
  filter(Drug %in% c("DMSO(-)","MS177")) %>% filter(Collection != "Immune") %>% filter(CellType == "T") %>%
  # filter(pathway != "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>%
  mutate(pathway = gsub("REACTOME_|HALLMARK_","",pathway),
         CellType = factor(CellType, levels = c("T","Myeloid")),
         Drug = factor(Drug, levels = c("DMSO(-)","MS177"), labels = c("Activation","MS177"))) %>%
  group_by(pathway) %>% mutate(sum = sum(-log10(padj))) %>%
  ggplot(aes(y = reorder(pathway, sum), x = -log10(padj), fill = Type)) + geom_col(position = 'dodge', show.legend = F) + facet_wrap(Drug~paste(CellType,Type), axes = 'all', axis.labels = 'margins', nrow = 1) + 
  theme_DC + labs(y = "Reactome Set", title = "MS177-Treated T Cells")
# save_plot(g, "./Plots/MS177_Responsive_GeneSets.pdf", w = 6, h = 1.5, show.legend = 'right')


## NF-kB Footprinting

kB <- motifPositions[getFeatures(proj_4, select = paste(c("REL","NFKB","HIVEP"), collapse="|"), useMatrix = "MotifMatrix") %>% grep(pattern = "z:", value = T) %>% gsub("z:","",.)]
kB <- Reduce(GenomicRanges::union, kB)

seKB <- getFootprints(
  ArchRProj = proj_4, 
  positions = GRangesList(kB = kB), 
  groupBy = "CellDrugDose2",
  flank = 2000,
  minCells = 10)

g <- ggFootprint(seFoot = seKB, name = "kB", flank = 250, 
                 subset = "Drug == 'MS177' & CellType != 'B'",
                 flankNorm = 50, normMethod = "Subtract", smoothWindow = 25, dupControl = T,
                 groupings = c("CellType","Drug","Dose"), returnDF = T)
g <- ggplot(g, aes(x, mean, color = Dose)) + geom_line() + 
  facet_grid(~CellType, scales = 'free_y') + scale_color_manual(values = col_dose2) + 
  xlim(c(-500,500)) + labs(x = "Distance From Motif Center", y = "Bias-Corrected Insertion Profile", title = "kB Sites - MS177") + theme_DC
save_plot(g, "./Plots/Footprints_kBv3.pdf", w = 3.5, h = 2, show.legend = 'right')  



# Histone Accessibility ---------------------------------------------------

g <- seu[,seu$Dose %in% c("DMSO(-)","DMSO(+)","1µM")] %>%
  CovPlot("chr6-26050000-26300000", group.by = "Drug", up = 0, down = 0, col_drug, window = 1000, ratio = c(6,1.5,0.5))
save_plot(g, "./Plots/CovPlot_HIST1.pdf", 6, 5, show.legend = 'none')

g <- seu[,seu$CellType_Major == "T" & seu$Dose %in% c("DMSO(-)","DMSO(+)","1µM")] %>%
  CovPlot("chr6-26050000-26300000", group.by = "Drug", up = 0, down = 0, col_drug, window = 1000, ratio = c(6,1.5,0.5))
save_plot(g, "./Plots/CovPlot_HIST1_T.pdf", 6, 5, show.legend = 'none')

g <- seu[,seu$CellType_Major == "Myeloid" & seu$Dose %in% c("DMSO(-)","DMSO(+)","1µM")] %>% # & seu$Drug %ni% c("MS177","EPZ6438","GNE-781")
  CovPlot("chr6-26050000-26300000", group.by = "Drug", up = 0, down = 0, col_drug, window = 1000, ratio = c(6,1.5,0.5))
save_plot(g, "./Plots/CovPlot_HIST1_M.pdf", 6, 5, show.legend = 'none')

    
## Histone genes by subtype
lm_score %>% .[c("p","coef")] %>%
  lapply(function(x) {
    x %>% ungroup %>% {.[,c("CellType","Drug",grep("^(HIST)|^H(1|2|3|4)(-|A|B|C)",colnames(scorecor),value=T))]}
  }) %>%
  { .$p[,-1:-2] <- -log10(.$p[,-1:-2]) * sign(.$coef[,-1:-2]); .$p } %>%
  pivot_longer(cols = 3:ncol(.), names_to = "Gene", values_to = "CorrP") %>%
  filter(CellType != "B") %>%
  mutate(
    CorrP = case_when(abs(CorrP) > -log10(0.01) ~ CorrP, TRUE ~ NA),
    # CorrP = case_when(CorrP < 0.01 ~ CorrP, TRUE ~ NA),
    Drug = factor(Drug, levels = levels(df$Drug) %>% rev),
    Gene = factor(Gene, 
                  getGeneAnnotation(proj_4)$genes %>% subset(symbol %in% grep("^(HIST)|^H(1|2|3|4)(-|A|B|C)",symbol, value = T)) %>% .$symbol,
                  getGeneAnnotation(proj_4)$genes %>% subset(symbol %in% grep("^(HIST)|^H(1|2|3|4)(-|A|B|C)",symbol, value = T)) %>% .$symbol %>% fixGene(new)),
    Type = gsub("(C|-|W).*","",Gene),
  ) %>%
  {
    ggplot(.,aes(x = Gene, y = Drug, fill = CorrP)) + geom_tile() + 
      centered_gradient("pals::ocean.curl", midquant = centerQuant(.$CorrP), direction = 1) +
      theme_DC + theme(axis.text.x = element_text(angle=-90, hjust = 0, vjust = 0.5), strip.clip = 'off') + 
      facet_grid(CellType~Type, scales = 'free_x', space = 'free_x') # axes = 'all_x', axis.labels = "margins")
  }

## Histone genes by locus
g <- lm_score %>%
  # g <- lm_rna %>%
  .[c("p","coef")] %>%
  lapply(function(x) {
    x %>% ungroup %>% {.[,c("CellType","Drug",grep("^(HIST)|^H(1|2|3|4)(-|A|B|C)",colnames(scorecor),value=T))]}
  }) %>%
  { .$p[,-1:-2] <- -log10(.$p[,-1:-2]) * sign(.$coef[,-1:-2]); .$p } %>%
  pivot_longer(cols = 3:ncol(.), names_to = "Gene", values_to = "CorrP") %>%
  filter(CellType != "B") %>%
  mutate(
    CorrP = case_when(abs(CorrP) > -log10(0.01) ~ CorrP, TRUE ~ NA),
    Drug = factor(Drug, levels = levels(df$Drug) %>% rev),
    Chr = getGeneAnnotation(proj_4)$genes %>% {.[match(Gene, getGeneAnnotation(proj_4)$genes$symbol)]} %>% seqnames %>% as.vector %>% gsub("chr","",.),
    Chr = factor(Chr, c(1:22,"X","Y")),
    Gene = factor(Gene, 
                  getGeneAnnotation(proj_4)$genes %>% subset(symbol %in% grep("^(HIST)|^H(1|2|3|4)(-|A|B|C)",symbol, value = T)) %>% .$symbol,
                  getGeneAnnotation(proj_4)$genes %>% subset(symbol %in% grep("^(HIST)|^H(1|2|3|4)(-|A|B|C)",symbol, value = T)) %>% .$symbol %>% fixGene(new)),
  ) %>%
  {
    ggplot(.,aes(x = Gene, y = Drug, fill = CorrP)) + geom_tile() + 
      centered_gradient("pals::ocean.curl", midquant = centerQuant(.$CorrP), direction = 1) +
      theme_DC + theme(axis.text.x = element_text(angle=-90, hjust = 0, vjust = 0.5), strip.clip = 'off') + 
      facet_grid(CellType~Chr, scales = 'free_x', space = 'free_x') # axes = 'all_x', axis.labels = "margins")
  }
save_plot(g, "./Plots/RepTest_Score_Hist.pdf", w = 8, h = 2.5, show.legend = 'right')
# save_plot(g, "./Plots/RepTest_RNA_Hist.pdf", w = 8, h = 2.5, show.legend = 'right')



proj_4 <- addFeatureCounts(proj_4, 
                           c(GRanges("chr6", IRanges(start = 26050000, end = 26300000)),
                             GRanges("chr7", IRanges(start = 76455305, end = 76455491))),
                           name = "HIST1")


plotEmbedding(proj_4, name = "log10(HIST1Ratio)", plotAs = 'points') + theme_DC + theme_UMAP
plotEmbedding(proj_4, name = "log10(HIST1Counts)", plotAs = 'points') + theme_DC + theme_UMAP

proj_4@cellColData %>% as.data.frame %>% 
  filter(CellType_Major != "B") %>%
  mutate(Drug = factor(Drug, levels(df$Drug))) %>% 
  group_by(CellType_Major, Drug, Dose, Replicate) %>% summarize(HIST1Ratio = mean(HIST1Ratio)) %>%
  ggplot(aes(x = Drug, y = HIST1Ratio, fill = Drug, alpha = Dose)) + geom_boxplot(outliers = F) + facet_wrap(CellType_Major~.) + scale_fill_manual(values = col_drug) + theme_DC

proj_4@cellColData %>% as.data.frame %>% 
  filter(CellType_Major != "B") %>%
  mutate(Drug = factor(Drug, levels(df$Drug))) %>% 
  group_by(CellType_Major, Drug, Dose, Replicate) %>% summarize(HIST1Counts = mean(HIST1Counts)) %>%
  ggplot(aes(x = Drug, y = HIST1Counts, fill = Drug, alpha = Dose)) + geom_boxplot(outliers = F) + facet_wrap(CellType_Major~.) + scale_fill_manual(values = col_drug) + theme_DC



### Export R Objects --------------------------------------------------------

dir.create("Export")

sessionOut("Export/DrugScreen_Packages.tsv")

## All significant markers
markers_lm_sig2 <- markers_lm_sig %>% select(-FDR)
saveRDS(markers_lm_sig2, file = "Export/DrugScreen_SigMarkers.rds")
write.table(markers_lm_sig2, file = "Export/DrugScreen_SigMarkers.csv", quote = F, sep = ',', row.names = F)

## RNA marker GSEA
saveRDS(gsea_lm_rna, file = "Export/DrugScreen_GSEA_RNA.rds")

## GeneScore marker GSEA
saveRDS(gsea_lm_score, file = "Export/DrugScreen_GSEA_GeneScores.rds")

## Single-cell metadata
saveRDS(df, file = "Export/DrugScreen_SingleCellMetadata.rds")



