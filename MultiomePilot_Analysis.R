### Analysis of Multiome Pilot w/ Katya ###

### Set Up Environment ---------------------------------------------------------

library(ArchR)
library(deMULTIplex2)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)

# '%ni%' <- Negate('%in%')

source("/Volumes/DannySSD/MULTI_ATAC/CustomFunctions.R")

addArchRThreads(threads = 1) 
addArchRGenome("mm10")

wd <- "/Volumes/DannySSD/MULTI_ATAC/mtDNA_Katya/MultiomePilot1/ArchR"
dir.create(wd)
setwd(wd)

dir.create("Plots")

### Create Arrow Files ---------------------------------------------------------

# User cellranger-arc output files
inputFiles <- c(
  Pilot = "/Volumes/DannySSD/MULTI_ATAC/mtDNA_Katya/MultiomePilot1/outs_comb/atac_fragments.tsv.gz"
)

ArrowFiles <- vector()

for (i in 1:length(inputFiles)) {
  ArrowFiles[i] <- createArrowFiles(
    inputFiles = inputFiles[i],
    sampleNames = names(inputFiles)[i],
    minTSS = 0, 
    minFrags = 100, 
    addTileMat = TRUE,
    addGeneScoreMat = TRUE, force = T)
}

ArrowFiles <- c("Pilot.arrow")

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
  input = c("/Volumes/DannySSD/MULTI_ATAC/mtDNA_Katya/MultiomePilot1/outs_comb/raw_feature_bc_matrix.h5"),
  names = names(inputFiles) 
)

proj_1 <- addGeneExpressionMatrix(input = proj_1, seRNA = GEX, force = TRUE)


### Filter By ATAC & GEX QC ----------------------------------------------------

## nFrags, TSS enrichment, and nUMI
## Use permissive/generous filtering because of epigenetic drugs used in the experiment

g <- ggplot(as.data.frame(proj_1@cellColData), aes(x = log10(nFrags), y = TSSEnrichment, color = TSSEnrichment > (-13*log10(nFrags) + 45) & Gex_nUMI > 1000)) + 
  ggrastr::rasterise(geom_point(size = 1, alpha = 0.3, stroke = 0), dpi = 600, scale = 1) +
  geom_abline(slope = -13, intercept = 45) + geom_vline(xintercept = 2) + geom_hline(yintercept = 2) + 
  labs(color = "PassFilter", title = "TSSEnrichment > (-13*log10(nFrags) + 45) & Gex_nUMI > 1000") + 
  theme_DC +
  scale_color_manual(values = c("black","dodgerblue")) +
  facet_grid(~Sample) +
  ylim(c(0,quantile(as.data.frame(proj_1@cellColData)$TSSEnrichment, 0.999)))

save_plot(g, "./Plots/proj1_PassFilter.pdf", w = 2.5, h = 2.5, show.legend = "none")

g <- ggplot(as.data.frame(proj_1@cellColData), aes(x = log10(nFrags), y = log10(Gex_nUMI), color = TSSEnrichment > (-13*log10(nFrags) + 45) & Gex_nUMI > 1000)) + 
  ggrastr::rasterise(geom_point(size = 1, alpha = 0.3, stroke = 0), dpi = 600, scale = 1) +
  # geom_abline(slope = -13, intercept = 45) + geom_vline(xintercept = 2) + geom_hline(yintercept = 2) + 
  labs(color = "PassFilter", title = "TSSEnrichment > (-13*log10(nFrags) + 45) & Gex_nUMI > 1000") + 
  theme_DC +
  scale_color_manual(values = c("black","dodgerblue")) +
  facet_grid(~Sample) 


cells <- proj_1$cellNames[proj_1$TSSEnrichment > (-13*log10(proj_1$nFrags) + 45) &
                            proj_1$Gex_nUMI > 1000]

proj_2 <- subsetArchRProject(proj_1,
                             cells = cells,
                             outputDirectory = "proj_2",
                             dropCells = F, 
                             force = T)
# proj_2 <- readRDS("proj_2/Save-ArchR-Project.rds")

## Compute Dim Reduction and UMAP embedding with proper scaling

proj_2 <- dimReduce(proj_2, scale = "Col", default = "NoScale", dims = 2:30)

## Plot Embedding
p1 <- plotEmbedding(proj_2, name = "Clusters", embedding = "UMAP", size = 1, labelAsFactors=F, labelMeans=T) + theme_void() + theme(legend.position = 'none') + ggtitle("Clusters")
p2 <- plotEmbedding(proj_2, name = "log10(nFrags)", embedding = "UMAP", size = 1, labelAsFactors=F, labelMeans=T, plotAs="points") + theme_void() + ggtitle("ATAC nFrags")
p3 <- plotEmbedding(proj_2, name = "log10(Gex_nUMI)", embedding = "UMAP", size = 1, labelAsFactors=F, labelMeans=T, plotAs="points") + theme_void() + ggtitle("GEX nUMI")
plot_grid(p1,p2,p3, nrow = 1)
# pdfPlot(g, width = 12, height = 4, "./Plots/proj2_MetadataUMAP.pdf")


### Add CellRanger Metadata ----------------------------------------------------

csv <- read.csv("/Volumes/DannySSD/MULTI_ATAC/mtDNA_Katya/MultiomePilot1/outs_comb/per_barcode_metrics.csv")
rownames(csv) <- paste("Pilot#", csv$barcode, sep = "")
csv <- csv[proj_2$cellNames,]
csv$Library <- gsub("#.*", "", rownames(csv))

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


### Sample Classification --------------------------------------------------

dir.create(paste(wd, "SampleClassification",sep="/"))

### Create read and barcode tables ###
readTable_MULTIseq <- readTags(dir = "/Volumes/DannySSD/MULTI_ATAC/mtDNA_Katya/MultiomePilot1/MULTI-seq/",
                               name = "MULTIseq", 
                               barcode.type = "MULTIseq", 
                               assay = "Multiome")

readTable_MULTIATAC <- readTags(dir = "/Volumes/DannySSD/MULTI_ATAC/mtDNA_Katya/MultiomePilot1/MULTI-ATAC/",
                                name = "MULTI_Multiome_Cyc7", 
                                barcode.type = "MULTI-ATAC", 
                                assay = "Multiome")

saveRDS(readTable_MULTIseq, "SampleClassification/readTable_MULTIseq.rds")
saveRDS(readTable_MULTIATAC, "SampleClassification/readTable_MULTIATAC.rds")

barTable_MULTIseq <- alignTags(read_table = readTable_MULTIseq,
                               tag.ref = multiseq_oligos[4:5],
                               filter.cells = gsub("-1","",proj_2$CB_CellRangerGEX))
barTable_MULTIATAC <- alignTags(read_table = readTable_MULTIATAC,
                                tag.ref = multiseq_oligos[6:7],
                                filter.cells = revComplement(gsub("-1","",proj_2$CB_CellRangerATAC), cbLength = 16))

saveRDS(barTable_MULTIseq, "SampleClassification/barTable_MULTIseq.rds")
saveRDS(barTable_MULTIATAC, "SampleClassification/barTable_MULTIATAC.rds")

rm(readTable_MULTIseq, readTable_MULTIATAC)
gc()

barTable_MULTIseq <- readRDS("SampleClassification/barTable_MULTIseq.rds")
barTable_MULTIATAC <- readRDS("SampleClassification/barTable_MULTIATAC.rds")

barTable <- cbind(barTable_MULTIseq, barTable_MULTIATAC)

# p1 <- ggplot(barTable_MULTIseq, aes(x = Bar1, y = Bar2)) +
p1 <- ggplot(barTable_MULTIseq, aes(x = A4, y = A5)) +
  geom_point(color = "dark grey") + 
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  ggtitle(label = "MULTIseq")
p2 <- ggplot(barTable_MULTIATAC, aes(x = A6, y = A7)) +
  geom_point(color = "dark grey") + 
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  ggtitle(label = "MULTIATAC")

ggExtra::ggMarginal(p1, type = "histogram", bins = 50, color = "white")
ggExtra::ggMarginal(p2, type = "histogram", bins = 50, color = "white")


### Demultiplexing

res <- demultiplexTags(barTable, 
                       plot.diagnostics = T, 
                       plot.path = "SampleClassification/",
                       plot.name = "Combined")

res_MULTIseq <- demultiplexTags(barTable_MULTIseq, 
                         plot.diagnostics = T, 
                         plot.path = "SampleClassification/",
                         plot.name = "MULTIseq")

res_MULTIATAC <- demultiplexTags(barTable_MULTIATAC, 
                                plot.diagnostics = T, 
                                plot.path = "SampleClassification/",
                                plot.name = "MULTIATAC")

saveRDS(res, "SampleClassification/res.rds")
saveRDS(res_MULTIseq, "SampleClassification/res_MULTIseq.rds")
saveRDS(res_MULTIATAC, "SampleClassification/res_MULTIATAC.rds")

calls <- res$final_assign
saveRDS(calls, "SampleClassification/calls.rds")

tagCallHeatmap(barTable_MULTIseq, res_MULTIseq$final_assign)
tagCallHeatmap(barTable_MULTIATAC, res_MULTIATAC$final_assign)
g <- tagCallHeatmap(cbind(barTable_MULTIseq, barTable_MULTIATAC), res$final_assign) + scale_fill_viridis_c(option="D")
# g[[2]][[2]] <- NULL
save_plot(g, "./Plots/res_Heatmap.pdf", w = 3, h = 3, show.legend = "none")

classDF <- cbind(barTable, res$assign_table, res$umap)
classDF$nUMI <- rowSums(classDF[,1:4])


g <- classDF %>% 
  {
    ggplot(., aes(x = UMAP_1, y = UMAP_2, color = droplet_type, size = droplet_type)) +
      scale_size_manual(values = c(0.5,0.5,1)) + 
      ggrastr::rasterise(geom_point(alpha = 0.6, stroke = 0), dpi = 600, scale = 1) +
      scale_color_manual(values = c("#AD72D6FF","#5AAE61FF","grey10")) +
      theme_DC + theme_UMAP +
      labs(color = "Droplet Type", x = "UMAP1", y = "UMAP2", title = "deMULTIplex2 Classifications") + 
      ggrepel::geom_label_repel(
        data = . %>% group_by(droplet_type, final_assign) %>% summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)), #%>% filter(droplet_type == "singlet"),
        aes(label = final_assign), 
        label.padding = 0.08, max.overlaps = 20, color = 'black', alpha = 0.8, size = 1)
  }
save_plot(g, "./Plots/res_UMAP.pdf", w = 2, h = 2, show.legend = 'none')

g <- classDF %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = final_assign, size = droplet_type, alpha = droplet_type)) +
  scale_size_manual(values = c(0.5,0.5,1)) + scale_alpha_manual(values = c(0.2,0.2,0.5)) + 
  ggrastr::rasterise(geom_point(stroke = 0), dpi = 600, scale = 1) +
  scale_color_manual(values = c(colorRampPalette(c("navyblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1", "maroon"))(4),"black","grey")) + 
  theme_DC + theme_UMAP +
  labs(x = "UMAP1", y = "UMAP2", title = "deMULTIplex2 Classifications")
save_plot(g, "./Plots/res_UMAP_Rainbow.pdf", w = 2, h = 2, show.legend = 'none')


classDF %>% 
  ggplot(aes(x = nUMI, fill = final_assign)) + 
  geom_histogram(bins = 100) + facet_wrap(~final_assign, ncol = 1, scales = 'free_y') + scale_x_log10()


### Merge Classifications ------------------------------------------------------

proj_2$MULTI <- calls[as.vector(gsub("-1","",proj_2$CB_CellRangerGEX))]
proj_2$MULTI[is.na(proj_2$MULTI)] <- "negative"

proj_2$tagUMI <- cbind(barTable_MULTIseq, barTable_MULTIATAC)[as.vector(gsub("-1","",proj_2$CB_CellRangerGEX)),] %>% Matrix::rowSums()

plotEmbedding(proj_2, name = "log10(tagUMI)", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("log10(tagUMI)")
plotEmbedding(proj_2, name = "MULTI", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("Sample")

rm(list = ls(pattern = "barTable"))
rm(list = ls(pattern = "res_"))

## Breakdown Classifications into Singlets, Multiplets, Negatives
cbind(Singlets = sum(calls %ni% c("multiplet","negative")),
      Multiplets = sum(calls == "multiplet"),
      Negatives = sum(calls == "negative"))

proj_2$MULTItype <- factor(proj_2$MULTI, levels = c(names(bar.ref.rev),"multiplet","negative"), labels = c(rep("Singlet",96),"Multiplet","Negative"))

plotEmbedding(proj_2, name = "MULTItype", embedding = "UMAP", plotAs = 'points', size = 1) + theme_void() + ggtitle("MULTItype")

## Labeling Type
proj_2$Label <- factor(proj_2$MULTI,
                       levels = c("A4","A5","A6","A7","multiplet","negative"),
                       labels = c("MULTIseq","MULTIseq","MULTI-ATAC","MULTI-ATAC","multiplet","unlabeled"))

### Analysis ----------------------------------------------------------------

proj_2 <- addIterativeLSI(
  ArchRProj = proj_2,
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix",
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA",
  scaleDims = F,
  corCutOff = 1,
  force = T
)

proj_2 <- addCombinedDims(proj_2, reducedDims = c("IterativeLSI", "LSI_RNA"), name = "Combined", corCutOff = 1, scaleDims = F)

colnames(proj_2@reducedDims$Combined$matRD) <- paste0("LSI",1:60)

proj_2 <- addUMAP(proj_2, reducedDims = "Combined", dimsToUse = 2:60, name = "UMAP_Combined", scaleDims = F, force = TRUE)
proj_2 <- addClusters(proj_2, reducedDims = "Combined", dimsToUse = 2:60, name = "Clusters_Combined", scaleDims = F, force = TRUE)
proj_2 <- addImputeWeights(proj_2, reducedDims = "Combined", dimsToUse = 2:60, scaleDims = F)


plotEmbedding(proj_2, embedding = "UMAP_Combined", name = "MULTI", plotAs = 'points', size = 1) + theme_void() + ggtitle("Sample")


### ATAC Peak Calling ----------------------------------------------------------

pathToMacs2 <- "/Users/dannyconrad/opt/miniconda3/bin/macs2"

proj_2 <- addGroupCoverages(proj_2, groupBy = "Clusters", force = T)
proj_2 <- addReproduciblePeakSet(
  ArchRProj = proj_2, 
  groupBy = "Clusters", 
  pathToMacs2 = pathToMacs2, 
  force = T)
proj_2 <- addPeakMatrix(proj_2, force = T)


g <- plot_grid(plotlist = lapply(c("Clusters_Combined", "MULTItype", "Label", "log10(nFrags)", "TSSEnrichment", "FRIP", "ATAC_PercentMito", "log10(Gex_nUMI)", "log10(tagUMI)"), function(x) {
  g <- plotEmbedding(proj_2, name = x, embedding = "UMAP_Combined", plotAs = 'points', size = 0.75, labelAsFactors=F, labelMeans=T) + labs(x = "UMAP1", y = "UMAP2") +
    theme_DC + theme_UMAP + ggtitle(x)
  if (x %in% c("Clusters_Combined","MULTItype","Label")) { g <- g + theme(legend.position = 'none') }
  g
}), ncol = 3) 
save_plot(g, "./Plots/proj2_QC_Multiome.pdf", w = 6, h = 6, show.legend = "none")



### Remove Doublets and Low-Quality -----------------------------------------


g <- ggplot(as.data.frame(proj_2@cellColData), aes(x = log10(nFrags), y = TSSEnrichment, color = MULTI != "multiplet" & Clusters_Combined %ni% c("C4","C2"))) + 
  ggrastr::rasterise(geom_point(size = 1, alpha = 0.3, stroke = 0), dpi = 600, scale = 1) +
  labs(color = "PassFilter", title = "TSSEnrichment > (-13*log10(nFrags) + 45) & Gex_nUMI > 1000") + 
  theme_DC +
  scale_color_manual(values = c("black","dodgerblue"))

save_plot(g, "./Plots/proj2_PassFilter.pdf", w = 2.5, h = 2.5, show.legend = "none")

g <- ggplot(as.data.frame(proj_2@cellColData), aes(x = log10(nFrags), y = log10(Gex_nUMI), color = MULTItype)) + 
  ggrastr::rasterise(geom_point(size = 1, alpha = 0.3, stroke = 0), dpi = 600, scale = 1) +
  labs(title = "Barcode") + 
  theme_DC
save_plot(g, "./Plots/proj2_QC_Barcode.pdf", w = 2.75, h = 2.5, show.legend = "right")

ggplot(as.data.frame(proj_2@cellColData), aes(x = log10(nFrags), y = log10(Gex_nUMI), color = Clusters_Combined == "C2")) + 
  ggrastr::rasterise(geom_point(size = 1, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
  labs(title = "Barcode") + 
  theme_DC
ggplot(as.data.frame(proj_2@cellColData), aes(x = log10(nFrags), y = log10(Gex_nUMI), color = Clusters_Combined == "C4")) + 
  ggrastr::rasterise(geom_point(size = 1, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
  labs(title = "Barcode") + 
  theme_DC

ggplot(as.data.frame(proj_2@cellColData), aes(x = log10(nFrags), y = log10(Gex_nUMI), color = FRIP)) + 
  ggrastr::rasterise(geom_point(size = 1, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
  labs(title = "Barcode") + scale_color_viridis_c(option = "H") +
  theme_DC
ggplot(as.data.frame(proj_2@cellColData), aes(x = log10(nFrags), y = FRIP, color = log10(Gex_nUMI))) + 
  ggrastr::rasterise(geom_point(size = 1, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
  labs(title = "Barcode") + scale_color_viridis_c(option = "H") +
  theme_DC


cells <- proj_2$cellNames[proj_2$MULTI != "multiplet" & 
                            proj_2$Clusters_Combined %ni% c("C4","C2") &
                            proj_2$FRIP > 0.5 & 
                            proj_2$nFrags > 1000 & 
                            proj_2$Gex_nUMI > 2000]

proj_3 <- subsetArchRProject(proj_2,
                             cells = cells,
                             outputDirectory = "proj_3",
                             dropCells = F, 
                             force = T)

proj_3 <- dimReduce(proj_3, scale = "Col", default = "NoScale", dims = 2:30)

proj_3 <- addIterativeLSI(
  ArchRProj = proj_3,
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix",
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA",
  scaleDims = F,
  corCutOff = 1,
  force = T
)

proj_3 <- addCombinedDims(proj_3, reducedDims = c("IterativeLSI", "LSI_RNA"), name = "Combined", corCutOff = 1, scaleDims = F)

colnames(proj_3@reducedDims$Combined$matRD) <- paste0("LSI",1:60)

proj_3 <- addUMAP(proj_3, reducedDims = "Combined", dimsToUse = 2:60, name = "UMAP_Combined", scaleDims = F, force = TRUE)
proj_3 <- addClusters(proj_3, reducedDims = "Combined", dimsToUse = 2:60, name = "Clusters_Combined", scaleDims = F, force = TRUE)
proj_3 <- addImputeWeights(proj_3, reducedDims = "Combined", dimsToUse = 2:60, scaleDims = F)

proj_3 <- addUMAP(proj_3, reducedDims = "LSI_RNA", dimsToUse = 1:30, name = "UMAP_RNA", scaleDims = F, force = TRUE)
# proj_3 <- addUMAP(proj_3, reducedDims = "LSI_RNA", dimsToUse = 1:15, name = "UMAP_RNA", scaleDims = F, force = TRUE)
proj_3 <- addClusters(proj_3, reducedDims = "LSI_RNA", dimsToUse = 1:30, name = "Clusters_RNA", scaleDims = F, force = TRUE)

proj_3@reducedDims$LSI_RNA$matSVD %>%
  apply(2, var) %>% plot


g <- plot_grid(plotlist = lapply(c("Clusters_Combined", "MULTItype", "Label", "log10(nFrags)", "TSSEnrichment", "FRIP", "ATAC_PercentMito", "log10(Gex_nUMI)", "log10(tagUMI)"), function(x) {
  g <- plotEmbedding(proj_3, name = x, embedding = "UMAP_Combined", plotAs = 'points', size = 0.75, labelAsFactors=F, labelMeans=T) + labs(x = "UMAP1", y = "UMAP2") +
    theme_DC + theme_UMAP + ggtitle(x)
  if (x %in% c("Clusters_Combined","MULTItype","Label")) { g <- g + theme(legend.position = 'none') }
  g
}), ncol = 3) 
save_plot(g, "./Plots/proj3_QC_Multiome.pdf", w = 6, h = 6, show.legend = "none")

g <- plot_grid(plotlist = lapply(c("Clusters_RNA", "MULTItype", "Label", "log10(nFrags)", "TSSEnrichment", "FRIP", "ATAC_PercentMito", "log10(Gex_nUMI)", "log10(tagUMI)"), function(x) {
  g <- plotEmbedding(proj_3, name = x, embedding = "UMAP_RNA", plotAs = 'points', size = 0.75, labelAsFactors=F, labelMeans=T) + labs(x = "UMAP1", y = "UMAP2") +
    theme_DC + theme_UMAP + ggtitle(x)
  if (x %in% c("Clusters_RNA","MULTItype","Label")) { g <- g + theme(legend.position = 'none') }
  g
}), ncol = 3) 
save_plot(g, "./Plots/proj3_QC_RNA.pdf", w = 6, h = 6, show.legend = "none")

g <- plot_grid(plotlist = lapply(c("Clusters", "MULTItype", "Label", "log10(nFrags)", "TSSEnrichment", "FRIP", "ATAC_PercentMito", "log10(Gex_nUMI)", "log10(tagUMI)"), function(x) {
  g <- plotEmbedding(proj_3, name = x, embedding = "UMAP_NoScale", plotAs = 'points', size = 0.75, labelAsFactors=F, labelMeans=T) + labs(x = "UMAP1", y = "UMAP2") +
    theme_DC + theme_UMAP + ggtitle(x)
  if (x %in% c("Clusters","MULTItype","Label")) { g <- g + theme(legend.position = 'none') }
  g
}), ncol = 3) 
save_plot(g, "./Plots/proj3_QC_ATAC.pdf", w = 6, h = 6, show.legend = "none")


proj_3@cellColData %>% as.data.frame %>%
  ggplot(aes(x = MULTI, y = log10(nFrags), fill = Label)) + geom_violin()
proj_3@cellColData %>% as.data.frame %>%
  ggplot(aes(x = MULTI, y = log10(Gex_nUMI), fill = Label)) + geom_boxplot()
proj_3@cellColData %>% as.data.frame %>%
  ggplot(aes(x = MULTI, y = TSSEnrichment, fill = Label)) + geom_boxplot()

proj_3@cellColData %>% as.data.frame %>%
  ggplot(aes(x = log10(nFrags), y = log10(Gex_nUMI), color = Label)) + geom_point() + theme_DC + 
  geom_abline(slope = 0.8, intercept = 1.2, color = 'black')
  
plotEmbedding(proj_3, name = "log10(Gex_nUMI) > 0.8 * log10(nFrags) + 1.2", embedding = "UMAP_Combined", plotAs = 'points', size = 0.75, labelAsFactors=F, labelMeans=T) + labs(x = "UMAP1", y = "UMAP2") +
  theme_DC + theme_UMAP

plotEmbedding(proj_3, name = "FRIP", embedding = "UMAP_Combined", plotAs = 'points', size = 0.75, labelAsFactors=F, labelMeans=T) + labs(x = "UMAP1", y = "UMAP2") +
  theme_DC + theme_UMAP


proj_3@cellColData %>% as.data.frame %>%
  group_by(MULTI, Label)
  
plot_grid(
  proj_3@cellColData %>% as.data.frame %>%
    ggplot(aes(x = log10(nFrags), y = MULTI, fill = Label)) + ggridges::geom_density_ridges() + theme_DC + scale_fill_brewer(palette = "Set2"),
  proj_3@cellColData %>% as.data.frame %>%
    ggplot(aes(x = TSSEnrichment, y = MULTI, fill = Label)) + ggridges::geom_density_ridges() + theme_DC + scale_fill_brewer(palette = "Set2"),
  proj_3@cellColData %>% as.data.frame %>%
    ggplot(aes(x = FRIP, y = MULTI, fill = Label)) + ggridges::geom_density_ridges() + theme_DC + scale_fill_brewer(palette = "Set2"),
  ncol = 1
)

g <- plot_grid(
  # proj_3@cellColData %>% as.data.frame %>%
  #   ggplot(aes(x = log10(tagUMI), y = Label, fill = Label)) + ggridges::geom_density_ridges() + theme_DC + scale_fill_brewer(palette = "Set2"),
  proj_3@cellColData %>% as.data.frame %>%
    ggplot(aes(x = log10(nFrags), y = Label, fill = Label)) + ggridges::geom_density_ridges() + theme_DC + scale_fill_brewer(palette = "Set2"),
  proj_3@cellColData %>% as.data.frame %>%
    ggplot(aes(x = log10(Gex_nUMI), y = Label, fill = Label)) + ggridges::geom_density_ridges() + theme_DC + scale_fill_brewer(palette = "Set2"),
  proj_3@cellColData %>% as.data.frame %>%
    ggplot(aes(x = TSSEnrichment, y = Label, fill = Label)) + ggridges::geom_density_ridges() + theme_DC + scale_fill_brewer(palette = "Set2"),
  proj_3@cellColData %>% as.data.frame %>%
    ggplot(aes(x = FRIP, y = Label, fill = Label)) + ggridges::geom_density_ridges() + theme_DC + scale_fill_brewer(palette = "Set2"),
  ncol = 1
)
save_plot(g, "./Plots/CompareLabeling.pdf", w = 3, h = 5)

g <- proj_3@cellColData %>% as.data.frame %>% select(Label, tagUMI, nFrags, Gex_nUMI, TSSEnrichment, FRIP) %>%
  select(-tagUMI) %>%
  tidyr::pivot_longer(cols = -1, names_to = "Metric", values_to = "Value") %>%
  mutate(Value = case_when(Metric %in% c("nFrags","Gex_nUMI") ~ log10(Value),
                           TRUE ~ Value),
         Metric = case_when(Metric %in% c("nFrags","Gex_nUMI") ~ paste("log10(",Metric,")",sep=""),
                            TRUE ~ Metric)) %>%
  mutate(Label = factor(Label, levels = c("MULTI-ATAC","MULTIseq","unlabeled"), labels = c("MULTI-ATAC","MULTI-seq","Unlabeled"))) %>%
  ggplot(aes(x = Value, y = Label, fill = Label)) + ggridges::geom_density_ridges() + 
  theme_DC + scale_fill_brewer(palette = "Set1") + facet_wrap(~Metric, ncol = 1, scales = 'free_x') +
  theme(axis.title = element_blank())
save_plot(g, "./Plots/CompareLabeling_v2.pdf", w = 2.5, h = 3, show.legend = 'right')

proj_3@cellColData %>% as.data.frame %>% select(Label, tagUMI, nFrags, Gex_nUMI, TSSEnrichment, FRIP) %>%
  select(-tagUMI) %>%
  tidyr::pivot_longer(cols = -1, names_to = "Metric", values_to = "Value") %>%
  mutate(Value = case_when(Metric %in% c("nFrags","Gex_nUMI") ~ log10(Value),
                           TRUE ~ Value),
         Metric = case_when(Metric %in% c("nFrags","Gex_nUMI") ~ paste("log10(",Metric,")",sep=""),
                            TRUE ~ Metric)) %>%
  ggplot(aes(x = Value, y = Label, fill = Label)) + 
  ggridges::geom_density_ridges() +
  # geom_boxplot() +
  theme_DC + scale_fill_brewer(palette = "Set2") + facet_wrap(~Metric, ncol = 1, scales = 'free_x') +
  theme(axis.title = element_blank()) + #coord_flip() +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "unlabeled", paired = F, )
  # stat_compare_means(label = "p.signif", method = "anova")

markers_atac <- getMarkerFeatures(proj_3, groupBy = "Clusters", useMatrix = "PeakMatrix")
markers_score <- getMarkerFeatures(proj_3, groupBy = "Clusters", useMatrix = "GeneScoreMatrix")
markers_rna <- getMarkerFeatures(proj_3, groupBy = "Clusters", useMatrix = "GeneExpressionMatrix")

getMarkers(markers_atac) %>% do.call(rbind,.)
getMarkers(markers_score) %>% do.call(rbind,.) %>% as.data.frame %>% arrange(log10(FDR)) %>% head
getMarkers(markers_rna) %>% do.call(rbind,.) %>% as.data.frame %>% arrange(log10(FDR))

plotEmbedding(proj_3, embedding = "UMAP_Combined")
plotEmbedding(proj_3, embedding = "UMAP_RNA")
plotEmbedding(proj_3, embedding = "UMAP_NoScale")

plotEmbedding(proj_3, embedding = "UMAP_Combined", colorBy = "GeneExpressionMatrix", name = "Tmed4", plotAs = 'points', size = 1)
plotEmbedding(proj_3, embedding = "UMAP_Combined", colorBy = "GeneScoreMatrix", name = "Tmed4", plotAs = 'points', size = 1)

plotEmbedding(proj_3, embedding = "UMAP_Combined", colorBy = "GeneExpressionMatrix", name = "Tmed4", plotAs = 'points', size = 1)
plotEmbedding(proj_3, embedding = "UMAP_Combined", colorBy = "GeneScoreMatrix", name = "Tmed4", plotAs = 'points', size = 1)

genes <- c("Serpina3k","Cyp2e1","Ass1", "Alb", "Oat")


cv <- c("Axin2", "Cyp1a2", "Gstm3", "Psmd4", "Glul", "Cyp2e1")
nm <- c("Hamp", "Igfbp2", "Cyp8b1", "Mup3")
pv <- c("Arg1", "Pck1", "C2", "Sdhd", "Ass1", "Asl", "Alb", "Cyp2f2")

g <- lapply(cv, function(x) {
  plotEmbedding(proj_3, embedding = "UMAP_RNA", colorBy = "GeneExpressionMatrix", name = x, plotAs = 'points', size = 1) + theme_UMAP + theme_DC
}) %>% plot_grid(plotlist = ., ncol = 2)
save_plot(g, "./Plots/MarkerGenes_CentralVein.pdf", w = 4, h = 6)

g <- lapply(pv, function(x) {
  plotEmbedding(proj_3, embedding = "UMAP_RNA", colorBy = "GeneExpressionMatrix", name = x, plotAs = 'points', size = 1) + theme_UMAP + theme_DC
}) %>% plot_grid(plotlist = ., ncol = 2)
save_plot(g, "./Plots/MarkerGenes_PortalVein.pdf", w = 4, h = 6)

g <- lapply(nm, function(x) {
  plotEmbedding(proj_3, embedding = "UMAP_RNA", colorBy = "GeneExpressionMatrix", name = x, plotAs = 'points', size = 1) + theme_UMAP + theme_DC
}) %>% plot_grid(plotlist = ., ncol = 2)
save_plot(g, "./Plots/MarkerGenes_PortalVein.pdf", w = 4, h = 6)


lapply(genes, function(x) {
  plotEmbedding(proj_3, embedding = "UMAP_RNA", colorBy = "GeneExpressionMatrix", name = x, plotAs = 'points', size = 1) + theme_UMAP + theme_DC
}) %>% plot_grid(plotlist = ., ncol = 3)


proj_3 <- addMotifAnnotations(ArchRProj = proj_3, motifSet = "cisbp", annoName = "Motif", force = T) ## cis-bp v1
proj_3 <- addBgdPeaks(proj_3)
proj_3 <- addDeviationsMatrix(
  ArchRProj = proj_3,
  peakAnnotation = "Motif", force = T
)
getVarDeviations(proj_3, plot = TRUE, name = "MotifMatrix")

plotEmbedding(proj_3, embedding = "UMAP_NoScale", colorBy = "MotifMatrix", name = "z:Ctcf_146", plotAs = 'points', size = 1)

getFeatures(proj_3, useMatrix = 'MotifMatrix', select = "irf")



proj_3$Zonation <- "Portal Vein"
proj_3$Zonation[proj_3$Clusters_RNA %in% c("C1","C2","C3")] <- "Central Vein"

plotEmbedding(proj_3, embedding = "UMAP_Combined", name = "Zonation", plotAs = 'points', size = 1)
plotEmbedding(proj_3, embedding = "UMAP_NoScale", name = "Zonation", plotAs = 'points', size = 1)



g <- proj_3@cellColData %>% data.frame(UMAP1 = proj_3@embeddings$UMAP_RNA$df$`LSI_RNA#UMAP_Dimension_1`,
                                       UMAP2 = proj_3@embeddings$UMAP_RNA$df$`LSI_RNA#UMAP_Dimension_2`,
                                       .) %>%
  mutate(Label = factor(Label, levels = c("MULTI-ATAC","MULTIseq","unlabeled"), labels = c("MULTI-ATAC","MULTI-seq","Unlabeled"))) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Label)) + geom_point(size = 0.75, stroke = 0) +
  scale_color_brewer(palette = "Set1") + guides(colour = guide_legend(override.aes = list(size=3))) + labs(title = "Mouse Hepatocytes") + 
  theme_DC + theme_UMAP
save_plot(g, "./Plots/proj3_Label.pdf", w = 2.5, h = 2, show.legend = 'right')
g <- proj_3@cellColData %>% data.frame(UMAP1 = proj_3@embeddings$UMAP_RNA$df$`LSI_RNA#UMAP_Dimension_1`,
                                       UMAP2 = proj_3@embeddings$UMAP_RNA$df$`LSI_RNA#UMAP_Dimension_2`,
                                       .) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Zonation)) + geom_point(size = 0.75, stroke = 0) +
  scale_color_brewer(palette = "Dark2") + guides(colour = guide_legend(override.aes = list(size=3))) + labs(title = "Mouse Hepatocytes") + 
  theme_DC + theme_UMAP
save_plot(g, "./Plots/proj3_Zonation.pdf", w = 2.5, h = 2, show.legend = 'right')

g <- proj_3@cellColData %>% as.data.frame %>%
  ggplot(aes(x = MULTI, fill = Zonation)) + geom_bar(position = 'fill')




### Export R Objects --------------------------------------------------------

dir.create("Export")

sessionOut("Export/MultiomePilot_Packages.tsv")

## Single-cell metadata
saveRDS(proj_3@cellColData, file = "Export/MultiomePilot_SingleCellMetadata.rds")

