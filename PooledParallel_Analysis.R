### Finalized Analysis of MULTI-ATAC Batch Effect Experiment ###

### Set Up Environment ------------------------------------------------------

library(ArchR)
library(dplyr)
library(cowplot)
library(pheatmap)
library(BSgenome.Hsapiens.UCSC.hg38)
library(deMULTIplex2)
library(ggpubr)
'%ni%' <- Negate('%in%')

source("/Volumes/DannySSD/MULTI_ATAC/CustomFunctions.R")

addArchRThreads(threads = 1)
addArchRGenome("hg38")

wd <- "/Volumes/DannySSD/MULTI_ATAC/20220805_BatchEffectExperiment5/ArchR"
dir.create(wd)
setwd(wd)

dir.create("./Plots")


### Create Arrow Files ------------------------------------------------------

inputFiles <- c(
  Para = "/Volumes/DannySSD/MULTI_ATAC/20220805_BatchEffectExperiment5/Para/fragments.tsv.gz",
  Pool = "/Volumes/DannySSD/MULTI_ATAC/20220805_BatchEffectExperiment5/Pool/fragments.tsv.gz"
)

ArrowFiles <- vector()

for (i in 1:length(inputFiles)) {
  ArrowFiles[i] <- createArrowFiles(
    inputFiles = inputFiles[i],
    sampleNames = names(inputFiles)[i],
    minTSS = 1, 
    minFrags = 100,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE, 
    force = T)
}

ArrowFiles <- c("Para.arrow", "Pool.arrow")


### Create ArchR Project ----------------------------------------------------

proj_1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "proj_1",
  copyArrows = TRUE
)

thresh.TSS <- 5
thresh.nFrag <- 200

ggplot(as.data.frame(proj_1@cellColData[proj_1$nFrags > thresh.nFrag & proj_1$TSSEnrichment > thresh.TSS,]), aes(x = log10(nFrags), y = TSSEnrichment)) + 
  geom_density_2d_filled(bins = 10) + 
  geom_vline(xintercept = log10(thresh.nFrag)) + 
  geom_hline(yintercept = thresh.TSS) + 
  facet_grid(~Sample) + 
  theme_minimal() + 
  theme(legend.position = "none")
ggplot(as.data.frame(proj_1@cellColData[proj_1$nFrags > thresh.nFrag & proj_1$TSSEnrichment > thresh.TSS,]), aes(x = log10(nFrags), y = TSSEnrichment)) + 
  geom_point() + 
  geom_vline(xintercept = log10(thresh.nFrag)) + 
  geom_hline(yintercept = thresh.TSS) + 
  facet_grid(~Sample) + 
  theme_minimal() + 
  theme(legend.position = "none")

proj_1 <- proj_1[proj_1$nFrags > thresh.nFrag & proj_1$TSSEnrichment > thresh.TSS]

proj_1 <- addIterativeLSI(
  ArchRProj = proj_1, 
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "IterativeLSI", 
  corCutOff = 1,
  scaleDims = F,
  force = T
)

# including Dim1 useful for doublet identification
proj_1 <- addClusters(proj_1, dimsToUse = 1:30, scaleDims = F, corCutOff = 1, force = T)
proj_1 <- addUMAP(proj_1, dimsToUse = 1:30, scaleDims = F, corCutOff = 1, force = T)

plot_grid(plotlist = 
            lapply(c("Sample", "Clusters", "log10(nFrags)", "TSSEnrichment"), function(x) {
              plotEmbedding(proj_1, name = x, plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle(x) 
            }),
          ncol = 2)



### Add CellRanger Metadata -------------------------------------------------

proj_1$CB_CellRanger <- gsub(".*#", "", proj_1$cellNames)
proj_1$CB <- gsub("-.*", "", proj_1$CB_CellRanger)
proj_1$CB_ReverseComp <- revComplement(proj_1$CB)

csv_para <- read.csv("/Volumes/DannySSD/MULTI_ATAC/20220805_BatchEffectExperiment5/Para/singlecell.csv")[-1,]
csv_pool <- read.csv("/Volumes/DannySSD/MULTI_ATAC/20220805_BatchEffectExperiment5/Pool/singlecell.csv")[-1,]

rownames(csv_para) <- paste("Para#", csv_para$barcode, sep="")
rownames(csv_pool) <- paste("Pool#", csv_pool$barcode, sep="")

csv_para_cells <- csv_para[csv_para$is__cell_barcode == 1, c("barcode", "is__cell_barcode")]
csv_pool_cells <- csv_pool[csv_pool$is__cell_barcode == 1, c("barcode", "is__cell_barcode")]

proj_1$CellRanger_Cell <- FALSE
proj_1$CellRanger_Cell[proj_1$Sample == "Para" & proj_1$CB_CellRanger %in% csv_para_cells$barcode] <- TRUE
proj_1$CellRanger_Cell[proj_1$Sample == "Pool" & proj_1$CB_CellRanger %in% csv_pool_cells$barcode] <- TRUE



### deMULTIplex2 Classification ---------------------------------------------

dir.create(paste(wd, "SampleClassification",sep="/"))

data(multiseq_oligos)
bar.ref <- multiseq_oligos[8:16] 

cellIDs_Para <- proj_1$CB_ReverseComp[proj_1$Sample == "Para"]
cellIDs_Pool <- proj_1$CB_ReverseComp[proj_1$Sample == "Pool"]


# Parallel
readTable_Para <- readTags(dir = "/Volumes/DannySSD/MULTI_ATAC/20220805_BatchEffectExperiment5/MULTI/",
                           name = "MULTI_Para", 
                           barcode.type = "MULTI-ATAC", 
                           assay = "ATAC",
                           filter.cells = cellIDs_Para)
saveRDS(readTable_Para, "SampleClassification/readTable_Para.rds")

barTable_Para <- alignTags(read_table = readTable_Para, 
                           tag.ref = bar.ref)
saveRDS(barTable_Para, "SampleClassification/barTable_Para.rds")
rm(readTable_Para)
gc()

# Pooled
readTable_Pool <- readTags(dir = "/Volumes/DannySSD/MULTI_ATAC/20220805_BatchEffectExperiment5/MULTI/",
                           name = "MULTI_Pool", 
                           barcode.type = "MULTI-ATAC", 
                           assay = "ATAC",
                           filter.cells = cellIDs_Pool)
saveRDS(readTable_Pool, "SampleClassification/readTable_Pool.rds")

barTable_Pool <- alignTags(read_table = readTable_Pool, 
                           tag.ref = bar.ref)
saveRDS(barTable_Pool, "SampleClassification/barTable_Pool.rds")
rm(readTable_Pool)
gc()

# readTable_Para <- readRDS("SampleClassification/readTable_Para.rds")
# readTable_Pool <- readRDS("SampleClassification/readTable_Pool.rds")
barTable_Para <- readRDS("SampleClassification/barTable_Para.rds")
barTable_Pool <- readRDS("SampleClassification/barTable_Pool.rds")

# barTable_Para <- barTable_Para[Matrix::rowSums(barTable_Para) >= 10,]
# barTable_Pool <- barTable_Pool[Matrix::rowSums(barTable_Pool) >= 10,]

## Visualize barcode distributions

tagHist(barTable_Para, scale_y_sqrt = T)
tagHist(barTable_Pool, scale_y_sqrt = T)

## Run Classification

calls_Para <- demultiplexTags(barTable_Para, 
                              plot.diagnostics = T, 
                              plot.path = "SampleClassification/", 
                              plot.name = "Para",
                              init.cos.cut = 0.5)

calls_Pool <- demultiplexTags(barTable_Pool, 
                              plot.diagnostics = T, 
                              plot.path = "SampleClassification/", 
                              plot.name = "Pool",
                              init.cos.cut = 0.5)

saveRDS(calls_Para, "./SampleClassification/calls_Para.rds")
saveRDS(calls_Pool, "./SampleClassification/calls_Pool.rds")

calls_Para <- readRDS("./SampleClassification/calls_Para.rds")
calls_Pool <- readRDS("./SampleClassification/calls_Pool.rds")

tagCallHeatmap(barTable_Para[names(calls_Para$final_assign),], calls_Para$final_assign)
tagCallHeatmap(barTable_Pool[names(calls_Pool$final_assign),], calls_Pool$final_assign)

calls <- c(calls_Para$final_assign, calls_Pool$final_assign)
names(calls) <- c(paste("Para#", revComplement(names(calls_Para$final_assign)), "-1", sep = ""),
                  paste("Pool#", revComplement(names(calls_Pool$final_assign)), "-1", sep = ""))


## Add classifications, sample info, & barcode data to ArchR project
proj_1$Calls <- calls[row.names(proj_1@cellColData)]

barTable <- rbind(barTable_Para, barTable_Pool)
row.names(barTable) <- c(paste("Para#", as.character(reverseComplement(DNAStringSet(row.names(barTable_Para)))), "-1", sep = ""),
                         paste("Pool#", as.character(reverseComplement(DNAStringSet(row.names(barTable_Pool)))), "-1", sep = ""))


proj_1$nUMI <- rowSums(barTable)[row.names(proj_1@cellColData)]

proj_1$Calls <- factor(proj_1$Calls,
                       levels = c(names(bar.ref),"multiplet", "negative"),
                       labels = c(paste("Bar", 1:9, sep = ""),"Doublet", "Negative"))
proj_1$Density <- factor(proj_1$Calls,
                         levels = c(paste("Bar", 1:9, sep = ""),"Doublet", "Negative"),
                         labels = c(rep(c("Low", "Medium", "High"), each = 3), "Doublet", "Negative"))

g <- plot_grid(plotlist =
                 lapply(c("Calls", "Density"), function(x) {
                   plotEmbedding(proj_1, name = x, plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle(x)
                 }),
               ncol = 2)


metadata <- c("Sample", "Clusters", "Calls", "CellRanger_Cell",
              "log10(nFrags)", "log10(nUMI)", "TSSEnrichment", "PromoterRatio")

g <- plot_grid(plotlist = 
                 lapply(metadata, function(x) {
                   plotEmbedding(proj_1[!is.na(proj_1$nUMI)], name = x, plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle(x) 
                 }),
               nrow = 2)



### Clean Up ArchR Project --------------------------------------------------

## Cluster Composition
plotEmbedding(proj_1, name = 'Calls %in% c("Doublet","Negative",NA)', size = 1) + theme_void()
plotEmbedding(proj_1, name = 'Clusters', size = 1) + theme_void() + ggtitle("Clusters") + theme(legend.position = "none")

ggplot(as.data.frame(proj_1@cellColData), aes(x = Clusters, fill = Calls == 'Doublet')) + geom_bar(position = 'fill')

## Remove clearly true doublet & true negative (debris/empty drop) cells/clusters
cells <- proj_1$cellNames[!proj_1$Calls %in% c("Doublet", "Negative", NA) & !proj_1$Clusters %in% c("C1","C5")]

proj_2 <- subsetArchRProject(proj_1,
                             cells = cells, 
                             outputDirectory = "proj_2",
                             dropCells = T, 
                             force = T)

proj_2 <- addIterativeLSI(
  ArchRProj = proj_2, 
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "IterativeLSI", 
  corCutOff = 1,
  scaleDims = F,
  force = T
)

proj_2 <- addClusters(proj_2, dimsToUse = 2:30, scaleDims = F, corCutOff = 1, force = T)
proj_2 <- addUMAP(proj_2, dimsToUse = 2:30, scaleDims = F, corCutOff = 1, force = T)
proj_2 <- addImputeWeights(proj_2, dimsToUse = 2:30, scaleDims = F, corCutOff = 1)

metadata <- c("Sample", "Clusters", "Calls", "CellRanger_Cell",
              "log10(nFrags)", "log10(nUMI)", "TSSEnrichment", "PromoterRatio")

proj_2$Transposition <- factor(proj_2$Sample, levels = c("Para","Pool"), labels = c("Parallel","Pooled"))

g <- plot_grid(plotlist = 
                 lapply(metadata, function(x) {
                   plotEmbedding(proj_2, name = x, plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle(x) 
                 }),
               nrow = 2)
pdfPlot(g, width = 12, height = 8, file = "./Plots/proj2_Metadata.pdf")

plot_grid(plotlist = 
            lapply("Density", function(x) {
              plotEmbedding(proj_2, name = x, plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle(x) 
            }),
          nrow = 1)



### Cell Type Annotation ----------------------------------------------------

plotEmbedding(proj_2, name = "Clusters", size = 1, labelAsFactors = F) + theme_void() + ggtitle("Clusters") + theme(legend.position = 'bottom')

plotEmbedding(proj_2,  colorBy = "GeneScoreMatrix", name = "NOTCH1", plotAs = "points") + theme_void()  # Jurkat Marker
plotEmbedding(proj_2,  colorBy = "GeneScoreMatrix", name = "TCF3", plotAs = "points") + theme_void()    # Jurkat Marker
plotEmbedding(proj_2,  colorBy = "GeneScoreMatrix", name = "GATA2", plotAs = "points") + theme_void()   # K562 Marker
plotEmbedding(proj_2,  colorBy = "GeneScoreMatrix", name = "ZFPM1", plotAs = "points") + theme_void()   # K562 Marker

proj_2$CellType <- "K562"
proj_2$CellType[proj_2@embeddings$UMAP$df$`IterativeLSI#UMAP_Dimension_1` > 0] <- "Jurkat"

g <- plotEmbedding(proj_2, name = "CellType", size = 1, labelAsFactors = F) + theme_void() + ggtitle("CellType") + theme(legend.position = 'bottom') + scale_color_manual(values = colors_cell)
pdfPlot(g, width = 5, height = 5, file = "./Plots/proj2_CellType.pdf")



### Plot Library Quality Stats ----------------------------------------------

## Color Palettes

colors_bar <- paletteer::paletteer_d("colorBlindness::SteppedSequential5Steps")[c(19:17,14:12,4:2)]
colors_bar %>% scales::show_col()
names(colors_bar) <- sort(unique(res$Calls))
colors_den <- colors_bar[c(2,5,8)]
names(colors_den) <- c("Low","Medium","High")
colors_cell <- paletteer::paletteer_d("colorBlindness::SteppedSequential5Steps")[c(8,23)]
colors_cell %>% scales::show_col()
names(colors_cell) <- c("Jurkat","K562")
c(colors_den,colors_cell) %>% scales::show_col()


# Cell recovery by Barcode/Library
ggplot(as.data.frame(proj_2@cellColData), aes(x = Calls, fill = Calls)) + 
  geom_bar(position = "stack") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(rows = "Sample") +
  scale_fill_manual(values = colors_bar)

# CellType proportions by Barcode/Library
ggplot(as.data.frame(proj_2@cellColData), aes(x = Density, fill = CellType)) + 
  geom_bar(position = "fill", alpha = 0.5) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(rows = "Sample") +
  scale_fill_manual(values = colors_cell)
ggplot(as.data.frame(proj_2@cellColData[proj_2$CellRanger_Cell == 1,]), aes(x = Density, fill = CellType)) + 
  geom_bar(position = "fill", alpha = 0.5) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(rows = "Sample") +
  scale_fill_manual(values = colors_cell)

# Ratio of cell types
as.data.frame(proj_2@cellColData) %>% group_by(Sample, Density) %>% 
  dplyr::summarize(Ratio = sum(CellType == "K562")/sum(CellType == "Jurkat")) %>%
  ggplot(aes(x = Density, y = Ratio, color = Density, fill = Density, group = Sample)) + 
  # geom_col(alpha = 0.5) +
  geom_col(alpha = 0.8) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(~Sample) + ylab("Ratio of K562::Jurkat Cells") +
  # scale_y_continuous(expand = c(0,0), limits = c(0, 1.4), breaks= seq(0, 1.4, by = .2)) + 
  scale_fill_manual(values = colors_den) +
  scale_color_manual(values = colors_den)

# Ratio of cell types
as.data.frame(proj_2@cellColData) %>% group_by(Sample, Density) %>% 
  dplyr::summarize(Jurkat = sum(CellType == "Jurkat")/sum(CellType == "K562"), K562 = 1) %>% reshape2::melt() %>%
  ggplot(aes(x = Density, y = value, color = Density, fill = Density, group = variable, linetype = rev(variable))) + 
  geom_line(color = "grey30", linewidth = 1) +
  geom_point(size = 3) +
  theme_bw() + theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(~Sample) + ylab("Relative Cell Type Proportions") +
  scale_color_manual(values = colors_den)
as.data.frame(proj_2@cellColData) %>% group_by(Sample, Density) %>% 
  dplyr::summarize(Jurkat = sum(CellType == "Jurkat")/sum(CellType == "K562"), K562 = 1) %>% 
  ggplot(aes(x = Density, y = Jurkat, color = Density, fill = Density, group = Sample, linetype = rev(Sample))) + 
  geom_point(size = 3) +
  geom_smooth(method = 'lm', se = F, color = "grey30", linewidth = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  ylab("Relative Proportion of Jurkat Cells") +
  scale_color_manual(values = colors_den)

g <- as.data.frame(proj_2@cellColData) %>% group_by(Sample, Density, CellType) %>% 
  summarize(Count = n()) %>% mutate(Freq = Count/sum(Count)) %>%
  ggplot(aes(x = Density, y = Freq, color = Density, group = CellType, linetype = CellType)) + 
  geom_smooth(mapping = aes(group = CellType), method = 'lm', se = F, linewidth = 1, color = "grey30") +
  geom_point(size = 3) +
  theme_classic(base_line_size=0.5/.pt, base_size = 7) + 
  facet_grid(~Sample) + ylab("Cell Type Proportions") +
  scale_color_manual(values = c(colors_den,colors_cell)) + scale_linetype_manual(values = c("solid","11"))
save_plot(g, "./Plots/proj2_ParaPool_CellTypeProportions_v1.pdf", w = 3, h = 2.5, show.legend = "right")

g <- as.data.frame(proj_2@cellColData) %>% group_by(Transposition, Density, Calls, CellType) %>% 
  summarize(Count = n()) %>% mutate(Freq = Count/sum(Count)) %>%
  mutate(CellType = factor(CellType, levels = c("K562","Jurkat"))) %>%
  ggplot(aes(x = Density, y = Freq, color = CellType, group = CellType)) + 
  geom_point(size = 2) + 
  geom_smooth(mapping = aes(group = CellType), method = 'lm', se = F, linewidth = 1, color = "grey30") +
  theme_DC + 
  theme(axis.text.x = element_text(color = colors_den)) +
  facet_grid(~Transposition) + labs(x = "Nuclei:Tn5 Ratio", y = "Cell Type Proportions") +
  scale_color_manual(values = c(colors_den,colors_cell)) + 
  scale_y_continuous(breaks=scales::pretty_breaks(n = 4)) +
  stat_cor(size = 2, p.accuracy = 0.01, 
           label.y = c(0.6,0.4),
           output.type = 'expression', show.legend = F)
save_plot(g, "./Plots/proj2_ParaPool_CellTypeProportions_v2.pdf", w = 3, h = 2, show.legend = "right")



### Fragments per Cell Type -------------------------------------------------

# K562 yields more fragments than Jurkat across each condition
as.data.frame(proj_2@cellColData) %>% ggplot(aes(x = Calls, y = nFrags, fill = CellType)) + 
  geom_boxplot(outliers = F, alpha = 0.5) + 
  scale_y_log10() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(rows = "Sample") + 
  scale_fill_manual(values = colors_cell)
as.data.frame(proj_2@cellColData) %>% ggplot(aes(x = Density, y = nFrags, fill = CellType)) + 
  geom_boxplot(outliers = F, alpha = 0.5) + 
  scale_y_log10() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(rows = "Sample") + 
  scale_fill_manual(values = colors_cell)

g <- as.data.frame(proj_2@cellColData) %>% 
  ggplot(aes(x = CellType, y = nFrags, fill = CellType)) + 
  geom_boxplot(outliers = F, alpha = 0.5) +
  # geom_violin(alpha = 0.5) +
  scale_y_log10() +
  theme_DC + 
  scale_fill_manual(values = colors_cell) + 
  stat_compare_means(method = 't.test', label = "p.signif", label.x.npc = "center", hjust = 0.5, vjust = 1)
save_plot(g, "./Plots/proj2_CellType_nFrags.pdf", w = 1.5, h = 2, show.legend = "none")

proj_2@cellColData %>% as.data.frame %>% 
  group_by(Sample,CellType) %>% summarize(nFrags = mean(nFrags)) %>%
  group_by(Sample) %>% summarize(Ratio = nFrags[CellType == "Jurkat"]/nFrags[CellType == "K562"])

proj_2@cellColData %>% as.data.frame %>% 
  group_by(CellType) %>% summarize(nFrags = mean(nFrags)) %>%
  summarize(Ratio = nFrags[CellType == "Jurkat"]/nFrags[CellType == "K562"])



### Fragments per Transposition Density -------------------------------------

as.data.frame(proj_2@cellColData) %>% group_by(Sample,Calls) %>% 
  dplyr::summarise(mean = mean(nFrags), sd = sd(nFrags), median = median(nFrags), geom = EnvStats::geoMean(nFrags), geosd = EnvStats::geoSD(nFrags)) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = Sample, y = geom, fill = Calls)) + geom_col(position = "dodge", color = 'white') + 
  geom_errorbar(aes(ymax = geom * geosd, ymin = geom / geosd), position = position_dodge(0.9), width = 0.5, alpha = 0.4) + 
  theme_bw() + scale_fill_manual(values = colors_bar)

as.data.frame(proj_2@cellColData) %>% group_by(Sample,Calls) %>% #mutate(nFrags = log10(nFrags)) %>%
  dplyr::summarise(mean = mean(nFrags), sd = sd(nFrags), median = median(nFrags), geom = EnvStats::geoMean(nFrags), geosd = EnvStats::geoSD(nFrags)) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = Sample, y = geom, color = Calls)) + geom_point(position = position_dodge(0.9), size = 3) + 
  geom_errorbar(aes(ymax = geom * geosd, ymin = geom / geosd), position = position_dodge(0.9), width = 0.5, alpha = 1, linewidth = 1) + 
  theme_bw() + scale_color_manual(values = colors_bar) + scale_y_log10()

as.data.frame(proj_2@cellColData) %>% group_by(Sample,Calls) %>% mutate(nFrags = log10(nFrags)) %>%
  dplyr::summarise(mean = mean(nFrags), sd = sd(nFrags), median = median(nFrags), geom = EnvStats::geoMean(nFrags), geosd = EnvStats::geoSD(nFrags)) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = Sample, y = mean, color = Calls)) + geom_point(position = position_dodge(0.9), size = 3) + 
  geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd), position = position_dodge(0.9), width = 0.5, alpha = 1, linewidth = 1) + 
  theme_bw() + scale_color_manual(values = colors_bar)

as.data.frame(proj_2@cellColData) %>% group_by(Sample,Calls) %>% 
  dplyr::summarise(mean = mean(nFrags), sd = sd(nFrags), median = median(nFrags), geom = EnvStats::geoMean(nFrags), geosd = EnvStats::geoSD(nFrags)) %>% 
  as.data.frame() %>% 
  # ggplot(aes(x = Sample, y = mean, fill = Calls)) + geom_col(position = "dodge", color = 'white') + 
  ggplot(aes(x = Sample, y = geom, fill = Calls)) + geom_col(position = "dodge", color = 'white') + 
  # geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd), position = position_dodge(0.9), width = 0.5, alpha = 0.4) + 
  theme_bw() + scale_fill_manual(values = colors_bar) 


as.data.frame(proj_2@cellColData) %>% ggplot(aes(x = Density, y = nFrags, fill = Calls)) + 
  geom_jitter(mapping = aes(color = Calls), size = 1, alpha = 0.5, position = position_jitterdodge(dodge.width = 1,jitter.width = 2)) + 
  geom_boxplot(outliers = F, linewidth = 0.75, width = 1, alpha = 0.5) +
  scale_y_log10() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(~Sample) + 
  scale_color_manual(values = colors_bar) +
  scale_fill_manual(values = colors_bar) 

g <- as.data.frame(proj_2@cellColData) %>% ggplot(aes(x = Density, y = nFrags, fill = Calls)) + 
  geom_boxplot(outliers = F, linewidth = 0.5, width = 1, alpha = 1, coef = 0.5) +
  theme_DC + theme(axis.text.x = element_text(color = colors_den)) +
  facet_grid(~Transposition) + labs(x = "Nuclei:Tn5 Ratio") +
  scale_fill_manual(values = colors_bar) +
  stat_compare_means(aes(group = Density), label.y = 5800, label.x = 2, hjust = 0.5, size = 2)
save_plot(g, "./Plots/proj2_Density_nFrags.pdf", w = 3, h = 2, show.legend = "none")



### Other Stats -------------------------------------------------------------

# TSSEnrichment
as.data.frame(proj_2@cellColData) %>% ggplot(aes(x = Calls, y = TSSEnrichment, fill = Calls)) + 
  geom_boxplot(outliers = F) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(rows = "Sample") +
  scale_fill_manual(values = colors_bar)

# FRIP
as.data.frame(proj_2@cellColData) %>% ggplot(aes(x = Calls, y = FRIP, fill = Calls)) + 
  geom_boxplot(outliers = F) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(rows = "Sample") +
  scale_fill_manual(values = colors_bar)

## CellRanger Cell Calling
as.data.frame(proj_2@cellColData) %>% ggplot(aes(x = log10(nFrags), y = TSSEnrichment, color = CellRanger_Cell)) + 
  geom_point() + 
  facet_grid(~Sample) + 
  theme_minimal() + 
  theme(legend.position = "none")
as.data.frame(proj_2@cellColData) %>% ggplot(aes(x = Calls, fill = Calls, alpha = CellRanger_Cell)) + 
  geom_bar(position = 'fill') + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(rows = "Sample") +
  scale_fill_manual(values = colors_bar) + 
  scale_alpha_manual(values = c(0.6,1))

as.data.frame(proj_2@cellColData) %>% ggplot(aes(x = Calls, y = nFrags, fill = CellRanger_Cell)) + 
  geom_boxplot(outliers = F) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(rows = "Sample") +
  scale_y_log10()
as.data.frame(proj_2@cellColData) %>% ggplot(aes(x = Calls, y = TSSEnrichment, fill = CellRanger_Cell)) + 
  geom_boxplot(outliers = F) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(rows = "Sample")
as.data.frame(proj_2@cellColData) %>% ggplot(aes(x = Calls, y = FRIP, fill = CellRanger_Cell)) + 
  geom_boxplot(outliers = F) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  facet_grid(rows = "Sample")


csv_para %>% filter(passed_filters > 10) %>% 
  ggplot(aes(x = peak_region_fragments/passed_filters, fill = as.character(is__cell_barcode))) + 
  geom_histogram(bins = 50) + facet_wrap(~is__cell_barcode, scales = 'free_y') + theme_bw()
csv_pool %>% filter(passed_filters > 10) %>% 
  ggplot(aes(x = peak_region_fragments/passed_filters, fill = as.character(is__cell_barcode))) + 
  geom_histogram(bins = 50) + facet_wrap(~is__cell_barcode, scales = 'free_y') + theme_bw()




### Separate Parallel/Pooled Projects ---------------------------------------

metadata <- c("Clusters", "CellType", "Density",
              "log10(nFrags)", "log10(nUMI)", "TSSEnrichment")

## Parallel
cells <- proj_2$cellNames[proj_2$Sample == "Para"]

proj_para <- subsetArchRProject(proj_2,
                                cells = cells, 
                                outputDirectory = "proj_para",
                                dropCells = T, force = T)

proj_para <- dimReduce(proj_para, default = 'NoScale', dims = 2:30)

proj_para <- addImputeWeights(proj_para, dimsToUse = 2:30, scaleDims = F, corCutOff = 1)



g <- plot_grid(plotlist = 
                 lapply(metadata, function(x) {
                   plotEmbedding(proj_para, name = x, plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle(x) 
                 }),
               nrow = 2)
pdfPlot(g, width = 12, height = 8, file = "./Plots/projPara_Metadata.pdf")

plot_grid(plotlist = 
            lapply(metadata[c(3,4,6)], function(x) {
              plotEmbedding(proj_para, name = x, plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle(x) 
            }),
          nrow = 1)

## Pooled
cells <- proj_2$cellNames[proj_2$Sample == "Pool"]

proj_pool <- subsetArchRProject(proj_2,
                                cells = cells, 
                                outputDirectory = "proj_pool",
                                dropCells = T, force = T)

proj_pool <- dimReduce(proj_pool, default = 'NoScale', dims = 2:30)

proj_pool <- addImputeWeights(proj_pool, dimsToUse = 2:30, scaleDims = F, corCutOff = 1)

g <- plot_grid(plotlist = 
                 lapply(metadata, function(x) {
                   plotEmbedding(proj_pool, name = x, plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle(x) 
                 }),
               nrow = 2)
pdfPlot(g, width = 12, height = 8, file = "./Plots/projPool_Metadata.pdf")

plot_grid(plotlist = 
            lapply(metadata[c(3,4,6)], function(x) {
              plotEmbedding(proj_pool, name = x, plotAs = 'points', size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle(x) 
            }),
          nrow = 1)

### How does parallel tagmentation affect data interpretation?

## By CellType
proj_para <- addGroupCoverages(ArchRProj = proj_para, groupBy = "CellType")
proj_para <- addReproduciblePeakSet(
  ArchRProj = proj_para, 
  groupBy = "CellType", 
  pathToMacs2 = pathToMacs2, 
  force = T)
proj_para <- addPeakMatrix(proj_para)

proj_pool <- addGroupCoverages(ArchRProj = proj_pool, groupBy = "CellType")
proj_pool <- addReproduciblePeakSet(
  ArchRProj = proj_pool, 
  groupBy = "CellType", 
  pathToMacs2 = pathToMacs2,
  force = T)
proj_pool <- addPeakMatrix(proj_pool)

## By Density
proj_para <- addGroupCoverages(ArchRProj = proj_para, groupBy = "Density", force = T)
proj_para <- addReproduciblePeakSet(
  ArchRProj = proj_para, 
  groupBy = "Density", 
  pathToMacs2 = pathToMacs2, 
  force = T)
proj_para <- addPeakMatrix(proj_para)

proj_pool <- addGroupCoverages(ArchRProj = proj_pool, groupBy = "Density", force = T)
proj_pool <- addReproduciblePeakSet(
  ArchRProj = proj_pool, 
  groupBy = "Density", 
  pathToMacs2 = pathToMacs2,
  force = T)
proj_pool <- addPeakMatrix(proj_pool)

## By Replicate
proj_para <- addGroupCoverages(ArchRProj = proj_para, groupBy = "Calls", force = T)
proj_para <- addReproduciblePeakSet(
  ArchRProj = proj_para, 
  groupBy = "Calls", 
  pathToMacs2 = pathToMacs2, 
  force = T)
proj_para <- addPeakMatrix(proj_para)

proj_pool <- addGroupCoverages(ArchRProj = proj_pool, groupBy = "Calls", force = T,)
proj_pool <- addReproduciblePeakSet(
  ArchRProj = proj_pool, 
  groupBy = "Calls", 
  pathToMacs2 = pathToMacs2,
  force = T)
proj_pool <- addPeakMatrix(proj_pool)


## By Replicate
proj_para <- addGroupCoverages(ArchRProj = proj_para, groupBy = "Calls", force = T)
proj_para <- addReproduciblePeakSet(
  ArchRProj = proj_para, reproducibility = "1",
  groupBy = "Calls", 
  pathToMacs2 = pathToMacs2, 
  force = T)
proj_para <- addPeakMatrix(proj_para)

proj_pool <- addGroupCoverages(ArchRProj = proj_pool, groupBy = "Calls", force = T,)
proj_pool <- addReproduciblePeakSet(
  ArchRProj = proj_pool, reproducibility = "1",
  groupBy = "Calls", 
  pathToMacs2 = pathToMacs2,
  force = T)
proj_pool <- addPeakMatrix(proj_pool)





### SVD/LSI Correlation -----------------------------------------------------

proj_k562 <- subsetArchRProject(proj_2, 
                                cells = proj_2$cellNames[proj_2$CellType == "K562"],
                                dropCells = T,
                                outputDirectory = 'proj_k562',
                                force = T)
proj_k562 <- dimReduce(proj_k562, default = 'NoScale', dims = 2:30)

proj_jurkat <- subsetArchRProject(proj_2, 
                                  cells = proj_2$cellNames[proj_2$CellType == "Jurkat"],
                                  dropCells = T,
                                  outputDirectory = 'proj_jurkat',
                                  force = T)
proj_jurkat <- dimReduce(proj_jurkat, default = 'NoScale', dims = 2:30)


samples <- paste(rep(c("L","M","H"),each = 3), rep(1:3,3), sep="")
# samples <- paste("Bar",1:9,sep="")

para_svd <- lapply(c(proj_k562, proj_jurkat), function(x) {
  # SVD <- x@reducedDims$IterativeLSI@listData$matSVD
  SVD <- x@reducedDims$IterativeLSI@listData$matSVD[,2:30]
  meanSVD <- lapply(1:9, function(y) {
    SVD[x$Sample == "Para" & x$Calls == paste("Bar", y, sep=""),] %>% colMeans()
  }) 
  meanSVD <- do.call(cbind, meanSVD) %>% as.data.frame()
  # colnames(meanSVD) <- paste("Bar", 1:9, sep="")
  colnames(meanSVD) <- samples
  meanSVD
})
para_svd_cor <- lapply(para_svd, cor, method = "spearman")
# para_svd_cor <- lapply(para_svd, cor, method = "pearson")

pool_svd <- lapply(c(proj_k562, proj_jurkat), function(x) {
  # SVD <- x@reducedDims$IterativeLSI@listData$matSVD
  SVD <- x@reducedDims$IterativeLSI@listData$matSVD[,2:30]
  meanSVD <- lapply(1:9, function(y) {
    SVD[x$Sample == "Pool" & x$Calls == paste("Bar", y, sep=""),] %>% colMeans()
  }) 
  meanSVD <- do.call(cbind, meanSVD) %>% as.data.frame()
  # colnames(meanSVD) <- paste("Bar", 1:9, sep="")
  colnames(meanSVD) <- samples
  meanSVD
})
pool_svd_cor <- lapply(pool_svd, cor, method = "spearman")
# pool_svd_cor <- lapply(pool_svd, cor, method = "pearson")


heatmapcolors <- paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 100) %>% rev()

svd_cor <- c(para_svd_cor, pool_svd_cor)

plot_grid(plotlist = lapply(svd_cor, function(x) {
  # p <- pheatmap::pheatmap(x, breaks = seq(min(do.call(rbind, svd_cor)), 1, length = length(heatmapcolors) + 1), color = heatmapcolors, treeheight_row = 0)
  p <- pheatmap::pheatmap(x, breaks = seq(-1, 1, length = length(heatmapcolors) + 1), color = heatmapcolors, treeheight_row = 0)
  # p <- pheatmap::pheatmap(x, breaks = seq(min(do.call(rbind, svd_cor), na.rm=T), max(do.call(rbind, svd_cor), na.rm=T), length = length(heatmapcolors) + 1), color = heatmapcolors, treeheight_row = 0, na_col = NA)
  p$gtable
}), nrow = 2)

plot_grid(plotlist = lapply(svd_cor, function(x) {
  p <- pheatmap::pheatmap(x, breaks = seq(-1, 1, length = length(heatmapcolors) + 1), color = heatmapcolors, treeheight_row = 0,
                          annotation_col = data.frame(Calls = samples, Density = rep(c("Low","Medium","High"), each = 3)) %>% tibble::column_to_rownames(var = "Calls"),
                          annotation_colors = list(Density = colors_den), annotation_legend = F)
  p$gtable
}), nrow = 2)


# Jurkat (Main Figure)
pheatmap::pheatmap(svd_cor[[2]], breaks = seq(-1, 1, length = length(heatmapcolors) + 1), color = heatmapcolors, treeheight_row = 0, treeheight_col = 20, border_color = NA, cutree_cols = 3, cutree_rows = 3,
                   annotation_col = data.frame(Calls = samples, Ratio = rep(c("Low","Medium","High"), each = 3)) %>% tibble::column_to_rownames(var = "Calls"),
                   annotation_colors = list(Ratio = colors_den), annotation_legend = F, 
                   main = "Parallel - Jurkat", filename = "./Plots/SVDCor_Jurkat_Parallel.pdf", width = 3, height = 3)
pheatmap::pheatmap(svd_cor[[4]], breaks = seq(-1, 1, length = length(heatmapcolors) + 1), color = heatmapcolors, treeheight_row = 0, treeheight_col = 20, border_color = NA, cutree_cols = 3, cutree_rows = 3,
                   annotation_col = data.frame(Calls = samples, Ratio = rep(c("Low","Medium","High"), each = 3)) %>% tibble::column_to_rownames(var = "Calls"),
                   annotation_colors = list(Ratio = colors_den), annotation_legend = F, 
                   main = "Pooled - Jurkat", filename = "./Plots/SVDCor_Jurkat_Pooled.pdf", width = 3, height = 3)

# K562 (SI Figure)
pheatmap::pheatmap(svd_cor[[1]], breaks = seq(-1, 1, length = length(heatmapcolors) + 1), color = heatmapcolors, treeheight_row = 0, treeheight_col = 20, border_color = NA, cutree_cols = 3, cutree_rows = 3,
                   annotation_col = data.frame(Calls = samples, Ratio = rep(c("Low","Medium","High"), each = 3)) %>% tibble::column_to_rownames(var = "Calls"),
                   annotation_colors = list(Ratio = colors_den), annotation_legend = F, 
                   main = "Parallel - K562", filename = "./Plots/SVDCor_K562_Parallel.pdf", width = 3, height = 3)
pheatmap::pheatmap(svd_cor[[3]], breaks = seq(-1, 1, length = length(heatmapcolors) + 1), color = heatmapcolors, treeheight_row = 0, treeheight_col = 20, border_color = NA, cutree_cols = 3, cutree_rows = 3,
                   annotation_col = data.frame(Calls = samples, Ratio = rep(c("Low","Medium","High"), each = 3)) %>% tibble::column_to_rownames(var = "Calls"),
                   annotation_colors = list(Ratio = colors_den), annotation_legend = F, 
                   main = "Pooled - K562", filename = "./Plots/SVDCor_K562_Pooled.pdf", width = 3, height = 3)



### Clustering Differences --------------------------------------------------

plotEmbedding(proj_para, name = "Clusters", size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle("Clusters - Parallel")
plotEmbedding(proj_pool, name = "Clusters", size = 1, labelAsFactors = F, labelMeans = F) + theme_void() + ggtitle("Clusters - Pooled")

# proj_para$Calls <- as.character(proj_para$Calls)
# proj_pool$Calls <- as.character(proj_pool$Calls)

proj_para@cellColData[, c("Calls","Clusters")] %>% as.data.frame() %>% table() %>% pheatmap(scale = "column", cluster_rows = F, cluster_cols = F)
proj_pool@cellColData[, c("Calls","Clusters")] %>% as.data.frame() %>% table() %>% pheatmap(scale = "column", cluster_rows = F, cluster_cols = F)

proj_para@cellColData[, c("Calls","Clusters")] %>% as.data.frame() %>% table() %>% 
  pheatmap(scale = "column", cluster_rows = F, cluster_cols = F,
           annotation_row = data.frame(Calls = paste("Bar",1:9,sep=""), Density = rep(c("Low","Medium","High"), each = 3)) %>% tibble::column_to_rownames(var = "Calls"),
           annotation_colors = list(Density = colors_den), annotation_legend = F, 
           filename = "./Plots/projPara_ClusterHeatmap.pdf", height = 3, width = 3)
proj_pool@cellColData[, c("Calls","Clusters")] %>% as.data.frame() %>% table() %>% 
  pheatmap(scale = "column", cluster_rows = F, cluster_cols = F,
           annotation_row = data.frame(Calls = paste("Bar",1:9,sep=""), Density = rep(c("Low","Medium","High"), each = 3)) %>% tibble::column_to_rownames(var = "Calls"),
           annotation_colors = list(Density = colors_den), annotation_legend = F, 
           filename = "./Plots/projPool_ClusterHeatmap.pdf", height = 3, width = 3)

(proj_para@cellColData[, c("Density","Clusters")] %>% as.data.frame() %>% table())[1:3,] %>% pheatmap(scale = "column", cluster_rows = F, cluster_cols = F)
(proj_pool@cellColData[, c("Density","Clusters")] %>% as.data.frame() %>% table())[1:3,] %>% pheatmap(scale = "column", cluster_rows = F, cluster_cols = F)




############################################
### Local Inverse Simpson's Index (LISI) ###
############################################

## Assess degree of mixing of cells across a categorical variable

library(lisi)

res_para <- compute_lisi(proj_para@reducedDims$IterativeLSI$matSVD, as.data.frame(proj_para@cellColData), c("Density","Calls"), perplexity = 30)
res_pool <- compute_lisi(proj_pool@reducedDims$IterativeLSI$matSVD, as.data.frame(proj_pool@cellColData), c("Density","Calls"), perplexity = 30)
# res_para_h <- compute_lisi(proj_para@reducedDims$Harmony$matDR, as.data.frame(proj_para@cellColData), c("Density","Calls"))
colnames(res_para) <- paste("LISI_",colnames(res_para),sep="")
colnames(res_pool) <- paste("LISI_",colnames(res_pool),sep="")
# colnames(res_para_h) <- paste("LISI_",colnames(res_para_h),"_Harmony",sep="")

res <- rbind(cbind(as.data.frame(proj_para@cellColData), res_para)[,-34],
             cbind(as.data.frame(proj_pool@cellColData), res_pool))
ggplot(res, aes(x = Sample, y = LISI_Density, fill = Density)) + geom_boxplot(outliers = F) + facet_wrap(~CellType) + theme_minimal()
ggplot(res, aes(x = Sample, y = LISI_Calls, fill = Density)) + geom_boxplot(outliers = F) + facet_wrap(~CellType) + theme_minimal()
ggplot(res, aes(x = Sample, y = LISI_Calls, fill = Calls)) + geom_boxplot(outliers = F) + facet_wrap(~CellType) + theme_minimal()

# res <- cbind(as.data.frame(proj_para@cellColData), res_para, res_para_h)
# ggplot(res, aes(x = Sample, y = LISI_Calls, fill = Density)) + geom_boxplot() + facet_wrap(~CellType)
# ggplot(res, aes(x = Sample, y = LISI_Calls_Harmony, fill = Density)) + geom_boxplot() + facet_wrap(~CellType)

plotEmbedding(proj_para, embedding = "UMAP", name = "log10(nFrags)", size = 1, labelAsFactors = F, labelMeans = F) + theme_classic() + ggtitle("Density")

plotEmbedding(proj_para, embedding = "UMAP", name = "Density", size = 1, labelAsFactors = F, labelMeans = F) + theme_classic() + ggtitle("Density")
plotEmbedding(proj_para, embedding = "UMAP_Harmony", name = "Density", size = 1, labelAsFactors = F, labelMeans = F) + theme_classic() + ggtitle("Density")

plotEmbedding(proj_pool, embedding = "UMAP", name = "Density", size = 1, labelAsFactors = F, labelMeans = F) + theme_classic() + ggtitle("Density")





plot_grid(pheatmap::pheatmap(cor(as.data.frame(proj_para@cellColData[,c("nFrags", "TSSEnrichment", "NucleosomeRatio", "PromoterRatio", "nUMI", "FRIP")])),
                             color = paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100), scale = "none", main = "Parallel",
                             border_color = NA, breaks = seq(-1,1,0.02))$gtable,
          pheatmap::pheatmap(cor(as.data.frame(proj_pool@cellColData[,c("nFrags", "TSSEnrichment", "NucleosomeRatio", "PromoterRatio", "nUMI", "FRIP")])),
                             color = paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100), scale = "none", main = "Pooled",
                             border_color = NA, breaks = seq(-1,1,0.02))$gtable,
          ncol = 2)

plot_grid(pheatmap::pheatmap(cor(proj_para$nFrags, as.data.frame(proj_para@cellColData[,c("TSSEnrichment", "NucleosomeRatio", "PromoterRatio", "nUMI", "FRIP")])),
                             color = paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100), scale = "none", main = "Parallel",
                             border_color = NA, breaks = seq(-1,1,0.02), cluster_rows = F, cluster_cols = F)$gtable,
          pheatmap::pheatmap(cor(proj_pool$nFrags, as.data.frame(proj_pool@cellColData[,c("TSSEnrichment", "NucleosomeRatio", "PromoterRatio", "nUMI", "FRIP")])),
                             color = paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100), scale = "none", main = "Pooled",
                             border_color = NA, breaks = seq(-1,1,0.02), cluster_rows = F, cluster_cols = F)$gtable,
          nrow = 2)


proj_para@cellColData %>% as.data.frame() %>% group_by(Calls, Density) %>% dplyr::summarise(nFrags = median(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags, color= Density)) + geom_point() + 
  theme_bw() + labs(y = "Mean nFrags", x = "nCells") + ggtitle("9 Parallel Transposition Reactions") + ylim(c(500,1800)) + 
  scale_color_manual(values = colors_den)

proj_pool@cellColData %>% as.data.frame() %>% group_by(Calls, Density) %>% dplyr::summarise(nFrags = median(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags, color= Density)) + geom_point() + 
  theme_bw() + labs(y = "Mean nFrags", x = "nCells") + ggtitle("1 Pooled Transposition Reaction") + ylim(c(500,1800)) + 
  scale_color_manual(values = colors_den)

proj_para@cellColData %>% as.data.frame() %>% group_by(Density, Calls) %>% summarize(Count = n()) %>% ungroup %>% mutate(Rank = rank(Count, ties.method = "first")) %>%
  ggplot(., aes(x = Rank, y = Count, fill = Density)) + geom_col() + 
  theme_bw() + labs(x = "Rank", y = "nCells") + ggtitle("9 Parallel Transposition Reactions") + 
  scale_fill_manual(values = colors_den)
proj_pool@cellColData %>% as.data.frame() %>% group_by(Density, Calls) %>% summarize(Count = n()) %>% ungroup %>% mutate(Rank = rank(Count, ties.method = "first")) %>%
  ggplot(., aes(x = Rank, y = Count, fill = Density)) + geom_col() + 
  theme_bw()+ labs(x = "Rank", y = "nCells") + ggtitle("1 Pooled Transposition Reaction") + 
  scale_fill_manual(values = colors_den)



### Export R Objects --------------------------------------------------------

dir.create("Export")

sessionOut("Export/PooledParallel_Packages.tsv")

## Single-cell metadata
saveRDS(proj_2@cellColData, file = "Export/PooledParallel_SingleCellMetadata.rds")

