### Finalized Analysis of Batch Effects in Public Datasets ###

### Set Up Environment ------------------------------------------------------

library(ArchR)
library(Seurat)
library(Signac)
library(ggplot2)
library(ggpubr)
library(deMULTIplex2)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(cowplot)
library(lisi)
'%ni%' <- Negate('%in%')

source("/Volumes/DannySSD/MULTI_ATAC/CustomFunctions.R")

addArchRThreads(threads = 1) 

wd <- "/Volumes/DannySSD/MULTI_ATAC/CompetingMethods"
setwd(wd) 

dir.create("./Plots")

##########################################################################################################################################################

### Pre-Process Datasets ----------------------------------------------------

### SNuBar-ATAC -------------------------------------------------------------

## Simple oligonucleotide-based multiplexing of single-cell chromatin accessibility
## https://www.sciencedirect.com/science/article/pii/S1097276521007954
## Barcode oligos introduced into each transposition reaction
## Barcodes carried into 10x Genomics GEMs by hybridizing to free ends of transposed gDNA fragments

## 96-plex Cell Line Experiment
setwd("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/SNuBar")
addArchRGenome("hg19")

# snubar_csv <- read.csv("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/SNuBar/GSE162798_RAW/GSM4960035_Snubar-96plex.cell_metainfo.csv.gz")[-1,]
# snubar_csv <- snubar_csv[snubar_csv$cell_id != "None",]
# hist(snubar_csv$passed_filters %>% log10, breaks = 100)
# hist(proj_snu$nFrags %>% log10, breaks = 100)

snubar_arrow <- createArrowFiles(
  inputFiles = "/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/SNuBar/GSE162798_RAW/GSM4960035_Snubar-96plex.fragments.tsv.gz",
  sampleNames = "96plex",
  minTSS = 2,
  minFrags = 100,
  addTileMat = TRUE,
  addGeneScoreMat = F) 
snubar_arrow <- "96plex.arrow"

proj_snu <- ArchRProject(
  ArrowFiles = snubar_arrow, 
  outputDirectory = "proj_snu",
  copyArrows = F
)

proj_snu@cellColData %>% as.data.frame() %>%
  ggplot(aes(x = log10(nFrags), y = TSSEnrichment)) + 
  # geom_density2d() +
  geom_point(size = 0.1, alpha = 0.2) +
  geom_hline(yintercept = 6) + geom_vline(xintercept = log10(1000))

proj_snu <- proj_snu[proj_snu$nFrags > 1000]

proj_snu <- addIterativeLSI(proj_snu, force = T)
proj_snu <- addClusters(proj_snu, force = T)
proj_snu <- addUMAP(proj_snu, force = T)

plotEmbedding(proj_snu, name = "Clusters")
plotEmbedding(proj_snu, name = "log10(nFrags)", plotAs = 'points', size = 1)
plotEmbedding(proj_snu, name = "TSSEnrichment", plotAs = 'points', size = 1)

proj_snu$CB_CellRanger <- gsub(".*#", "", proj_snu$cellNames)
proj_snu$CB <- gsub("-.*", "", proj_snu$CB_CellRanger)

snubar_mtx <- read.table("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/SNuBar/GSE162798_RAW/GSM4960035_Snubar-96plex.snubar_count_matrix.mtx.gz")
snubar_mtx <- t(snubar_mtx[-97,])
tagHist(snubar_mtx) + theme(legend.position = 'none') + scale_y_sqrt()

snubar_res <- demultiplexTags(snubar_mtx,
                              plot.path = "/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/SNuBar",
                              plot.name = "96plex",
                              plot.diagnostics = T)
names(snubar_res$final_assign) %in% proj_snu$CB %>% table()

cbind(snubar_res$umap,snubar_res$assign_table) %>%
  {
    ggplot(., aes(x = UMAP_1, y = UMAP_2, color = final_assign)) + 
      # geom_point(data = g[g$droplet_type != 'singlet',], size = 0.5, alpha = 0) + geom_point(data = g[g$droplet_type == 'singlet',], size = 0.5, alpha = 0.5) +
      geom_point(data = filter(., droplet_type != 'singlet'), size = 0.3, alpha = 0.2) + 
      geom_point(data = filter(., droplet_type == 'singlet'), size = 1, alpha = 0.5) +
      scale_color_manual(values = c(colorRampPalette(c("navyblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1", "maroon"))(96),"black","grey")) + 
      theme_classic() + theme(legend.position = "none") + 
      labs(x = "UMAP1", y = "UMAP2", title = "deMULTIplex2 Classifications")
  }

proj_snu$nUMI <- rowSums(snubar_mtx)[proj_snu$CB]
proj_snu$MultipletType <- snubar_res$assign_table[proj_snu$CB, "barcode_count"]
proj_snu$Tn5 <- snubar_res$final_assign[proj_snu$CB]
proj_snu$Tn5 <- gsub("-.*","",proj_snu$Tn5)
proj_snu$CallType <- factor(proj_snu$Tn5, levels = c(paste("s",1:96,sep=""),"multiplet","negative"), labels = c(rep("singlet",96),"multiplet","negative"))
plotEmbedding(proj_snu, name = "Tn5", plotAs = 'points', size = 1) + theme(legend.position = "none")
plotEmbedding(proj_snu, name = "CallType", plotAs = 'points', size = 1)
plotEmbedding(proj_snu[!is.na(proj_snu$nUMI)], name = "log10(nUMI)", plotAs = 'points', size = 1) + theme(legend.position = "none")
plotEmbedding(proj_snu[!is.na(proj_snu$nUMI)], name = "MultipletType", plotAs = 'points', size = 1) + theme(legend.position = "none")

# Orthogonal Doublet Calling
proj_snu <- addDoubletScores(proj_snu)
plotEmbedding(proj_snu, name = "DoubletEnrichment")
plotEmbedding(proj_snu, name = "CallType")

# Filter project
proj_snu2 <- subsetArchRProject(proj_snu,
                                outputDirectory = "./proj_snu2",
                                cells = proj_snu$cellNames[proj_snu$Tn5 %ni% c("multiplet","negative",NA)], 
                                dropCells = F,
                                force = T)

# proj_snu2 <- addIterativeLSI(proj_snu2, force = T)
# proj_snu2 <- addClusters(proj_snu2, force = T)
# proj_snu2 <- addUMAP(proj_snu2, force = T)

proj_snu2 <- dimReduce(proj_snu2, what = c("iLSI","UMAP"), scale = c("Row","Col"), default = "Col", dims = 2:30)
proj_snu2 <- dimDefault(proj_snu2, what = c("iLSI","UMAP"), default = "NoScale")

plotEmbedding(proj_snu2, name = "Clusters")
plotEmbedding(proj_snu2, name = "Tn5") + theme(legend.position = 'none')
plotEmbedding(proj_snu2, name = "log10(nFrags)", plotAs = 'points', size = 1)
plotEmbedding(proj_snu2, name = "nFrags > 10000", plotAs = 'points', size = 1)
plotEmbedding(proj_snu2, name = "TSSEnrichment", plotAs = 'points', size = 1)
plotEmbedding(proj_snu2, name = "log10(nUMI)", plotAs = 'points', size = 1)
plotEmbedding(proj_snu2, name = "DoubletEnrichment", plotAs = 'points', size = 1)

g <- pheatmap::pheatmap(table(proj_snu2@cellColData[,c("Clusters","Tn5")]),scale = "row")

proj_snu2$CellType <- NA
proj_snu2$CellType[proj_snu2$Tn5 %in% g$tree_col$labels[g$tree_col$order][1:32]] <- "MDA-MB-231"
proj_snu2$CellType[proj_snu2$Tn5 %in% g$tree_col$labels[g$tree_col$order][33:64]] <- "K562"
proj_snu2$CellType[proj_snu2$Tn5 %in% g$tree_col$labels[g$tree_col$order][65:96]] <- "MDA-MB-436"

plotEmbedding(proj_snu2, name = "CellType", plotAs = 'points', size = 1)
proj_snu2@cellColData %>% as.data.frame() %>% ggplot(aes(x = log10(nFrags), y = TSSEnrichment, color = CellType)) + geom_point()
proj_snu2@cellColData %>% as.data.frame() %>% ggplot(aes(x = log10(nFrags), y = log10(nUMI), color = CellType)) + geom_point() + facet_wrap(~CellType)
proj_snu2@cellColData %>% as.data.frame() %>% ggplot(aes(x = log2(nUMI), fill = CellType)) + geom_histogram(bins = 40) + facet_wrap(~CellType) + theme_classic()
proj_snu2@cellColData %>% as.data.frame() %>% group_by(CellType) %>% summarize(nUMI = mean(nUMI), Count = n())


plot_grid(plotEmbedding(proj_snu2, embedding = "UMAP_NoScale", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_snu2, embedding = "UMAP_NoScale", name = "log10(nUMI)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_snu2, embedding = "UMAP_NoScale", name = "CellType", labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_snu2, embedding = "UMAP_Row", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_snu2, embedding = "UMAP_Row", name = "log10(nUMI)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_snu2, embedding = "UMAP_Row", name = "CellType", labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_snu2, embedding = "UMAP_Col", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_snu2, embedding = "UMAP_Col", name = "log10(nUMI)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_snu2, embedding = "UMAP_Col", name = "CellType", labelAsFactors = F, size = 1) + theme_void(),
          nrow = 3)

proj_snu2@cellColData %>% as.data.frame() %>% group_by(Tn5, CellType) %>% 
  summarize(Count = n(), nFrags = mean(nFrags), nUMI = mean(nUMI)) %>% 
  mutate(Tn5 = as.numeric(gsub("s","",Tn5))) %>%
  # mutate(Row = ceiling(Tn5 / 12), Col = (as.numeric(Tn5)-1) %% 12) %>%
  mutate(Row = (as.numeric(Tn5)-1) %% 8, Col = ceiling((as.numeric(Tn5)) / 8)) %>%
  # ggplot(aes(x = Col, y = -Row, fill = CellType, alpha = log10(Count), label = Tn5)) + geom_tile() + geom_text()
  # ggplot(aes(x = Col, y = -Row, fill = CellType, alpha = log10(nFrags), label = Tn5)) + geom_tile() + geom_text()
  ggplot(aes(x = Col, y = -Row, fill = CellType, alpha = log10(nUMI), label = Tn5)) + geom_tile() + geom_text()


proj_snu2@cellColData %>% as.data.frame() %>% group_by(Tn5, CellType) %>% summarize(nFrags = mean(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags, color = CellType)) + geom_point() + 
  theme_bw() + labs(y = "Mean nFrags", x = "nCells") + ggtitle("96 Transposition Reactions")
proj_snu2@cellColData %>% as.data.frame() %>% group_by(Tn5, CellType) %>% summarize(nFrags = mean(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags, color = CellType, label = Tn5)) + geom_text() + 
  theme_bw() + labs(y = "Mean nFrags", x = "nCells") + ggtitle("96 Transposition Reactions")
proj_snu2@cellColData %>% as.data.frame() %>% group_by(Tn5, CellType) %>% summarize(nFrags = mean(nFrags), Count = n(), nUMI = mean(nUMI)) %>% 
  ggplot(., aes(x = Count, y = nFrags, color = CellType, size = log10(nUMI))) + geom_point() + 
  theme_bw() + labs(y = "Mean nFrags", x = "nCells") + ggtitle("96 Transposition Reactions")
proj_snu2@cellColData %>% as.data.frame() %>% group_by(Tn5, CellType) %>% summarize(nFrags = mean(nFrags), Count = n(), nUMI = mean(nUMI)) %>% 
  ggplot(., aes(x = Count, y = nFrags, color = CellType, size = log10(nUMI))) + geom_point() + 
  theme_bw() + labs(y = "Mean nFrags", x = "nCells") + ggtitle("96 Transposition Reactions")
proj_snu2@cellColData %>% as.data.frame() %>% group_by(Tn5, CellType) %>% summarize(nFrags = mean(nFrags), Count = n(), nUMI = mean(nUMI)) %>% 
  ggplot(., aes(x = Count, y = log10(nUMI), color = CellType, size = log10(nUMI), label = Tn5)) + geom_text() + 
  theme_bw() + labs(y = "Mean nUMI", x = "nCells") + ggtitle("96 Transposition Reactions")


# filterDoublets(proj_snu2) %>% .@cellColData %>%
proj_snu2@cellColData %>%
  as.data.frame() %>% 
  # filter(Tn5 != "s94") %>%
  # filter(nFrags < 10000) %>%
  group_by(Tn5, CellType) %>% summarize(nFrags = mean(nFrags), nUMI = mean(nUMI), Count = n()) %>% 
  # filter(nUMI > 400) %>%
  ggplot(., aes(x = Count, y = log10(nFrags), color= CellType)) + geom_point() + geom_smooth(method = "lm", se=F) +
  theme_classic() + labs(y = "Mean nFrags", x = "nCells") #+ ggtitle("96 Transposition Reactions")

proj_snu2@cellColData %>% as.data.frame() %>% 
  ggplot(., aes(x = Tn5, y = nFrags, color= CellType)) + geom_boxplot(outliers = F) + scale_y_log10() + 
  theme_bw() + facet_wrap(CellType~., scales = "free", ncol = 1) + labs(y = "nFrags", x = "Tn5 Reaction") + ggtitle("96 Transposition Reactions")
proj_snu2@cellColData %>% as.data.frame() %>%
  ggplot(aes(x = Tn5, y = log10(nFrags), fill = CellType)) + geom_violin() + 
  theme_bw() + facet_wrap(~CellType, scales='free', ncol = 1)


proj_snu2@cellColData %>% as.data.frame() %>% group_by(Tn5, CellType) %>% summarize(Count = n()) %>% group_by(CellType) %>% mutate(Rank = rank(-Count, ties.method = "first")) %>%
  ggplot(., aes(x = Rank, y = Count, fill = CellType)) + geom_col()  + facet_wrap(.~CellType) +
  theme_bw() + labs(x = "Rank", y = "nCells") + ggtitle("96 Transposition Reactions")

proj_snu2@cellColData %>% as.data.frame() %>% group_by(Tn5, CellType) %>% summarize(nFrags = mean(nFrags), Count = n()) %>% group_by(CellType) %>% 
  summarize(Cor = cor(nFrags,Count,method = "spearman"),
            Slope = (lm(nFrags~Count) %>% coef())[2],
            Robust = (MASS::rlm(nFrags~Count) %>% coef())[2])
proj_snu2@cellColData[proj_snu2$Tn5 %ni% c("s47","s94","s95"),] %>% as.data.frame() %>% group_by(Tn5, CellType) %>% summarize(nFrags = mean(nFrags), Count = n()) %>% group_by(CellType) %>% 
  summarize(Cor = cor(nFrags,Count,method = "spearman"),
            Slope = (lm(nFrags~Count) %>% coef())[2])


## HBCA (dissected human breast tissue) Datasets

breast_arrow1 <- createArrowFiles(
  inputFiles = "/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/SNuBar/GSE162798_RAW/GSM4960040_Snubar-HBCA1.fragments.tsv.gz",
  sampleNames = "HBCA1",
  minTSS = 2,
  minFrags = 100,
  addTileMat = TRUE,
  addGeneScoreMat = FALSE) 
breast_arrow2 <- createArrowFiles(
  inputFiles = "/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/SNuBar/GSE162798_RAW/GSM4960041_Snubar-HBCA2.fragments.tsv.gz",
  sampleNames = "HBCA2",
  minTSS = 2,
  minFrags = 100,
  addTileMat = TRUE,
  addGeneScoreMat = FALSE) 

proj_breast <- ArchRProject(
  ArrowFiles = list.files(pattern = "HBCA.*arrow"), 
  outputDirectory = "proj_breast",
  copyArrows = F
)
proj_breast <- addIterativeLSI(proj_breast, force = T)
proj_breast <- addClusters(proj_breast, force = T)
proj_breast <- addUMAP(proj_breast, force = T)

plot_grid(plotEmbedding(proj_breast, name = "Sample", size = 1, labelAsFactors = F, labelMeans = F) + theme_void(),
          plotEmbedding(proj_breast, name = "Clusters", size = 1, labelAsFactors = F) + theme_void(),
          plotEmbedding(proj_breast, name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          nrow = 1)

proj_breast$CB_CellRanger <- gsub(".*#", "", proj_breast$cellNames)
proj_breast$CB <- gsub("-.*", "", proj_breast$CB_CellRanger)

breastbar_mtx1 <- read.table("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/SNuBar/GSE162798_RAW/GSM4960040_Snubar-HBCA1.snubar_count_matrix.mtx.gz")
breastbar_mtx1 <- t(breastbar_mtx1[-33,])
colnames(breastbar_mtx1) <- gsub("-.*","",colnames(breastbar_mtx1))
# tagHist(breastbar_mtx1) + theme(legend.position = 'none') + scale_y_sqrt()
breastbar_res1 <- demultiplexTags(breastbar_mtx1,
                                  plot.path = "/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/SNuBar",
                                  plot.name = "HBCA1",
                                  plot.diagnostics = T)

breastbar_mtx2 <- read.table("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/SNuBar/GSE162798_RAW/GSM4960041_Snubar-HBCA2.snubar_count_matrix.mtx.gz")
breastbar_mtx2 <- t(breastbar_mtx2[-33,])
colnames(breastbar_mtx2) <- gsub("-.*","",colnames(breastbar_mtx2))
# tagHist(breastbar_mtx2) + theme(legend.position = 'none') + scale_y_sqrt()
breastbar_res2 <- demultiplexTags(breastbar_mtx2,
                                  plot.path = "/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/SNuBar",
                                  plot.name = "HBCA2",
                                  plot.diagnostics = T)

names(breastbar_res1$final_assign) %in% proj_breast$CB %>% table()
names(breastbar_res2$final_assign) %in% proj_breast$CB %>% table()

proj_breast$Tn5 <- NA
proj_breast$Tn5[proj_breast$Sample == "HBCA1"] <- breastbar_res1$final_assign[proj_breast$CB[proj_breast$Sample == "HBCA1"]]
proj_breast$Tn5[proj_breast$Sample == "HBCA2"] <- breastbar_res2$final_assign[proj_breast$CB[proj_breast$Sample == "HBCA2"]]

plotEmbedding(proj_breast, name = "Tn5", plotAs = 'points', size = 1) + theme(legend.position = "none")


proj_breast2 <- subsetArchRProject(proj_breast,
                                   outputDirectory = "./proj_breast2",
                                   cells = proj_breast$cellNames[proj_breast$Tn5 %ni% c("multiplet","negative",NA)], 
                                   dropCells = F,
                                   force = T)

# proj_breast2 <- addIterativeLSI(proj_breast2, force = T)
# proj_breast2 <- addClusters(proj_breast2, force = T)
# proj_breast2 <- addUMAP(proj_breast2, force = T)

proj_breast2 <- dimReduce(proj_breast2, what = c("iLSI","UMAP"), scale = c("Row","Col"), default = "Col", dims = 2:30)
proj_breast2 <- dimDefault(proj_breast2, what = c("iLSI","UMAP"), default = "NoScale")

plot_grid(plotEmbedding(proj_breast2, embedding = "UMAP_NoScale", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_breast2, embedding = "UMAP_NoScale", name = "Region", labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_breast2, embedding = "UMAP_NoScale", name = "CellType", labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_breast2, embedding = "UMAP_Row", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_breast2, embedding = "UMAP_Row", name = "Region", labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_breast2, embedding = "UMAP_Row", name = "CellType", labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_breast2, embedding = "UMAP_Col", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_breast2, embedding = "UMAP_Col", name = "Region", labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_breast2, embedding = "UMAP_Col", name = "CellType", labelAsFactors = F, size = 1) + theme_void(),
          nrow = 3)


plot_grid(plotEmbedding(proj_breast2, name = "Sample", size = 1, labelAsFactors = F, labelMeans = F) + theme_void(),
          plotEmbedding(proj_breast2, name = "Clusters", size = 1, labelAsFactors = F) + theme_void(),
          plotEmbedding(proj_breast2, name = "Tn5", size = 1, labelAsFactors = F) + theme_void(),
          nrow = 1)

# proj_breast2@cellColData %>% as.data.frame() %>% group_by(Tn5, CellType) %>% dplyr::summarise(nFrags = mean(nFrags), Count = n()) %>%
#   ggplot(., aes(x = Count, y = nFrags, color= CellType)) + geom_point() + scale_x_log10() + scale_y_log10() + theme_bw()
proj_breast2@cellColData %>% as.data.frame() %>% group_by(Tn5) %>% dplyr::summarise(nFrags = mean(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags)) + geom_point() + 
  # scale_x_log10() + scale_y_log10() + 
  theme_bw()+ labs(y = "Mean nFrags", x = "nCells") + ggtitle("32 Transposition Reactions")

proj_breast2@cellColData %>% as.data.frame() %>% group_by(Tn5) %>% dplyr::summarize(nFrags = mean(nFrags), Count = n()) %>% mutate(Rank = dense_rank(-Count)) %>%
  ggplot(., aes(x = Rank, y = Count)) + geom_col() + 
  theme_bw()+ labs(x = "Rank", y = "nCells") + ggtitle("32 Transposition Reactions")

proj_breast2@cellColData %>% as.data.frame() %>% group_by(Tn5) %>% dplyr::summarise(nFrags = mean(nFrags), Count = n()) %>% 
  dplyr::summarise(Cor = cor(nFrags,Count,method = "spearman"))

proj_breast2 <- addGeneScoreMatrix(proj_breast2)
proj_breast2 <- addImputeWeights(proj_breast2)

plot_grid(plotlist = lapply(c(LS = "PROM1",
                              LHr = "ANKRD30A",
                              BS = "KRT17",
                              Adp = "ACACB",
                              `T` = "CD3G",
                              Myl = "HLA-DRA",
                              Endo = "VWF",
                              Fb = "COL1A2"), function(x) {
                                plotEmbedding(proj_breast2, colorBy = "GeneScoreMatrix", name = x, plotAs = 'points', size = 1) + theme_void() + ggtitle(x)
                              }), nrow = 2)

proj_breast2$CellType <- NA
proj_breast2$CellType[proj_breast2$Clusters %in% c("C3","C4","C5")] <- "LS"
proj_breast2$CellType[proj_breast2$Clusters %in% c("C2")] <- "LHr"
proj_breast2$CellType[proj_breast2$Clusters %in% c("C1")] <- "Bs"
proj_breast2$CellType[proj_breast2$Clusters %in% c("C9")] <- "Adp"
proj_breast2$CellType[proj_breast2$Clusters %in% c("C7")] <- "T"
proj_breast2$CellType[proj_breast2$Clusters %in% c("C6")] <- "Myl"
proj_breast2$CellType[proj_breast2$Clusters %in% c("C8")] <- "Endo"
proj_breast2$CellType[proj_breast2$Clusters %in% c("C10")] <- "Fb"

plotEmbedding(proj_breast2, name = "CellType", labelAsFactors = F) + theme_void() + ggtitle("Cell Types")

proj_breast2$Region <- "Epithelial"
proj_breast2$Region[proj_breast2$Tn5 %in% c("s14","s15","s16","s18","s19","s20","s26","s27","s31","s32")] <- "Adipose"

plotEmbedding(proj_breast2, name = "Region", labelAsFactors = F) + theme_void() + ggtitle("Cell Types")


proj_breast2@cellColData %>% as.data.frame() %>% group_by(Tn5) %>% dplyr::summarise(nFrags = mean(nFrags), Count = n(), percA = sum(CellType == "Adp") / n(), percE = sum(CellType %in% c("LS","LHr","Bs")) / n()) %>% 
  ggplot(., aes(x = Count, y = nFrags, color = percE)) + geom_point() + 
  theme_bw()+ labs(y = "mean nFrags", x = "nCells") + ggtitle("32 Transposition Reactions")
proj_breast2@cellColData %>% as.data.frame() %>% group_by(Tn5, Region) %>% dplyr::summarise(nFrags = mean(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags, color = Region)) + geom_point() + 
  theme_bw()+ labs(y = "mean nFrags", x = "nCells") + ggtitle("32 Transposition Reactions")
proj_breast2@cellColData %>% as.data.frame() %>% group_by(Tn5, Region) %>% dplyr::summarise(nFrags = mean(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags, color = Region)) + geom_point() + geom_smooth(method = "lm", se = F) + facet_wrap(.~Region, scales = "free") +
  theme_bw()+ labs(y = "mean nFrags", x = "nCells") + ggtitle("32 Transposition Reactions")

proj_breast2@cellColData %>% as.data.frame() %>% 
  # filter(CellType != "Adp") %>%
  # group_by(Tn5, Region, Sample) %>% 
  group_by(Tn5, Region) %>% 
  summarise(Count = n(), nFrags = mean(nFrags)) %>%
  ggplot(., aes(x = Count, y = nFrags, color = Region)) + #facet_wrap(~Sample) +
  # ggplot(., aes(x = Count, y = nFrags)) + 
  geom_point() + geom_smooth(method = "lm", se = F) + 
  theme_classic()+ labs(y = "mean nFrags", x = "nCells") #+ ggtitle("32 Transposition Reactions")

proj_breast2@cellColData %>% as.data.frame() %>% 
  # filter(CellType != "Adp") %>%
  # group_by(Tn5, Region, Sample) %>% 
  group_by(Tn5, Region) %>% 
  summarise(Count = n(), nFrags = median(nFrags)) %>%
  ggplot(., aes(x = Count, y = nFrags, color = Region)) + 
  geom_point() + geom_smooth(method = "lm", se = F) + 
  theme_classic()+ labs(y = "nFrags", x = "Count") + 
  ggtitle("SNuBar_HBCA - Standard")

proj_breast2@cellColData %>% as.data.frame() %>% group_by(Tn5) %>% mutate(Count = n()) %>% group_by(Tn5, CellType, Count) %>% dplyr::summarise(nFrags = mean(nFrags)) %>% group_by(CellType) %>%
  dplyr::summarise(Cor = cor(nFrags,Count,method = "spearman"))
proj_breast2@cellColData %>% as.data.frame() %>% group_by(Tn5) %>% mutate(Count = n()) %>% group_by(Tn5, Region, Count) %>% dplyr::summarise(nFrags = mean(nFrags)) %>% group_by(Region) %>%
  dplyr::summarise(Cor = cor(nFrags,Count,method = "spearman"))


### sciPlex-ATAC-seq --------------------------------------------------------

## High-capacity sample multiplexing for single cell chromatin accessibility profiling (2023)
## https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-023-09832-1

## 2-Level (sciPlex-ATAC) Chemical Perturbation (SciChem) Experiment
# Hash Label -> Pool/Split -> Barcoded Hash Capture + Barcoded Tn5 -> Pool/Split -> Indexed PCR
# A549 cells treated 24h with 8 concentrations of 4 different compounds in triplicate

sc2 <- read.table("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/sciPlex-ATAC-seq/GSE178953_RAW/GSM5401969_sciPlex2_cellMetaData.txt.gz", "\t", header = T, comment.char = "")
sc2_hash <- read.table("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/sciPlex-ATAC-seq/GSE178953_RAW/GSM5401969_sciPlex2_hashCellAssignments.txt.gz", "\t", header = T, comment.char = "")

ref <- readxl::read_excel("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/sciPlex-ATAC-seq/media-3.xlsx")
ref$Tn5 <- ref$`Sequence (5'-3')` %>% gsub(".*CCACGC|.*CCTGTCC","",.) %>% gsub("GCGATC.*|CACCGT.*","",.)
ref <- ref[1:20,]

# breakdown of cell barcode structure (Tn5 Oligos in Supp Table 4)
sc2$barcode <- gsub(".*#","",sc2$cell)
lapply(2:36, function(x) { data.frame(End = x, Unique = length(unique(str_sub(sc2$barcode, 1, x)))) } ) %>% do.call(rbind,.) %>% 
  ggplot(., aes(x = End, y = Unique)) + geom_line() + geom_point() + geom_vline(xintercept = c(8, 18, 28, 36)) + scale_y_log10()
table(str_sub(sc2$barcode, 1, 8)) %>% length() # Tn5 i7 Oligo, 8 nt, 12 unique (i.e. revComp(ATTACTCG))
table(str_sub(sc2$barcode, 9, 18)) %>% length() # PCR i7 Oligo, 10 nt, 12 unique
table(str_sub(sc2$barcode, 19, 28)) %>% length() # PCR i5 Oligo, 10 nt, 56 unique
table(str_sub(sc2$barcode, 29, 36)) %>% length() # Tn5 i5 Oligo, 8 nt, 8 unique (i.e. TATAGCCT)

str_sub(sc2$barcode, 1, 8) %>% unique %in% revComplement(ref$Tn5, cbLength = 8)
str_sub(sc2$barcode, 29,36) %>% unique %in% ref$Tn5

table(paste(str_sub(sc2$barcode, 1, 8), str_sub(sc2$barcode, 29, 36))) %>% length() # 87 Tn5 barcode combos in dataset
table(paste(str_sub(sc2$barcode, 9, 18), str_sub(sc2$barcode, 19, 28))) %>% length() # 471 PCR index combos in dataset

# breakdown of wellID (e.g. "A04_F04_RT_29")
sc2$wellID <- sc2_hash$wellID[match(sc2$barcode, sc2_hash$cell)]
strsplit(sc2_hash$wellID, split="_") %>% do.call(rbind,.) %>% apply(., 2, function(x) length(table(x)) )
strsplit(sc2_hash$wellID, split="_") %>% do.call(rbind,.) %>% apply(., 2, function(x) table(x) )

# first two indices of wellID refer to the PCR indices used
# RT_* refers to the Tn5 well
table(data.frame(CellBarcode = paste(str_sub(sc2$barcode, 1, 8),str_sub(sc2$barcode, 29, 36)),
                 WellID = str_sub(sc2$wellID, 9, 13))) %>% pheatmap::pheatmap(cluster_rows = T, cluster_cols = T, show_rownames = F)
table(data.frame(CellBarcode = paste(str_sub(sc2$barcode, 9, 18),str_sub(sc2$barcode, 19, 28)),
                 WellID = str_sub(sc2$wellID, 1, 7))) %>% pheatmap::pheatmap(cluster_rows = T, cluster_cols = T, show_rownames = F)

sc2$Sample <- paste(sc2$treatment, sc2$dose, sc2$replicate, sep = "_")
sc2$Tn5 <- paste(str_sub(sc2$barcode, 1, 8),str_sub(sc2$barcode, 29, 36), sep = "_")
sc2$PCR <- paste(str_sub(sc2$barcode, 9, 18),str_sub(sc2$barcode, 19, 28), sep = "_")

sc2$Relative_dose <- factor(sc2$Relative_dose, levels = sort(unique(sc2$Relative_dose), decreasing = F))


# Plot sample metadata by treatment & relative dose
sc2 %>% group_by(treatment, dose, replicate, Relative_dose, vehicle) %>% summarize(nFrags = mean(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Relative_dose, y = Count, fill = Relative_dose)) + geom_col() +
  facet_grid(~treatment) + scale_fill_viridis_d(option = "B") + theme_bw()
ggplot(sc2, aes(x = Relative_dose, y = log10(nFrags), fill = Relative_dose)) + geom_boxplot(outliers = F) + 
  facet_grid(~treatment) + scale_fill_viridis_d(option = "B") + theme_bw()
ggplot(sc2, aes(x = Relative_dose, y = FRIP, fill = Relative_dose)) + geom_boxplot(outliers = F) + 
  facet_grid(~treatment) + scale_fill_viridis_d(option = "B") + theme_bw()
ggplot(sc2, aes(x = Relative_dose, y = TSSEnrichment, fill = Relative_dose)) + geom_boxplot(outliers = F) + 
  facet_grid(~treatment) + scale_fill_viridis_d(option = "B") + theme_bw()


# Plot Counts & mean nFrags by Tn5, PCR, or Sample
sc2 %>% group_by(Tn5) %>% summarize(nFrags = mean(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags)) + geom_point() + geom_smooth(method = 'lm', se = F) + 
  theme_classic()
sc2 %>% group_by(PCR) %>% summarize(nFrags = mean(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags)) + geom_point() + geom_smooth(method = 'lm', se = F) +
  theme_minimal()
sc2 %>% group_by(treatment, Sample) %>% summarize(nFrags = mean(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags, color = treatment)) + geom_point(size = 3) + geom_smooth(method = 'lm', se = F) +
  theme_minimal()
sc2 %>% group_by(Relative_dose, Sample, replicate) %>% summarize(nFrags = mean(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags, color = Relative_dose)) + geom_point(size = 3) +
  theme_minimal() + scale_color_manual(values = viridis::inferno(n = 9)[1:8])

# Correlation testing and linear regression modeling
sc2 %>% group_by(Tn5) %>% summarize(nFrags = mean(nFrags), Count = n()) %>% summarize(p.value = cor.test(nFrags, Count, alternative = "t", method = "pearson", conf.level = 0.95)$p.value,
                                                                                      cor = cor.test(nFrags, Count, alternative = "t", method = "pearson", conf.level = 0.95)$estimate,
                                                                                      slope = (lm(nFrags~Count) %>% coef())[2],
                                                                                      Robust = (MASS::rlm(nFrags~Count) %>% coef())[2])
sc2 %>% group_by(PCR) %>% summarize(nFrags = mean(nFrags), Count = n()) %>% summarize(p.value = cor.test(nFrags, Count, alternative = "t", method = "pearson", conf.level = 0.95)$p.value,
                                                                                      cor = cor.test(nFrags, Count, alternative = "t", method = "pearson", conf.level = 0.95)$estimate,
                                                                                      slope = (lm(nFrags~Count) %>% coef())[2],
                                                                                      Robust = (MASS::rlm(nFrags~Count) %>% coef())[2])
sc2 %>% group_by(treatment,Sample) %>% summarize(nFrags = mean(nFrags), Count = n()) %>% summarize(p.value = cor.test(nFrags, Count, alternative = "t", method = "pearson", conf.level = 0.95)$p.value,
                                                                                                   cor = cor.test(nFrags, Count, alternative = "t", method = "pearson", conf.level = 0.95)$estimate,
                                                                                                   slope = (lm(nFrags~Count) %>% coef())[2],
                                                                                                   Robust = (MASS::rlm(nFrags~Count) %>% coef())[2])

## Analysis of SciChem dataset with ArchR
wd <- "/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/sciPlex-ATAC-seq/ArchR"
dir.create(wd)
setwd(wd)

# Create Arrow Files
# inputFiles <- list.files(path = "../GSE178953_RAW/", pattern = ".*MLR.*.fragments.txt.gz", full.names = T, recursive = T)
inputFiles <- list.files(path = "../GSE178953_RAW/", pattern = ".*SC.*.fragments.txt.gz$", full.names = T, recursive = T)
names(inputFiles) <- gsub(".frag.*","",gsub("^.*/.*/.*?_","",inputFiles))

addArchRGenome("hg19")

ArrowFiles <- vector()
for (i in 1:length(inputFiles)) {
  arrow <- createArrowFiles(
    inputFiles = inputFiles[i],
    sampleNames = names(inputFiles)[i],
    minTSS = 0,
    minFrags = 100,
    addTileMat = F,
    addGeneScoreMat = F)
  ArrowFiles[i] <- arrow
  rm(arrow)
}
ArrowFiles <- "SC2.arrow"

# Create ArchR Project
proj_plex <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "proj_plex",
  copyArrows = F
)
proj_plex <- addTileMatrix(proj_plex, force = T)
proj_plex <- addIterativeLSI(proj_plex, force = T)
proj_plex <- addClusters(proj_plex, force = T)
proj_plex <- addUMAP(proj_plex, force = T)

plot_grid(plotEmbedding(proj_plex, labelMeans = T, labelAsFactors = F) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_plex, name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_plex, name = "TSSEnrichment", plotAs = 'points', size = 1) + theme_void(),
          nrow = 1)

proj_plex$barcode <- gsub(".*#","",proj_plex$cellNames)

proj_plex2 <- subsetArchRProject(proj_plex, 
                                 cells = toupper(sc2$cell), 
                                 outputDirectory = "proj_plex2",
                                 dropCells = F, force = T)

# proj_plex2 <- addIterativeLSI(proj_plex2, force = T)
# proj_plex2 <- addClusters(proj_plex2, dimsToUse = 2:30, force = T)
# proj_plex2 <- addUMAP(proj_plex2, dimsToUse = 2:30, force = T)

proj_plex2 <- dimReduce(proj_plex2, what = c("iLSI","UMAP"), scale = c("Row","Col"), default = "Col", dims = 2:30)
proj_plex2 <- dimDefault(proj_plex2, what = c("iLSI","UMAP"), default = "NoScale")

proj_plex2@cellColData <- left_join(proj_plex2@cellColData %>% as.data.frame %>% rownames_to_column(), sc2, by = "barcode", suffix = c("","_old")) %>% column_to_rownames() %>% DataFrame

plot_grid(plotEmbedding(proj_plex2, name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_plex2, name = "TSSEnrichment", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_plex2, name = "FRIP", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_plex2, name = "Clusters", labelMeans = T, labelAsFactors = F, size = 1) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_plex2, name = "treatment", labelMeans = T, labelAsFactors = F, size = 1) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_plex2, name = "Relative_dose", labelMeans = F, labelAsFactors = F, size = 1) + theme_void() + scale_color_viridis_d(option = "B"),
          nrow = 2)

plot_grid(plotEmbedding(proj_plex2, embedding = "UMAP_NoScale", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_plex2, embedding = "UMAP_NoScale",name = "treatment", labelMeans = T, labelAsFactors = F, size = 1) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_plex2, embedding = "UMAP_NoScale",name = "Relative_dose", labelMeans = F, labelAsFactors = F, size = 1) + theme_void() + scale_color_viridis_d(option = "B"),
          plotEmbedding(proj_plex2, embedding = "UMAP_Row", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_plex2, embedding = "UMAP_Row",name = "treatment", labelMeans = T, labelAsFactors = F, size = 1) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_plex2, embedding = "UMAP_Row",name = "Relative_dose", labelMeans = F, labelAsFactors = F, size = 1) + theme_void() + scale_color_viridis_d(option = "B"),
          plotEmbedding(proj_plex2, embedding = "UMAP_Col", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_plex2, embedding = "UMAP_Col",name = "treatment", labelMeans = T, labelAsFactors = F, size = 1) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_plex2, embedding = "UMAP_Col",name = "Relative_dose", labelMeans = F, labelAsFactors = F, size = 1) + theme_void() + scale_color_viridis_d(option = "B"),
          nrow = 3)


### sci-ATAC-seq3 -----------------------------------------------------------

## A human cell atlas of fetal chromatin accessibility (2020)
## https://www.science.org/doi/10.1126/science.aba7612
## Samples from 23 fetuses 89-125 days gestational age
## Processed in 3 batches, each with a control mix of sentinel human fetal brain tissues and a mouse suspension cell line
## Per batch, 200,000 nuclei from each of 24 samples where tagmented in 4 replicate rxns (50k each), then replicates were pooled and split into 16 wells for ligation
## Each entry in table represents the aggregation of 4 replicate Tn5 reactions (also 16 ligation reactions)
## Number of cells reflects how many were recovered after full sci-ATAC-seq3 protocol, but presumably is correlated with original loading in well

sci3_meta <- readxl::read_excel("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/sci-ATAC-seq3/aba7612_domcke_table-s1.xlsx", skip = 1) %>% as.data.frame()
sci3_meta <- sci3_meta[1:72,]
sci3_meta <- sci3_meta[sci3_meta$`Included_in_downstream_analysis (yes/no)` == "Y",]

ggplot(sci3_meta, aes(x = Number_of_cells, y = Median_total_fragments)) + geom_point() + 
  scale_x_log10() + scale_y_log10() +
  facet_grid(~Batch) + theme_bw()+ labs(y = "Median nFrags", x = "nCells") + ggtitle("288 Transposition Reactions",subtitle = "3 Batches of 24 Samples * 4 Replicates (Aggregated)")

sci3_meta %>% group_by(Batch) %>% summarize(cor = cor(Median_total_fragments, Number_of_cells, method = "spearman"))
sci3_meta %>% group_by(Batch) %>% mutate(Rank = rank(-Number_of_cells, ties.method = "first")) %>%
  ggplot(., aes(x = Rank, y = Number_of_cells, fill = Batch)) + geom_col() + facet_wrap(.~Batch) + 
  theme_bw() + labs(x = "Rank", y = "nCells") + ggtitle("288 Transposition Reactions")

## Load and concatentate metadata from all Seurat Objects
sci3_fet <- list.files(path = "/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/sci-ATAC-seq3/GSE149683_RAW", full.names = T)

sci3_fet <- lapply(sci3_fet, function(x) {
  seu <- readRDS(x)
  seu@meta.data
}) %>% do.call(rbind, .)

ggplot(sci3_fet, aes(x = tissue_umap_1, y = tissue_umap_2, color = sample_name)) + geom_point() + facet_wrap(.~tissue, scales = "free") + theme_minimal() + theme(legend.position = "none")
ggplot(sci3_fet, aes(x = tissue_umap_1, y = tissue_umap_2, color = log10(total_deduplicated))) + geom_point() + facet_wrap(.~tissue, scales = "free") + theme_minimal() + scale_color_viridis_c()

## Figure out and break down sci-ATAC-seq3 cell barcode into constituent parts to identify well-specific barcode
# from oligos table-s7:
# N5 oligo barcode = 10nt
# N7 oligo barcode = 10nt
# P5 barcode = 10nt
# P7 barcode = 10nt

lapply(2:40, function(x) { data.frame(End = x, Unique = length(unique(str_sub(sci3_fet$cell, 1, x)))) } ) %>% do.call(rbind,.) %>% 
  ggplot(., aes(x = End, y = Unique)) + geom_line() + geom_point() + geom_vline(xintercept = c(10, 20, 30, 40)) + scale_y_log10()

exl <- readxl::read_xlsx("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/sci-ATAC-seq3/aba7612_domcke_table-s7.xlsx",sheet = 3) ## N5 Barcode, Plate 1/4
exl <- readxl::read_xlsx("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/sci-ATAC-seq3/aba7612_domcke_table-s7.xlsx",sheet = 7) ## N7 Barcode, Plate 1/4
exl <- readxl::read_xlsx("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/sci-ATAC-seq3/aba7612_domcke_table-s7.xlsx",sheet = 11) ## P5 Barcode, Plate 1/4
exl <- readxl::read_xlsx("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/sci-ATAC-seq3/aba7612_domcke_table-s7.xlsx",sheet = 12) ## P7 Barcode, Plate 1/4

revComplement(toupper(exl$Barcode), cbLength = 10) %in% substr(sci3_fet$cell, 1, 10) ## N7 Barcode (2nd Barcode), 384 unique
revComplement(toupper(exl$Barcode), cbLength = 10) %in% substr(sci3_fet$cell, 11, 20) ## P7 Barcode (3rd Barcode), 96 unique
toupper(exl$Barcode) %in% substr(sci3_fet$cell, 21, 30) ## P5 Barcode (3rd Barcode), 96 unique
toupper(exl$Barcode) %in% substr(sci3_fet$cell, 31, 40) ## N5 Barcode (1st Barcode), 384 unique

substr(sci3_fet$cell, 1, 10) %>% table() %>% length()
substr(sci3_fet$cell, 11, 20) %>% table() %>% length()
substr(sci3_fet$cell, 21, 30) %>% table() %>% length()
substr(sci3_fet$cell, 31, 40) %>% table() %>% length()

revComplement(unique(substr(sci3_fet$cell, 11, 20)), cbLength = 10) %in% toupper(readxl::read_xlsx("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/sci-ATAC-seq3/aba7612_domcke_table-s7.xlsx",sheet = 12)$Barcode)
unique(substr(sci3_fet$cell, 21, 30)) %in% toupper(readxl::read_xlsx("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/sci-ATAC-seq3/aba7612_domcke_table-s7.xlsx",sheet = 11)$Barcode)

sci3_fet$Lig1 <- str_sub(sci3_fet$cell, 31, 40)
sci3_fet$Lig2 <- str_sub(sci3_fet$cell, 1, 10)
sci3_fet$PCR <- paste(str_sub(sci3_fet$cell, 11, 20), str_sub(sci3_fet$cell, 21, 30), sep = "_")

table(sci3_fet$Lig1) %>% sort(decreasing = T) %>% head()
sci3_fet[sci3_fet$Lig1 == "TGGATCAGGC", c("batch","sample_name")] %>% table()
sci3_fet[sci3_fet$Lig1 == "ACCTAGGAGA", c("batch","sample_name")] %>% table()

sci3_fet %>% group_by(sample_name, batch) %>% summarize(nFrags = mean(total_deduplicated), Count = n()) %>%
  ggplot(., aes(x = Count, y = nFrags, color = batch)) + geom_point() + geom_smooth(method = "lm", se = F) +
  # scale_x_log10() + scale_y_log10() + 
  # facet_wrap(~batch, scales = 'free') + 
  theme_bw() + labs(y = "mean nFrags", x = "nCells") + ggtitle("Sample Aggregates")
sci3_fet %>% group_by(Lig1, batch) %>% summarize(nFrags = mean(total_deduplicated), Count = n()) %>%
  ggplot(., aes(x = Count, y = nFrags, color = batch)) + geom_point() + geom_smooth(method = "lm", se = F) +
  scale_x_log10() + scale_y_log10() +
  # facet_wrap(~batch, scales = 'free') + 
  theme_bw() + labs(y = "mean nFrags", x = "nCells") + ggtitle("1st-Round Ligation Reactions")
sci3_fet %>% group_by(Lig2, batch) %>% summarize(nFrags = mean(total_deduplicated), Count = n()) %>%
  ggplot(., aes(x = Count, y = nFrags, color = batch)) + geom_point() + geom_smooth(method = "lm") +
  scale_x_log10() + scale_y_log10() +
  # facet_wrap(~batch, scales = 'free') + 
  theme_bw() + labs(y = "mean nFrags", x = "nCells") + ggtitle("2nd-Round Ligation Reactions")
sci3_fet %>% group_by(PCR, batch) %>% summarize(nFrags = mean(total_deduplicated), Count = n()) %>%
  ggplot(., aes(x = Count, y = nFrags, color = batch)) + geom_point() + geom_smooth(method = "lm") +
  scale_x_log10() + scale_y_log10() +
  # facet_grid(~batch) + 
  theme_bw() + labs(y = "mean nFrags", x = "nCells") + ggtitle("3rd-Round PCR Reactions")

goodPCR <- sci3_fet %>% group_by(PCR, batch) %>% summarize(Count = n())
ggplot(goodPCR, aes(x = Count, fill = batch)) + geom_histogram() + facet_wrap(~batch, scales = "free") + geom_vline(xintercept = 100)
goodPCR <- goodPCR[goodPCR$Count >= 100,]

ggplot(sci3_fet %>% head(n = 200000), aes(x = tissue_umap_1, y = tissue_umap_2, color = PCR %in% goodPCR$PCR)) + geom_point(size = 0.3) + facet_wrap(.~tissue, scales = "free") + theme_minimal() + theme(legend.position = "none")
ggplot(sci3_fet, aes(x = tissue, y = log10(total_deduplicated), fill = PCR %in% goodPCR$PCR)) + geom_boxplot(outliers=F) + theme_minimal()

sci3_fet_sub <- sci3_fet[paste(sci3_fet$PCR, sci3_fet$batch) %in% paste(goodPCR$PCR, goodPCR$batch),]

sci3_fet_sub %>% group_by(batch) %>% summarize(Unique = length(unique(PCR)))

sci3_fet_sub %>% group_by(sample_name, batch) %>% summarize(nFrags = mean(total_deduplicated), Count = n()) %>%
  ggplot(., aes(x = Count, y = nFrags, color = batch)) + geom_point() + geom_smooth(method = "lm", se = F) +
  # scale_x_log10() + scale_y_log10() + 
  # facet_wrap(~batch, scales = 'free') + 
  theme_bw() + labs(y = "mean nFrags", x = "nCells") + ggtitle("Sample Aggregates")
sci3_fet_sub %>% group_by(Lig1, batch) %>% summarize(nFrags = mean(total_deduplicated), Count = n()) %>%
  ggplot(., aes(x = Count, y = nFrags, color = batch)) + geom_point() + geom_smooth(method = "lm", se = F) +
  scale_x_log10() + scale_y_log10() +
  # facet_wrap(~batch, scales = 'free') + 
  theme_bw() + labs(y = "mean nFrags", x = "nCells") + ggtitle("1st-Round Ligation Reactions")
sci3_fet_sub %>% group_by(Lig2, batch) %>% summarize(nFrags = mean(total_deduplicated), Count = n()) %>%
  ggplot(., aes(x = Count, y = nFrags, color = batch)) + geom_point() + geom_smooth(method = "lm") +
  scale_x_log10() + scale_y_log10() +
  # facet_wrap(~batch, scales = 'free') + 
  theme_bw() + labs(y = "mean nFrags", x = "nCells") + ggtitle("2nd-Round Ligation Reactions")
sci3_fet_sub %>% group_by(PCR, batch) %>% summarize(nFrags = mean(total_deduplicated), Count = n()) %>%
  ggplot(., aes(x = Count, y = nFrags, color = batch)) + geom_point() + geom_smooth(method = "lm") +
  scale_x_log10() + scale_y_log10() +
  # facet_grid(~batch) + 
  theme_bw() + labs(y = "mean nFrags", x = "nCells") + ggtitle("3rd-Round PCR Reactions")

## Each set of 4 replicates was pooled and split across 16 first-round ligation reactions, so individual Tn5 wells cannot be discerned
sci3_fet %>% group_by(Lig1, batch) %>% summarize(nFrags = mean(total_deduplicated), Count = n()) %>%
  ggplot(., aes(x = Count, y = nFrags)) + geom_point() + geom_smooth(method = "lm") +
  scale_x_log10() + scale_y_log10() + facet_grid(~batch) + theme_bw() + 
  labs(y = "mean nFrags", x = "nCells") + ggtitle("Individual 1st-Round Ligation Reactions",subtitle = "3 Batches of 24 Samples * 4 Tn5 Replicates (Pooled & Split Across 16 Ligation Rxns)")
sci3_fet %>% group_by(Lig1, batch, sample_name) %>% summarize(nFrags = mean(total_deduplicated), Count = n()) %>%
  ggplot(., aes(x = Count, y = nFrags, color = sample_name)) + geom_point() + geom_smooth(method = "lm",se = F) +
  # scale_x_log10() + scale_y_log10() + 
  facet_wrap(~batch, scales = 'free') + theme_bw() + theme(legend.position = 'none') +
  labs(y = "mean nFrags", x = "nCells") + ggtitle("Individual Ligation Reactions",subtitle = "3 Batches of 24 Samples * 4 Replicates (Pooled & Split Across 16 Ligation Rxns)")

## Group by sample instead - aggregate value of 4 transposition replicates and 16 first-round ligations
sci3_fet %>% group_by(batch, sample_name) %>% summarize(nFrags = mean(total_deduplicated), Count = n()) %>%
  ggplot(., aes(x = Count, y = nFrags, color = batch)) + geom_point() + geom_smooth(method = "lm", se = F) +
  # scale_x_log10() + scale_y_log10() + 
  facet_wrap(~batch, scales = "free") + theme_bw() + 
  labs(y = "mean nFrags", x = "nCells") + ggtitle("288 Transposition Reactions",subtitle = "3 Batches of 24 Samples * 4 Replicates (Aggregated)")
sci3_fet %>% group_by(batch, sample_name) %>% summarize(nFrags = mean(total_deduplicated), Count = n()) %>%
  summarize(Cor = cor(nFrags,Count,method = "spearman"),
            Slope = (lm(nFrags~Count) %>% coef())[2],
            Robust = (MASS::rlm(nFrags~Count) %>% coef())[2])
sci3_fet %>% group_by(batch, sample_name) %>% summarize(nFrags = mean(total_deduplicated), Count = n()) %>%
  ggplot(., aes(x = Count, y = nFrags, color = batch)) + geom_point() + geom_smooth(method = "lm", se = F) +
  # scale_x_log10() + scale_y_log10() + 
  # facet_wrap(~batch, scales = "free") + 
  theme_classic() + labs(y = "mean nFrags", x = "nCells")


## The continuum of Drosophila embryonic development at single-cell resolution
## https://www.science.org/doi/10.1126/science.abn5800

sci3_dros <- readRDS("/Volumes/DannySSD/MULTI_ATAC/OtherDatasets/DrosophilaDev/ATAC_xx_xy_annotated.Rds")
sci3_dros_meta <- readRDS("/Volumes/DannySSD/MULTI_ATAC/OtherDatasets/DrosophilaDev/atac_meta.rds")

sci3_dros$time %>% table()
head(sci3_dros$cell)

sci3_dros <- column_to_rownames(sci3_dros, "cell")[sci3_dros_meta$cell,] %>% rownames_to_column("cell")

sci3_dros$CellType <- sci3_dros_meta[sci3_dros$cell,"refined_annotation"]
table(sci3_dros[,c("time","CellType")]) %>% heatmap()

sci3_dros[1:10000,] %>% {lapply(2:40, function(x) { data.frame(End = x, Unique = length(unique(str_sub(.$cell, 1, x)))) } )} %>% 
  do.call(rbind,.) %>% 
  ggplot(., aes(x = End, y = Unique)) + geom_line() + geom_point() + geom_vline(xintercept = c(10, 20, 30, 40)) + scale_y_log10()

str_sub(sci3_dros$cell, 1, 10) %>% table() %>% length() # 384 detected
str_sub(sci3_dros$cell, 11, 20) %>% table() %>% length() # 96 detected
str_sub(sci3_dros$cell, 21, 30) %>% table() %>% length() # 64 detected
str_sub(sci3_dros$cell, 31, 40) %>% table() %>% length() # 352 detected

sci3_dros$Lig1 <- str_sub(sci3_dros$cell, 31, 40)
sci3_dros$Lig2 <- str_sub(sci3_dros$cell, 1, 10)
sci3_dros$PCR <- paste(str_sub(sci3_dros$cell, 11, 20), str_sub(sci3_dros$cell, 21, 30), sep = "_")

sci3_dros$Batch <- gsub("_.*","",sci3_dros$time)
sci3_dros$Timepoint <- sci3_dros$time %>% gsub("exp._","",.) %>% gsub("_.*","",.)

table(sci3_dros[,c("Lig1","time")]) %>% heatmap()
table(sci3_dros[,c("Lig2","time")]) %>% heatmap()
table(sci3_dros[,c("PCR","time")]) %>% heatmap()

table(sci3_dros[,c("Lig1","Timepoint")]) %>% heatmap()
table(sci3_dros[,c("Lig2","Timepoint")]) %>% heatmap()
table(sci3_dros[,c("PCR","Timepoint")]) %>% heatmap()

# looks like each timepoint/sample was divided into 44 unique ligation reactions
sci3_dros %>% group_by(Batch, Timepoint) %>% summarize(nLig = length(unique(Lig1)))

sci3_dros %>% group_by(Batch, Timepoint) %>% summarize(Count = n(), nFrags = median(total_reads_in_peaks)) %>%
  ggplot(aes(x = Count, y = nFrags, color = Batch)) + geom_point() + geom_smooth(method = 'lm', se = F) + 
  scale_x_log10() +
  theme_classic() + ggtitle("sci-ATAC-seq3 - Drosophila Development")
sci3_dros %>% group_by(Batch, Timepoint) %>% summarize(Count = n(), nFrags = median(total_reads_in_peaks)) %>%
  ggplot(aes(x = Count, y = nFrags)) + geom_point() + geom_smooth(method = 'lm', se = F) + 
  scale_x_log10() +
  theme_classic() + ggtitle("sci-ATAC-seq3 - Drosophila Development")



### sci-ATAC-seq ------------------------------------------------------------

## A single-cell atlas of chromatin accessibility in the human genome (2021)
## https://www.cell.com/cell/fulltext/S0092-8674(21)01279-4
## Each tissue sample was processed in its own 96-well plate (2000 nuclei loaded per well), in batches of 2-4 samples/plates

atlas <- readxl::read_excel("/Volumes/DannySSD/MULTI_ATAC/OtherDatasets/Atlas/mmc2.xlsx", skip = 0) %>% as.data.frame()
atlas <- atlas[1:92,]
ggplot(atlas, aes(x = `Nuclei passing QC`, y = `Median fragments per nucleus for nuclei passing QC`)) + geom_point() + geom_smooth(method = "lm", se=F) + theme_classic()
ggplot(atlas, aes(x = `Nuclei passing QC`, y = `Percent doublets (%) per dataset`)) + geom_point() + scale_x_log10()
ggplot(atlas, aes(x = `Nuclei passing QC`, y = `Median % reads duplicated`)) + geom_point() + scale_x_log10()

sci_atlas <- read.table("/Volumes/DannySSD/MULTI_ATAC/OtherDatasets/Atlas/GSE184462_metadata.tsv", sep = "\t", header = T)
sci_atlas$barcode <- gsub(".*\\+", "", sci_atlas$cellID)

sci_atlas %>% group_by(sample) %>% summarize(barLength = mean(nchar(barcode))) %>% data.frame()

sci_atlas$barcodeLength <- nchar(sci_atlas$barcode)

plot_grid(plotlist = lapply(unique(sci_atlas$barcodeLength), function(L) { 
  tmp <- sci_atlas[sci_atlas$barcodeLength == L,]
  tmp <- tmp[sample(nrow(tmp),size = min(nrow(tmp),10000)),]
  lapply(2:L, function(x) { data.frame(End = x, Unique = length(unique(str_sub(tmp$barcode, 1, x)))) } ) %>% do.call(rbind,.) %>% 
    ggplot(., aes(x = End, y = Unique)) + geom_line() + geom_point() + geom_vline(xintercept = c(10, 20, 30, 40)) + 
    # scale_y_log10() + 
    ggtitle(L)
}))

# these data are from the Domcke sci-ATAC-seq3 fetal chromatin accessibility paper
table(str_sub(sci_atlas[sci_atlas$barcodeLength == 40,]$barcode, 1, 10)) %>% length() # Tn5 i7 Oligo, 10 nt, 384 unique
table(str_sub(sci_atlas[sci_atlas$barcodeLength == 40,]$barcode, 11, 20)) %>% length() # PCR i7 Oligo, 10 nt, 96 unique
table(str_sub(sci_atlas[sci_atlas$barcodeLength == 40,]$barcode, 21, 30)) %>% length() # Tn5 i5 Oligo, 10 nt, 96 unique
table(str_sub(sci_atlas[sci_atlas$barcodeLength == 40,]$barcode, 31, 40)) %>% length() # PCR i5 Oligo, 8 nt, 384 unique

# other barcode lenghts are from other papers

# these are the data generated in this paper
table(str_sub(sci_atlas[sci_atlas$barcodeLength == 22,]$barcode, 1, 10)) %>% length() # Tn5 i7 Oligo, 10 nt, 384 unique
table(str_sub(sci_atlas[sci_atlas$barcodeLength == 22,]$barcode, 11, 22)) %>% length() # PCR i5 Oligo, 12 nt, 768 unique

sci_atlas <- sci_atlas[sci_atlas$barcodeLength == 22,]
sci_atlas$Tn5 <- str_sub(sci_atlas$barcode, 1, 10)

# sci_atlas <- sci_atlas[sci_atlas$barcodeLength == 40,]
# sci_atlas$Tn5 <- str_sub(sci_atlas$barcode, 31, 40)

sci_atlas %>% group_by(sample, Tn5) %>% summarize(nFrags = mean(logUMI), Count = n()) %>% 
  ggplot(aes(x = Count, y = nFrags, color = sample)) +  
  geom_point() + geom_smooth(method = 'lm', se=F) + #facet_wrap(~sample) +
  theme_classic() + theme(legend.position = 'none')
sci_atlas %>% group_by(sample, Tn5) %>% summarize(nFrags = mean(logUMI), Count = n()) %>% 
  ggplot(aes(x = Count, y = nFrags, color = sample)) +  
  geom_point() + geom_smooth(method = 'lm', se=F) + facet_wrap(~sample) +
  theme_minimal() + theme(legend.position = 'none')

sci_atlas %>% group_by(sample) %>% summarize(Tn5 = length(unique(Tn5)))

table(sci_atlas[,c("tissue", "Tn5")]) %>% heatmap() 

sci_atlas %>% group_by(tissue, Tn5) %>% summarize(Count = n()) %>% summarize(var = var(Count)) %>% arrange(-var)
sci_atlas %>% group_by(tissue, Tn5) %>% summarize(Count = n()) %>% summarize(sd = sd(Count)) %>% arrange(-sd)



### dsci-ATAC-seq -----------------------------------------------------------

## Droplet-based combinatorial indexing for massive-scale single-cell chromatin accessibility (2019)
## https://www.nature.com/articles/s41587-019-0147-6
## BMMC dataset: stimulated cells from 2 donors
## 96 barcoded Tn5 complexes; 8k cells, 20µL reactions (4µL indexed Tn5)
## Pooled tagmented nuclei split across 16 ddSeq lanes for library prep with varying input (20k, 40k, or 80k cells)

dsci <- rbind(data.frame(read.table("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/dsciATAC-seq/Baseline-Cells.tsv", sep = "\t"), Condition = "Baseline"),
              data.frame(read.table("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/dsciATAC-seq/Stimulation-Cells.tsv", sep = "\t"), Condition = "Stimulated"))
dsci$Tn5 <- gsub(".*Tn5-", "", dsci$V1)
dsci$Tn5 <- gsub("_.*", "", dsci$Tn5)
rownames(dsci) <- dsci$V1

tmp <- readxl::read_excel("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/dsciATAC-seq/41587_2019_147_MOESM8_ESM.xlsx", sheet = 3, skip = 11)
dsci$Tn5 <- factor(dsci$Tn5, levels = tmp$Barcode, labels = tmp$Oligo)

dsci %>% separate(V1,into=as.character(1:7),sep="_|-| |\\.") %>% {apply(.[,1:7], 2, table)}
dsci %>% separate(V1,into=as.character(1:7),sep="_|-| |\\.") %>% select(Tn5, Cell = 6, Rep = 7,Condition) %>% 
  group_by(Tn5, Condition) %>% summarize(Count = n()) %>% ggplot(aes(x = as.numeric(Tn5), y = Condition, fill = Count)) + geom_tile()
# group_by(Tn5, Rep) %>% summarize(Count = n()) %>% ggplot(aes(x = Tn5, y = Rep, fill = Count)) + geom_tile()
dsci %>% group_by(Tn5, Condition) %>% summarize(Count = n()) %>% 
  mutate(Row = ceiling((as.numeric(Tn5)) / 12), Col = (as.numeric(Tn5)-1) %% 12) %>%
  # mutate(Row = (as.numeric(Tn5)-1) %% 8, Col = ceiling((as.numeric(Tn5)) / 8)) %>%
  ggplot(aes(x = Col, y = Row, fill = log10(Count), alpha = Condition)) + geom_tile() + 
  scale_alpha_manual(values = c(0.6,1))# + scale_fill_viridis_c(option = "A")

dsci$Donor <- factor(as.numeric(dsci$Tn5), levels = 1:96, labels = rep(c("A","B"), each = 48))

# ArchR Analysis
addArchRGenome("hg19")

wd <- "/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/dsciATAC-seq"
setwd(wd)

# Create Arrow Files
inputFiles <- list.files(path = "GSE123581_RAW", pattern = ".fragments.tsv.gz$", full.names = T)
names(inputFiles) <- paste("Channel",1:16,sep="")

ArrowFiles <- vector()
for (i in 1:length(inputFiles)) {
  arrow <- createArrowFiles(
    inputFiles = inputFiles[i],
    sampleNames = names(inputFiles)[i],
    minTSS = 0,
    minFrags = 0,
    addTileMat = TRUE,
    addGeneScoreMat = F)
  ArrowFiles[i] <- arrow
  rm(arrow)
}

# Create ArchR Project
proj_1 <- ArchRProject(
  ArrowFiles = list.files(pattern = ".arrow"), 
  outputDirectory = "proj_1",
  copyArrows = T
)

proj_1$Lane <- gsub("#.*","",proj_1$cellNames) %>% gsub("Channel","",.)
proj_1$Cell <- gsub(".*#","",proj_1$cellNames)
proj_1$Tn5 <- gsub(".*Tn5-", "", proj_1$Cell) %>% gsub("_.*", "", .)
table(proj_1$Tn5)

table(proj_1$Cell %in% dsci$V1)
table(dsci$V1 %in% proj_1$Cell)

as.data.frame(proj_1@cellColData[proj_1$Cell %in% dsci$V1,]) %>% group_by(Tn5) %>% dplyr::summarise(nFrags = mean(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags)) + geom_point() + scale_x_log10() + scale_y_log10()

thresh.TSS <- 4
thresh.nFrag <- 100

ggplot(as.data.frame(proj_1@cellColData[proj_1$nFrags > thresh.nFrag & proj_1$TSSEnrichment > thresh.TSS,]), aes(x = log10(nFrags), y = TSSEnrichment)) + 
  # geom_density_2d_filled(bins = 30) + 
  geom_point(size = 0.1, alpha = 0.1) + 
  geom_vline(xintercept = log10(thresh.nFrag)) + 
  geom_hline(yintercept = thresh.TSS) + 
  facet_wrap(~Sample) + 
  theme_minimal() + 
  theme(legend.position = "none")

proj_2 <- subsetArchRProject(proj_1,
                             outputDirectory = "./proj_2",
                             cells = proj_1$cellNames[proj_1$Cell %in% dsci$V1 & proj_1$nFrags > 100], 
                             dropCells = F,
                             force = T)

# proj_2 <- addIterativeLSI(proj_2, force = T)
# proj_2 <- addClusters(proj_2, dimsToUse = 2:30, force = T)
# proj_2 <- addUMAP(proj_2, dimsToUse = 2:30, force = T)

proj_2 <- dimReduce(proj_2, what = c("iLSI","UMAP"), scale = c("Row","Col"), default = "Col", dims = 2:30)
proj_2 <- dimDefault(proj_2, what = c("iLSI","UMAP"), default = "NoScale")

proj_2$CellType <- dsci$V2[match(proj_2$Cell, dsci$V1)]
proj_2$Condition <- dsci$Condition[match(proj_2$Cell, dsci$V1)]
proj_2$Tn5 <- dsci$Tn5[match(proj_2$Cell, dsci$V1)]
proj_2$Donor <- dsci$Donor[match(proj_2$Cell, dsci$V1)]
# table(proj_2$Lane) %>% plot()
proj_2$Load <- factor(as.numeric(proj_2$Lane), levels = 1:16, labels = c(rep("20k",4),rep("40k",8),rep("80k",4)))


plot_grid(plotEmbedding(proj_2, name = "Load", labelAsFactors = F, labelMeans = F) + theme_void() + theme(legend.position = 'none') + scale_color_viridis_d(),
          plotEmbedding(proj_2, name = "Donor", labelAsFactors = F, labelMeans = F) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_2, name = "Condition", labelAsFactors = F, labelMeans = F) + theme_void(),
          plotEmbedding(proj_2, name = "CellType", labelAsFactors = F) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_2, name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_2, name = "TSSEnrichment", plotAs = 'points', size = 1) + theme_void(),
          nrow = 2, byrow = F)

plot_grid(plotEmbedding(proj_2, embedding = "UMAP_NoScale", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_2, embedding = "UMAP_NoScale", name = "Donor", labelMeans = T, labelAsFactors = F, size = 1) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_2, embedding = "UMAP_NoScale", name = "Condition", labelMeans = F, labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_2, embedding = "UMAP_NoScale", name = "CellType", labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_2, embedding = "UMAP_Row", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_2, embedding = "UMAP_Row", name = "Donor", labelMeans = T, labelAsFactors = F, size = 1) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_2, embedding = "UMAP_Row", name = "Condition", labelMeans = F, labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_2, embedding = "UMAP_Row", name = "CellType", labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_2, embedding = "UMAP_Col", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_2, embedding = "UMAP_Col", name = "Donor", labelMeans = T, labelAsFactors = F, size = 1) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_2, embedding = "UMAP_Col", name = "Condition", labelMeans = F, labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_2, embedding = "UMAP_Col", name = "CellType", labelAsFactors = F, size = 1) + theme_void(),
          nrow = 3)



as.data.frame(proj_2@cellColData) %>% group_by(Tn5) %>% dplyr::summarise(nFrags = mean(nFrags), Count = n()) %>% ggplot(., aes(x = Count, y = nFrags)) + geom_point() + theme_bw()
as.data.frame(proj_2@cellColData) %>% group_by(Tn5, Lane = as.numeric(Lane),Load) %>% dplyr::summarise(nFrags = mean(nFrags), Count = n()) %>% ggplot(., aes(x = Count, y = nFrags, color = Load)) + geom_point() + facet_wrap(~Lane) + theme_bw()
as.data.frame(proj_2@cellColData) %>% group_by(Tn5, Load) %>% dplyr::summarise(nFrags = mean(nFrags), Count = n()) %>% ggplot(., aes(x = Count, y = nFrags, color = Load)) + geom_point() + facet_wrap(~Load) + theme_bw()
as.data.frame(proj_2@cellColData) %>% group_by(Tn5, Load) %>% dplyr::summarise(nFrags = mean(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags, color = Load)) + geom_point() + geom_smooth(method = "lm", se = F) + theme_classic()
as.data.frame(proj_2@cellColData) %>% group_by(Tn5, Load, Condition, Donor) %>% dplyr::summarise(nFrags = mean(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags, color = Load, group = paste(Load,Condition,Donor))) + geom_point() + geom_smooth(method = "lm", se = F) + theme_bw()

as.data.frame(proj_2@cellColData) %>% group_by(Lane = as.numeric(Lane),Load) %>% dplyr::summarise(nFrags = mean(nFrags), Count = n()) %>% ggplot(., aes(x = Count, y = nFrags, color = Load)) + geom_point() + theme_bw()
as.data.frame(proj_2@cellColData) %>% group_by(Lane = as.numeric(Lane),Load) %>% dplyr::summarise(Collisions = sum(CellType == "Collision"), Count = n()) %>% ggplot(., aes(x = Lane, y = Collisions/Count, color = Load)) + geom_point() + theme_bw()



as.data.frame(proj_2@cellColData) %>% group_by(Tn5, Condition) %>% summarize(nFrags = median(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags, color = Condition)) + geom_point() + geom_smooth(method = 'lm', se = F) +
  theme_classic() + ggtitle("dsci-ATAC-seq_BMMC - Barcoded")

as.data.frame(proj_2@cellColData) %>% group_by(Tn5, Condition, Donor) %>% summarize(nFrags = mean(nFrags), Count = n()) %>% 
  ggplot(., aes(x = Count, y = nFrags, color = Condition)) + geom_point() + geom_smooth(method = 'lm', se = F) +
  theme_minimal() + facet_wrap(~Donor)



### txci-ATAC-seq -----------------------------------------------------------

## txci-ATAC-seq: a massive-scale single-cell technique to profile chromatin accessibility
## https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03150-1
## 96-well plate pre-loaded with 5µL 500nM indexed Tn5
## To each well, added 5k nuclei in 20µL 1.25X Tagment DNA buffer
## All wells pools and concentrated before loading 75k into 10x lane
## To each well, added 20k nuclei (7µL nuclei + 13µL TBS)

# addArchRGenome("hg38")
addArchRGenome("mm10")

wd <- "/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/txci-ATAC-seq"
setwd(wd)

# Create Arrow Files
inputFiles <- list.files(path = getwd(), pattern = "mm10.*.fragments.tsv.gz$", full.names = T)
names(inputFiles) <- c("mm10_liver", "mm10_lung")

ArrowFiles <- vector()

for (i in 1:length(inputFiles)) {
  arrow <- createArrowFiles(
    inputFiles = inputFiles[i],
    sampleNames = names(inputFiles)[i],
    minTSS = 0,
    minFrags = 0,
    addTileMat = TRUE,
    addGeneScoreMat = F)
  ArrowFiles[i] <- arrow
  rm(arrow)
}

ArrowFiles <- list.files(pattern = "mm.*arrow")

# Create ArchR Project
proj_tx <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "proj_tx",
  copyArrows = F
)

proj_tx <- addIterativeLSI(proj_tx, force = T)
proj_tx <- addClusters(proj_tx, force = T)
proj_tx <- addUMAP(proj_tx, force = T)

proj_tx$Cell <- gsub(".*#","",proj_tx$cellNames)

nt <- consensusMatrix(DNAStringSet(proj_tx$Cell), as.prob = T)[1:4,] %>% t() %>% as.data.frame()
nt$Position <- 1:nrow(nt)
nt <- tidyr::pivot_longer(nt, cols = 1:4, names_to = "Base", values_to = "Frequency")
ggplot(nt, aes(x = Position, y = Frequency, color = Base)) + geom_line(linewidth = 1) + theme_minimal()

proj_tx$Tn5 <- stringr::str_sub(proj_tx$cellNames, -8, -1)

plot_grid(plotEmbedding(proj_tx) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_tx, name = "Tn5") + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_tx, name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_tx, name = "TSSEnrichment", plotAs = 'points', size = 1) + theme_void(),
          nrow = 2)

## mapping from wells of mixed mouse lung/liver + human lung nuclei
# mapping <- read.table("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/txci-ATAC-seq/GSE231708_standard_txci.native_tissue_mixed_species_mapping_report_100k.txt.gz", header = T)
# mapping <- read.table("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/txci-ATAC-seq/GSE231708_standard_txci.native_tissue_mixed_species_mapping_report_200k.txt.gz", header = T)
# mapping$Tn5 <- gsub(".*_N...", "", mapping$barcode)
# 
# mapping %>% filter(total_unique_reads > 300) %>% ggplot(aes(x = hg38_unique_reads, y = mm10_unique_reads)) + geom_point() + facet_wrap(~sample) + scale_x_log10() + scale_y_log10()
# mapping %>% filter(total_unique_reads > 300) %>% ggplot(aes(x = hg38_unique_reads/total_unique_reads)) + geom_density() + facet_wrap(~sample) #+ scale_x_log10() + scale_y_log10()
# mapping %>% filter(total_unique_reads > 300) %>% group_by(Tn5) %>% mutate(Count = n(), hg38 = mean(hg38_unique_reads/total_unique_reads)) %>% ggplot(aes(x = Count, y = hg38)) + geom_point() + geom_smooth(method = 'lm')
# mapping %>% filter(total_unique_reads > 300) %>% group_by(Tn5, Geno = hg38_unique_reads/total_unique_reads > .5) %>% mutate(Count = n(), hg38 = mean(hg38_unique_reads/total_unique_reads)) %>% 
#   # filter(Count > 300) %>% 
#   ggplot(aes(x = Count, y = hg38, color = Geno)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~Geno, scales = "free")

tx_lungH <- read.table("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/txci-ATAC-seq/GSE231708_standard_txci.native_hg38_lung_metadata.txt.gz", header = T)
tx_lungM <- read.table("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/txci-ATAC-seq/GSE231708_standard_txci.native_mm10_lung_metadata.txt.gz", header = T)
tx_liverM <- read.table("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/txci-ATAC-seq/GSE231708_standard_txci.native_mm10_liver_metadata.txt.gz", header = T)
tx <- rbind(tx_lungH %>% select(-peaks_snn_res.0.3) %>% mutate(Replicate = Tissue),
            tx_lungM %>% select(-peaks_snn_res.0.8),
            tx_liverM %>% select(-peaks_snn_res.0.2))

ref <- readxl::read_excel("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/txci-ATAC-seq/13059_2023_3150_MOESM5_ESM.xlsx", sheet = 2)

lapply(2:32, function(x) { data.frame(End = x, Unique = length(unique(str_sub(tx$Barcodes, 1, x)))) } ) %>% do.call(rbind,.) %>% 
  ggplot(., aes(x = End, y = Unique)) + geom_line() + geom_point() + geom_vline(xintercept = c(8, 24)) + scale_y_log10()
# this plot reveals something about 10x barcode synthesis, seems broken up into 3 parts, two 7-base oligos merged with a set of linkers

table(stringr::str_sub(tx$Barcodes,1,8)) %>% length() # i7 representing library (100k vs 200k input)
table(stringr::str_sub(tx$Barcodes,9,24)) %>% length() # 10X cell barcodes
table(stringr::str_sub(tx$Barcodes,25,32)) %>% length() # Tn5 barcode

tx$Tn5 <- stringr::str_sub(tx$Barcodes,25,32)

table(tx$Tn5 %in% ref$Tn5_bc)

tx$Tn5ID <- factor(tx$Tn5, levels = ref$Tn5_bc, labels = ref$Well_id)
tx$Row <- factor(tx$Tn5, levels = ref$Tn5_bc, labels = ref$Row_id)
tx$Col <- factor(tx$Tn5, levels = ref$Tn5_bc, labels = ref$Column_id)

tx %>% group_by(Tn5) %>% mutate(Count = n()) %>% 
  group_by(Tn5, Count, Barnyard, Tissue, Source = paste(Tissue,Replicate), Col, Row) %>% summarize(nFrags = mean(Complexity)) %>%
  filter(Tissue != "Human_lung") %>%
  ggplot(aes(x = Col, y = Row, fill = Source, color = Barnyard)) + 
      geom_dotplot(binwidth = 0.7, alpha = 0.8, stroke = 3) + scale_y_discrete(limits = rev) + scale_color_manual(values = c(NA, "black"), na.value = NA)

tx %>% group_by(Tn5) %>% mutate(Count = n()) %>% 
  group_by(Tn5, Count, Tissue = paste(Tissue,Replicate), NucInput) %>% summarize(nFrags = mean(Complexity)) %>%
  ggplot(., aes(x = Count, y = log10(nFrags), color = Tissue, shape = NucInput, group = paste(NucInput,Tissue))) +
  geom_point() + theme_classic() + geom_smooth(method = 'lm', se = F)

tx %>% group_by(Tn5) %>% mutate(Count = n()) %>% 
  group_by(Tn5, Count, Tissue = paste(Tissue,Replicate)) %>% summarize(nFrags = mean(Complexity)) %>%
  ggplot(., aes(x = Count, y = log10(nFrags), color = Tissue, group = Tissue)) +
  geom_point() + theme_classic() + geom_smooth(method = 'lm', se = F)

tx2 <- tx %>% filter(Barnyard != "True") %>% 
  mutate(Cell = case_when(
    Tissue == "liver" ~ paste("mm10_liver",Barcodes,sep="#"),
    Tissue == "lung" ~ paste("mm10_lung",Barcodes,sep="#"))
  )
rownames(tx2) <- tx2$Cell
  

proj_tx2 <- subsetArchRProject(proj_tx,
                               outputDirectory = "./proj_tx2",
                               cells = tx2$Cell,
                               dropCells = F,
                               force = T)

proj_tx2 <- addDoubletScores(
  input = proj_tx2,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

proj_tx3 <- filterDoublets(proj_tx2)

# proj_tx3 <- addIterativeLSI(proj_tx3, force = T)
# proj_tx3 <- addClusters(proj_tx3, force = T)
# proj_tx3 <- addUMAP(proj_tx3, force = T)

proj_tx3 <- dimReduce(proj_tx3, what = c("iLSI","UMAP"), scale = c("Row","Col"), default = "NoScale", dims = 2:30)
# proj_tx3 <- dimDefault(proj_tx3, what = c("iLSI","UMAP"), default = "NoScale")

proj_tx3@embeddings$UMAP_Split <- rbind(addUMAP(proj_tx3[proj_tx3$Tissue == "liver"], reducedDims = "IterativeLSI_NoScale", name = "UMAP_Split", dimsToUse = 2:30, corCutOff = 1, scale = F, force = TRUE)@embeddings$UMAP_Split,
                                        addUMAP(proj_tx3[proj_tx3$Tissue == "lung"], reducedDims = "IterativeLSI_NoScale", name = "UMAP_Split", dimsToUse = 2:30, corCutOff = 1, scale = F, force = TRUE)@embeddings$UMAP_Split)

proj_tx3@embeddings$UMAP_Split <- proj_tx3@embeddings$UMAP
proj_tx3@embeddings$UMAP_Split$df <- rbind(dimReduce(proj_tx3[proj_tx3$Tissue == "liver"], what = c("iLSI","UMAP"), scale = NULL, default = "NoScale", dims = 2:30)@embeddings$UMAP$df,
                                           dimReduce(proj_tx3[proj_tx3$Tissue == "lung"], what = c("iLSI","UMAP"), scale = NULL, default = "NoScale", dims = 2:30)@embeddings$UMAP$df)[proj_tx3$cellNames,]

plotEmbedding(proj_tx3[proj_tx3$Tissue == "liver"], embedding = "UMAP_Split", name = "CellType") + theme_DC
plotEmbedding(proj_tx3[proj_tx3$Tissue == "lung"], embedding = "UMAP_Split", name = "CellType") + theme_DC


proj_tx3@cellColData <- left_join(proj_tx3@cellColData %>% as.data.frame %>% rownames_to_column(), tx2 %>% mutate(Cell = Barcodes, Sample = paste("mm10",Tissue,sep="_")), by = c("Tn5", "Cell", "Sample")) %>% column_to_rownames() %>% DataFrame

plot_grid(plotEmbedding(proj_tx3, name = "NucInput") + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_tx3, name = "Tissue") + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_tx3, name = "Replicate") + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_tx3, name = "CellType", labelAsFactors = F) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_tx3, name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_tx3, name = "TSSEnrichment", plotAs = 'points', size = 1) + theme_void(),
          nrow = 2)

plot_grid(plotEmbedding(proj_tx3, embedding = "UMAP_NoScale", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_tx3, embedding = "UMAP_NoScale", name = "Replicate", labelMeans = T, labelAsFactors = F, size = 1) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_tx3, embedding = "UMAP_NoScale", name = "Tissue", labelMeans = F, labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_tx3, embedding = "UMAP_NoScale", name = "CellType", labelMeans = F, labelAsFactors = F, size = 1) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_tx3, embedding = "UMAP_Row", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_tx3, embedding = "UMAP_Row", name = "Replicate", labelMeans = T, labelAsFactors = F, size = 1) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_tx3, embedding = "UMAP_Row", name = "Tissue", labelMeans = F, labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_tx3, embedding = "UMAP_Row", name = "CellType", labelMeans = F, labelAsFactors = F, size = 1) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_tx3, embedding = "UMAP_Col", name = "log10(nFrags)", plotAs = 'points', size = 1) + theme_void(),
          plotEmbedding(proj_tx3, embedding = "UMAP_Col", name = "Replicate", labelMeans = T, labelAsFactors = F, size = 1) + theme_void() + theme(legend.position = 'none'),
          plotEmbedding(proj_tx3, embedding = "UMAP_Col", name = "Tissue", labelMeans = F, labelAsFactors = F, size = 1) + theme_void(),
          plotEmbedding(proj_tx3, embedding = "UMAP_Col", name = "CellType", labelMeans = F, labelAsFactors = F, size = 1) + theme_void() + theme(legend.position = 'none'),
          nrow = 3)

as.data.frame(proj_tx3@cellColData) %>%
  group_by(Tn5, Replicate, NucInput, Tissue) %>%
  summarize(nFrags = mean(nFrags), Count = n()) %>%
  ggplot(., aes(x = Count, y = nFrags, color = Tissue, shape = NucInput, group = paste(NucInput,Tissue,Replicate))) +
  geom_point() + theme_classic() + geom_smooth(method = 'lm', se = F)

as.data.frame(proj_tx3@cellColData) %>% group_by(Tn5) %>% mutate(Count = n()) %>% group_by(Tissue, Tn5, CellType, Count) %>% summarize(nFrags = mean(nFrags)) %>%
  ggplot(., aes(x = Count, y = nFrags, color = CellType)) + geom_point() + geom_smooth(method = "lm", se = F) +
  theme_minimal() + facet_wrap(~Tissue)

as.data.frame(proj_tx3@cellColData) %>% group_by(Tissue,Replicate, NucInput, Barnyard, Tn5) %>% summarize(nFrags = mean(nFrags), Count = n()) %>% group_by(Tissue, Replicate, NucInput) %>%
  summarize(Cor = cor(nFrags,Count,method = "spearman"),
            Slope = (lm(nFrags~Count) %>% coef())[2])



### scifi-ATAC-seq ----------------------------------------------------------

## scifi‐ATAC‐seq: massive‐scale single‐cell chromatin accessibility sequencing using combinatorial fluidic indexing (2024)
## Basically txci-ATAC-seq, but with double-indexed Tn5
## Experiments done with maize nuclei from different genotypes
## Pre-indexed with 96 different Tn5 complexes (combo of 8 A and 12 B barcodes) and then superloaded 300k nuclei into 10x lane (100k-200k used in pilot experiments with just 2 genotypes)

scifi <- readxl::read_excel("/Volumes/DannySSD/MULTI_ATAC/CompetingMethods/scifi-ATAC_seq/13059_2024_3235_MOESM2_ESM.xlsx", sheet=4, skip = 1)
ggplot(scifi, aes(x=umap1, y=umap2, color = Cell_type)) + geom_point(size = 0.1)
ggplot(scifi, aes(x=umap1, y=umap2, color = genotype)) + geom_point(size = 0.1)

scifi$Tn5 <- paste(scifi$tn5A_bc, scifi$tn5B_bc)
scifi$Cell <- gsub(".*:", "", scifi$cellID) %>% gsub("-.*","",.)

lapply(2:26, function(x) { data.frame(End = x, Unique = length(unique(str_sub(scifi$Cell, 1, x)))) } ) %>% do.call(rbind,.) %>% 
  ggplot(., aes(x = End, y = Unique)) + geom_line() + geom_point() + geom_vline(xintercept = c(16, 21)) + scale_y_log10()

scifi %>% group_by(genotype) %>% summarize(nTn5 = length(unique(Tn5)))
scifi %>% group_by(library) %>% summarize(nTn5 = length(unique(Tn5)))

scifi %>% group_by(Tn5, genotype) %>% summarize(Count = n(), nFrags = median(total)) %>%
  ggplot(aes(x = Count, y = nFrags, color = genotype)) + geom_point() + geom_smooth(method = 'lm', se=F) +
  theme_classic() + ggtitle("scifi-ATAC-seq - 7 maize genotypes")



### EasySci-ATAC ------------------------------------------------------------

## A global view of aging and Alzheimer’s pathogenesis-associated cell population dynamics and molecular signatures in human and mouse brains
## https://www.nature.com/articles/s41588-023-01572-y
## Similar to sci-ATAC-seq

easy <- read.csv("/Volumes/DannySSD/MULTI_ATAC/OtherDatasets/EasySci-ATAC/GSE212606_RAW/GSM6538357_ATAC_cell_annotation.csv.gz")

table(gsub(".*\\.","",easy$sample) %>% str_sub(., 1, 10)) %>% length()
table(gsub(".*\\.","",easy$sample) %>% str_sub(., 11, 20)) %>% length()

easy$BC1 <- gsub(".*\\.","",easy$sample) %>% str_sub(., 1, 10)
easy$BC2 <- gsub(".*\\.","",easy$sample) %>% str_sub(., 11, 20)
easy$Type <- gsub("_.*","",easy$Conditions)
easy$Replicate <- gsub("^.*?_","",easy$Conditions)

easy %>% group_by(BC2, Type) %>% summarize(Count = n(), nFrags = median(reads_in_peaks)) %>%
  ggplot(aes(x = Count, y = nFrags, color = Type)) +
  geom_point() + geom_smooth(method = 'lm', se= F)

easy %>% group_by(BC2, Type) %>% summarize(Count = n(), nFrags = median(total_reads)) %>%
  ggplot(aes(x = Count, y = nFrags, color = Type)) +
  geom_point() + geom_smooth(method = 'lm', se= F)

easy %>% group_by(BC2, Conditions) %>% summarize(Count = n(), nFrags = median(total_reads)) %>%
  ggplot(aes(x = Count, y = nFrags, color = Conditions)) +
  geom_point() + geom_smooth(method = 'lm', se= F)

easy %>% group_by(BC2, Type) %>% summarize(Count = n(), nFrags = mean(total_reads)) %>% lm(log10(nFrags)~log10(Count), .) %>% coef

easy %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = Main_cluster_name)) + geom_point() + theme(legend.position = 'none')
easy %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = Type)) + geom_point()

easy %>% 
  group_by(BC2, Type, Replicate) %>% mutate(Tn5Count = n()) %>%           # count nCells per Tn5 reaction per batch/dataset
  group_by(Type, Replicate) %>% mutate(Tn5Rank = rank(Tn5Count),          # rank Tn5 reaction counts per batch/dataset
                                      Bin = factor(Tn5Rank, levels = sort(unique(Tn5Rank)), labels = ntile(sort(unique(Tn5Rank)), nbin)),
                                      Random = sample(Bin, size = n(), replace=F)) %>%
# binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(CellType = Main_cluster_name, Group = paste(Type,Replicate), Bin) %>% summarize(n = n()) %>% group_by(Group, Bin) %>% mutate(freq = n / sum(n)) %>%
  # group_by(CellType, Group) %>% mutate(Diff = freq[Bin == 3] - freq[Bin == 1]) %>% ungroup %>%
  # group_by(Group) %>% filter(Diff %in% c(min(Diff),max(Diff))) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = CellType, y = freq, fill = Bin)) + geom_col(position = "dodge") + facet_wrap(~Group, scales = "free") +
  labs(subtitle = paste("EASY -", group), y = "Cell Type Frequency", x = "Cell Type", fill = NULL) + scale_fill_manual(values = colors_bins) + theme_DC +
  theme(axis.text.x = element_text(angle = -25, hjust = 0))



### 10x scATAC-seq ----------------------------------------------------------

## Single-cell epigenomics reveals mechanisms of human cortical development
## https://www.nature.com/articles/s41586-021-03209-8
## Lots of individual 10x lanes

tmp <- readxl::read_excel("/Volumes/DannySSD/MULTI_ATAC/OtherDatasets/CorticalDev/41586_2021_3209_MOESM3_ESM.xlsx", skip = 0)
tmp <- c(colnames(tmp)[1:10],tmp[1,11:19]) %>% unname %>% unlist
tenx <- readxl::read_excel("/Volumes/DannySSD/MULTI_ATAC/OtherDatasets/CorticalDev/41586_2021_3209_MOESM3_ESM.xlsx", skip = 1)
colnames(tenx) <- tmp

tenx %>% ggplot(aes(x = `Cell Ranger Cell Count`, y = `Median Fragments per cell`, color = Region)) + 
  geom_point() + geom_smooth(method = 'lm', se=F) #+ scale_y_log10()
tenx %>% ggplot(aes(x = `Cell Ranger Cell Count`, y = `Median Fragments per cell`, color = `Fresh or Frozen?`)) + 
  geom_point() + geom_smooth(method = 'lm', se=F) #+ scale_y_log10()
tenx %>% ggplot(aes(x = `Cell Ranger Cell Count`, y = `Median Fragments per cell`, color = `Specimen of Origin`)) + 
  geom_point() + geom_smooth(method = 'lm', se=F) #+ scale_y_log10()

tenx %>% filter(Region != "GE") %>% 
  ggplot(aes(x = `Cell Ranger Cell Count`, y = `Median Fragments per cell`, color = paste(Region, `Fresh or Frozen?`))) + 
  geom_point() + geom_smooth(method = 'lm', se=F) + scale_y_log10() + theme_classic() + labs(color = "SampleType")

tenx %>% filter(Region != "GE") %>% 
  ggplot(aes(x = `Cell Ranger Cell Count`, y = `Median Fragments per cell`, color = Region)) + 
  geom_point() + geom_smooth(method = 'lm', se=F) + scale_y_log10() + theme_classic() + labs(color = "SampleType")
tenx %>% filter(Region != "GE") %>% 
  ggplot(aes(x = `Cell Ranger Cell Count`, y = `Median Fragments per cell`, color = `Fresh or Frozen?`)) + 
  geom_point() + geom_smooth(method = 'lm', se=F) + scale_y_log10() + theme_classic() + labs(color = "SampleType")
tenx %>% filter(Region != "GE") %>% 
  ggplot(aes(x = `Cell Ranger Cell Count`, y = `Median Fragments per cell`, color = paste(`Fresh or Frozen?`,Region))) + 
  geom_point() + geom_smooth(method = 'lm', se=F) + scale_y_log10() + theme_classic() + labs(color = "SampleType")
tenx %>%# filter(Region != "GE") %>% 
  ggplot(aes(x = `Cell Ranger Cell Count`, y = `Median Fragments per cell`, color = `Fresh or Frozen?`)) + 
  geom_point() + geom_smooth(method = 'lm', se=F) + scale_y_log10() + theme_classic() + labs(color = "SampleType")

tenx %>% lm(log10(`Median Fragments per cell`)~log10(`Cell Ranger Cell Count`), .) %>% coef


tenx %>% select(`sample ID`, `Specimen of Origin`, Region, `Fresh or Frozen?`, Assay, `Number of Cells`) %>% print(n = 50)
tenx %>% select(`sample ID`, `Specimen of Origin`, Region, `Fresh or Frozen?`, Assay, `Number of Cells`) %>% mutate(Date = gsub(".*_", "", `Specimen of Origin`)) %>% .$Date %>% table

max(tenx$`Cell Ranger Cell Count`) / min(tenx$`Cell Ranger Cell Count`)
max(tenx$`Number of Cells`) / min(tenx$`Number of Cells`)



### Spear-ATAC --------------------------------------------------------------

## High-throughput single-cell chromatin accessibility CRISPR screens enable unbiased identification of regulatory networks in cancer
## https://www.nature.com/articles/s41467-021-23213-w
## CRISPR screens using 10x Genomics platform

tmp <- list.files("../OtherDatasets/Spear-ATAC", pattern = ".singlecell.csv", full.names = T)

spear <- lapply(tmp, function(x) {
  df <- data.frame(Sample = x %>% gsub(".*/","",.) %>% gsub(".singlecell.csv","",.),
                   read.csv(x)[-1,c("barcode","passed_filters","is__cell_barcode")])
  df <- df[df$is__cell_barcode == 1,]
}) %>% do.call(rbind,.)

spear$CellType <- gsub("-.*","",spear$Sample)
spear$Experiment <- gsub("^.*?-","",spear$Sample) %>% gsub("-.*","",.)

spear %>% group_by(Experiment) %>% summarize(Lanes = length(unique(Sample)))

spear %>% group_by(Sample) %>% summarize(Count = n(), nFrags = mean(passed_filters)) %>%
  ggplot(aes(x = Count, y = nFrags)) + geom_point() + geom_smooth(method = 'lm', se = F)
spear %>% group_by(Sample, CellType, Experiment) %>% summarize(Count = n(), nFrags = mean(passed_filters)) %>%
  ggplot(aes(x = Count, y = nFrags, shape = Experiment, color = CellType, group = CellType)) + geom_point() + geom_smooth(method = 'lm', se = F) + 
  theme_classic()
spear %>% ggplot(aes(x = log10(passed_filters), color = CellType)) + geom_density() + theme_classic()
spear %>% ggplot(aes(x = log10(passed_filters), color = Experiment)) + geom_density() + theme_classic() + facet_wrap(CellType~.)
spear %>% filter(Experiment == "TimeCourse") %>% ggplot(aes(x = Sample, y = log10(passed_filters))) + geom_violin()

spear %>% 
  # filter(Experiment %in% c("TimeCourse","LargeScreen")) %>%
  # filter(passed_filters > 1000) %>%
  group_by(Sample, CellType, Experiment) %>% summarize(Count = n(), nFrags = mean(passed_filters)) %>%
  ggplot(aes(x = Count, y = nFrags, shape = Experiment, color = CellType, group = CellType)) + geom_point() + geom_smooth(method = 'lm', se = F) + 
  theme_classic()

spear %>% 
  filter(Experiment %in% c("TimeCourse","LargeScreen")) %>%
  # filter(passed_filters > 1000) %>%
  group_by(CellType, Experiment) %>% summarize(Count = n(), nFrags = mean(passed_filters)) %>%
  ggplot(aes(x = Count, y = nFrags, shape = Experiment, color = CellType, group = CellType)) + geom_point() + geom_smooth(method = 'lm', se = F) + 
  theme_classic()

spear %>% 
  filter(Experiment %in% c("TimeCourse","LargeScreen")) %>%
  # filter(passed_filters > 1000) %>%
  group_by(Sample, Experiment, CellType) %>% summarize(Count = n(), nFrags = mean(passed_filters)) %>%
  ggplot(aes(x = Count, y = nFrags, color = CellType, group = CellType)) + geom_point() + geom_smooth(method = 'lm', se = F) + 
  theme_classic()


##########################################################################################################################################################

### Combined Meta-Analysis --------------------------------------------------

comb <- rbind(sci_atlas %>% 
                select(Tn5, CellType = cell.type, nFrags = logUMI, Group = tissue) %>% 
                mutate(nFrags = 10^nFrags, Tn5Type = "Indexed", UMAP1 = NA, UMAP2 = NA, Dataset = "sci-ATAC-seq_HumanAtlas"),
              easy %>% 
                select(Tn5 = BC2, CellType = Main_cluster_name, nFrags = total_reads, Group = Type, UMAP1 = UMAP_1, UMAP2 = UMAP_2) %>% 
                mutate(Tn5Type = "Indexed", Dataset = "EasySci-ATAC_MouseBrain"),
              proj_plex2@cellColData %>% as.data.frame %>%
                mutate(UMAP1 = proj_plex2@embeddings$UMAP$df[,1], UMAP2 = proj_plex2@embeddings$UMAP$df[,2]) %>%
                select(Tn5, CellType = cell_type, nFrags = nFrags, Group = cell_type, UMAP1, UMAP2) %>% 
                mutate(Tn5Type = "Indexed", Dataset = "sciPlex-ATAC-seq2_DrugScreen"),
              scifi %>%
                select(Tn5, CellType = Cell_type, nFrags = unique, Group = genotype, UMAP1 = umap1, UMAP2 = umap2) %>% 
                mutate(Tn5Type = "Indexed", Dataset = "scifi-ATAC_Maize"),
              proj_2@cellColData %>% as.data.frame %>% 
                mutate(Group = paste(Condition,Donor,Load), UMAP1 = proj_2@embeddings$UMAP$df[,1], UMAP2 = proj_2@embeddings$UMAP$df[,2]) %>% 
                filter(Load == "40k") %>%
                select(Tn5, CellType, nFrags, Group, UMAP1, UMAP2) %>% 
                mutate(Tn5Type = "Indexed", Dataset = "dsci-ATAC-seq_BMMC"),
              proj_tx3@cellColData %>% as.data.frame %>% 
                # mutate(Group = paste(Tissue,Replicate,NucInput), UMAP1 = proj_tx2@embeddings$UMAP$df[,1], UMAP2 = proj_tx2@embeddings$UMAP$df[,2]) %>%
                # mutate(Group = paste(Tissue,Replicate,NucInput), UMAP1 = proj_tx2@embeddings$UMAP_Split$df[,1], UMAP2 = proj_tx2@embeddings$UMAP_Split$df[,2]) %>%
                mutate(Group = paste(Tissue,Replicate,NucInput), UMAP1 = proj_tx3@embeddings$UMAP_Split$df[,1], UMAP2 = proj_tx3@embeddings$UMAP_Split$df[,2]) %>%
                filter(NucInput == "100k") %>%
                select(Tn5, CellType, nFrags, Group, UMAP1, UMAP2) %>% 
                mutate(Tn5Type = "Indexed", Dataset = "txci-ATAC-seq_Barnyard"),
              sci3_fet %>% 
                select(Tn5 = sample_name, CellType = cell_type, nFrags = total_deduplicated, Group = batch, UMAP1 = tissue_umap_1, UMAP2 = tissue_umap_2) %>% 
                mutate(Tn5Type = "Standard", Dataset = "sci-ATAC-seq3_FetalAtlas"),
              sci3_dros %>% 
                select(Tn5 = time, CellType, nFrags = total_reads_in_peaks, Group = Batch) %>% 
                mutate(Tn5Type = "Standard", UMAP1 = NA, UMAP2 = NA, Dataset = "sci-ATAC-seq3_DrosophilaDev"),
              spear %>% filter(Experiment %in% c("LargeScreen","TimeCourse")) %>% 
                select(Tn5 = Sample, CellType, nFrags = passed_filters) %>% 
                mutate(Tn5Type = "Standard", Group = CellType, UMAP1 = NA, UMAP2 = NA, Dataset = "Spear-ATAC_CellLines"),
              proj_snu2@cellColData %>% as.data.frame() %>% 
                mutate(UMAP1 = proj_snu2@embeddings$UMAP$df[,1], UMAP2 = proj_snu2@embeddings$UMAP$df[,2]) %>%
                filter(Tn5 %ni% "s94") %>% 
                select(Tn5, CellType, nFrags, UMAP1, UMAP2) %>% 
                mutate(Group = CellType, Tn5Type = "Standard", Dataset = "SNuBar_CellLines"),
              proj_breast2@cellColData %>% as.data.frame() %>% 
                mutate(UMAP1 = proj_breast2@embeddings$UMAP$df[,1], UMAP2 = proj_breast2@embeddings$UMAP$df[,2]) %>%
                select(Tn5, CellType, nFrags, Group = Region, UMAP1, UMAP2) %>% 
                mutate(Tn5Type = "Standard", Dataset = "SNuBar_HBCA"))

comb2 <- comb %>% 
  group_by(Dataset, Tn5Type, Group, Tn5) %>% mutate(Count = n()) %>% # count cells per Tn5 reaction
  group_by(Dataset, Tn5Type, Group, Count, Tn5) %>% summarize(nFrags = median(nFrags)) # get median fragment counts per Tn5 reaction by biological group

comb2 <- rbind(comb2,
               tenx %>% filter(Region != "GE") %>% mutate(Group = paste(Region, `Fresh or Frozen?`)) %>%
               # tenx %>% mutate(Group = `Fresh or Frozen?`) %>%
                 select(Tn5 = `sample ID`, Count = `Cell Ranger Cell Count`, nFrags = `Median Fragments per cell`, Group) %>% 
                 mutate(Tn5Type = "Standard", Dataset = "10xATACv1_CorticalDev"))


colors_datasets <- pals::cols25(length(unique(comb2$Dataset)))
names(colors_datasets) <- c(grep("SNu", unique(comb2$Dataset), value = T),
                            grep("Spear|10x", unique(comb2$Dataset), value = T),
                            grep("sci-ATAC-seq3_", unique(comb2$Dataset), value = T),
                            grep("^sci-ATAC-seq_|Easy", unique(comb2$Dataset), value = T),
                            grep("dsci|txci|scifi", unique(comb2$Dataset), value = T),
                            grep("sciPlex|sciMAP", unique(comb2$Dataset), value = T))

colors_short <- colors_datasets
names(colors_short) <- c("SNU_A","SNU_B","SPEAR","10X","SCI3_A","SCI3_B","EASY","SCI","DSCI","SCIFI","TXCI","PLEX")

comb$Dataset <- factor(comb$Dataset,
                       levels = names(colors_datasets))
comb2$Dataset <- factor(comb2$Dataset,
                        levels = names(colors_datasets))
comb$Short <- factor(comb$Dataset, levels = names(colors_datasets), 
                     labels = names(colors_short))
comb2$Short <- factor(comb2$Dataset, levels = names(colors_datasets), 
                      labels = names(colors_short))
comb$Idx <- factor(comb$Dataset, levels = names(colors_datasets), labels = 1:12)
comb2$Idx <- factor(comb2$Dataset, levels = names(colors_datasets), labels = 1:12)

comb$Tn5Type <- factor(comb$Tn5Type, levels = c("Standard", "Indexed"))
comb2$Tn5Type <- factor(comb2$Tn5Type, levels = c("Standard", "Indexed"))

# comb2 %>%
#   ggplot(aes(x = log10(Count), y = log10(nFrags), group = paste(Dataset,Group), color = Dataset)) +
#   geom_smooth(method = 'lm', se = F) + theme_bw() +
#   geom_point(size = 0.3, alpha = 0.5) +
#   facet_wrap(~Tn5Type, scales = 'free') +
#   scale_color_manual(values = colors_datasets) +
#   labs(color = "Dataset")

g <- comb2 %>% 
  group_by(Dataset, Group, Tn5Type) %>% summarize(Slope = coef(lm(log10(nFrags)~log10(Count)))[2], Count = mean(Count)) %>%
  # group_by(Dataset, Tn5Type, Group) %>% summarize(Slope = coef(lm(nFrags~Count))[2], Count = mean(Count)) %>%
  group_by(Dataset, Tn5Type) %>% summarize(Slope = mean(Slope)) %>% 
  ggplot(aes(x = Tn5Type, y = Slope, fill = factor(Dataset, levels = Dataset))) +
  geom_col(position = "dodge") + geom_hline(yintercept = 0, linetype = "dashed") + theme_classic() +
  scale_fill_manual(values = colors_datasets) + 
  labs(fill = "Dataset", x = "Transposome", y = "Mean Slope") + theme_DC
save_plot(g, "Plots/SlopeByTn5.pdf", w=3, h=2, x.angle=0, h.just=0.5, xlab=NA, ylab=NA, show.legend="right")

g <- comb2 %>% 
  group_by(Dataset, Group, Tn5Type) %>% summarize(Slope = coef(lm(log10(nFrags)~log10(Count)))[2], Count = mean(Count)) %>%
  group_by(Dataset, Tn5Type) %>% summarize(SD = sd(Slope), Slope = mean(Slope)) %>%
  ggplot(aes(x = Slope, y = Dataset, fill = factor(Dataset, levels = names(colors_datasets)))) +
  geom_col(position = "dodge", alpha = 0.8) + geom_vline(xintercept = 0, linetype = "dashed") + theme_classic() +
  geom_errorbar(mapping = aes(xmin = Slope - SD, xmax = Slope + SD, color = Dataset), width = 0.4, size = 0.5, show.legend = F) + 
  scale_fill_manual(values = colors_datasets) + scale_color_manual(values = colors_datasets) + facet_wrap(~Tn5Type, ncol = 1, scales='free_y', labeller = labeller(Tn5Type = c(Standard = "Standard Tn5", Indexed = "Indexed Tn5"))) + scale_y_discrete(limits = rev) +
  labs(fill = "Dataset", x = "Mean Slope", y = "Dataset") + theme_DC + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
save_plot(g, "Plots/SlopeByTn5_v2.pdf", w=2, h=2, x.angle=0, h.just=0.5, xlab=NA, ylab=NA, show.legend="none")


## Range of Sample Sizes

g <- comb2 %>%
  ggplot(aes(x = Count, y = Dataset, color = Dataset)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.5) + 
  theme_linedraw() + theme(legend.position = 'none') +
  scale_color_manual(values = colors_datasets) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) + 
  scale_y_discrete(limits = rev) +
  coord_cartesian(clip="off") +
  annotation_logticks(sides = "b", outside = T,
                      short = unit(0.05, "cm"),
                      mid = unit(0.075, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2) +
  theme_DC
save_plot(g, "Plots/CountRanges.pdf", w=4, h=2, x.angle=0, h.just=0.5, show.legend="none")

g <- comb2 %>%
  ggplot(aes(x = Count, y = Short, color = Dataset)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.5) + 
  theme_linedraw() + theme(legend.position = 'none') +
  scale_color_manual(values = colors_datasets) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) + 
  scale_y_discrete(limits = rev) +
  coord_cartesian(clip="off") +
  annotation_logticks(sides = "b", outside = T,
                      short = unit(0.05, "cm"),
                      mid = unit(0.075, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2) +
  theme_DC + labs(x = "Nuclei per Transposition", y = "Dataset")
# save_plot(g, "Plots/CountRanges_v2.pdf", w=4, h=2, x.angle=0, h.just=0.5, show.legend="none")
save_plot(g, "Plots/CountRanges_v2.pdf", w=3, h=2, x.angle=0, h.just=0.5, show.legend="none")

# comb2 %>%
#   ggplot(aes(x = Count, y = Dataset, color = Dataset)) +
#   geom_jitter(height = 0.2) +
#   scale_x_log10() + theme_linedraw() + theme(legend.position = 'none') +
#   scale_color_manual(values = colors_datasets)

# comb2 %>%
#   ggplot(aes(x = Count, y = Dataset, color = Dataset)) + 
#   geom_boxplot(outliers = F, linewidth = 1) + 
#   geom_jitter(size = 1, alpha = 0.5) +
#   scale_x_log10() + theme_linedraw() + theme(legend.position = 'none') +
#   scale_color_manual(values = colors_datasets)

# comb2 %>% group_by(Dataset) %>% summarize(CV = sd(Count)/mean(Count)) %>%
#   ggplot(aes(x = CV, y = Dataset, fill = Dataset, color = Dataset)) + 
#   geom_col() + theme_linedraw() + theme(legend.position = 'none') +
#   scale_color_manual(values = colors_datasets) + scale_fill_manual(values = colors_datasets)

# comb2 %>% group_by(Dataset) %>% 
#   summarize(FoldRange = (quantile(Count, 0.99)/ quantile(Count, 0.01)) %>% round(digits = 0)) %>%
#   ggplot(aes(x = FoldRange, y = Dataset, fill = Dataset, color = Dataset)) + 
#   geom_col() + theme_linedraw() + theme(legend.position = 'none') +
#   scale_color_manual(values = colors_datasets) + scale_fill_manual(values = colors_datasets)

comb2 %>% group_by(Dataset) %>% 
  summarize(Min = min(Count),
            Max = max(Count),
            `01` = quantile(Count, 0.01) %>% round,
            `99` = quantile(Count, 0.99) %>% round,
            Range_q1q99 = paste(`01`,`99`,sep="-"),
            FoldRange_q1q99 = (`99`/`01`) %>% round(digits = 0),
            nTn5 = length(unique(paste(Group,Tn5))),
            # nBatch = length(unique(Batch))
            ) %>% .[length(unique(comb2$Dataset)):1,] %>% View()

comb2 %>% group_by(Dataset) %>% 
  summarize(Min = min(Count),
            Max = max(Count),
            `01` = quantile(Count, 0.01) %>% round,
            `99` = quantile(Count, 0.99) %>% round,
            Min2 = Count %>% sort %>% head(n = 2) %>% tail(n = 1),
            Max2 = Count %>% sort %>% tail(n = 2) %>% head(n = 1),
            FoldRange_q1q99 = (`99`/`01`) %>% round(digits = 0),
            nTn5 = length(unique(paste(Group,Tn5))),
  ) #%>% .[length(unique(comb2$Dataset)):1,] %>% View()

comb2 %>% group_by(Dataset) %>% 
  mutate(Min = min(Count),
            Max = max(Count),
            `01` = quantile(Count, 0.01), # %>% round,
            `99` = quantile(Count, 0.99), # %>% round,
            FoldRange = (`99`/`01`) %>% round(digits = 0),
            nTn5 = length(unique(paste(Group,Tn5)))) %>% 
  group_by(Dataset, Min, Max, `01`, `99`, FoldRange, nTn5) %>%
  filter(Count > `01`, Count < `99`) %>% 
  summarize(FoldRange2 = (max(Count)/min(Count)) %>% round(digits = 0)) %>% View()
  # .[length(unique(comb2$Dataset)):1,] %>% View()

# fetal chromatin sci-ATAC-seq3 dataset: each tissue sample transposed in 4 rxns and pooled for ligation 1
# drosophila dev sci-ATAC-seq3 dataset: each timepoint transposed in 11 rxns and pooled for ligation 1
# sciPlex-ATAC-seq dataset generated from 96 pooled transposiition reactions

##########################################################################################################################################################

### Quantifying Batch Effects -----------------------------------------------

## LISI Batch Mixing

lisi_umap <- plotLISI(binned %>% filter(!is.na(UMAP1)), "UMAP", returnDF = T) %>% filter(Dataset %ni% c("SNuBar_HBCA", "sci-ATAC-seq3_FetalAtlas"))
lisi_lsi <- rbind(data.frame(plotLISI(proj_2, "LSI", "Tn5", c("Condition","Donor","Load"), bins = 3, returnDF = T), Dataset = "dsci-ATAC-seq_BMMC"),
                  data.frame(plotLISI(proj_tx2, "LSI", "Tn5", c("Tissue","Replicate","NucInput"), bins = 3, returnDF = T), Dataset = "txci-ATAC-seq_Barnyard"),
                  data.frame(plotLISI(proj_snu2, "LSI", "Tn5", c("CellType"), bins = 3, returnDF = T), Dataset = "SNuBar_CellLines"),
                  data.frame(plotLISI(proj_plex2, "LSI", "Tn5", c("cell_type"), bins = 3, returnDF = T), Dataset = "sciPlex-ATAC-seq2_DrugScreen"))
lisi_lsi_no1 <- rbind(data.frame(plotLISI(proj_2, "LSI", "Tn5", c("Condition","Donor","Load"), bins = 3, dims = 2:30, returnDF = T), Dataset = "dsci-ATAC-seq_BMMC"),
                      data.frame(plotLISI(proj_tx2, "LSI", "Tn5", c("Tissue","Replicate","NucInput"), bins = 3, dims = 2:30, returnDF = T), Dataset = "txci-ATAC-seq_Barnyard"),
                      data.frame(plotLISI(proj_snu2, "LSI", "Tn5", c("CellType"), bins = 3, dims = 2:30, returnDF = T), Dataset = "SNuBar_CellLines"),
                      data.frame(plotLISI(proj_plex2, "LSI", "Tn5", c("cell_type"), bins = 3, returnDF = T), Dataset = "sciPlex-ATAC-seq2_DrugScreen"))

lisi_umap$Dataset <- factor(lisi_umap$Dataset, levels = names(colors_datasets))
lisi_lsi$Dataset <- factor(lisi_lsi$Dataset, levels = names(colors_datasets))
lisi_lsi_no1$Dataset <- factor(lisi_lsi_no1$Dataset, levels = names(colors_datasets))

g <- ggplot(lisi_lsi %>% reshape2::melt(), aes(x = Group, y = value, fill = Dataset, linetype = variable)) +
  geom_boxplot(outliers = F, linewidth = 0.5) +
  theme_classic() + labs(y = "LISI", x = NULL, fill = "Dataset", linetype = "Grouping") +
  scale_fill_manual(values = colors_datasets) +
  scale_linetype_manual(values = c("solid","11")) +
  guides(fill = 'none') +
  facet_grid(~Dataset, scales = "free_x", space = 'free_x') +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 2, vjust = 0) + theme_DC + theme(strip.clip = "off", strip.text = element_blank())
save_plot(g, "./Plots/LISI.pdf", w = 6, h = 3, x.angle = 90, h.just = 1, show.legend = 'bottom')

g <- ggplot(lisi_lsi_no1 %>% reshape2::melt(), aes(x = Group, y = value, fill = Dataset, linetype = variable)) +
  geom_boxplot(outliers = F, linewidth = 0.5) +
  theme_classic() + labs(y = "LISI", x = NULL, fill = "Dataset", linetype = "Grouping") +
  scale_fill_manual(values = colors_datasets) +
  scale_linetype_manual(values = c("solid","11")) +
  facet_grid(~Dataset, scales = "free_x", space = 'free_x') +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 2, vjust = 0) + theme_DC + theme(strip.clip = "off", strip.text = element_blank())
save_plot(g, "./Plots/LISI_NoDim1.pdf", w = 6, h = 3, x.angle = 90, h.just = 1, show.legend = 'bottom')


# g <- 
lisi_lsi %>% reshape2::melt() %>%
  mutate(Short = factor(Dataset, levels = names(colors_datasets), labels = names(colors_short)),
         variable = factor(variable, levels = c("Bin","Random"), labels = c("Binned", "Shuffled"))) %>%
  group_by(Dataset, Short, variable, Group) %>% summarize(value = mean(value)) %>%
  # group_by(Dataset, Short, variable) %>% summarize(sd = sd(value), value = mean(value)) %>%
  ggplot(aes(x = variable, y = value, fill = Dataset)) +
  geom_boxplot(outliers = F, linewidth = 0.5) +
  # geom_jitter() +
  # geom_col() + geom_errorbar(aes(ymin = value-sd, ymax = value+sd)) +
  theme_classic() + labs(y = "LISI", x = NULL, fill = "Dataset") +
  scale_fill_manual(values = colors_datasets) +
  scale_linetype_manual(values = c("solid","11")) +
  guides(fill = 'none') +
  facet_grid(~Short, scales = "free_x", space = 'free_x') +
  # stat_compare_means(method = "wilcox.test", label = "p.signif", size = 2, vjust = 0) + theme_DC + theme(strip.clip = "off", strip.text = element_blank())
  stat_compare_means(aes(group = variable), method = "t.test", paired = T, ref.group = "Shuffled", label = "p.signif", size = 2, vjust = 0) + theme_DC #+ theme(strip.clip = "off", strip.text = element_blank())


# Ratio Between Shuffled & Binned LISI Scores
g <- lisi_umap %>% group_by(Dataset, Group) %>% summarize(Bin = mean(Bin), Random = mean(Random)) %>% mutate(Ratio = Bin/Random) %>%
  ggplot(aes(x = Dataset, y = Ratio, color = Dataset)) + 
  geom_boxplot(outliers = F) + geom_point() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  labs(y = "Relative Mixing") + scale_color_manual(values = colors_datasets) + theme_DC

g <- lisi_lsi %>% group_by(Dataset, Group) %>% summarize(Bin = mean(Bin), Random = mean(Random)) %>% mutate(Ratio = Bin/Random) %>%
  # filter(Dataset != "sciPlex-ATAC-seq2_DrugScreen") %>%
  mutate(Short = factor(Dataset, levels = names(colors_datasets), labels = names(colors_short))) %>%
  ggplot(aes(x = Short, y = Ratio, color = Dataset)) + 
  geom_boxplot(outliers = F) + 
  ggbeeswarm::geom_beeswarm() +
  # geom_jitter(width = 0.1) +
  facet_grid(~Short, scales = 'free_x') +
  labs(y = "Shuffled LISI / Binned LISI", x = "Dataset", subtitle = "Batch Mixing (LSI1-30)") + scale_color_manual(values = colors_datasets) + 
  theme_DC + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.clip = 'off')
save_plot(g, "./Plots/LISI_Ratio.pdf", w = 1.75, h = 2, show.legend = 'none')

g <- lisi_lsi_no1 %>% group_by(Dataset, Group) %>% summarize(Bin = mean(Bin), Random = mean(Random)) %>% mutate(Ratio = Bin/Random) %>%
 # filter(Dataset != "sciPlex-ATAC-seq2_DrugScreen") %>%
  mutate(Short = factor(Dataset, levels = names(colors_datasets), labels = names(colors_short))) %>%
  ggplot(aes(x = Short, y = Ratio, color = Dataset)) + 
  geom_boxplot(outliers = F) + 
  ggbeeswarm::geom_beeswarm() +
  # geom_jitter(width = 0.1) +
  facet_grid(~Short, scales = 'free_x') +
  labs(y = "Shuffled LISI / Binned LISI", x = "Dataset", subtitle = "Batch Mixing (LSI2-30)") + scale_color_manual(values = colors_datasets) + 
  theme_DC + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.clip = 'off')
save_plot(g, "./Plots/LISI_Ratio_NoDim1.pdf", w = 1.75, h = 2, show.legend = 'none')

# LISI Mean Values
# g <- lisi_umap %>% reshape2::melt() %>% group_by(Dataset, Group, variable) %>% summarize(value = mean(value)) %>%
g <- lisi_lsi %>% reshape2::melt() %>% group_by(Dataset, Group, variable) %>% summarize(value = mean(value)) %>%
  # ggplot(aes(x = variable, y = value, color = Dataset, linetype = variable)) + 
  ggplot(aes(x = Dataset, y = value, color = Dataset, linetype = variable)) + 
  geom_boxplot(outliers = F) + scale_linetype_manual(values = c("solid","11")) +
  # geom_point(mapping = aes(group = variable), position = position_dodge(width = 0.75)) +
  labs(x = NULL, y = "LISI") + scale_color_manual(values = colors_datasets) +
  facet_grid(~Dataset, scales = 'free_x',margins = F) +
  theme(legend.position = 'none') + theme_DC + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 2, vjust = 0)
save_plot(g, "./Plots/LISI_Means.pdf", w = 2.5, h = 2.5, x.angle = 25, h.just = 1, show.legend = 'none')


### Binned Analysis ---------------------------------------------------------

### Cell Type Fold Changes

## Pseudo-code
# 1. select datasets where individual, heterogeneous tissue samples were sampled many times for transposition
# 2. tabulate how many cells per Tn5 reaction, accounting for experimental batch
# 3. rank Tn5 reaction sizes (cell counts), accounting for experimental batch
# 4. assign cells into a Tn5 reaction size tercile, as well as a random assignment
# 5. keep celltypes for each biogroup that represent some minimum proportion of total population
# 6. calculate fold changes for top 10 cell types between Bin 3 (high Tn5 Count) and Bin 1 (low Tn5 Count)
# 7. calculate fold changes for top 10 cell types between cells randomly assigned 3 and cells randomly assigned 1

nbin <- 3
freq <- 0.05
# freq <- 0.1

colors_bins <- viridis::plasma(n = 11)[c(2,6,10)] %>% rev

set.seed(1)
binned <- comb %>% 
  group_by(Dataset, Tn5, Group) %>% mutate(Tn5Count = n()) %>%           # count nCells per Tn5 reaction per batch/dataset
  group_by(Dataset, Group) %>% mutate(Tn5Rank = rank(Tn5Count),          # rank Tn5 reaction counts per batch/dataset
                                             Bin = factor(Tn5Rank, levels = sort(unique(Tn5Rank)), labels = ntile(sort(unique(Tn5Rank)), nbin)),
                                             Random = sample(Bin, size = n(), replace=F))  # scramble Bin assignments (to get same counts per bin)

binned <- binned %>% 
  group_by(Dataset, Group, CellType) %>% mutate(CellTypeN = n()) %>%          # count nCells per cell type per biogroup/batch/dataset
  group_by(Dataset, Group) %>% mutate(CellTypeFreq = CellTypeN/n(),           # calculate frequency of each cell type per biogroup/dataset
                                      Keep = case_when(                       # determine which cell types to keep for fold-change analysis
                                        CellTypeFreq >= freq ~ "Yes",
                                        CellTypeFreq < freq ~ "No"
                                      )) 

# calculate fold-changes of each cell type between top bin and bottom bin using binned and scrambled assignments
celltypefreq <- binned %>% 
  filter(Dataset %in% c("scifi-ATAC_Maize",
                        "dsci-ATAC-seq_BMMC",
                        "txci-ATAC-seq_Barnyard",
                        "sci-ATAC-seq_HumanAtlas",
                        "EasySci-ATAC_MouseBrain"))

# Top vs Bottom Bin
# celltypefreq <- rbind(data.frame(celltypefreq %>% 
#                                    group_by(Dataset, Group, CellType) %>% mutate(nFragsMean = mean(nFrags)) %>%
#                                    group_by(Dataset, Group, Bin, Keep, CellType, nFragsMean) %>% 
#                                    summarize(n = n()) %>% group_by(Dataset, Group, Bin) %>% mutate(Freq = n / sum(n)) %>%
#                                    group_by(Dataset, Group, Keep, CellType, nFragsMean) %>% summarize(FoldChange = Freq[Bin %in% nbin] / Freq[Bin %in% 1]) %>%
#                                    group_by(Dataset, Group) %>% mutate(nFragsRank = rank(nFragsMean)),
#                                  Comparison = "Binned"),
#                       data.frame(celltypefreq %>% 
#                                    group_by(Dataset, Group, CellType) %>% mutate(nFragsMean = mean(nFrags)) %>%
#                                    group_by(Dataset, Group, Random, Keep, CellType, nFragsMean) %>% 
#                                    summarize(n = n()) %>% group_by(Dataset, Group, Random) %>% mutate(Freq = n / sum(n)) %>%
#                                    group_by(Dataset, Group, Keep, CellType, nFragsMean) %>% summarize(FoldChange = Freq[Random %in% nbin] / Freq[Random %in% 1]) %>%
#                                    group_by(Dataset, Group) %>% mutate(nFragsRank = rank(nFragsMean)),
#                                  Comparison = "Shuffled"))

# Bottom vs Top Bin
celltypefreq <- rbind(data.frame(celltypefreq %>% 
                                   group_by(Dataset, Group, CellType) %>% mutate(nFragsMean = mean(nFrags)) %>%
                                   group_by(Dataset, Short, Group, Bin, Keep, CellType, nFragsMean) %>% 
                                   summarize(n = n()) %>% group_by(Dataset, Group, Bin) %>% mutate(Freq = n / sum(n)) %>%
                                   group_by(Dataset, Short, Group, Keep, CellType, nFragsMean) %>% summarize(FoldChange = Freq[Bin %in% 1] / Freq[Bin %in% nbin]) %>%
                                   group_by(Dataset, Group) %>% mutate(nFragsRank = rank(nFragsMean)),
                                 Comparison = "Binned"),
                      data.frame(celltypefreq %>% 
                                   group_by(Dataset, Group, CellType) %>% mutate(nFragsMean = mean(nFrags)) %>%
                                   group_by(Dataset, Short, Group, Random, Keep, CellType, nFragsMean) %>% 
                                   summarize(n = n()) %>% group_by(Dataset, Group, Random) %>% mutate(Freq = n / sum(n)) %>%
                                   group_by(Dataset, Short, Group, Keep, CellType, nFragsMean) %>% summarize(FoldChange = Freq[Random %in% 1] / Freq[Random %in% nbin]) %>%
                                   group_by(Dataset, Group) %>% mutate(nFragsRank = rank(nFragsMean)),
                                 Comparison = "Shuffled"))


# plot each celltype per dataset/group
celltypefreq %>% filter(Comparison == "Binned") %>%
  filter(Keep == "Yes") %>%
  # filter(Dataset == "dsci-ATAC-seq_BMMC") %>%
  filter(Dataset == "txci-ATAC-seq_Barnyard") %>%
  # filter(Dataset == "scifi-ATAC_Maize") %>%
  # filter(Dataset == "EasySci-ATAC_MouseBrain") %>%
  mutate(FoldChange = log2(FoldChange)) %>%
  # group_by(Dataset, Group, Comparison) %>% mutate(nFragsRank = nFragsRank / max(nFragsRank) * 100) %>%
  ggplot(aes(x = nFragsRank, y = FoldChange, fill = CellType)) + 
  geom_col() + geom_hline(yintercept = 0, linetype = '22') + 
  facet_wrap(~Group, scales = 'free', nrow = 2) + theme_DC +
  labs(x = "Cell Types Ranked by nFrags", y = "Log2 Fold Change (Top vs Bottom Bin)")

# plot each grouping per dataset
celltypefreq %>%
  filter(Keep == "Yes") %>%
  mutate(FoldChange = log2(FoldChange)) %>%
  group_by(Dataset, Group, Comparison) %>% mutate(nFragsRank = nFragsRank / max(nFragsRank) * 100) %>%
  ggplot(aes(x = nFragsRank, y = FoldChange, color = Dataset, group = Group)) + 
  geom_smooth(method = 'lm', se = F, size = 1.5, alpha = 0.5) + geom_hline(yintercept = 0, linetype = '22') + 
  facet_wrap(Dataset~Comparison, scales = 'free_x', nrow = 5) + theme_DC + theme(legend.position = 'none') + 
  scale_color_manual(values = colors_datasets) + scale_fill_manual(values = colors_datasets) +
  labs(x = "Cell Types Ranked by nFrags", y = "Log2 Fold Change (Top vs Bottom Bin)")

# aggregate groupings per dataset
g <- celltypefreq %>% 
  filter(Keep == "Yes") %>%
  mutate(FoldChange = log2(FoldChange)) %>%
  group_by(Dataset, Group, Comparison) %>% mutate(nFragsRank = nFragsRank / max(nFragsRank) * 100) %>%
  ggplot(aes(x = nFragsRank, y = FoldChange, color = Short, fill = Short)) + 
  facet_wrap(~Comparison) + theme_DC + theme(axis.ticks.x = element_blank(), legend.background = element_blank()) +
  # geom_point(stroke = 1, alpha = 0.2) +
  geom_smooth(method = 'lm', se = F, size = 1, alpha = 0) +
  geom_smooth(method = 'lm', se = T, size = 1, alpha = 0.2, show.legend = F) + geom_hline(yintercept = 0, linetype = '22') + 
  scale_color_manual(values = colors_short) + scale_fill_manual(values = colors_short) + scale_x_continuous(breaks = c(15,90), labels = c("Fewer", "More")) + 
  # labs(x = "Cell Types Ordered by nFrags", y = "Log2 Fold Change (Top vs Bottom Bin)") 
  labs(x = "Cell Types Ordered by nFrags", y = "Cell Type Proportion Log2FC", color = "", fill = "")
save_plot(g, "./Plots/CellTypeProportionFC_nFragsRanked.pdf", w = 2, h = 2, show.legend = c(0.75,0.3))

# violin plots
g <- ggplot(filter(celltypefreq, Comparison == "Binned", Keep == "Yes"), aes(x = log2(FoldChange), y = factor(Dataset, levels = names(colors_datasets)), color = Dataset, fill = Dataset)) +
  # geom_point(alpha = 0.5) +
  ggbeeswarm::geom_beeswarm(size = 1, alpha = 0.5, corral = "random", corral.width = 0.6, stroke = 0.3) +
  geom_vline(xintercept = 0, linetype = '11', color = 'grey60', size = 1) +
  geom_violin(trim = T, size = 0.75, alpha = 0.3, scale = "width") +
  geom_violin(data = filter(celltypefreq, Comparison == "Shuffled", Keep == "Yes"), size = 0.75, color = "black", fill = "black", alpha = 0.3, trim = T, scale = 'width') +
  theme_linedraw() +
  scale_color_manual(values = colors_datasets) + scale_fill_manual(values = colors_datasets) +
  scale_x_continuous(limits=ggpmisc::symmetric_limits) +
  scale_y_discrete(limits = rev) +
  # labs(x = "Cell Type Proportion Log2FC", y = "Dataset", title = "Cell Type Proportion Fold Change (Top vs Bottom Bin)") +
  labs(x = "Cell Type Proportion Log2FC", y = "Dataset") + #ggtitle("Cell Type Proportion Fold Change (Bottom vs Top Bin)")
  ggrepel::geom_text_repel(
    data = celltypefreq %>% filter(Comparison == "Binned", Keep == "Yes") %>% group_by(Dataset) %>% slice_max(order_by = -log2(FoldChange), n = 1),
    aes(x = log2(FoldChange), y = Dataset, label = CellType), min.segment.length = 0,
    # size = 2.5, alpha = 0.7, nudge_x = -0.6, nudge_y = -0.4, color = "black") +
    nudge_x = -0.75, nudge_y = -0.5,
    size = 2, color = "black", bg.color = "white") +
  ggrepel::geom_text_repel(
    data = celltypefreq %>% filter(Comparison == "Binned", Keep == "Yes") %>% group_by(Dataset) %>% slice_max(order_by = log2(FoldChange), n = 1),
    aes(x = log2(FoldChange), y = Dataset, label = CellType), min.segment.length = 0,
    # size = 2.5, alpha = 0.7, nudge_x = 1, nudge_y = 0.4, color = "black") +
    nudge_x = 0.75, nudge_y = 0.5,
    size = 2, color = "black", bg.color = "white") +
  # geom_label(data = dplyr::select(celltypefreq, Dataset, Short, Comparison),
  #            # mapping = aes(y = Dataset, fill = Dataset, label = gsub("_.*","",Dataset)),
  #            mapping = aes(y = Dataset, fill = Dataset, label = Short),
  #            # x = -Inf + 1, stat = "unique", size = 3, color = 'black', nudge_y = 0.25, nudge_x = 0.1, hjust=-0.01) +
  #            x = -Inf + 1, stat = "unique", size = 2, color = 'black', nudge_y = 0, nudge_x = 1, vjust=1, angle = 90) +
  # annotate(geom = "label", vjust = -0.5, label = "Shuffled", x = 0, y = -0.25, color = 'black', fill = 'grey80', alpha = 1, size = 3) +
  # annotate(geom = "label", vjust = -0.5, label = "Shuffled", x = 0, y = -0.25, color = 'black', fill = 'grey80', alpha = 1, size = 2) +
  # geom_label(data = data.frame(Short = unique(celltypefreq$Short), Label = c(NA,NA,NA,NA,"Shuffled")), mapping = aes(label = Label),  x = 0, y = 0.5, color = 'black', fill = 'grey80', alpha = 1, size = 2) + 
  # geom_label(data = data.frame(Short = unique(celltypefreq$Short), Label = c("Shuffled",NA,NA,NA,NA)), mapping = aes(label = Label),  x = 0, y = 1.75, color = 'black', fill = 'grey80', alpha = 1, size = 2) + 
  facet_grid(rows = "Short", scales = 'free', space = 'free_y', switch = "y") + 
  theme_DC + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                   axis.text.y = element_blank(),legend.position = "none",
                   panel.spacing = unit(0, "lines")) + coord_cartesian(clip = 'off')
# save_plot(g, "./Plots/CellTypeProportionFC.pdf", w = 4, h = 3, show.legend = "none")
save_plot(g, "./Plots/CellTypeProportionFC.pdf", w = 3, h = 2, show.legend = "none")



### scifi-ATAC-seq

ds <- "scifi-ATAC_Maize"
group <- "M162W"
# group <- "B97"
# group <- "B73"
# group <- "Ky21"

## Tn5 x nFrags
# binned %>% filter(Dataset == ds) %>%
g1 <- binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(Group, Tn5, Bin) %>% summarize(Tn5Count = mean(Tn5Count), nFrags = median(nFrags)) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = Tn5Count, y = nFrags, color = Bin)) + 
  geom_point(size = 2) + geom_smooth(method = 'lm', se = F, color = "black") + facet_wrap(~Group, nrow = 2, scales = "free") + 
  scale_color_manual(values = colors_bins) + labs(title = "Individual Tn5 Reactions", y = "Median nFrags", x = "Cell Count") + theme_DC
save_plot(g1, "./Plots/scifiATAC_Bins.pdf", w = 2.5, h = 3, show.legend = "right")

## Cell Type Frequencies by Bin
# binned %>% filter(Dataset == ds) %>%
g2 <- binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(CellType, Group, Bin) %>% summarize(n = n()) %>% group_by(Group, Bin) %>% mutate(freq = n / sum(n)) %>%
  group_by(CellType, Group) %>% mutate(Diff = freq[Bin == 3] - freq[Bin == 1]) %>% ungroup %>%
  group_by(Group) %>% filter(Diff %in% c(min(Diff),max(Diff))) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = CellType, y = freq, fill = Bin)) + geom_col(position = "dodge") + facet_wrap(~Group, nrow = 1, scales = "free") +
  labs(title = "Cell Type Frequencies", y = "Frequency", x = "Cell Type") + scale_fill_manual(values = colors_bins) + theme_DC
save_plot(g2, "./Plots/scifiATAC_CellTypeFreq.pdf", w = 2.5, h = 3, show.legend = "none")

g2 <-
  binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(CellType, Group, Bin) %>% summarize(n = n()) %>% group_by(Group, Bin) %>% mutate(total = sum(n), freq = n / total) %>%
  group_by(CellType, Group) %>% mutate(Diff = freq[Bin == 3] - freq[Bin == 1]) %>% ungroup %>%
  group_by(Group) %>% filter(Diff %in% c(min(Diff),max(Diff))) %>% 
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>% ungroup %>%
  {
    res1 <- pairwise.prop.test(as.matrix(select(., n, total))[1:3,]) %>% broom::tidy()
    res2 <- pairwise.prop.test(as.matrix(select(., n, total))[4:6,]) %>% broom::tidy()
    df <- cbind(Bin = as.numeric(.$Bin), CellType = .$CellType, freq = .$freq, rbind(res1,res2), Type = c(1,1,1,-1,-1,-1))
    df$p.value[6] <- df$p.value[6] + 0.000001
    ggplot(., aes(x = as.numeric(Bin), y = freq)) +
      geom_col(position = "dodge", mapping = aes(fill = Bin)) +
      labs(subtitle = paste("SCIFI -", group), y = "Cell Type Frequency", x = NULL, fill = "Bin") + scale_fill_manual(values = colors_bins) + theme_DC +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), strip.clip = 'off') +
      facet_wrap(~gsub("_pr","\npr",CellType) %>% gsub("_"," ",.), switch = "x") + coord_cartesian(clip = 'off') +
      geom_signif(data = df,
                  mapping = aes(x = as.numeric(Bin),
                                y = freq,
                                xmin = as.numeric(group1),
                                xmax = as.numeric(group2),
                                y_position = (0.3 + -Type * 0.175) + ((as.numeric(Bin)-0.5) * 0.035) * Type,
                                annotations = format(p.value, digits = 2)),
                  textsize = 2,
                  manual = T)
  }
save_plot(g2, "./Plots/scifiATAC_CellTypeFreq_v2.pdf", w = 3, h = 2, show.legend = "right")

## UMAP
g3 <- binned %>% filter(Dataset == ds) %>% filter(Group == group) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Bin)) + 
  ggrastr::rasterise(geom_point(size = 0.5, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
  labs(title = "Tn5 Reaction Bin") + facet_wrap(Group~.) + scale_color_manual(values = colors_bins) + 
  theme_DC + theme(legend.position = 'none', axis.ticks = element_blank(), axis.text = element_blank())
save_plot(g3, "./Plots/scifiATAC_BinsUMAP.pdf", w = 2.5, h = 2.5, show.legend = "none")

g4 <- binned %>% filter(Dataset == ds) %>%
  {
    ggplot(.,aes(x = UMAP1, y = UMAP2, color = CellType)) + 
      ggrastr::rasterise(geom_point(size = 0.5, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
      ggrepel::geom_label_repel(data = summarize(group_by(., CellType), UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), CellTypeFreq = mean(CellTypeFreq)),
                                mapping = aes(label = CellType, alpha = CellTypeFreq > freq),
                                size = 2, label.padding = 0.1, max.overlaps = 20, show.legend = F) +
      labs(title = "Cell Type") + scale_alpha_manual(values = c(0.2, 0.9)) + scale_color_manual(values = paletteer::paletteer_d("ggsci::category20_d3")) +
      # labs(title = "Cell Type") + scale_alpha_manual(values = c(0.2, 0.9)) + scale_color_manual(values = colorRampPalette(c("navyblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1", "maroon"))(length(unique(.$CellType)))) +
      theme_DC + theme(legend.position = 'none', axis.ticks = element_blank(), axis.text = element_blank())
  }
save_plot(g4, "./Plots/scifiATAC_CellTypeUMAP.pdf", w = 2.5, h = 2.5, show.legend = "none")

g <- plot_grid(g1, g3, g2, g4, nrow = 2, rel_widths = c(3,4,3,4))
save_plot(g, "./Plots/scifiATAC_Summary.pdf", w = 5, h = 5, show.legend = "none")

g <- plot_grid(g2 + theme(axis.text.x = element_text(angle = -25, hjust = 0)), 
               g3, g4, nrow = 1, rel_widths = c(3.5,4,4))
save_plot(g, "./Plots/scifiATAC_Summary_v2.pdf", w = 6, h = 2, show.legend = "none")



### Binned: dsci-ATAC-seq ---------------------------------------------------

ds <- "dsci-ATAC-seq_BMMC"
# group <- "Stimulated A 20k"
# group <- "Stimulated A 40k"
# group <- "Baseline B 20k"
# group <- "Baseline B 40k"
group <- "Baseline B 80k"

## Tn5 x nFrags
# binned %>% filter(Dataset == ds) %>%
g1 <- binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(Group, Tn5, Bin) %>% summarize(Tn5Count = mean(Tn5Count), nFrags = median(nFrags)) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = Tn5Count, y = nFrags, color = Bin)) + 
  geom_point(size = 2) + geom_smooth(method = 'lm', se = F, color = "black") + facet_wrap(~Group, nrow = 2, scales = "free") + 
  scale_color_manual(values = colors_bins) + labs(title = "Individual Tn5 Reactions", y = "Median nFrags", x = "Cell Count") + theme_DC
save_plot(g1, "./Plots/dsciATAC_Bins.pdf", w = 2.5, h = 3, show.legend = "right")

## Cell Type Frequencies by Bin
# binned %>% filter(Dataset == ds) %>%
g2 <- binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(CellType, Group, Bin) %>% summarize(n = n()) %>% group_by(Group, Bin) %>% mutate(freq = n / sum(n)) %>%
  group_by(CellType, Group) %>% mutate(Diff = freq[Bin == 3] - freq[Bin == 1]) %>% ungroup %>%
  group_by(Group) %>% filter(Diff %in% c(min(Diff),max(Diff))) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = CellType, y = freq, fill = Bin)) + geom_col(position = "dodge") + facet_wrap(~Group, nrow = 1, scales = "free") +
  labs(title = "Cell Type Frequencies", y = "Frequency", x = "Cell Type") + scale_fill_manual(values = colors_bins) + theme_DC
save_plot(g2, "./Plots/dsciATAC_CellTypeFreq.pdf", w = 2.5, h = 3, show.legend = "none")


g2 <- binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(CellType, Group, Bin) %>% summarize(n = n()) %>% group_by(Group, Bin) %>% mutate(total = sum(n), freq = n / total) %>%
  group_by(CellType, Group) %>% mutate(Diff = freq[Bin == 3] - freq[Bin == 1]) %>% ungroup %>%
  group_by(Group) %>% filter(Diff %in% c(min(Diff),max(Diff))) %>% 
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>% ungroup %>%
  {
    res1 <- pairwise.prop.test(as.matrix(select(., n, total))[1:3,]) %>% broom::tidy()
    res2 <- pairwise.prop.test(as.matrix(select(., n, total))[4:6,]) %>% broom::tidy()
    df <- cbind(Bin = as.numeric(.$Bin), CellType = .$CellType, freq = .$freq, rbind(res1,res2), Type = c(1,1,1,-1,-1,-1))
    ggplot(., aes(x = as.numeric(Bin), y = freq)) + 
      geom_col(position = "dodge", mapping = aes(fill = Bin)) +
      labs(subtitle = paste("DSCI -", group), y = "Cell Type Frequency", x = NULL, fill = "Bin") + scale_fill_manual(values = colors_bins) + theme_DC +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), strip.clip = 'off') +
      facet_wrap(~gsub("$","\n",CellType), switch = 'x') + coord_cartesian(clip = 'off') +
      geom_signif(data = df, 
                  mapping = aes(x = as.numeric(Bin), 
                                y = freq, 
                                xmin = as.numeric(group1), 
                                xmax = as.numeric(group2),
                                y_position = 0.28 + ((as.numeric(Bin)-0.5) * 0.04) * Type,
                                annotations = format(p.value, digits = 2)), 
                  textsize = 2,
                  manual = T)
  }
save_plot(g2, "./Plots/dsciATAC_CellTypeFreq_v2.pdf", w = 3, h = 2, show.legend = "right")

## UMAP
g3 <- binned %>% filter(Dataset == ds) %>% filter(Group == group) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Bin)) + 
  ggrastr::rasterise(geom_point(size = 0.5, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
  labs(title = "Tn5 Reaction Bin") + facet_wrap(Group~.) + scale_color_manual(values = colors_bins) + 
  theme_DC + theme(legend.position = 'none', axis.ticks = element_blank(), axis.text = element_blank())
save_plot(g3, "./Plots/dsciATAC_BinsUMAP.pdf", w = 2.5, h = 2.5, show.legend = "none")

g4 <- binned %>% filter(Dataset == ds) %>%
  {
    ggplot(.,aes(x = UMAP1, y = UMAP2, color = CellType)) + 
      ggrastr::rasterise(geom_point(size = 0.5, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
      ggrepel::geom_label_repel(data = summarize(group_by(., CellType), UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), CellTypeFreq = mean(CellTypeFreq)),
                                mapping = aes(label = CellType, alpha = CellTypeFreq > freq),
                                size = 2, label.padding = 0.1, max.overlaps = 20, show.legend = F) +
      labs(title = "Cell Type") + scale_alpha_manual(values = c(0.2, 0.9)) + scale_color_manual(values = paletteer::paletteer_d("ggsci::category20_d3")) +
      # labs(title = "Cell Type") + scale_alpha_manual(values = c(0.2, 0.9)) + scale_color_manual(values = colorRampPalette(c("navyblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1", "maroon"))(length(unique(.$CellType)))) +
      theme_DC + theme(legend.position = 'none', axis.ticks = element_blank(), axis.text = element_blank())
  }
save_plot(g4, "./Plots/dsciATAC_CellTypeUMAP.pdf", w = 2.5, h = 2.5, show.legend = "none")

g <- plot_grid(g1, g3, g2, g4, nrow = 2, rel_widths = c(3,4,3,4))
save_plot(g, "./Plots/dsciATAC_Summary.pdf", w = 5, h = 5, show.legend = "none")

g <- plot_grid(g2 + theme(axis.text.x = element_text(angle = -25, hjust = 0)), 
               g3, g4, nrow = 1, rel_widths = c(3.5,4,4))
save_plot(g, "./Plots/dsciATAC_Summary_v2.pdf", w = 6, h = 2, show.legend = "none")



### Binned: dsci-ATAC-seq Donor Differences ---------------------------------

g <- proj_2@cellColData %>% as.data.frame %>% 
  filter(Load == "40k") %>%
  group_by(Load, Condition, Donor, Tn5) %>% summarize(Tn5Count = n(), nFrags = median(nFrags)) %>% 
  # filter(Tn5Count > 500) %>% filter(Condition == "Baseline") %>%
  ggplot(aes(x = Donor, y = nFrags, color = Donor)) + 
  geom_boxplot(outliers = F) + 
  # facet_wrap(~Load) +
  facet_wrap(~Condition) +
  # stat_compare_means(method = "wilcox.test", label = "p.signif", size = 2, show.legend = F) +
  theme_DC + scale_color_manual(values = RColorBrewer::brewer.pal(2, "Set1"))
save_plot(g, "./Plots/dsciATAC_Donors_nFrags.pdf", w = 2, h = 2, show.legend = "none")

# proj_2@cellColData %>% as.data.frame %>% filter(Load == "40k") %>%
#   mutate(Freq = Count/sum(Count)) %>%
#   ggplot(aes(x = paste(Condition, Donor), y = nFrags, fill = paste(Condition, Donor))) + geom_boxplot(outliers = F) + facet_wrap(~CellType) + scale_y_log10()

tmp <- binned %>% filter(Dataset == "dsci-ATAC-seq_BMMC", Group %in% grep(".*40k",Group,value =T))
g <- ggplot(mapping = aes(x = Donor, y = freq, fill = Bin)) + 
  geom_col(data = tmp %>% group_by(CellType, Group, Bin) %>% summarize(n = n()) %>% group_by(Group, Bin) %>% mutate(freq = n / sum(n)) %>%
             group_by(CellType, Group) %>% mutate(Diff = freq[Bin == max(as.numeric(Bin))] - freq[Bin == min(as.numeric(Bin))]) %>% ungroup %>%
             filter(CellType %in% c("CD4-like")) %>%
             mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
             separate(Group, into = c("Condition","Donor","Load")),
           mapping = aes(x = Donor, y = freq, fill = Bin),
           position = "dodge") + 
  geom_col(data = tmp %>% group_by(CellType, Group) %>% summarize(n = n()) %>% group_by(Group) %>% mutate(freq = n / sum(n)) %>%
             filter(CellType %in% c("CD4-like")) %>%
             separate(Group, into = c("Condition","Donor","Load")),
           mapping = aes(x = Donor, y = freq),
           color = 'black', linetype = "22", fill = NA) + 
  facet_wrap(~Condition) + scale_fill_manual(values = colors_bins) + ggtitle("CD4-like Cell Frequency") + theme_DC
save_plot(g, "./Plots/dsciATAC_Donors_CD4.pdf", w = 3, h = 2, show.legend = "right")



### Binned: txci-ATAC-seq ---------------------------------------------------

ds <- "txci-ATAC-seq_Barnyard"
group <- "liver m7 100k"
# group <- "lung m1 100k"

## Tn5 x nFrags
# binned %>% filter(Dataset == ds) %>%
g1 <- binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(Group, Tn5, Bin) %>% summarize(Tn5Count = mean(Tn5Count), nFrags = median(nFrags)) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = Tn5Count, y = nFrags, color = Bin)) + 
  geom_point(size = 2) + geom_smooth(method = 'lm', se = F, color = "black") + facet_wrap(~Group, nrow = 2, scales = "free") + 
  scale_color_manual(values = colors_bins) + labs(title = "Individual Tn5 Reactions", y = "Median nFrags", x = "Nuclei Count", color = NULL) + theme_DC
save_plot(g1, "./Plots/txciATAC_Bins.pdf", w = 2.5, h = 3, show.legend = "right")

g1 <- binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(Group, Tn5, Bin) %>% summarize(Tn5Count = mean(Tn5Count), nFrags = median(nFrags)) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = Tn5Count, y = nFrags, color = Bin)) + 
  geom_point(size = 2) + geom_smooth(method = 'lm', se = F, color = "black") + #facet_wrap(~Group, nrow = 2, scales = "free") + 
  scale_color_manual(values = colors_bins) + theme_DC +
  # labs(title = "Individual Tn5 Reactions", subtitle = paste("TXCI (",group,")",sep=""), y = "Median nFrags", x = "Nuclei Count", color = NULL)
  labs(subtitle = paste("TXCI -",group), y = "Median nFrags", x = "Nuclei Count", color = NULL)
  # labs(title = "TXCI", subtitle = group, y = "Median nFrags", x = "Nuclei Count", color = NULL)
save_plot(g1, "./Plots/txciATAC_Bins_v2.pdf", w = 1.5, h = 2, show.legend = c(0.35,0.875))

## Cell Type Frequencies by Bin
# binned %>% filter(Dataset == ds) %>%
# g2 <- binned %>% filter(Dataset == ds, Group == group) %>%
#   group_by(CellType, Group, Bin) %>% summarize(n = n()) %>% group_by(Group, Bin) %>% mutate(freq = n / sum(n)) %>%
#   group_by(CellType, Group) %>% mutate(Diff = freq[Bin == 3] - freq[Bin == 1]) %>% ungroup %>%
#   group_by(Group) %>% filter(Diff %in% c(min(Diff),max(Diff))) %>%
#   mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
#   ggplot(aes(x = CellType, y = freq, fill = Bin)) + geom_col(position = "dodge") + facet_wrap(~Group, nrow = 1, scales = "free") +
#   labs(title = "Cell Type Frequencies", y = "Frequency", x = "Cell Type") + scale_fill_manual(values = colors_bins) + theme_DC
# save_plot(g2, "./Plots/txciATAC_CellTypeFreq.pdf", w = 2.5, h = 3, show.legend = "none")
g2 <- binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(CellType, Group, Bin) %>% summarize(n = n()) %>% group_by(Group, Bin) %>% mutate(freq = n / sum(n)) %>%
  group_by(CellType, Group) %>% mutate(Diff = freq[Bin == 3] - freq[Bin == 1]) %>% ungroup %>%
  group_by(Group) %>% filter(Diff %in% c(min(Diff),max(Diff))) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = CellType, y = freq, fill = Bin)) + geom_col(position = "dodge") + #facet_wrap(~Group, nrow = 1, scales = "free") +
  labs(subtitle = paste("TXCI -", group), y = "Cell Type Frequency", x = "Cell Type", fill = NULL) + scale_fill_manual(values = colors_bins) + theme_DC
save_plot(g2, "./Plots/txciATAC_CellTypeFreq.pdf", w = 1.5, h = 2, show.legend = c(0.7,0.875))

g2 <- binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(CellType, Group, Bin) %>% summarize(n = n()) %>% group_by(Group, Bin) %>% mutate(total = sum(n), freq = n / total) %>%
  group_by(CellType, Group) %>% mutate(Diff = freq[Bin == 3] - freq[Bin == 1]) %>% ungroup %>%
  group_by(Group) %>% filter(Diff %in% c(min(Diff),max(Diff))) %>% 
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>% ungroup %>%
  {
    res1 <- pairwise.prop.test(as.matrix(select(., n, total))[1:3,]) %>% broom::tidy()
    res2 <- pairwise.prop.test(as.matrix(select(., n, total))[4:6,]) %>% broom::tidy()
    df <- cbind(Bin = as.numeric(.$Bin), CellType = .$CellType, freq = .$freq, rbind(res1,res2), Type = c(1,1,1,-1,-1,-1))
    ggplot(., aes(x = as.numeric(Bin), y = freq)) + 
      geom_col(position = "dodge", mapping = aes(fill = Bin)) +
      labs(subtitle = paste("TXCI -", group), y = "Cell Type Frequency", x = NULL, fill = NULL) + scale_fill_manual(values = colors_bins) + theme_DC +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), strip.clip = 'off') +
      # facet_wrap(~gsub(" n","\nn",CellType)) + coord_cartesian(clip = 'off') +
      facet_wrap(~CellType, switch = "x") + coord_cartesian(clip = 'off') +
      geom_signif(data = df, 
                  mapping = aes(x = as.numeric(Bin), 
                                y = freq, 
                                xmin = as.numeric(group1), 
                                xmax = as.numeric(group2),
                                y_position = (0.5 + Type * 0.325) + ((as.numeric(Bin)-0.5) * 0.06) * -Type,
                                annotations = format(p.value, digits = 2)), 
                  textsize = 2,
                  manual = T)
  }
save_plot(g2, "./Plots/txciATAC_CellTypeFreq_v2.pdf", w = 1.5, h = 2, show.legend = c(0.75,0.875))

## UMAP
g3 <- binned %>% filter(Dataset == ds) %>% filter(Group == group) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Bin)) + 
  ggrastr::rasterise(geom_point(size = 0.5, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
  labs(title = "Tn5 Reaction Bin") + facet_wrap(Group~.) + scale_color_manual(values = colors_bins) + 
  theme_DC + theme(legend.position = 'none', axis.ticks = element_blank(), axis.text = element_blank())
save_plot(g3, "./Plots/txciATAC_BinsUMAP.pdf", w = 2.5, h = 2.5, show.legend = "none")

# g3 <- 
  binned %>% filter(Dataset == ds) %>% filter(Group == group) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Bin)) + 
  geom_point(inherit.aes = F, mapping = aes(x = UMAP1, y = UMAP2), color = 'grey50', size = 0.5, stroke = 0) + 
  # ggrastr::rasterise(geom_point(size = 0.5, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
  geom_density_2d(bins = 10) +
  labs(title = "Tn5 Reaction Bin") + facet_wrap(Bin~.) + scale_color_manual(values = colors_bins) + 
  theme_DC + theme(legend.position = 'none', axis.ticks = element_blank(), axis.text = element_blank())
# save_plot(g3, "./Plots/txciATAC_BinsUMAP_v2.pdf", w = 2.5, h = 2.5, show.legend = "none")
  
binned %>% filter(Dataset == ds) %>% filter(Group == group) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  mutate(Cluster = factoextra::hkmeans(data.frame(UMAP1,UMAP2), k = 5)$cluster) %>%
  # ggplot(aes(x = UMAP1, y = UMAP2, color = as.character(Cluster))) + 
  ggplot(aes(x = UMAP1, y = UMAP2, color = Bin)) + 
  geom_point(inherit.aes = F, mapping = aes(x = UMAP1, y = UMAP2), color = 'grey50', size = 0.5, stroke = 0) + 
  ggrastr::rasterise(geom_point(size = 0.5, alpha = 1, stroke = 0, color = "lightgrey"), dpi = 600, scale = 1) +
  # ggrastr::rasterise(geom_point(size = 0.5, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
  # geom_density_2d(bins = 20) +
  geom_path(data = . %>% group_by(Cluster,Bin) %>% summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), mapping = aes(group = Cluster), color = 'black') + 
  geom_point(data = . %>% group_by(Cluster,Bin) %>% summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), size = 3) +
  # geom_point(data = . %>% group_by(Bin) %>% summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), size = 3) +
  labs(title = "Tn5 Reaction Bin") + #facet_wrap(Bin~.) + 
  scale_color_manual(values = colors_bins) +
  theme_DC + theme(legend.position = 'none', axis.ticks = element_blank(), axis.text = element_blank())

g4 <- binned %>% filter(Dataset == ds) %>% filter(Group %in% grep(str_sub(group, 1, 4), Group, value = T)) %>%
  {
    ggplot(.,aes(x = UMAP1, y = UMAP2, color = CellType)) + 
      ggrastr::rasterise(geom_point(size = 0.5, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
      ggrepel::geom_label_repel(data = summarize(group_by(., CellType), UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), CellTypeFreq = mean(CellTypeFreq)),
                                mapping = aes(label = CellType, alpha = CellTypeFreq > freq),
                                size = 3, label.padding = 0.1, max.overlaps = 20, show.legend = F) +
      labs(title = "Cell Type") + scale_alpha_manual(values = c(0.4, 0.9)) + scale_color_manual(values = paletteer::paletteer_d("ggsci::category20_d3")) +
      # labs(title = "Cell Type") + scale_alpha_manual(values = c(0.2, 0.9)) + scale_color_manual(values = colorRampPalette(c("navyblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1", "maroon"))(length(unique(.$CellType)))) +
      theme_DC + theme(legend.position = 'none', axis.ticks = element_blank(), axis.text = element_blank())
  }
save_plot(g4, "./Plots/txciATAC_CellTypeUMAP.pdf", w = 2.5, h = 2.5, show.legend = "none")

g <- plot_grid(g1, g3, g2, g4, nrow = 2, rel_widths = c(3,4,3,4))
save_plot(g, "./Plots/txciATAC_Summary.pdf", w = 5, h = 5, show.legend = "none")

g <- plot_grid(g2, g3, g4, nrow = 1, rel_widths = c(3.5,4,4))
save_plot(g, "./Plots/txciATAC_Summary_v2.pdf", w = 6, h = 2, show.legend = "none")



### Binned: EasySci-ATAC ----------------------------------------------------

ds <- "EasySci-ATAC_MouseBrain"
# group <- "Young"
group <- "Aged"
# group <- "WT"
# group <- "5xFAD"

## Tn5 x nFrags
# binned %>% filter(Dataset == ds) %>%
g1 <- binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(Group, Tn5, Bin) %>% summarize(Tn5Count = mean(Tn5Count), nFrags = median(nFrags)) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = Tn5Count, y = nFrags, color = Bin)) + 
  geom_point(size = 2) + geom_smooth(method = 'lm', se = F, color = "black") + facet_wrap(~Group, nrow = 2, scales = "free") + 
  scale_color_manual(values = colors_bins) + labs(title = "Individual Tn5 Reactions", y = "Median nFrags", x = "Cell Count") + theme_DC
save_plot(g1, "./Plots/EasySciATAC_Bins.pdf", w = 2.5, h = 3, show.legend = "right")

## Cell Type Frequencies by Bin
# binned %>% filter(Dataset == ds) %>%
g2 <- binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(CellType, Group, Bin) %>% summarize(n = n()) %>% group_by(Group, Bin) %>% mutate(freq = n / sum(n)) %>%
  group_by(CellType, Group) %>% mutate(Diff = freq[Bin == 3] - freq[Bin == 1]) %>% ungroup %>%
  group_by(Group) %>% filter(Diff %in% c(min(Diff),max(Diff))) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = CellType, y = freq, fill = Bin)) + geom_col(position = "dodge") + facet_wrap(~Group, nrow = 1, scales = "free") +
  labs(title = "Cell Type Frequencies", y = "Frequency", x = "Cell Type") + scale_fill_manual(values = colors_bins) + theme_DC
save_plot(g2, "./Plots/EasySciATAC_CellTypeFreq.pdf", w = 2.5, h = 3, show.legend = "none")

g2 <- binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(CellType, Group, Bin) %>% summarize(n = n()) %>% group_by(Group, Bin) %>% mutate(freq = n / sum(n)) %>%
  group_by(CellType, Group) %>% mutate(Diff = freq[Bin == 3] - freq[Bin == 1]) %>% ungroup %>%
  group_by(Group) %>% filter(Diff %in% c(min(Diff),max(Diff))) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = CellType, y = freq, fill = Bin)) + geom_col(position = "dodge") + #facet_wrap(~Group, nrow = 1, scales = "free") +
  labs(subtitle = paste("EASY -", group), y = "Cell Type Frequency", x = "Cell Type", fill = NULL) + scale_fill_manual(values = colors_bins) + theme_DC +
  theme(axis.text.x = element_text(angle = -25, hjust = 0))
save_plot(g2, "./Plots/EasySciATAC_CellTypeFreq_v2.pdf", w = 2, h = 2, show.legend = "right")

g2 <- binned %>% filter(Dataset == ds, Group == group) %>%
  group_by(CellType, Group, Bin) %>% summarize(n = n()) %>% group_by(Group, Bin) %>% mutate(total = sum(n), freq = n / total) %>%
  group_by(CellType, Group) %>% mutate(Diff = freq[Bin == 3] - freq[Bin == 1]) %>% ungroup %>%
  group_by(Group) %>% filter(Diff %in% c(min(Diff),max(Diff))) %>% 
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>% ungroup %>%
  {
    res1 <- pairwise.prop.test(as.matrix(select(., n, total))[1:3,]) %>% broom::tidy()
    res2 <- pairwise.prop.test(as.matrix(select(., n, total))[4:6,]) %>% broom::tidy()
    df <- cbind(Bin = as.numeric(.$Bin), CellType = .$CellType, freq = .$freq, rbind(res1,res2), Type = c(1,1,1,-1,-1,-1))
    ggplot(., aes(x = as.numeric(Bin), y = freq)) + 
      geom_col(position = "dodge", mapping = aes(fill = Bin)) +
      labs(subtitle = paste("EASY -", group), y = "Cell Type Frequency", x = NULL, fill = "Bin") + scale_fill_manual(values = colors_bins) + theme_DC +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), strip.clip = 'off') +
      facet_wrap(~gsub(" n","\nn",CellType), switch = "x") + coord_cartesian(clip = 'off') +
      geom_signif(data = df, 
                  mapping = aes(x = as.numeric(Bin), 
                                y = freq, 
                                xmin = as.numeric(group1), 
                                xmax = as.numeric(group2),
                                y_position = 0.25 + ((as.numeric(Bin)-0.5) * 0.03) * Type,
                                annotations = format(p.value, digits = 2)), 
                  textsize = 2,
                  manual = T)
  }
save_plot(g2, "./Plots/EasySciATAC_CellTypeFreq_v2.pdf", w = 3, h = 2, show.legend = "right")
  

## UMAP
g3 <- binned %>% filter(Dataset == ds) %>% filter(Group == group) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Bin)) + 
  ggrastr::rasterise(geom_point(size = 0.5, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
  labs(title = "Tn5 Reaction Bin") + facet_wrap(Group~.) + scale_color_manual(values = colors_bins) + 
  theme_DC + theme(legend.position = 'none', axis.ticks = element_blank(), axis.text = element_blank())
save_plot(g3, "./Plots/EasySciATAC_BinsUMAP.pdf", w = 2.5, h = 2.5, show.legend = "none")

g4 <- binned %>% filter(Dataset == ds) %>%
  {
    ggplot(.,aes(x = UMAP1, y = UMAP2, color = CellType)) + 
      ggrastr::rasterise(geom_point(size = 0.5, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
      ggrepel::geom_label_repel(data = summarize(group_by(., CellType), UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), CellTypeFreq = mean(CellTypeFreq)),
                                mapping = aes(label = CellType, alpha = CellTypeFreq > freq),
                                size = 1.5, label.padding = 0.1, max.overlaps = 20, show.legend = F) +
      # labs(title = "Cell Type") + scale_color_manual(values = paletteer::paletteer_d("ggsci::category20_d3")) +
      labs(title = "Cell Type") + scale_alpha_manual(values = c(0.2, 0.9)) +scale_color_manual(values = colorRampPalette(c("navyblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1", "maroon"))(length(unique(.$CellType)))) +
      theme_DC + theme(legend.position = 'none', axis.ticks = element_blank(), axis.text = element_blank())
  }
save_plot(g4, "./Plots/EasySciATAC_CellTypeUMAP.pdf", w = 2.5, h = 2.5, show.legend = "none")

g <- plot_grid(g1, g3, g2, g4, nrow = 2, rel_widths = c(3,4,3,4))
save_plot(g, "./Plots/EasySciATAC_Summary.pdf", w = 5, h = 5, show.legend = "none")

g <- plot_grid(g2 + theme(axis.text.x = element_text(angle = -25, hjust = 0)), 
               g3, g4, nrow = 1, rel_widths = c(3.5,4,4))
save_plot(g, "./Plots/EasySciATAC_Summary_v2.pdf", w = 6, h = 2, show.legend = "none")




### Binned: Combined --------------------------------------------------------

## Tn5 x nFrags
# ds <- c("txci-ATAC-seq_Barnyard", "dsci-ATAC-seq_BMMC", "scifi-ATAC_Maize", "EasySci-ATAC_MouseBrain")
# group <- c("liver m7 100k", "Baseline B 80k", "M162W", "Aged", "Epithelial","batch_1","A549","MDA-MB-231","K562","LungMap_D122")
# group <- c("liver m7 100k", "Baseline B 80k", "M162W", "Aged", "Epithelial","batch_1","A549","exp1","K562","muscle_SM-C1PWV","Cortex Fresh")
group <- c("K562","Epithelial","K562","Cortex Fresh","exp1","batch_1","Aged","muscle_SM-C1PWV","Baseline B 80k","M162W","liver m7 100k","A549")

g <- rbind(binned,
           tenx %>% filter(Region != "GE") %>% mutate(Group = paste(Region, `Fresh or Frozen?`)) %>%
             select(Tn5 = `sample ID`, Tn5Count = `Cell Ranger Cell Count`, nFrags = `Median Fragments per cell`, Group) %>% 
             mutate(Tn5Type = "Standard", Dataset = "10xATACv1_CorticalDev") %>%
             group_by(Dataset, Group) %>% mutate(Tn5Rank = rank(Tn5Count), Bin = factor(Tn5Rank, levels = sort(unique(Tn5Rank)), labels = ntile(sort(unique(Tn5Rank)), nbin)))) %>%
  filter(Group %in% group) %>%
  group_by(Dataset, Group, Tn5Type, Tn5, Bin) %>% summarize(Tn5Count = mean(Tn5Count), nFrags = median(nFrags)) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  mutate(Order = factor(paste(Dataset, Group, sep="\n"), levels = paste(names(colors_datasets), group, sep="\n")),
         Tn5Type = factor(Tn5Type, levels = c("Standard","Barcoded"))) %>% 
  ggplot(aes(x = log10(Tn5Count), y = log10(nFrags), color = Bin, group = paste(Dataset, Group))) + 
  geom_point(size = 2) + geom_smooth(method = 'lm', se = F, color = "black") + 
  facet_wrap(Tn5Type~Order, nrow = 2, scales = "free") +
  scale_color_manual(values = colors_bins) + labs(title = "Individual Tn5 Reactions", y = "Median nFrags (log10)", x = "Cell Count (log10)") + 
  theme_DC + coord_cartesian(clip="off") + scale_x_continuous(breaks=scales::pretty_breaks(n = 3)) + scale_y_continuous(breaks=scales::pretty_breaks(n = 3))
save_plot(g, "./Plots/Datasets_Bins.pdf", w = 10, h = 5, show.legend = "right")

g <- rbind(binned,
      tenx %>% filter(Region != "GE") %>% mutate(Group = paste(Region, `Fresh or Frozen?`)) %>%
        select(Tn5 = `sample ID`, Tn5Count = `Cell Ranger Cell Count`, nFrags = `Median Fragments per cell`, Group) %>% 
        mutate(Tn5Type = "Standard", Dataset = "10xATACv1_CorticalDev", Short = "10X") %>%
        group_by(Dataset, Group) %>% mutate(Tn5Rank = rank(Tn5Count), Bin = factor(Tn5Rank, levels = sort(unique(Tn5Rank)), labels = ntile(sort(unique(Tn5Rank)), nbin)))) %>%
  filter(Group %in% group) %>%
  group_by(Dataset, Short, Group, Tn5Type, Tn5, Bin) %>% summarize(Tn5Count = mean(Tn5Count), nFrags = median(nFrags)) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  mutate(Order = factor(paste(Dataset, Group, sep="\n"), levels = paste(names(colors_datasets), group, sep="\n")),
         Short = factor(Short, levels = names(colors_short)),
         Tn5Type = factor(Tn5Type, levels = c("Standard","Indexed"))) %>% 
  # ggplot(aes(x = Bin, y = log10(nFrags), color = Bin)) + 
  ggplot(aes(x = log10(Tn5Count), y = log10(nFrags), color = Bin, group = paste(Dataset, Group))) + 
  # geom_boxplot(size = 1) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_point(size = 1, alpha = 0.75) + geom_smooth(method = 'lm', se = F, color = "black") + scale_x_continuous(breaks=scales::pretty_breaks(n = 3)) + 
  facet_wrap(Short~Group, nrow = 2, scales = "free") +
  scale_color_manual(values = colors_bins) + labs(title = "Example Samples or Cell Types", y = "Median nFrags (log10)", x = "Nuclei Count (log10)") + 
  theme_DC +  coord_cartesian(clip="off") + 
  scale_y_continuous(breaks=scales::pretty_breaks(n = 3), expand = expansion(mult = c(.05, 0.15))) +
  stat_cor(aes(group = Dataset), size = 1.75, p.accuracy = 0.05, show.legend = F, label.y = Inf, vjust = 1, label.x.npc = "center", hjust = 0.5)
  # stat_compare_means(method = 'anova', label = 'p.signif')
save_plot(g, "./Plots/Datasets_Bins_v2.pdf", w = 5.5, h = 3, show.legend = "right")



### Scatterpie UMAPs --------------------------------------------------------

binned %>%
  filter(Dataset == "EasySci-ATAC_MouseBrain") %>% filter(Group == "Young") %>%
  # filter(Dataset == "scifi-ATAC_Maize") %>% filter(Group == "M162W") %>%
  # filter(Dataset == "txci-ATAC-seq_Barnyard") %>% filter(Group == "liver m7 200k") %>%
  # filter(Dataset == "dsci-ATAC-seq_BMMC") %>% filter(Group == "Baseline B 40k", CellType != "Collision") %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  group_by(CellType, Group) %>%
  {
    ggplot(., aes(x = UMAP1, y = UMAP2, color = CellType)) +
      ggrastr::rasterise(geom_point(size = 0.4, alpha = 0.5, stroke = 0, show.legend = F), dpi = 300, scale = 1) +
      ggrepel::geom_label_repel(data = summarize(group_by(., Group, CellType), UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)),
                                mapping = aes(label = CellType),
                                size = 3, alpha = 0.7, max.overlaps = 5, nudge_x = 1, show.legend = F) +
      scatterpie::geom_scatterpie(data = mutate(group_by(pivot_wider(summarize(group_by(mutate(., UMAP1 = median(UMAP1) - 0.5, UMAP2 = median(UMAP2)),
                                                                                        Group, CellType, Bin, UMAP1, UMAP2),
                                                                               Count = n()),
                                                                     names_from = "Bin", values_from = "Count", values_fill = 0), Group),
                                                # `Bottom 33%` = `Bottom 33%`/max(`Bottom 33%`), `Middle 33%` = `Middle 33%`/max(`Middle 33%`), `Top 33%` = `Top 33%`/max(`Top 33%`)),
                                                `Bottom 33%` = `Bottom 33%`/sum(`Bottom 33%`), `Middle 33%` = `Middle 33%`/sum(`Middle 33%`), `Top 33%` = `Top 33%`/sum(`Top 33%`)),
                                  mapping = aes(x = UMAP1, y = UMAP2, group = CellType),
                                  cols = c("Bottom 33%","Middle 33%","Top 33%"), pie_scale = 1.5, alpha = 0.8) +
      facet_wrap(~Group) +
      labs(title = "Binned by Tn5Count") + theme_DC +
      coord_equal() +
      scale_fill_manual(values = colors_bins)
  }




### Significantly Correlated LSI Dims ---------------------------------------

## Plot number of significanly correlated dimensions per dataset
threshCor <- 0.5
threshP <- 0.05

list <- list(LSICor(proj_plex2, metric = "nFrags", groups = c("treatment"), dims = 1:30, absValue = T, threshCor = threshCor, threshP = threshP),
             LSICor(proj_2, metric = "nFrags", groups = c("CellType"), dims = 1:30, absValue = T, threshCor = threshCor, threshP = threshP),
             LSICor(proj_tx3, metric = "nFrags", groups = c("CellType"), dims = 1:30, absValue = T, threshCor = threshCor, threshP = threshP),
             LSICor(proj_snu2, metric = "nFrags", groups = c("CellType"), dims = 1:30, absValue = T, threshCor = threshCor, threshP = threshP),
             LSICor(proj_breast2, metric = "nFrags", groups = c("CellType"), dims = 1:30, absValue = T, threshCor = threshCor, threshP = threshP))

list <- lapply(1:length(list), function(x) { 
  tmp <- list[[x]]
  colnames(tmp)[2] <- "group" 
  tmp$Dataset <- c("sciPlex-ATAC-seq2_DrugScreen", "dsci-ATAC-seq_BMMC", "txci-ATAC-seq_Barnyard", "SNuBar_CellLines", "SNuBar_HBCA")[x]
  tmp$GroupID <- factor(tmp$group, levels = unique(tmp$group), labels = 1:length(unique(tmp$group)))
  tmp
})

df <- do.call(rbind, list)

df %>% group_by(Dataset, GroupID) %>% 
  # summarize(nSig = sum(Cor > 0.4 & Cor.p < 0.01)) %>%
  summarize(nSig = sum(Cor > threshCor & Cor.p < threshP)) %>%
  ggplot(., aes(x = Dataset, y = nSig, color = Dataset)) + 
  geom_boxplot(outlier.colour = NA) + geom_jitter(width = 0.1, height = 0.1) +
  theme_classic() + 
  # theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  theme(axis.text.x = element_blank()) + 
  scale_color_manual(valu es = colors_datasets) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  labs(y = "Number of Significantly Correlated Dimensions")

g <- df %>% group_by(Dataset) %>% 
  filter(Cor > threshCor & Cor.p < threshP) %>% 
  summarize(nSig = length(unique(Dim))) %>%
  mutate(Dataset = factor(Dataset, levels = names(colors_datasets), labels = names(colors_short))) %>%
  ggplot(., aes(y = factor(Dataset, levels = Dataset[order(nSig)]), x = nSig, fill = Dataset)) + 
  geom_col() +
  scale_fill_manual(values = colors_short) + 
  ggtitle("LSI vs. Depth Correlation ") +
  labs(x = "Significantly Correlated\nLSI Components", y = "Dataset") + theme_DC
save_plot(g, "Plots/LSI_CorrelatedDims.pdf", w = 2, h = 2, show.legend="none")

df %>% group_by(Dataset) %>% 
  filter(Cor > threshCor & Cor.p < threshP) %>% 
  summarize(nSig = length(unique(Dim)),
            nGroups = length(unique(group))) %>%
  ggplot(., aes(x = nSig, y = nGroups, color = Dataset)) + 
  geom_point() +
  scale_color_manual(values = colors_datasets) + 
  ggtitle("LSI Dimensions Significantly Correlated w/ nFrags") +
  labs(x = "Significantly Correlated LSI Components", y = "Dataset") + theme_DC


# dsci-ATAC-seq_BMMC
LSICor(proj_2, metric = "nFrags", groups = c("CellType"), plotBiPlot = T, absValue = T, threshCor = threshCor, threshP = threshP)

g <- cbind(proj_2@reducedDims$IterativeLSI$matSVD, proj_2@cellColData) %>% as.data.frame %>%
  pivot_longer(cols = c(paste("LSI",1:30,sep="")), names_to = "Dim", values_to = "Value") %>%
  filter(Dim %in% c("LSI1","LSI2","LSI3")) %>%
  ggplot(aes(x = log10(nFrags), y = Value, color = CellType)) +
  ggrastr::rasterise(geom_point(size = 0.5, alpha = 0.5, stroke = 0), dpi = 600, scale = 1) +
  # geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
  paletteer::scale_color_paletteer_d("ggsci::category20_d3") + ggtitle("DSCI") + 
  facet_wrap(~Dim, scales = 'free', nrow = 1) + theme_DC + guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1)))
save_plot(g, "Plots/LSI_ExampleDims_DSCI.pdf", w = 4.5, h = 1.5, show.legend="right")

# SNuBar_CellLines
LSICor(proj_snu2, metric = "nFrags", groups = c("CellType"), plotBiPlot = T, absValue = T, threshCor = threshCor, threshP = threshP)

g <- cbind(proj_snu2@reducedDims$IterativeLSI$matSVD, proj_snu2@cellColData) %>% as.data.frame %>%
  pivot_longer(cols = c(paste("LSI",1:30,sep="")), names_to = "Dim", values_to = "Value") %>%
  filter(Dim %in% c("LSI1","LSI2","LSI3")) %>%
  ggplot(aes(x = log10(nFrags), y = Value, color = CellType)) + 
  ggrastr::rasterise(geom_point(size = 0.5, alpha = 0.5, stroke = 0), dpi = 600, scale = 1) +
  # geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
  scale_color_brewer(palette = 'Set1') + ggtitle("SNU_A") + 
  facet_wrap(~Dim, scales = 'free', nrow = 1) + theme_DC + guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
save_plot(g, "Plots/LSI_ExampleDims_SNU_A.pdf", w = 4, h = 1.5, show.legend="right")




### Binned: UMAP Batch Effects ----------------------------------------------

binned %>% filter(Dataset == ds) %>% filter(Group == group) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  # mutate(Cluster = factoextra::hkmeans(data.frame(UMAP1,UMAP2), k = 5)$cluster) %>%
  # filter(UMAP1 > coord[[1]][1] & UMAP1 < coord[[1]][2]) %>%
  # filter(UMAP2 > coord[[2]][1] & UMAP2 < coord[[2]][2]) %>%
  # filter(CellType %in% celltypes) %>%
  mutate(Cluster = CellType) %>%
  # ggplot(aes(x = UMAP1, y = UMAP2, color = as.character(Cluster))) +
  ggplot(aes(x = UMAP1, y = UMAP2, color = Bin, fill = Bin)) +
  # ggplot(aes(x = UMAP1, y = UMAP2, color = log10(nFrags), fill = Bin)) +
  ggrastr::rasterise(geom_point(size = 1, alpha = 1, stroke = 0), dpi = 600, scale = 1) +
  # ggrastr::rasterise(geom_point(size = 1, alpha = 0.5, color = 'grey80', stroke = 0), dpi = 600, scale = 1) + geom_density_2d(bins = 15, size = 1) + facet_grid(~Bin) +
  # geom_path(data = . %>% group_by(Cluster,Bin) %>% summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), mapping = aes(group = Cluster), color = 'black') +
  # geom_point(data = . %>% group_by(Cluster,Bin) %>% summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), mapping = aes(shape = as.character(Cluster)), size = 4, color = 'black', stroke = 1) +
  geom_path(data = . %>% group_by(Cluster,Bin) %>% summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), mapping = aes(group = Cluster), color = 'black') +
  geom_point(data = . %>% group_by(Cluster,Bin) %>% summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), shape = 21, size = 2, color = 'black', stroke = 1) +
  # geom_path(data = . %>% group_by(Bin) %>% summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), color = 'black') +
  # geom_point(data = . %>% group_by(Bin) %>% summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), size = 4, color = 'black', stroke = 1) +
  labs(title = "Tn5 Reaction Bin") + 
  scale_color_manual(values = colors_bins) + scale_fill_manual(values = colors_bins) + #scale_shape_manual(values = c(21,22,24)) +
  theme_DC + theme_UMAP

binned %>% filter(Dataset == ds) %>% filter(Group == group) %>%
  mutate(Bin = factor(Bin, levels = 1:3, labels = c("Bottom 33%","Middle 33%","Top 33%"))) %>%
  # mutate(Cluster = factoextra::hkmeans(data.frame(UMAP1,UMAP2), k = 5)$cluster) %>%
  # filter(UMAP1 > coord[[1]][1] & UMAP1 < coord[[1]][2]) %>%
  # filter(UMAP2 > coord[[2]][1] & UMAP2 < coord[[2]][2]) %>%
  # filter(CellType %in% celltypes) %>%
  mutate(Cluster = CellType) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = as.character(Cluster))) +
  # ggplot(aes(x = UMAP1, y = UMAP2, color = Cluster, fill = Bin)) +
  # ggrastr::rasterise(geom_point(size = 1, alpha = 1, stroke = 0, color = 'grey80'),  dpi = 600, scale = 1) +
  ggrastr::rasterise(geom_point(size = 1, alpha = 0.5, stroke = 0), dpi = 600, scale = 1) +
  geom_label(. %>% group_by(Cluster) %>% summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), mapping = aes(label = Cluster), size = 4, alpha = 0.75) + 
  # geom_path(data = . %>% group_by(Cluster,Bin) %>% summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), mapping = aes(group = Cluster), color = 'black') +
  # geom_point(data = . %>% group_by(Cluster,Bin,CellType) %>% summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), mapping = aes(shape = CellType), size = 4, color = 'black', stroke = 1) +
  # geom_path(data = . %>% group_by(Bin) %>% summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), color = 'black') +
  # geom_point(data = . %>% group_by(Bin) %>% summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), size = 3) +
  labs(title = "Tn5 Reaction Bin") +
  # scale_fill_manual(values = colors_bins) + scale_shape_manual(values = c(21,22,24)) +
  theme_DC + theme_UMAP + theme(legend.position = 'none')


# binned %>% filter(Dataset == ds) %>% filter(Group == group) %>% .$CellType %>% unique %>% length

ds <- "txci-ATAC-seq_Barnyard"
group <- "liver m7 100k" # 7 celltypes
coord <- list(c(-3,10),c(-5,6))
celltypes <- c("Hep","HPC/Cho","Lym")
# group <- "lung m1 100k" # 24 celltypes

ds <- "EasySci-ATAC_MouseBrain" # 31 celltypes
# group <- "Young"
group <- "Aged"
coord <- list(c(-3,10),c(-5,6))
celltypes <- c("Hep","HPC/Cho","Lym")
# group <- "WT"
# group <- "5xFAD"

ds <- "scifi-ATAC_Maize" # 11 celltypes
group <- "M162W" 
# group <- "B97"
# group <- "B73"
# group <- "Ky21"

ds <- "dsci-ATAC-seq_BMMC"
# group <- "Stimulated A 20k"
# group <- "Stimulated A 40k"
# group <- "Baseline B 20k"
# group <- "Baseline B 40k"
group <- "Baseline B 80k"



### Binned: JSD Analysis ----------------------------------------------------

binned %>% filter(Dataset == ds) %>% filter(Group == group) %>% group_by(Bin) %>%
  ungroup %>% select(UMAP1, UMAP2, Group = Bin) %>%
  # ungroup %>% select(UMAP1, UMAP2, Group = Random) %>%
  {
    x <- range(.$UMAP1)
    y <- range(.$UMAP2)
    tmp <- MASS::kde2d(filter(., Group == 1)$UMAP1, filter(., Group == 1)$UMAP2, n = 500, lims = c(x, y))$z %>% as.vector
    tmp2 <- MASS::kde2d(filter(., Group == 3)$UMAP1, filter(., Group == 3)$UMAP2, n = 500, lims = c(x, y))$z %>% as.vector
    
    philentropy::JSD(rbind(tmp, tmp2))
  }



jsd <- lapply(unique(binned %>% filter(!is.na(UMAP1)) %>% .$Dataset), function(d) {
  binned %>% filter(Dataset == d) %>% group_by(Dataset, Group) %>%
    summarize(BinJSD = 
                {
                  x <- range(UMAP1)
                  y <- range(UMAP2)
                  a <- filter(., Bin == 1)
                  b <- filter(., Bin == 3)
                  tmp1 <- MASS::kde2d(a$UMAP1, a$UMAP2, n = 500, lims = c(x, y))$z %>% as.vector
                  tmp2 <- MASS::kde2d(b$UMAP1, b$UMAP2, n = 500, lims = c(x, y))$z %>% as.vector
                  philentropy::JSD(rbind(tmp1, tmp2))
                },
              RandomJSD = 
                {
                  x <- range(UMAP1)
                  y <- range(UMAP2)
                  a <- filter(., Random == 1)
                  b <- filter(., Random == 3)
                  tmp1 <- MASS::kde2d(a$UMAP1, a$UMAP2, n = 500, lims = c(x, y))$z %>% as.vector
                  tmp2 <- MASS::kde2d(b$UMAP1, b$UMAP2, n = 500, lims = c(x, y))$z %>% as.vector
                  philentropy::JSD(rbind(tmp1, tmp2))
                }
    )
})

jsd2 <- binned %>% filter(!is.na(UMAP1)) %>%
  # filter(Dataset == "EasySci-ATAC_MouseBrain", Group %in% c("5xFAD","Aged")) %>% 
  group_by(Dataset, Group) %>%
  summarize(BinJSD = philentropy::JSD(rbind(as.vector(MASS::kde2d(UMAP1[which(Bin == 1)], UMAP2[which(Bin == 1)], n = 500, lims = c(range(UMAP1), range(UMAP2)))$z),
                                            as.vector(MASS::kde2d(UMAP1[which(Bin == 3)], UMAP2[which(Bin == 3)], n = 500, lims = c(range(UMAP1), range(UMAP2)))$z))),
            RandomJSD = philentropy::JSD(rbind(as.vector(MASS::kde2d(UMAP1[which(Random == 1)], UMAP2[which(Random == 1)], n = 500, lims = c(range(UMAP1), range(UMAP2)))$z),
                                               as.vector(MASS::kde2d(UMAP1[which(Random == 3)], UMAP2[which(Random == 3)], n = 500, lims = c(range(UMAP1), range(UMAP2)))$z))))



# jsd %>% do.call(rbind,.) %>%
jsd2 %>%
  pivot_longer(cols = 3:4,names_to = "Grouping", values_to = "JSD") %>%
  mutate(Grouping = gsub("JSD","",Grouping)) %>%
  # ggplot(aes(x = Dataset, y = 1/JSD, fill = Dataset, linetype = Grouping)) + geom_boxplot() + 
  ggplot(aes(x = Dataset, y = JSD, fill = Dataset, linetype = Grouping)) + 
  # geom_boxplot() + 
  geom_jitter() +
  scale_fill_manual(values = colors_datasets) + theme_DC + 
  stat_compare_means(method = "t.test", paired = T, ref.group = ".all.", label = "p.signif", size = 2, vjust = 0)

jsd2 %>%
  pivot_longer(cols = 3:4,names_to = "Grouping", values_to = "JSD") %>%
  mutate(Grouping = gsub("JSD","",Grouping)) %>%
  # ggplot(aes(x = Dataset, y = 1/JSD, fill = Dataset, linetype = Grouping)) + geom_boxplot() + 
  ggplot(aes(x = Grouping, y = JSD, color = Dataset, linetype = Grouping)) + 
  # geom_boxplot() + 
  geom_jitter(width = 0.2) +
  facet_wrap(~Dataset, nrow = 1) + 
  scale_color_manual(values = colors_datasets) + theme_DC + 
  stat_compare_means(method = "t.test", paired = T, ref.group = "Random", label = "p.signif", size = 2, label.x = 1.5, mapping = aes(group = Grouping))

# jsd %>% do.call(rbind,.) %>%
jsd2 %>%
  # group_by(Dataset) %>% mutate(BinJSD = (1/BinJSD) / mean(1/RandomJSD),
  #                              RandomJSD = (1/RandomJSD) / mean(1/RandomJSD)) %>%
  group_by(Dataset) %>% mutate(BinJSD = BinJSD / mean(RandomJSD),
                               RandomJSD = RandomJSD / mean(RandomJSD)) %>%
  pivot_longer(cols = 3:4,names_to = "Grouping", values_to = "JSD") %>%
  mutate(Grouping = gsub("JSD","",Grouping)) %>%
  ggplot(aes(x = Dataset, y = JSD, fill = Dataset, linetype = Grouping)) + geom_boxplot() + 
  scale_fill_manual(values = colors_datasets) + theme_DC

# jsd %>% do.call(rbind,.) %>%
jsd2 %>%
  group_by(Dataset) %>% summarize(SD = sd(RandomJSD/BinJSD), Mean = mean(RandomJSD/BinJSD)) %>%
  ggplot(aes(x = Dataset, y = Mean, fill = Dataset)) + 
  geom_col() + geom_errorbar(mapping = aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.5) + 
  scale_fill_manual(values = colors_datasets) + theme_DC




### Export R Objects --------------------------------------------------------

dir.create("Export")

sessionOut("Export/PublishedDatasets_Packages.tsv")

## Master dataset (single-cell entries for 11 datasets)
saveRDS(binned, file = "Export/PublishedDatasets_Master.rds")

## Transposition reaction summary (Tn5 reaction entries for 12 datasets)
saveRDS(comb2, file = "Export/PublishedDatasets_Tn5Summary.rds")

## Cell Type Frequency Fold-Change Between Bins
saveRDS(celltypefreq, file = "Export/PublishedDatasets_CellTypeFoldChange.rds")


