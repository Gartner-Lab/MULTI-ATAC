### Custom Functions for MULTI-ATAC Analysis

sessionOut <- function(file) {
  r <- sessionInfo()$R.version$version.string
  tmp <- lapply(seq_along(sessionInfo()[["otherPkgs"]]), function(x) {
    paste(names(sessionInfo()[["otherPkgs"]])[x], 
          sessionInfo()[["otherPkgs"]][[x]]$Version,
          sep = "_")
  }) %>% unlist %>% sort 
  tmp <- c(r, tmp) %>% paste(collapse = "\n")
  write.table(tmp, file = file, quote = F, sep = '\t', col.names = F, row.names = F)
}
# sessionOut("~/Desktop/session.tsv")

# Aesthetics --------------------------------------------------------------

theme_DC <- theme_classic(base_line_size=0.5/.pt, base_size = 7) +
  theme(legend.text = element_text(colour="black", size=5), 
        legend.key.height = unit(3,"mm"),
        legend.key.width = unit(3,"mm"),
        strip.placement = "outside", 
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

theme_UMAP <- theme(axis.ticks = element_blank(), 
                    axis.text = element_blank())


centered_gradient <- function(pal = "ggthemes::Classic Orange-White-Blue", midquant = 0.5, direction = 1) {
  cols <- paletteer::paletteer_c(palette = pal, direction = direction, 255)
  na <- cols[128]
  if (midquant > 0.5) {
    cols <- cols[c(1:127,rep(128,3),head(129:255, n = ceiling(127 * midquant/(1-midquant))))]
  } else if (midquant < 0.5) {
    cols <- cols[c(tail(1:127, n = floor(127 * midquant/(1-midquant))),rep(128,3),129:255)] 
  }
  ggplot2::scale_fill_gradientn(colours = cols, na.value = na)
}

centerQuant <- function(x, midpoint = 0) {
  x <- x[!is.na(x)]
  (midpoint - min(x)) / (max(x) - min(x))
}


# Data Wrangling ----------------------------------------------------------

colZscores <- function (m = NULL) { z <- sweep(t(t(m) - colMeans(m)), 2, matrixStats::colSds(m), `/`) }

smooth_it <- function(x, y, n = 100, method = "natural") {
  t <- seq_along(x)
  new_t <- seq(min(t), max(t), length.out = n)
  new_x <- spline(t, x, xout = new_t, method = method)$y
  new_y <- spline(t, y, xout = new_t, method = method)$y
  data.frame(t = new_t, x = new_x, y = new_y)
}

dupControl <- function(df, controls = c("DMSO(-)","DMSO(+)","NegCtrl","PosCtrl","Resting Control","Stimulated Control"), drugs = drug[3:8]) {
  lapply(drugs, function(d) { df %>% filter(Drug %in% c(controls,d)) %>% mutate(Drug = d) }) %>% do.call(rbind,.) %>%
    mutate(Drug = factor(Drug, levels = levels(df$Drug)))
}

# old <- getGeneAnnotation(proj_4)$gene$symbol %>% unname
# old <- old[!is.na(old)]
# new <- GeneSymbolThesarus(old)
# new <- background_GeneSymbolThesaurus_results$new
# saveRDS(new, "./newGeneSymbols.rds")
fixGene <- function(old = NULL, new = new, named = F) {
  new2 <- new %>% subset(. %ni% getGeneAnnotation(proj_4)$gene$symbol)
  new2 <- new2 %>% subset(. %ni% grep("^MT-",.,value=T))
  new2 <- new2[names(new2) %in% old]
  if (named) { names(old) <- old }
  idx <- match(names(new2),old)
  old[idx] <- new2
  old
}

se <- function(proj, useMatrix = "PeakMatrix", groupBy = "CellDrugDose") {
  se <- getGroupSE(proj, useMatrix = useMatrix, groupBy = groupBy, scaleTo = NULL, divideN = T)
  mtx <- assay(se)
  if (useMatrix == "PeakMatrix") {
    rownames(mtx) <- rowData(se) %>% {paste(.$seqnames, .$start, .$end, sep="-")}
  } else {
    rownames(mtx) <- rowData(se)$name %>% fixGene(new)
  }
  mtx
}

# Analysis ----------------------------------------------------------------

corTest <- function(a, b, value = c("p","coef","pdir","comb"), test = c("cor","lm")) {
  if (test == "lm") {
    res <- summary(lm(b~a))
    p <- res$coefficients %>% {.[nrow(.),4]}
    coef <- res$coefficients %>% {.[nrow(.),1]}
    dir <- sign(coef)
  } else if (test == "cor") {
    res <- cor.test(a, b, alternative = "t", method = "pearson", conf.level = 0.95)
    p <- res$p.value
    coef <- res$estimate
    dir <- sign(coef)
  }
  # if (value == "p") { -log10(p) }
  if (value == "p") { p }
  else if (value == "coef") { coef }
  else if (value == "pdir") { -log10(p) * dir }
  else if (value == "comb") { -log10(p) * coef }
}

repTest <- function(assay = "RNA", bins = 50, repsumm = c("mean","median"), test = c("cor","lm"), sample = F, log = T) {
  
  if (assay == "ATAC") { 
    se <- getMatrixFromProject(proj_4, useMatrix = "PeakMatrix")
    rownames(se@assays@data$PeakMatrix) <- se@rowRanges %>% {paste(seqnames(.), start(.), end(.), sep = "-")}
    colnames(se@assays@data$PeakMatrix) <- se@colData %>% rownames
    mtx <- se@assays@data$PeakMatrix
  }
  if (assay == "RNA") { 
    mtx <- seu@assays$RNA$counts 
    colnames(mtx) <- gsub("_","#",colnames(mtx))
  }
  if (assay == "SCT") { 
    mtx <- seu@assays$SCT$counts 
    colnames(mtx) <- gsub("_","#",colnames(mtx))
  }
  if (assay == "Score") { 
    se <- getMatrixFromProject(proj_4, useMatrix = "GeneScoreMatrix")
    rownames(se@assays@data$GeneScoreMatrix) <- se@elementMetadata$name
    colnames(se@assays@data$GeneScoreMatrix) <- se@colData %>% rownames
    mtx <- se@assays@data$GeneScoreMatrix
  }
  if (assay == "Dev1") {
    se <- getMatrixFromProject(proj_4, useMatrix = "MotifMatrix", useSeqnames = "z")
    mtx <- se@assays@data$z
  }
  if (assay == "Dev2") {
    se <- getMatrixFromProject(proj_4, useMatrix = "MotifMatrix2", useSeqnames = "z")
    mtx <- se@assays@data$z
  }
  features <- rownames(mtx)
  features <- features[grep("^MT-", features, invert = T)]
  if (sample) { features <- head(features, n = 50); bins <- 2 }
  n <- 0
  
  repL <- summL <- summL2 <- log2fcL <- log2fcL2 <- pL <- coefL <- list()
  
  lapply(split(features, cut(seq_along(features), bins, labels = F)), function(g) {
    n <<- n + 1
    cat(n,"/",bins," bins\n",sep = "")
    g <- unlist(g)
    
    ## replicate mean values (average of cells)
    r <- cbind(proj_4@cellColData[colnames(mtx),c("CellType_Major","Drug","Dose","Replicate")],
               mtx[g,] %>% t) %>% as.data.frame %>%
      group_by(CellType_Major, Drug, Dose, Replicate) %>% summarize_all(mean)
    if (log) { r <- r %>% mutate_at(.vars = colnames(r)[-1:-4], .funs = "log1p") }
    
    ## summarize replicates
    if (repsumm == "mean") {
      m <- r %>% group_by(CellType_Major, Drug, Dose) %>% select(-Replicate) %>% summarize_all(mean)
      m.tmp <- m %>% dupControl(controls = "DMSO(+)", drugs = drug[c(1,3:8)]) %>% filter(Drug %ni% c("DMSO(+)")) %>% group_by(CellType_Major, Drug)
      ## max mean per drug
      m2 <- m.tmp %>% select(-Dose) %>% summarize_all(.funs = list(~ max(abs(.x))))
    } else if (repsumm == "median") { # less sensitive to single outliers
      m <- r %>% group_by(CellType_Major, Drug, Dose) %>% select(-Replicate) %>% summarize_all(median)
      m.tmp <- m %>% dupControl(controls = "DMSO(+)", drugs = drug[c(1,3:8)]) %>% filter(Drug %ni% c("DMSO(+)")) %>% group_by(CellType_Major, Drug)
      ## max median per drug
      m2 <- m.tmp %>% select(-Dose) %>% summarize_all(.funs = list(~ max(abs(.x))))
    }

    if (assay %in% c("Dev1","Dev2")) { l <- l2 <- data.frame() } else {
      ## log2fc compared to positive control
      l <- m.tmp %>% 
        mutate_at(.vars = colnames(m)[-1:-3], .funs = list(~ case_when(.x == 0 ~ 0, TRUE ~ log2((.x + 0.01)/(.x[Dose == "DMSO(+)"] + 0.01)))), .groups = "keep") %>%
        filter(!(Dose == "DMSO(+)" & Drug != "DMSO(-)")) %>% arrange(CellType_Major, Drug, Dose)
      ## max log2fc per drug
      l2 <- l %>% group_by(CellType_Major, Drug) %>% select(-Dose) %>% summarize_all(.funs = list(~ max(abs(.x))))
    }
    
    r.tmp <- r %>% 
      dupControl(controls = "DMSO(+)", drugs = drug[c(1,3:8)]) %>% filter(Drug %ni% c("DMSO(+)")) %>% mutate(Dose = as.numeric(Dose)) %>% group_by(CellType_Major, Drug) 
    ## regression/correlation p.value
    p <- r.tmp %>% summarize_at(.vars = colnames(r.tmp)[-1:-4],
                                .funs = list(~ case_when(mean(abs(.x)) == 0 ~ 1,
                                                         TRUE ~ suppressWarnings(corTest(Dose, .x, test = test, value = "p")))),
                                .groups = "keep")
    ## slope or r-value
    c <- r.tmp %>% summarize_at(.vars = colnames(r.tmp)[-1:-4],
                                .funs = list(~ case_when(mean(abs(.x)) == 0 ~ 0,
                                                         TRUE ~ suppressWarnings(corTest(Dose, .x, test = test, value = "coef")))),
                                .groups = "keep")
    
    repL <<- c(repL, list(r))
    summL <<- c(summL, list(m))
    summL2 <<- c(summL2, list(m2))
    log2fcL <<- c(log2fcL, list(l))
    log2fcL2 <<- c(log2fcL2, list(l2))
    pL <<- c(pL, list(p))
    coefL <<- c(coefL, list(c))
  })
  
  rep <- Reduce(join, repL)
  summ <- Reduce(join, summL)
  summ2 <- Reduce(join, summL2)
  log2fc <- Reduce(join, log2fcL)
  log2fc2 <- Reduce(join, log2fcL2)
  p <- Reduce(join, pL)
  coef <- Reduce(join, coefL)
  
  if (ncol(log2fc) > 0) {
    colnames(log2fc) <- c("CellType","Drug","Dose",features)
    colnames(log2fc2) <- c("CellType","Drug",features)
  }
  colnames(p) <- c("CellType","Drug",features)
  colnames(coef) <- c("CellType","Drug",features)
  colnames(summ) <- c("CellType","Drug","Dose",features)
  colnames(summ2) <- c("CellType","Drug",features)
  colnames(rep) <- c("CellType","Drug","Dose","Replicate",features)
  
  out <- list(p = p, coef = coef, log2fc = log2fc, log2fc.Max = log2fc2, summ = summ, summ.Max = summ2, rep = rep, summstat = repsumm)
  # names(out)[5:6] <- c(repsumm, paste(repsumm,"Max",sep="."))
  
  return(out)
}

dimReduce <- function(proj, what = c("iLSI","UMAP","Clustering"), scale = c("Row","Col"), default = "Col", dims = 2:30) {
  if (default %ni% scale) { default <- "NoScale" }
  if ("iLSI" %in% what) {
    proj <- addIterativeLSI(
      ArchRProj = proj,
      name = "IterativeLSI_NoScale",
      saveIterations = FALSE,
      corCutOff = 1, scaleDims = F,
      force = TRUE
    )
    if ("Row" %in% scale) { 
      proj@reducedDims$IterativeLSI_Row <- proj@reducedDims$IterativeLSI_NoScale
      proj@reducedDims$IterativeLSI_Row$matSVD <- ArchR:::.rowZscores(proj@reducedDims$IterativeLSI_Row$matSVD)
    }
    if ("Col" %in% scale) {
      proj@reducedDims$IterativeLSI_Col <- proj@reducedDims$IterativeLSI_NoScale
      proj@reducedDims$IterativeLSI_Col$matSVD <- colZscores(proj@reducedDims$IterativeLSI_Col$matSVD)
    }
    proj@reducedDims$IterativeLSI <- proj@reducedDims[[grep(paste("IterativeLSI",default,sep="_"), names(proj@reducedDims))]]
  }
  if ("UMAP" %in% what) {
    proj <- addUMAP(proj, reducedDims = "IterativeLSI_NoScale", name = "UMAP_NoScale", dimsToUse = dims, corCutOff = 1, scale = F, force = TRUE)
    if ("Row" %in% scale) {
      proj <- addUMAP(proj, reducedDims = "IterativeLSI_Row", name = "UMAP_Row", dimsToUse = dims, corCutOff = 1, scale = F, force = TRUE)
    }
    if ("Col" %in% scale) {
      proj <- addUMAP(proj, reducedDims = "IterativeLSI_Col", name = "UMAP_Col", dimsToUse = dims, corCutOff = 1, scale = F, force = TRUE)
    }
    proj@embeddings$UMAP <- proj@embeddings[[grep(paste("UMAP",default,sep="_"), names(proj@embeddings))]]
  }
  if ("Clustering" %in% what) {
    proj <- addClusters(proj, reducedDims = "IterativeLSI_NoScale", name = "Clusters_NoScale", dimsToUse = dims, corCutOff = 1, scale = F, force = TRUE)
    if ("Row" %in% scale) {
      proj <- addClusters(proj, reducedDims = "IterativeLSI_Row", name = "Clusters_Row", dimsToUse = dims, corCutOff = 1, scale = F, force = TRUE)
    }
    if ("Col" %in% scale) {
      proj <- addClusters(proj, reducedDims = "IterativeLSI_Col", name = "Clusters_Col", dimsToUse = dims, corCutOff = 1, scale = F, force = TRUE)
    }
    proj$Clusters <- proj@cellColData[[grep(paste("Clusters",default,sep="_"), names(proj@cellColData))]]
  }
  proj
}

dimDefault <- function(proj, what = c("iLSI","UMAP","Clustering"), default = "NoScale") {
  if ("iLSI" %in% what) {
    proj@reducedDims$IterativeLSI <- proj@reducedDims[[grep(paste("IterativeLSI",default,sep="_"), names(proj@reducedDims))]]
  }
  if ("UMAP" %in% what) {
    proj@embeddings$UMAP <- proj@embeddings[[grep(paste("UMAP",default,sep="_"), names(proj@embeddings))]]
  }
  if ("Clustering" %in% what) {
    proj$Clusters <- proj@cellColData[[grep(paste("Clusters",default,sep="_"), names(proj@cellColData))]]
  }
  proj
}

LSICor <- function(ArchRProject = NULL,
                   SeuratObject = NULL,
                   svdMatrix = NULL,
                   metaDF = NULL,
                   reducedDim = "IterativeLSI",
                   metric = c("nFrags"),
                   groups = c("Sample"),
                   dims = 1:30,
                   method = "pearson",
                   rankCor = F,
                   rankSlope = F,
                   absValue = F,
                   plotVolcano = F,
                   plotBiPlot = F,
                   threshCor = 0.4,
                   threshP = 0.01
) {
  require(dplyr)
  if (!is.null(ArchRProject)) {
    SVD <- ArchRProject@reducedDims[[reducedDim]]$matSVD
    meta <- ArchRProject@cellColData[,c(metric,groups)]
  } else if (!is.null(SeuratObject)) {
    SVD <- SeuratObject@reductions[[reducedDim]]@cell.embeddings
    meta <- SeuratObject@meta.data[,c(metric,groups)]
  } else {
    SVD <- svdMatrix
    meta <- metaDF
  }
  
  df <- lapply(dims, function(x) {
    tmp <- data.frame(Dim = x, SVD = SVD[,x], meta)
    df <- tmp %>% group_by(.dots = c("Dim",groups)) %>%
      summarize(Cor.p = cor.test(SVD, get(metric), alternative = "t", method = method, conf.level = 0.95)$p.value,
                Cor = cor.test(SVD, get(metric), alternative = "t", method = method, conf.level = 0.95)$estimate,
                Slope = (lm(SVD~get(metric)) %>% coef())[2], 
                .groups = "keep")
  }) %>% do.call(rbind,.)
  
  if (absValue) { 
    df$Cor <- abs(df$Cor) 
    df$Slope <- abs(df$Slope)
  }
  
  if (rankCor) {
    max <- df %>% group_by(Dim) %>% summarize(max = max(abs(Cor)))
    # df$Dim <- factor(df$Dim, levels = dims[order(df$Cor, decreasing = T)])
    df$Dim <- factor(df$Dim, levels = dims[order(max$max, decreasing = T)])
  } else if (rankSlope) {
    max <- df %>% group_by(Dim) %>% summarize(max = max(abs(Slope)))
    # df$Dim <- factor(df$Dim, levels = dims[order(df$Slope, decreasing = T)])
    df$Dim <- factor(df$Dim, levels = dims[order(max$max, decreasing = T)])
  } else { 
    df$Dim <- factor(df$Dim, levels = sort(dims, decreasing = F)) 
  }
  
  if (plotVolcano) {
    g <- ggplot(df, aes(x = Cor, y = -log10(Cor.p), label = Dim)) + 
      geom_point(data = df[abs(df$Cor) < threshCor | df$Cor.p > threshP,]) + 
      geom_text(data = df[abs(df$Cor) > threshCor & df$Cor.p < threshP,], vjust = "inward") +
      geom_vline(xintercept = c(-threshCor,threshCor), color = "blue", linetype = 'dashed') +
      geom_hline(yintercept = -log10(threshP), color = "red", linetype = 'dashed') +
      # facet_wrap(get(groups)~.) +
      facet_wrap(facets=head(groups,n=2)) + 
      theme_bw()
    g
  } else if (plotBiPlot) {
    sigdim <- unique(df$Dim[df$Cor >= threshCor & df$Cor.p < threshP]) %>% sort()
    lapply(sigdim, function(x) {
      tmp <- data.frame(Dim = x, SVD = SVD[,x], Metric = meta[,c(metric)], Color = meta[,c(groups)])
      siggroup <- df[df$Dim == x & df$Cor >= threshCor & df$Cor.p < threshP,] %>% as.data.frame %>% .[,2]
      tmp$SigGroup <- tmp$Color %in% siggroup
      g <- ggplot(tmp, aes(x = Metric, y = SVD, color = Color, alpha = SigGroup)) + geom_point(size = 0.1) + 
        xlab(metric) + ggtitle(paste("Dim",x)) + theme_bw() +
        scale_alpha_manual(values = c(`FALSE` = 0.2, `TRUE` = 0.6)) + 
        geom_smooth(data = tmp[tmp$SigGroup == T,], 
                    # method = "loess", formula = 'y ~ x',
                    method = "lm", formula = 'y ~ x',
                    se = F, show.legend = F)
      if (metric %in% c("nFrags","total_deduplicated","nCount_RNA")) { g <- g + scale_x_log10() }
      g
    }) %>% plot_grid(plotlist=.)
  } else { df }
}

drugScore <- function(se, 
                      markers,
                      modality) {
  
  # to account for gene names that were "fixed" in the process of marker calculation with repTest()
  rownames(se) <- fixGene(rownames(se), new)
  
  # prepare the metadata for each replicate and cell type in the summarizedExperiment object
  score <- data.frame(Total = se %>% colSums)
  score <- score %>% rownames_to_column(var = "Condition") %>% separate(Condition, into = c("Condition", "Replicate"), sep = "_") %>% separate(Condition, into = c("Cell_Type", "Drug", "Dose"), sep = " ",)
  score[is.na(score$Dose),]$Dose <- score[is.na(score$Dose),]$Drug
  score <- score %>% mutate(Drug = factor(Drug, levels = levels(df$Drug)), Dose = factor(Dose, levels = levels(df$Dose)))
  score <- score %>% filter(Cell_Type %ni% c("Other","B"))
  score <- score %>% arrange(Cell_Type,Drug,Dose)
  
  # select the markers (peaks, gene scores, genes) to be used
  markers <- markers %>% filter(MarkerType == modality)
  if (modality == "Genes") { 
    markers <- markers %>% filter(name %in% rownames(se))
  }
  
  # subset the activaton-associated markers for activation scoring
  activ <- markers %>% filter(Treatment == "Activation") %>% dplyr::select(CellType,name,Type)
  
  # iterate over each replicate/cell type combination
  score <- lapply(1:nrow(score), function(x) {
    df <- score[x,] %>% as_tibble
    c <- score$Cell_Type[x] %>% as.character
    d <- score$Drug[x] %>% as.character
    m <- score$Dose[x] %>% as.character
    r <- score$Replicate[x] %>% as.character
    
    # sum of values for activation-associated markers
    df$Activation <- se[activ %>% filter(CellType == c, name %in% rownames(se)) %>% .$name, colnames(se) %>% grep(c,.,value=T) %>% grep(d,.,value=T,fixed=T) %>% grep(m,.,value=T,fixed=T) %>% grep(paste("_",r,"$",sep=""),.,value=T,fixed=F)] %>% sum
    
    # sum of values for drug-associated markers (skip control replicates)
    if (df$Drug %in% c("DMSO(-)","DMSO(+)")) { 
      df$Drug_Markers <- df$Up_Markers <- df$Down_Markers <- df$Bgd_Markers <- list(c()) 
      df$Drug_Accessibility <- df$Up_Accessibility <- df$Down_Accessibility <- df$Bgd_Accessibility <- 1
      df$Drug_Baseline <- df$Up_Baseline <- df$Down_Baseline <- df$Bgd_Baseline <- 1
    }
    else {
      df <- mutate(df,
                   ## first generate marker list
                   # subset markers specific to each cell type and drug combo
                   # combine markers from paired PROTAC/inhibitors
                   # specifically exclude any that overlap with activation-associated markers
                   # separate markers by upregulated and downregulated
                   Drug_Markers = list(markers %>% filter(CellType == c, Target == unique(markers$Target[markers$Treatment == d]), name %ni% activ$name[CellType == c])),
                   Up_Markers = Drug_Markers[[1]] %>% filter(Type == "Up") %>% .$name %>% unique %>% list,
                   Down_Markers = Drug_Markers[[1]] %>% filter(Type == "Down") %>% .$name %>% unique %>% list,
                   Bgd_Markers = rownames(se)[rownames(se) %ni% markers$name] %>% unique %>% list,
                   ## then compute the aggregate value (accessibility or expression) of these marker lists for each replicate/cell type combination
                   Drug_Accessibility = se[Drug_Markers[[1]]$name, grep(paste(c, d, m, r), gsub("_"," ",colnames(se)))] %>% sum,
                   Up_Accessibility = se[Up_Markers[[1]], grep(paste(c, d, m, r), gsub("_"," ",colnames(se)))] %>% sum,
                   Down_Accessibility = se[Down_Markers[[1]], grep(paste(c, d, m, r), gsub("_"," ",colnames(se)))] %>% sum,
                   Bgd_Accessibility = se[Bgd_Markers[[1]], grep(paste(c, d, m, r), gsub("_"," ",colnames(se)))] %>% sum)
      
      df <- df %>% mutate(Drug_Markers = Drug_Markers[[1]] %>% .$name %>% unique %>% list,
                          Drug_Accessibility = se[Drug_Markers[[1]], grep(paste(c, d, m, r), gsub("_"," ",colnames(se)))] %>% sum)
      
      ## calculate the "baseline" for each marker list, i.e. the mean of the aggregated values for all 12 DMSO(+) control replicates
      df <- mutate(df,
                   Drug_Baseline = se[Drug_Markers[[1]], grep(paste(c, "DMSO.\\+"), colnames(se))] %>% colSums %>% mean,
                   Up_Baseline = se[Up_Markers[[1]], grep(paste(c, "DMSO.\\+"), colnames(se))] %>% colSums %>% mean,
                   Down_Baseline = se[Down_Markers[[1]], grep(paste(c, "DMSO.\\+"), colnames(se))] %>% colSums %>% mean,
                   Bgd_Baseline = se[Bgd_Markers[[1]], grep(paste(c, "DMSO.\\+"), colnames(se))] %>% colSums %>% mean)
    }
    df
  }) %>% do.call(rbind,.)
  
  ## replace the marker lists with the number of markers
  score <- score %>% mutate(Up_Markers = purrr::map_int(Up_Markers, length),
                            Down_Markers = purrr::map_int(Down_Markers, length),
                            Drug_Markers = purrr::map_int(Drug_Markers, length),
                            Bgd_Markers = purrr::map_int(Bgd_Markers, length))
  
  ## score relative activation by substracting the activation of DMSO(-) controls and then normalizing to DMSO(+) controls
  ## score up- and downregulated drug markers sets as the log2 fold-change relative to the pre-determined "baseline"/DMSO(+) value
  ## the aggregate drug score is the weighted average of absolute up and down scores, according to relative proportion of up vs down markers
  score <- score %>% group_by(Cell_Type) %>% 
    mutate(ActivationScore = (Activation - mean(Activation[Drug == "DMSO(-)"]))/(mean(Activation[Drug == "DMSO(+)"])-mean(Activation[Drug == "DMSO(-)"])),
           Up_Score = log2(Up_Accessibility / Up_Baseline),
           Down_Score = log2(Down_Accessibility / Down_Baseline),
           Drug_Score = case_when(Drug_Markers > 0 ~ abs(Up_Score) * Up_Markers/Drug_Markers + abs(Down_Score) * Down_Markers/Drug_Markers,
                                  Drug_Markers == 0 ~ 0),
           Bgd_Score = log2(Bgd_Accessibility / Bgd_Baseline)) 
  
  score <- score %>% dplyr::select(-grep("Baseline",colnames(.))) %>% pivot_longer(cols = grep(c("Bgd|Down|Up|Drug_"), colnames(.)), names_to = c("Type",".value"), names_sep="_")
  
  score$Score[!is.finite(score$Score)] <- 0
  
  score <- score %>% mutate(Drug = factor(Drug, levels = levels(df$Drug), labels = c("Resting Ctrl", "Stimulated Ctrl", levels(df$Drug)[3:8])), 
                            Dose = factor(Dose, levels = levels(df$Dose)),
                            Target = factor(Drug, c("Resting Ctrl", "Stimulated Ctrl", levels(df$Drug)[3:8]), labels = rep(c("Control","PRC2","SWI/SNF","p300/CBP"), each = 2)),
                            Type = factor(Type, levels = c("Up","Down","Bgd","Drug")),
                            Cell_Type = factor(Cell_Type, levels = c("T","Myeloid")),
                            DrugType = factor(Drug, c("Resting Ctrl", "Stimulated Ctrl", levels(df$Drug)[3:8]), labels = c("Vehicle","Vehicle",rep(c("PROTAC","Inhibitor"), 3))))
  
  
  score <- score %>% group_by(Cell_Type, Drug, Dose, Type, Target, DrugType, Markers) %>% 
    summarize(ActivationScore = mean(ActivationScore),
              Score = mean(Score))
  
  ## fit smooth curve through dose response trajectory
  spline <- score %>% filter(Dose %ni% c("DMSO(-)","DMSO(+)"), abs(Score) > 0) %>% group_by(Cell_Type, Drug, Type, Target, DrugType) %>% arrange(Dose) %>% group_split() %>% lapply(., function(x) { 
    spline <- smooth_it(x$ActivationScore, x$Score) %>% as.data.frame
    spline$Drug <- x$Drug[1]
    spline$Target <- x$Target[1]
    spline$Cell_Type <- x$Cell_Type[1]
    spline$Type <- x$Type[1]
    spline$DrugType <- x$DrugType[1]
    spline
  }) %>% do.call(rbind,.)
  
  g <- score %>% filter(Target != "Control", Type == "Drug") %>%
    {
      ggplot(., aes(x = ActivationScore, y = Score, color = Drug)) +
        geom_hline(yintercept = 0, color = "grey40", size = 1) +
        geom_point(x = 1, y = 0, color = 'grey65', size = 2.5) +
        geom_vline(xintercept = 0, size = 1, color="grey80") + geom_vline(xintercept = 1, size = 1, color="grey65") +
        geom_point(mapping = aes(size = Dose), alpha = 0.75, show.legend = T) + 
        geom_path(data = filter(spline, Cell_Type %in% .$Cell_Type, Type %in% .$Type),
                  mapping = aes(x = x, y = y, group = paste(Drug,Type), linetype = DrugType), size = 1, alpha = 0.75, show.legend = F) +
        geom_label(data = . %>% ungroup %>% dplyr::select(Cell_Type,Target,Markers,Type) %>% distinct,
                   mapping = aes(label = Markers), color = 'black', alpha = 0.5, size = 2, x = 0.5, y = Inf, vjust = 1) +
        scale_x_continuous(breaks = c(seq(0,1,length=3)), limits = c(0.49 - max(abs(.$ActivationScore-0.5)), 0.51 + max(abs(.$ActivationScore-0.5)))) +
        scale_y_continuous(breaks = scales::pretty_breaks(3)) +
        scale_size_manual(values = c(1,2,3)) + scale_linetype_manual(values = c("11","solid")) + scale_color_manual(values = col_drug) + 
        facet_grid(Cell_Type~Target, scales = 'free_y', drop = T) + coord_cartesian(clip = 'off') +
        theme_DC + theme(
          axis.line = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.key.size = unit(0.75, 'lines')) + 
        guides(colour = guide_legend(override.aes = list(size=3, shape = 16), byrow = T), linetype = 'none') +
        ggtitle("Drug Response") + labs(x="Normalized Activation Score", y = "Marker Gene Accessibility Score")
    }
  return(g)
}


# Plotting ----------------------------------------------------------------

save_plot <- function(g, file_name, w=2, h=1.5, x.angle=NA, h.just=NA, xlab=NA, ylab=NA, show.legend="none") {
  if(!(is.na(xlab)|is.na(ylab))) { g <- g + labs(x=xlab, y=ylab) }
  if(!(is.na(x.angle)|is.na(h.just))) { g <- g + theme(axis.text.x=element_text(angle=x.angle, hjust=h.just)) }
  g <- g + theme(legend.position = show.legend)
  ggsave(plot = g, filename = file_name, device = "pdf", width = w, height = h, units = "in")
}

plotHeatmap <- function(repCor, assay = c("Peaks","Scores","RNA","Motifs"), p = 0.01, coef = 0, log2fc = 0, cellType = c("T","Myeloid"), drugsSig = drug, drugsPlot = drug, cutCol = 8, pal = NULL, filename = NA, width = 8, height = 4, returnRes = F) {
  if (assay == "Peaks") { pal <- c("#FFC748","#425BB4") %>% rev }
  if (assay == "Scores") { pal <- c("#FFB148","#397AAC") %>% rev }
  if (assay == "RNA") { pal <- c("#FF9248","#30AA96") %>% rev }
  if (assay == "Motifs") { pal <- c("#D84347FF","#6B757DFF") %>% rev }
  start <- (col2rgb(pal[1]) %>% t * 0.5) %>% rgb(maxColorValue = 255)
  end <- (col2rgb(pal[2]) %>% t * 0.5) %>% rgb(maxColorValue = 255)
  start2 <- (col2rgb(pal[1]) %>% t * 0.75) %>% rgb(maxColorValue = 255)
  end2 <- (col2rgb(pal[2]) %>% t * 0.75) %>% rgb(maxColorValue = 255)
  pal <- colorRampPalette(c(start,start2,pal[1],"white",pal[2],end2,end))(100)
  if (assay == "Motifs") { log2 <- F } else { log2 <- T }
  
  heatmap <- repCor[c("p","coef","log2fc.Max","summ")] %>% {.[lapply(.,nrow) > 0]} %>%
    lapply(filter, CellType == cellType[1]) %>%
    { 
      tmp <- lapply(.,filter,Drug %in% drugsSig)
      if (assay == "Motifs") {
        idx <- which(colSums(tmp$p[,-1:-2] < p & abs(tmp$coef[,-1:-2]) > coef) >= 1)
      } else {
        idx <- which(colSums(tmp$p[,-1:-2] < p & abs(tmp$coef[,-1:-2]) > coef & tmp$log2fc.Max[,-1:-2] > log2fc) >= 1)
      }
      cat(length(idx), "Markers Found\n")
      tmp <- lapply(.,filter,Drug %in% drugsPlot)
      tmp$summ[,c(1:3, idx + 3)]
    } %>% unite("rowname", Drug, Dose, sep = ' ') %>% ungroup %>% select(-CellType) %>% 
    mutate(rowname = factor(rowname, 
                            levels = c(paste(c("DMSO(-)","DMSO(+)"),c("DMSO(-)","DMSO(+)")), paste(rep(levels(df$Drug)[3:8], each = 3), rep(levels(df$Dose)[3:5], 6))),
                            labels = c("NegCtrl","PosCtrl",paste(rep(levels(df$Drug)[3:8], each = 3), rep(levels(df$Dose)[3:5], 6))))) %>%
    arrange(rowname) %>% column_to_rownames() %>% 
    as.matrix
  
  cat("Clustering Markers...\n")
  hc <- hclust(dist(t(heatmap %>% colZscores)), method = "ward.D2")
  
  anno_row <- data.frame(rowname = proj_4$CellDrugDose2) %>% mutate(rowname = gsub("^T |^Myeloid |^B ", "", rowname)) %>% distinct %>% separate(col = "rowname", into = c("Drug","Dose"), sep = " ", remove = F) %>% mutate(Dose = case_when(is.na(Dose) ~ Drug, TRUE ~ Dose)) %>% column_to_rownames()
  
  if (cutCol > 0) {
    cut <- cutree(hc, k = cutCol)
    anno_col <- data.frame(Group = cut[hc$order]) %>% mutate(Group = as.character(Group))
  } else {
    cut <- NA
    cutCol <- NA
    anno_col <- NA
  }
  
  cat("Plotting Heatmap...\n")
  pheatmap(heatmap, scale = 'column', cluster_rows = F, show_colnames = F, border_color = NA, color = pal,
           cluster_cols = hc, cutree_col = cutCol, treeheight_col = 0, fontsize_row = 7,
           annotation_row = anno_row,
           annotation_col = anno_col, annotation_names_col = F,
           annotation_colors = list(Drug = col_drug2[-1], Dose = col_dose2, Group = col_clust), annotation_legend = F,
           filename = filename, width = width, height = height, main = paste(cellType[1], "Cells -", assay, "-", length(hc$order),"Markers"))
  if (returnRes) { list(heatmap = heatmap,
                        hclust = hc,
                        cut = cut) }
}

CovPlot <- function(object, region, group.by, up, down, colors, window = 100, assay = "ATAC", ratio = c(5,1,1)) {
  if (region %in% getGeneAnnotation(proj_4)$genes$symbol) {
    region <- getGeneAnnotation(proj_4)$genes %>% subset(symbol == region) %>% as.data.frame %>% {paste(.$seqnames, .$start, .$end, sep='-')}
  } 
  p1 <- CoveragePlot(object, region, group.by = group.by, extend.upstream = up, extend.downstream = down, assay = assay, window = window,
                     annotation = F, peaks = F, links = F, tile = F) + scale_fill_manual(values = colors)
  p2 <- AnnotationPlot(object, region, extend.upstream = up, extend.downstream = down)
  p3 <- PeakPlot(object, region, extend.upstream = up, extend.downstream = down, assay = assay)
  CombineTracks(list(p1,p2,p3), heights = ratio)
}

plotGeneSet <- function(sets, repCor, assay, log2fc = F, celltypes = c("T","Myeloid"), returnList = F) {
  label <- 
    if (log2fc & assay == "RNA") { "Log2FC Expression" } 
  else if (log2fc & assay == "Score") { "Log2FC Accessibility" }
  else if (!log2fc & assay == "RNA") { "Scaled Expression" }
  else if (!log2fc & assay == "Score") { "Scaled Accessibility" }
  g <- lapply(seq_along(sets), function(x) {
    title <- gsub(".*HALLMARK_|.*REACTOME_","",names(sets)[x])
    repCor$summ[,c("CellType","Drug","Dose", sets[x] %>% unlist %>% unique %>% subset(. %in% colnames(repCor$mean)))] %>%
      dupControl(controls = "DMSO(+)", drug[c(1,3:8)]) %>% filter(!(Dose == "DMSO(+)" & Drug == "DMSO(-)")) %>%
      filter(CellType %in% celltypes) %>% filter(Drug != "DMSO(-)") %>%
      mutate(Drug = factor(Drug, levels = levels(df$Drug)),
             CellType = factor(CellType, c("T","Myeloid","B"))) %>%
      {
        if (log2fc) {
          melt(.) %>% group_by(variable, CellType, Drug) %>% 
            mutate(fc = value/value[Dose == "DMSO(+)"]) %>% mutate(value = log2(fc))
        } else {
          melt(.) %>% group_by(variable, CellType, Drug) %>%
            mutate(fc = value/value[Dose == "DMSO(+)"]) %>% 
            group_by(CellType) %>% mutate(value = fc * mean(value[Dose == "DMSO(+)"]))
        }
      } %>%
      ggplot(aes(x = Dose, y = value, color = Drug)) + 
      geom_smooth(mapping = aes(group = Drug, fill = Drug), alpha = 0.2, method = 'lm', formula = y~poly(x,2)) +
      # geom_smooth(mapping = aes(group = Drug, fill = Drug), alpha = 0.2, method = 'lm', formula = y~x) +
      # geom_smooth(mapping = aes(group = Drug, fill = Drug), alpha = 0.2) +
      facet_wrap(CellType~.) + scale_color_manual(values=col_drug) + scale_fill_manual(values=col_drug) + labs(y = label, title = title) +
      theme_DC + theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_x_discrete(expand = expansion(add = 0.2))
  }) 
  if (returnList) { return(g) }
  else {
    plot_grid(plotlist = g, nrow = 1)
  }
} 

plotGeneSet2 <- function(sets, assay, log2fc = F, celltypes = c("T","Myeloid"), drugs = drug[3:8], facets = c("CellType","Type"), hline = NULL, nrow = 1) {
  label <- 
    if (log2fc) { "Log2FC " }
  else if (!log2fc) { "Scaled" }
  g <- lapply(seq_along(sets), function(x) {
    title <- gsub(".*HALLMARK_|.*REACTOME_","",names(sets)[x])
    
    rbind(data.frame(lm_score$summ[,c("CellType","Drug","Dose", sets[x] %>% unlist %>% unique %>% subset(. %in% colnames(lm_score$summ)))] %>%
                       dupControl(controls = "DMSO(+)", drug[c(1,3:8)]) %>% filter(!(Dose == "DMSO(+)" & Drug == "DMSO(-)")) %>%
                       filter(CellType %in% celltypes) %>% filter(Drug != "DMSO(-)") %>%
                       mutate(Drug = factor(Drug, levels = levels(df$Drug)),
                              CellType = factor(CellType, c("T","Myeloid","B"))) %>%
                       {
                         if (log2fc) {
                           melt(.) %>% group_by(variable, CellType, Drug) %>% 
                             mutate(fc = value/value[Dose == "DMSO(+)"]) %>% mutate(value = log2(fc))
                         } else {
                           melt(.) %>% group_by(variable, CellType, Drug) %>%
                             mutate(fc = value/value[Dose == "DMSO(+)"]) %>% 
                             group_by(CellType) %>% mutate(value = fc * mean(value[Dose == "DMSO(+)"]))
                         }
                       },
                     Type = "Score"),
          data.frame(lm_rna$summ[,c("CellType","Drug","Dose", sets[x] %>% unlist %>% unique %>% subset(. %in% colnames(lm_rna$summ)))] %>%
                       dupControl(controls = "DMSO(+)", drug[c(1,3:8)]) %>% filter(!(Dose == "DMSO(+)" & Drug == "DMSO(-)")) %>%
                       filter(CellType %in% celltypes) %>% filter(Drug != "DMSO(-)") %>%
                       mutate(Drug = factor(Drug, levels = levels(df$Drug)),
                              CellType = factor(CellType, c("T","Myeloid","B"))) %>%
                       {
                         if (log2fc) {
                           melt(.) %>% group_by(variable, CellType, Drug) %>% 
                             mutate(fc = value/value[Dose == "DMSO(+)"]) %>% mutate(value = log2(fc))
                         } else {
                           melt(.) %>% group_by(variable, CellType, Drug) %>%
                             mutate(fc = value/value[Dose == "DMSO(+)"]) %>% 
                             group_by(CellType) %>% mutate(value = fc * mean(value[Dose == "DMSO(+)"]))
                         }
                       },
                     Type = "RNA")) %>%
      filter(Drug %in% drugs) %>%
      ggplot(aes(x = Dose, y = value, color = Drug)) + 
      # geom_smooth(mapping = aes(group = Drug, fill = Drug), alpha = 0.2, method = 'lm', formula = y~poly(x,2)) + 
      # geom_smooth(mapping = aes(group = Drug, fill = Drug), alpha = 0.2, method = 'lm', formula = y~x) +
      # geom_smooth(mapping = aes(group = Drug, fill = Drug), alpha = 0.2) +
      geom_smooth(mapping = aes(group = Drug, fill = Drug), alpha = 0.2, span = 1.5) +
      geom_hline(yintercept = hline, color = 'grey40', linetype = '22') +
      facet_grid(as.formula(paste(facets[1], "~", facets[2])), scales = 'free_y', axes = 'all_x') + scale_color_manual(values=col_drug) + scale_fill_manual(values=col_drug) + labs(y = label, title = title) +
      theme_DC + theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_x_discrete(expand = expansion(add = 0.2))
  }) 
  plot_grid(plotlist = g, nrow = nrow)
} 


ggFootprint <- function (seFoot = NULL, name = NULL, subset = NULL,
                         smoothWindow = NULL, flank = NULL, flankNorm = NULL, normMethod = NULL, dupControl = NULL,
                         colorBy = NULL, groupings = NULL, facets = NULL, xlim = NULL, pal = NULL,
                         returnDF = F) 
{
  rowDF <- SummarizedExperiment::rowData(seFoot)
  footMat <- ArchR:::.getAssay(seFoot[BiocGenerics::which(rowDF[, 2] == "footprint"), ], name)
  biasMat <- ArchR:::.getAssay(seFoot[BiocGenerics::which(rowDF[, 2] == "bias"), ], name)
  footDF <- rowDF[BiocGenerics::which(rowDF[, 2] == "footprint"), ]
  biasDF <- rowDF[BiocGenerics::which(rowDF[, 2] == "bias"), ]
  footMat <- apply(footMat, 2, function(x) ArchR:::.centerRollMean(x, smoothWindow))
  biasMat <- apply(biasMat, 2, function(x) ArchR:::.centerRollMean(x, smoothWindow))
  idx <- which(abs(footDF$x) >= flank - flankNorm)
  footMat <- t(t(footMat)/colMeans(footMat[idx, , drop = FALSE]))
  biasMat <- t(t(biasMat)/colMeans(biasMat[idx, , drop = FALSE]))
  if (tolower(normMethod) == "none") {
    title <- ""
  }
  else if (tolower(normMethod) == "subtract") {
    title <- "Tn5 Bias Subtracted\n"
    footMat <- footMat - biasMat
  }
  else if (tolower(normMethod) == "divide") {
    title <- "Tn5 Bias Divided\n"
    footMat <- footMat/biasMat
  }
  else {
    stop("normMethod not recognized!")
  }
  footMatMean <- ArchR:::.groupMeans(footMat, SummarizedExperiment::colData(seFoot)$Group)
  footMatSd <- ArchR:::.groupSds(footMat, SummarizedExperiment::colData(seFoot)$Group)
  
  smoothFoot <- rowMaxs(apply(footMatMean, 2, function(x) ArchR:::.centerRollMean(x, 11)))
  
  plotIdx <- seq_len(nrow(footMatMean))
  plotFootDF <- lapply(seq_len(ncol(footMatMean)), function(x) {
    data.frame(x = footDF$x, mean = footMatMean[, x], sd = footMatSd[, x], group = colnames(footMatMean)[x])[plotIdx, , drop = FALSE]
  }) %>% Reduce("rbind", .)
  
  plotFootDF$group <- factor(paste0(plotFootDF$group), levels = gtools::mixedsort(unique(paste0(plotFootDF$group))))
  
  if (!is.null(groupings)) {
    tmp <- strsplit(as.character(plotFootDF$group), split = "_| ")
    # if (!is.null(groupings)) {
    for (x in 1:length(groupings)) {
      values <- unlist(lapply(tmp, function(y) { y[x] }))
      values <- data.frame(values)
      colnames(values) <- groupings[x]
      values[is.na(values)] <- unlist(lapply(tmp, function(y) { y[x-1] }))[is.na(values)]
      plotFootDF <- cbind(plotFootDF, values)
    }
  }
  
  if (!is.null(plotFootDF$Dose)) {
    plotFootDF$Dose <- factor(plotFootDF$Dose,
                              levels = c("DMSO(-)","DMSO(+)","NegCtrl","PosCtrl","10nM","100nM","1µM","HighDose"))
  }
  if (!is.null(plotFootDF$CellType)) {
    plotFootDF$CellType <- factor(plotFootDF$CellType,
                                  levels = c("B","T","Myeloid"))
  }
  if (!is.null(plotFootDF$Drug)) {
    plotFootDF$Drug <- factor(plotFootDF$Drug,
                              levels = c("DMSO(-)","DMSO(+)","NegCtrl","PosCtrl","MS177","EPZ6438","AU-15330","BRM014","dCBP-1","GNE-781"))
  }
  
  if (is.null(colorBy)) { colorBy <- "group" }
  if (is.null(pal)) { pal <- paletteDiscrete(values = levels(plotFootDF[,colorBy])) }
  if (is.null(xlim)) {
    xlim <- c(min(plotFootDF$x), max(plotFootDF$x))
  }
  
  if (!is.null(dupControl)) {
    plotFootDF <- plotFootDF %>% dupControl() #%>% filter(Drug != dupControl)
  }
  
  if (!is.null(subset)) {
    plotFootDF <- plotFootDF %>% 
      filter(eval(rlang::parse_expr(subset)))
  }
  
  if (returnDF) { plotFootDF }
  else {
    g <- ggplot(plotFootDF, aes(x = x, y = mean, color = .data[[colorBy]], group = group)) + 
      # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, #linetype = NA, 
      #                 fill = .data[[colorBy]]), alpha = 0.4) + 
      geom_line(linewidth = 0.5) + 
      scale_color_manual(values = pal) + 
      scale_fill_manual(values = pal) + 
      xlab("Distance to motif center (bp)") + 
      coord_cartesian(expand = FALSE,
                      ylim = c(quantile(plotFootDF$mean, 1e-04), 1.15 * quantile(smoothFoot, 0.999)),
                      xlim = xlim) +
      theme_minimal() +
      theme(legend.position = 'bottom') + 
      ggtitle(name) + 
      ylab(paste0(title, "Normalized Insertions")) 
    if (is.null(facets)) { g }
    else if (length(facets) == 1) { g + facet_grid(as.formula(paste(". ~", facets))) }
    else { g + facet_grid(as.formula(paste(facets[2], "~", facets[1]))) }
  }
}

plotLSIHeatmap <- function(proj = c("individual","full"), 
                           replicates = F, 
                           celltypes = c("T","Myeloid","B"),
                           drugs = levels(df$Drug),
                           doses = levels(df$Dose),
                           dims = 2:30, 
                           title = NULL,
                           fontsize = 7,
                           filename = NA,
                           width = NULL,
                           height = NULL) {
  
  list <- lapply(celltypes, function(y) {
    if (proj == "individual") {
      proj <- get(paste("proj",str_sub(y, 1, 1),sep="_"))
      SVD <- proj@reducedDims$IterativeLSI$matSVD[,dims]
    } else if (proj == "full") {
      proj <- proj_4
      SVD <- proj@reducedDims$IterativeLSI$matSVD[,dims]
    }
    if (replicates) { 
      heat_col <- paletteer::paletteer_c("grDevices::PRGn", 101, direction = -1)
      anno <- df %>% filter(Type == "Control") %>% mutate(rowname = paste(Drug,Replicate)) %>% dplyr::select(Drug, PlateSide, Plate, rowname) %>% distinct %>% remove_rownames() %>% column_to_rownames()
      condition <- unique(paste(df$DrugDose, df$Replicate)) 
      condition <- grep("DMSO", condition, value = T)
      tmp <- lapply(condition, function(x) {
        d <- str_split(x, " ",simplify = T)
        cells <- proj$cellNames[proj$Drug == d[,1] & proj$Replicate == d[,2] & proj$CellType_Major == y]
        df <- SVD[cells,] %>% colMeans() %>% data.frame()
        colnames(df) <- x
        df
      }) %>% do.call(cbind,.) %>% cor()
    } else {
      heat_col <- paletteer::paletteer_c("grDevices::RdBu", 101, direction = -1)
      anno <- df %>% dplyr::select(Drug, Dose, DrugDose) %>% distinct %>% remove_rownames() %>% column_to_rownames("DrugDose") 
      anno <- anno %>% filter(Drug %in% drugs, Dose %in% doses)
      # condition <- unique(df$DrugDose) 
      # tmp <- lapply(condition, function(x) {
      tmp <- lapply(rownames(anno), function(x) {
        d <- str_split(x, " ",simplify = T)
        if (length(d) == 1) { 
          cells <- proj$cellNames[proj$Drug == d[,1] & proj$CellType_Major == y]
        } else { 
          cells <- proj$cellNames[proj$Drug == d[,1] & proj$Dose == d[,2] & proj$CellType_Major == y]
        }
        df <- SVD[cells,] %>% colMeans() %>% data.frame()
        colnames(df) <- x
        df
      }) %>% do.call(cbind,.) %>% cor()
    }
    
    anno_col <- list(Drug = col_drug[names(col_drug) %in% anno$Drug], Dose = col_dose, PlateSide = col_side, Plate = c(A = "red",B = "blue"))
    
    pheatmap(tmp, main = title, breaks = seq(-1,1,.02), legend = T, fontsize = fontsize,
             color = heat_col, 
             cutree_rows = 2, cutree_cols = 2,
             border_color = NA, show_colnames = F, treeheight_col = 0, treeheight_row = 10,
             annotation_row = anno, annotation_colors = anno_col,
             filename = filename, width = width, height = height)
  })
}


plotLISI <- function(x, embedding = c("UMAP","LSI"), tn5 = "Tn5", grouping = "Group", bins = 3, dims = 1:30, returnDF = F) {
  
  if (embedding == "UMAP" & is.data.frame(x)) {
    meta <- x %>% select(Dataset,Group,UMAP1,UMAP2,Bin,Random) %>% group_by(Dataset, Group) %>% group_modify(~ compute_lisi(.x[,c("UMAP1","UMAP2")], .x, label_colnames = c("Bin","Random")))
    meta <- meta %>% select(Dataset,Group,Bin,Random) #%>% reshape2::melt()
  } else if (embedding == "LSI" & class(x) == "ArchRProject") {
    lsi <- x@reducedDims$IterativeLSI$matSVD[x$cellNames,dims]
    meta <- x@cellColData[x$cellNames,] %>% as.data.frame()
    
    meta <- meta %>% tibble::rownames_to_column(var = "row.name") %>% tidyr::unite("Group", grouping) %>% dplyr::select(Cell = row.name, Tn5 = tn5, Group) %>% 
      group_by(Tn5, Group) %>% mutate(Tn5Count = n()) %>%
      group_by(Group) %>% mutate(Tn5Rank = rank(Tn5Count),
                                 Bin = factor(Tn5Rank, levels = sort(unique(Tn5Rank)), labels = ntile(sort(unique(Tn5Rank)), bins)),
                                 Random = sample(Bin, size = n(), replace=F))
    
    meta <- meta %>% group_by(Group) %>% group_modify(~ compute_lisi(lsi[.x$Cell,], .x, label_colnames = c("Bin","Random")))
    
    meta <- meta %>% select(Group,Bin,Random) #%>% reshape2::melt()
  }
  
  if (returnDF) {
    meta
  } else {
    ggplot(meta, aes(x = variable, y = value, color = variable)) + geom_boxplot(outliers = F) + theme_bw() + labs(y = "LISI", x = "Cell Grouping")
  }
}





