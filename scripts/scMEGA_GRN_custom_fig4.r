library(ArchR)
library(Seurat)
library(Signac)
library(scMEGA)
library(harmony)
library(Nebulosa)
library(ggplot2)
library(dplyr)
library(TFBSTools)
library(igraph)
library(ggraph)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(igraph)
library(JASPAR2024)
library(JASPAR2020)
library(GenomeInfoDb)
library(TxDb.Athaliana.BioMart.plantsmart51)
library(org.At.tair.db)


setwd("/anvil/projects/x-mcb130189/trlee/ethylene/fig4")




###### redefine functions for Arabidopsis analysis


ggPoint_custom <- function(
  x = NULL, 
  y = NULL, 
  color = NULL, 
  discrete = TRUE, 
  discreteSet = "stallion",
  continuousSet = "solarExtra", 
  labelMeans = TRUE,  
  pal = NULL, 
  defaultColor = "lightGrey",
  highlightPoints = NULL,
  colorDensity = FALSE,
  size = 1, 
  xlim = NULL, 
  ylim = NULL, 
  extend = 0.05, 
  xlabel = "x", 
  ylabel = "y", 
  title = "", 
  randomize = FALSE, 
  seed = 1,
  colorTitle = NULL, 
  colorOrder = NULL, 
  colorLimits = NULL,
  alpha = 1, 
  baseSize = 10, 
  legendSize = 3,
  ratioYX = 1, 
  labelAsFactors = TRUE,
  fgColor = "black", 
  bgColor = "white", 
  bgWidth = 1,
  labelSize = 3,
  addFit = NULL, 
  rastr = FALSE, 
  dpi = 300,
  ...
  ){

  stopifnot(length(y) == length(x))
  if(length(x) < 5){
    stop("x must be at least length 5 to plot!")
  }

  if(randomize){
    set.seed(seed)
    idx <- sample(seq_along(x), length(x))
  }else{
    idx <- seq_along(x)
  }

  df <- data.frame(x = x, y = y)
  include <- which(is.finite(x) & is.finite(y))
  
  if(length(include) != length(x)){
    message("Some values are not finite! Excluding these points!")
    df <- df[include,]
    x <- x[include]
    y <- y[include]
    if(!is.null(color)){
      color <- color[include]
    }
  }

  if(is.null(xlim)){
      xlim <- range(df$x) %>% extendrange(f = extend)
  }

  if(is.null(ylim)){
      ylim <- range(df$y) %>% extendrange(f = extend)
  }

  ratioXY <- ratioYX * diff(xlim)/diff(ylim)

  #Plot

  if (is.null(color) & !colorDensity) {

    p <- ggplot(df[idx,], aes(x = x, y = y)) + 
        coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = F) + 
        xlab(xlabel) + ylab(ylabel) + 
        ggtitle(title) + 
        theme_ArchR(baseSize = baseSize)

    if(rastr){
      p <- p + .geom_point_rast2(
            size = size, raster.dpi = dpi, alpha = alpha, color = defaultColor)
      # if(!requireNamespace("ggrastr", quietly = TRUE)){
      #   message("ggrastr is not available for rastr of points, continuing without rastr!")
      #   p <- p + geom_point(size = size, alpha = alpha, color = defaultColor)
      # }else{
      #   .requirePackage("ggrastr")
      #   p <- p + geom_point_rast(
      #       size = size, raster.dpi = dpi, alpha = alpha, color = defaultColor)
      # }
    }else{
      p <- p + geom_point(size = size, alpha = alpha, color = defaultColor)
    }
      
  }else {

      if(colorDensity){
        
        discrete <- FALSE
        df <- .getDensity(x, y, n = 100, sample = NULL) #change
        df <- df[order(df$density), ,drop=FALSE]
        df$color <- df$density
        
        if(is.null(colorTitle)){
          colorTitle <- "density"
        }

      }else if(discrete){

        if(!is.null(highlightPoints)){
          if(length(highlightPoints) < length(color)){
            color[-highlightPoints] <- "Non.Highlighted"
            idx <- c(idx[-highlightPoints], idx[highlightPoints])
          }
        }
        color <- paste0(color)

        if(!is.null(colorOrder)){
          if(!all(color %in% colorOrder)){
            stop("Not all colors are in colorOrder!")
          }
        }else{
          colorOrder <- gtools::mixedsort(unique(color))
        }

        if(is.null(colorTitle)){
          colorTitle <- "color"
        }
        
        stopifnot(length(color) == nrow(df))
        df$color <- factor(color, levels = colorOrder)
        
        if(labelAsFactors){
          df$color <- factor(
            x = paste0(paste0(match(paste0(df$color), paste0(levels(df$color)))), "-", paste0(df$color)), 
            levels = paste0(seq_along(levels(df$color)), "-", levels(df$color))
          )
          if(!is.null(pal)){
            #print(pal)
            #print(paste0(levels(df$color))[match(names(pal), colorOrder)])
            names(pal) <- paste0(levels(df$color))[match(names(pal), colorOrder)]
          }
          colorOrder <- paste0(levels(df$color))
        }

      }else{
        stopifnot(length(color) == nrow(df))
        if(!is.null(highlightPoints)){
          if(length(highlightPoints) < length(color)){
            color[-highlightPoints] <- NA
            idx <- c(idx[-highlightPoints], idx[highlightPoints])
          }
        }
        if(!is.null(colorLimits)){
          color[color < min(colorLimits)] <- min(colorLimits)
          color[color > max(colorLimits)] <- max(colorLimits)
        }
        df$color <- color
      }

      p <- ggplot(df[idx,], aes(x = x, y = y, color = color)) +  
            coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = FALSE) + 
            xlab(xlabel) + ylab(ylabel) + 
            ggtitle(title) + theme_ArchR(baseSize = baseSize) +
            theme(legend.direction = "horizontal", legend.box.background = element_rect(color = NA)) +
            labs(color = colorTitle)

    if(rastr){
      
      p <- p + .geom_point_rast2(
            size = size, raster.dpi = dpi, alpha = alpha, 
            raster.width = min(par('fin')), 
            raster.height = (ratioYX * min(par('fin')))
          )
    
    }else{

        p <- p + geom_point(size = size, alpha = alpha)

    }

    if (discrete) {
        
        if (!is.null(pal)) {
            p <- p + scale_color_manual(values = pal)
        }else {
            pal <- paletteDiscrete(set = discreteSet, values = colorOrder)
            if(!is.null(highlightPoints)){
              pal[grep("Non.Highlighted", names(pal))] <- "lightgrey"
            }
            #print(pal)
            p <- p + scale_color_manual(values = pal) +
              .gg_guides(color = guide_legend(override.aes = list(size = legendSize, shape = 15)))
        }

        if (labelMeans) {
            
            dfMean <- split(df, df$color) %>% lapply(., function(x) {
              data.frame(x = median(x[, 1]), y = median(x[, 2]), color = x[1, 3])
            }) %>% Reduce("rbind", .)

            if(labelAsFactors){
              dfMean$label <- stringr::str_split(paste0(seq_len(nrow(dfMean))), pattern = "\\-", simplify=TRUE)[,1]
            }else{
              dfMean$label <- dfMean$color
            }
            dfMean$text <- stringr::str_split(dfMean$color, pattern = "-", simplify = TRUE)[,1]

            # make halo layers, similar to https://github.com/GuangchuangYu/shadowtext/blob/master/R/shadowtext-grob.R#L43
            theta <- seq(pi / 8, 2 * pi, length.out = 16)
            xo <- bgWidth * diff(range(df$x)) / 300
            yo <- bgWidth * diff(range(df$y)) / 300
            for (i in theta) {
              p <- p + 
                geom_text(data = dfMean, 
                    aes_q(
                      x = bquote(x + .(cos(i) * xo)),
                      y = bquote(y + .(sin(i) * yo)),
                      label = ~text
                    ),
                    size = labelSize,
                    color = bgColor
                )
            }

            if(is.null(fgColor)){
              p <- p + geom_text(data = dfMean, aes(x = x, y = y, color = color, label = label), size = labelSize, show.legend = FALSE)
            }else{
              p <- p + geom_text(data = dfMean, aes(x = x, y = y, label = label), color = fgColor, size = labelSize, show.legend = FALSE) 
            }

        }

    }else{

        if (!is.null(pal)) {
            if(!is.null(colorLimits)){
              p <- p + scale_colour_gradientn(colors = pal, limits=colorLimits, na.value = "lightgrey")
            }else{
              p <- p + scale_colour_gradientn(colors = pal, na.value = "lightgrey")
            }
        }else {
          if(!is.null(colorLimits)){
            p <- p + scale_colour_gradientn(colors = c("1" = "#ca6dca", "2"="#5826a3", "3"="#ffa22b", "4"="#b6a364", "5"="#f50501"), limits=colorLimits, na.value = "lightgrey")
          }else{
            p <- p + scale_colour_gradientn(colors = c("1" = "#ca6dca", "2"="#5826a3", "3"="#ffa22b", "4"="#b6a364", "5"="#f50501"), na.value = "lightgrey")
          }
        }
    }

  }

  if (!is.null(addFit)) {
      p <- p + geom_smooth(data = df, aes(color = NULL), method = addFit, color = "black") + 
        ggtitle(paste0(title, "\nPearson = ", round(cor(df$x, df$y), 3), "\nSpearman = ", round(cor(df$x, df$y, method = "spearman"), 3)))
  }

  p <- p + theme(legend.position = "bottom", legend.key = element_rect(size = 2))#, legend.spacing.x = unit(0.1, 'cm'), legend.spacing.y = unit(0.1, 'cm'))

  if(!is.null(ratioYX)){
    attr(p, "ratioYX") <- ratioYX
  }

  return(p)

}


########
########


TrajectoryPlot_custom <- function(object = NULL,
                           trajectory = "Trajectory",
                           reduction = NULL,
                           #name = "Trajectory",
                           size = 0.2,
                           rastr = FALSE,
                           quantCut = c(0.01, 0.99),
                           #quantHex = 0.5,
                           continuousSet = NULL,
                           discreteSet = NULL,
                           randomize = TRUE,
                           keepAxis = FALSE,
                           baseSize = 6,
                           addArrow = FALSE,
                           smoothWindow = 5) {
  dfT <- object@meta.data[, trajectory] %>%
    as.data.frame()

  rownames(dfT) <- colnames(object)
  idxRemove <- which(is.na(dfT[, 1]))
  df <- object@reductions[[reduction]]@cell.embeddings[, 1:2] %>%
    as.data.frame()

  dfT <- cbind(df, dfT[rownames(df),])
  colnames(dfT) <- c("DM1", "DM2", "PseudoTime")


  plotParams <- list()
  plotParams$x <- df[, 1]
  plotParams$y <- df[, 2]
  plotParams$title <- reduction
  plotParams$baseSize <- baseSize
  plotParams$color <- as.vector(dfT$PseudoTime)
  plotParams$discrete <- ArchR:::.isDiscrete(plotParams$color)
  plotParams$continuousSet <- "horizonExtra"
  plotParams$discreteSet <- "stallion"
  plotParams$title <- trajectory
  plotAs <- "points"
  #    paste(plotParams$title, " colored by\ncolData : ", name)
  #  if (is.null(plotAs)) {
  #    plotAs <- "hexplot"
  #  }

  plotParams$xlabel <- "DM1"
  plotParams$ylabel <- "DM2"

  if (!is.null(continuousSet)) {
    plotParams$continuousSet <- continuousSet
  }
  if (!is.null(continuousSet)) {
    plotParams$discreteSet <- discreteSet
  }
  plotParams$rastr <- rastr
  plotParams$size <- size
  plotParams$randomize <- randomize
  plotParams$color <- as.vector(plotParams$color)

  out <- do.call(ggPoint_custom, plotParams)
  out <- out +
    cowplot::theme_cowplot() +
    ggtitle(trajectory) +
    xlab(colnames(df)[1]) +
    ylab(colnames(df)[2]) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "cm")
    )


  dfT$value <- plotParams$color
  dfT <- dfT[order(dfT$PseudoTime), ]
  dfT <- dfT[!is.na(dfT$PseudoTime), ]

  colnames(dfT) <- c("x", "y", "PseudoTime", "value")
  dfT <- as.data.frame(dfT)

  if (addArrow) {
    dfArrow <-  split(dfT, floor(dfT$PseudoTime / 1.01)) %>%
      lapply(colMeans) %>%
      Reduce("rbind", .) %>%
      data.frame

    dfArrow$x <- ArchR:::.centerRollMean(dfArrow$x, smoothWindow)
    dfArrow$y <- ArchR:::.centerRollMean(dfArrow$y, smoothWindow)
    dfArrow <- rbind(dfArrow, dfT[nrow(dfT), , drop = FALSE])

    out <- out + geom_path(
      data = dfArrow,
      aes(x, y, color = NULL),
      size = 1,
      arrow = arrow(type = "open", length = unit(0.1, "inches"))
    )
  }

  return(out)
}



########
########


PeakToGene <- function(peak.mat,
                       gene.mat,
                       genome = "Ath",
                       max.dist = 10000,
                       method = "corrleation") {
  if (!(genome %in% c("hg19", "hg38", "mm9", "mm10", "Ath"))) {
    stop("Available genome are: hg19, hg38, mm9, mm10, Ath!")
  }

  if (genome == "Ath") {
    gene_anno <- geneAnnoAth
  } else if (genome == "hg38") {
    gene_anno <- geneAnnoHg38
  } else if (genome == "mm9") {
    gene_anno <- geneAnnoMm9
  } else if (genome == "mm10") {
    gene_anno <- geneAnnoMm10
  }

  ## create object for RNA data
  genes <- gene_anno$genes
  gene.use <-
    intersect(names(elementMetadata(genes)[, "symbol"]), rownames(gene.mat))

  genes <- genes[names(elementMetadata(genes)[, "symbol"]) %in% gene.use]
  gene.mat <- gene.mat[gene.use, ]

  gene_start <- ifelse(genes@strand == "+",
                       genes@ranges@start,
                       genes@ranges@start + genes@ranges@width - 1)

  genes <- GRanges(
    genes@seqnames,
    ranges = IRanges(gene_start,
                     width = 1),
    name = genes$symbol,
    gene_id = genes$gene_id,
    strand = genes@strand
  )

  seRNA <-
    SummarizedExperiment(assays = SimpleList(RNA = gene.mat),
                         rowRanges = genes)

  ## create object for ATAC data
  df_peak <-
    stringr::str_split_fixed(rownames(peak.mat), "-", 3)

  peakSet <- GRanges(df_peak[, 1],
                     IRanges(start = as.numeric(df_peak[, 2]),
                             end = as.numeric(df_peak[, 3])))

  seATAC <-
    SummarizedExperiment(assays = SimpleList(ATAC = peak.mat),
                         rowRanges = peakSet)

  ## find putative peak-to-gene
  o <-
    data.frame(findOverlaps(
      resize(seRNA, 2 * max.dist + 1, "center"),
      resize(rowRanges(seATAC), 1, "center"),
      ignore.strand = TRUE
    ))
  o$distance <- IRanges::distance(rowRanges(seRNA)[o[, 1]],
                                  rowRanges(seATAC)[o[, 2]])
  colnames(o) <- c("gene_idx", "peak_idx", "distance")

  df <- rowRanges(seATAC)[o$peak_idx, ]

  o$gene <- rowData(seRNA)[o$gene_idx, ]$name
  o$peak <- paste0(
    df@seqnames,
    "-",
    as.data.frame(df@ranges)$start,
    "-",
    as.data.frame(df@ranges)$end
  )


  ## compute correlation
  o$Correlation <- rowCorCpp(as.integer(o$peak_idx),
                             as.integer(o$gene_idx),
                             assay(seATAC),
                             assay(seRNA))

  ## compute p-value
  o$TStat <-
    (o$Correlation / sqrt((
      pmax(1 - o$Correlation ^ 2, 0.00000000000000001, na.rm = TRUE)
    ) / (ncol(seATAC) - 2))) #T-statistic P-value

  o$Pval <- 2 * pt(-abs(o$TStat), ncol(seATAC) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  o <- o[!is.na(o$FDR), ]

  return(o)

}


#######
#######

SelectTFs <- function(object,
                      tf.assay = "chromvar",
                      rna.assay = "RNA",
                      atac.assay= "ATAC",
                      trajectory.name = "Trajectory",
                      groupEvery = 1,
                      p.cutoff = 0.01,
                      cor.cutoff = 0.3,
                      return.heatmap = TRUE) {

    trajMM <- GetTrajectory(
    object,
    assay = tf.assay,
    trajectory.name=trajectory.name,
    groupEvery=groupEvery,
    slot = "data",
    smoothWindow = 7,
    log2Norm = FALSE
  )

  rownames(trajMM) <- object@assays[[atac.assay]]@motifs@motif.names

  trajRNA <- GetTrajectory(
    object,
    assay = rna.assay,
    trajectory.name=trajectory.name,
    groupEvery=groupEvery,
    slot = "data",
    smoothWindow = 7,
    log2Norm = TRUE
  )

rownames(trajRNA) = gene_name_ID$name[match(rownames(trajRNA), gene_name_ID$id)]


  df.cor <- GetCorrelation(trajMM, trajRNA)

  # we only select TFs that show significant correlation
  df.cor <- df.cor[df.cor$adj_p < p.cutoff & df.cor$correlation > cor.cutoff, ]

  matMM <- (TrajectoryHeatmap(
    trajMM,
    varCutOff = 0,
    pal = paletteContinuous(set = "solarExtra"),
    limits = c(-2, 2),
    name = "TF activity",
    returnMatrix = TRUE
  ))

  df_tf_time_point <- data.frame(tfs = rownames(matMM),
                                 time_point = seq(1, 100, length.out = nrow(matMM)))
  rownames(df_tf_time_point) <- df_tf_time_point$tfs

  df_tf_time_point <- df_tf_time_point[df.cor$tfs,]
  df.cor$time_point <- df_tf_time_point$time_point
  df.cor <- df.cor[order(df.cor$time_point),]

  trajMM <- trajMM[df.cor$tfs,]
  trajRNA <- trajRNA[df.cor$tfs,]


  if (return.heatmap) {
    ht <- (
      CorrelationHeatmap(
        trajectory1 = trajMM,
        trajectory2 = trajRNA,
        name1 = "TF activity",
        name2 = "Gene expression"
      )
    )

    res <- list("tfs" = df.cor, "heatmap" = ht)

  } else{
    res <- list("tfs" = df.cor)

  }

  return(res)

}



########
########


SelectGenes <- function(object,
                        atac.assay = "ATAC",
                        rna.assay = "RNA",
                        var.cutoff.gene = 0.9,
                        trajectory.name = "Trajectory",
                        distance.cutoff = 2000,
                        groupEvery = 1,
                        cor.cutoff = 0,
                        fdr.cutoff = 1e-04,
                        return.heatmap = TRUE,
                        labelTop1 = 10,
                        labelTop2 = 10,
                        genome = "Ath"
                        ) {
  trajRNA <- GetTrajectory(
    object,
    assay = rna.assay,
    trajectory.name = trajectory.name,
    groupEvery = groupEvery,
    slot = "data",
    smoothWindow = 7,
    log2Norm = FALSE
  )

  trajATAC <- GetTrajectory(
    object,
    assay = atac.assay,
    groupEvery = groupEvery,
    trajectory.name = trajectory.name,
    slot = "data",
    smoothWindow = 7,
    log2Norm = TRUE
  )

  # note here we only use the top 10% most variable genes
  groupMatRNA <- (
    TrajectoryHeatmap(
      trajRNA,
      varCutOff = var.cutoff.gene,
      pal = paletteContinuous(set = "horizonExtra"),
      limits = c(-2, 2),
      returnMatrix = TRUE
    )
  )

  groupMatATAC <- (
    TrajectoryHeatmap(
      trajATAC,
      varCutOff = 0,
      maxFeatures = nrow(trajATAC),
      pal = paletteContinuous(set = "solarExtra"),
      limits = c(-2, 2),
      name = "Chromatin accessibility",
      returnMatrix = TRUE
    )
  )

  rownames(trajRNA) = gene_name_ID$name[match(rownames(trajRNA), gene_name_ID$id)]

  message("Linking cis-regulatory elements to genes...")
  df.p2g <- PeakToGene(peak.mat = groupMatATAC,
                       gene.mat = groupMatRNA,
                       genome = genome)

  df.p2g <- df.p2g %>%
    subset(distance > distance.cutoff) %>%
    subset(Correlation > cor.cutoff & FDR < fdr.cutoff)

  trajATAC <- trajATAC[df.p2g$peak,]
  trajRNA <- trajRNA[df.p2g$gene,]

  if (return.heatmap) {
    ht <- (
      CorrelationHeatmap(
        trajectory1 = trajATAC,
        trajectory2 = trajRNA,
        name1 = "Chromatin accessibility",
        name2 = "Gene expression",
        labelTop1 = labelTop1,
        labelTop2 = labelTop2,
        labelRows1 = FALSE,
        labelRows2 = FALSE
      )
    )

    res <- list("p2g" = df.p2g, "heatmap" = ht)

  } else{
    res <- list("p2g" = df.p2g)
  }


  return(res)

}



#######
#######


SelectGenes_custom <- function(object,
                        atac.assay = "ATAC",
                        rna.assay = "RNA",
                        var.cutoff.gene = 0.9,
                        trajectory.name = "Trajectory",
                        distance.cutoff = 2000,
                        groupEvery = 1,
                        cor.cutoff = 0,
                        fdr.cutoff = 1e-04,
                        return.heatmap = TRUE,
                        labelTop1 = 10,
                        labelTop2 = 10,
                        genome = "Ath"
                        ) {
  trajRNA <- GetTrajectory(
    object,
    assay = rna.assay,
    trajectory.name = trajectory.name,
    groupEvery = groupEvery,
    slot = "data",
    smoothWindow = 7,
    log2Norm = FALSE
  )

  trajATAC <- GetTrajectory(
    object,
    assay = atac.assay,
    groupEvery = groupEvery,
    trajectory.name = trajectory.name,
    slot = "data",
    smoothWindow = 7,
    log2Norm = TRUE
  )

  # note here we only use the top 10% most variable genes
  groupMatRNA <- (
    TrajectoryHeatmap(
      trajRNA,
      varCutOff = var.cutoff.gene,
      pal = paletteContinuous(set = "horizonExtra"),
      limits = c(-2, 2),
      returnMatrix = TRUE
    )
  )

  groupMatATAC <- (
    TrajectoryHeatmap(
      trajATAC,
      varCutOff = 0,
      maxFeatures = nrow(trajATAC),
      pal = paletteContinuous(set = "solarExtra"),
      limits = c(-2, 2),
      name = "Chromatin accessibility",
      returnMatrix = TRUE
    )
  )

  rownames(trajRNA) = gene_name_ID$name[match(rownames(trajRNA), gene_name_ID$id)]

  message("Linking cis-regulatory elements to genes...")
  df.p2g <- PeakToGene(peak.mat = groupMatATAC,
                       gene.mat = groupMatRNA,
                       genome = genome)

  df.p2g <- df.p2g %>%
    subset(distance > distance.cutoff) %>%
    subset(Correlation > cor.cutoff & FDR < fdr.cutoff)

  trajATAC <- trajATAC[df.p2g$peak,]
  trajRNA <- trajRNA[df.p2g$gene,]

  if (return.heatmap) {
    ht <- (
      CorrelationHeatmap_custom(
        trajectory1 = trajATAC,
        trajectory2 = trajRNA,
        name1 = "Chromatin accessibility",
        name2 = "Gene expression",
        labelTop1 = labelTop1,
        labelTop2 = labelTop2,
        labelRows1 = FALSE,
        labelRows2 = FALSE
      )
    )

    res <- list("p2g" = df.p2g, "heatmap" = ht)

  } else{
    res <- list("p2g" = df.p2g)
  }


  return(res)

}

########
########


GetTFGeneCorrelation <- function(object,
                                 tf.use = NULL,
                                 gene.use = NULL,
                                 tf.assay = "chromvar",
                                 gene.assay = "RNA",
                                 atac.assay="ATAC",
                                 trajectory.name = "Trajectory",
                                 groupEvery=1) {
  ## get tf activity and gene expression along trajectory
  trajMM <- GetTrajectory(
    object,
    assay = tf.assay,
    slot = "data",
    trajectory.name = trajectory.name,
    groupEvery=groupEvery,
    smoothWindow = 7,
    log2Norm = FALSE
  )

  trajRNA <- GetTrajectory(
    object,
    assay = gene.assay,
    slot = "data",
    trajectory.name = trajectory.name,
    groupEvery=groupEvery,
    smoothWindow = 7,
    log2Norm = FALSE
  )

  rownames(trajRNA) = gene_name_ID$name[match(rownames(trajRNA), gene_name_ID$id)]


  rownames(trajMM) <- object@assays[[atac.assay]]@motifs@motif.names

  tf_activity <- (
    TrajectoryHeatmap(
      trajMM,
      varCutOff = 0,
      pal = paletteContinuous(set = "solarExtra"),
      limits = c(-2, 2),
      name = "TF activity",
      returnMatrix = TRUE
    )
  )

  gene_expression <- (
    TrajectoryHeatmap(
      trajRNA,
      varCutOff = 0.9,
      pal = paletteContinuous(set = "solarExtra"),
      limits = c(-2, 2),
      name = "Gene expression",
      returnMatrix = TRUE
    )
  )

  ## here we filter the TFs according to our correlation analysis
  if (!is.null(tf.use)) {
    tf_activity <- tf_activity[tf.use, ]
  }

  ## here we filter the genes by only considering genes that are linked to peaks
  if (!is.null(gene.use)) {
    sel_genes <- intersect(rownames(gene_expression), gene.use)

    ## subset the gene expression matrix
    gene_expression <- gene_expression[sel_genes, ]
  }

  ## compute the correlation of TF activity and gene expression along the trajectory
  ## df.cor -> gene by TF matrix
  df.cor <- t(cor(t(tf_activity), t(gene_expression))) %>%
    as.data.frame()

  if (!is.null(tf.use)) {
    df.cor <- df.cor[, tf.use]
  }

  df.cor$gene <- rownames(df.cor)

  df.cor <- df.cor %>%
    tidyr::pivot_longer(!gene, names_to = "tf", values_to = "correlation") %>%
    dplyr::select(c(tf, gene, correlation))

  df.cor$t_stat <-
    (df.cor$correlation / sqrt((
      pmax(1 - df.cor$correlation ^ 2, 0.00000000000000001, na.rm = TRUE)
    ) / (ncol(tf_activity) - 2))) #T-statistic P-value

  df.cor$p_value <-
    2 * pt(-abs(df.cor$t_stat), ncol(tf_activity) - 2)
  df.cor$fdr <- p.adjust(df.cor$p_value, method = "fdr")

  return(df.cor)

}

########
########


GRNPlot <- function(df.grn,
                    tfs.use = NULL,
                    show.tf.labels = TRUE,
                    tfs.timepoint = NULL,
                    genes.cluster = NULL,
                    genes.use = NULL,
                    genes.highlight = NULL,
                    cols.highlight = "#984ea3",
                    seed = 42,
                    plot.importance = TRUE,
                    min.importance = 2,
                    remove.isolated = FALSE) {
  if (is.null(tfs.timepoint)) {
    stop("Need time point for each TF!")
  }

  if (!is.null(tfs.use)){
    df.grn <- subset(df.grn, tf %in% tfs.use)

  }
  if (!is.null(genes.use)){
    df.grn <- subset(df.grn, gene %in% genes.use)
  }

  tf.list <- unique(df.grn$tf)
  gene.list <- setdiff(unique(df.grn$gene), tf.list)

  # create graph from data frame
  g <- igraph::graph_from_data_frame(df.grn, directed = TRUE)

  # remove the isolated if indicated
  if(remove.isolated){
    isolated <- which(degree(g)==0)
    g <- igraph::delete.vertices(g, isolated)
  }

  # compute pagerank and betweenness
  pagerank <- page_rank(g, weights = E(g)$weights)
  bet <-
    betweenness(g,
                weights = E(g)$weights,
                normalized = TRUE)

  df_measure <- data.frame(
    tf = V(g)$name,
    pagerank = pagerank$vector,
    betweenness = bet
  ) %>%
    subset(tf %in% df.grn$tf) %>%
    mutate(pagerank = scale(pagerank)[, 1]) %>%
    mutate(betweenness = scale(betweenness)[, 1])

  # compute importance only for TFs based on centrality and betweenness
  min.page <- min(df_measure$pagerank)
  min.bet <- min(df_measure$betweenness)
  df_measure$importance <-
    sqrt((df_measure$pagerank - min.page) ** 2 +
           (df_measure$betweenness - min.bet) ** 2)

  if (plot.importance) {
    p <- ggplot(data = df_measure) + aes(x = reorder(tf, -importance),
                                         y = importance) +
      geom_point() +
      xlab("TFs") + ylab("Importance") +
      cowplot::theme_cowplot() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))

    print(p)
  }

  df_measure_sub <- subset(df_measure, importance > 2)

  # assign size to each node
  # for TFs, the size is proportional to the importance
  tf_size <- df_measure$importance
  names(tf_size) <- df_measure$tf

  ## for genes, we use the minimum size of TFs
  gene_size <-
    rep(min(df_measure$importance), length(unique(df.grn$gene)))
  names(gene_size) <- gene.list
  v_size <- c(tf_size, gene_size)
  V(g)$size <- v_size[V(g)$name]

  # assign color to each node
  ## TFs are colored by pseudotime point
  cols.tf <- ArchR::paletteContinuous(set = "blueYellow",
                                      n = length(tfs.timepoint))
  names(cols.tf) <- names(tfs.timepoint)

  ## genes are colored based on the clustering
  if (is.null(genes.cluster)) {
    cols.gene <- rep("gray", length(gene.list))
    names(cols.gene) <- gene.list
  } else{
    genes.cluster <- genes.cluster %>%
      subset(gene %in% gene.list)

    cols <-
      ArchR::paletteDiscrete(values = as.character(genes.cluster$cluster))

    df.gene <- lapply(1:length(cols), function(x) {
      df <- subset(genes.cluster, cluster == x)
      df$color <- rep(cols[[x]], nrow(df))
      return(df)

    }) %>% Reduce(rbind, .)

    cols.gene <- df.gene$color
    names(cols.gene) <- df.gene$gene
  }
  v_color <- c(cols.tf, cols.gene)
  v_color <- v_color[V(g)$name]

  ## assign alpha
  tf_alpha <- rep(1, length(tf.list))
  gene_alpha <- rep(0.5, length(gene.list))
  names(tf_alpha) <- tf.list
  names(gene_alpha) <- gene.list
  v_alpha <- c(tf.list, gene.list)
  V(g)$alpha <- v_alpha[V(g)$name]

  # compute layout
  set.seed(seed)
  layout <- layout_with_fr(
    g,
    weights = E(g)$weights,
    dim = 2,
    niter = 1000
  )

  p <- ggraph(g, layout = layout) +
    geom_edge_link(edge_colour = "gray", edge_alpha = 0.25) +
    geom_node_point(aes(
      size = V(g)$size,
      color = as.factor(name),
      alpha = V(g)$alpha
    ),
    show.legend = FALSE) +
    scale_size(range = c(1, 10)) +
    scale_color_manual(values = v_color)

  if(show.tf.labels){
    p <- p +geom_node_label(
      aes(
        #filter = V(g)$name %in% df_measure_sub$tf,
        filter = V(g)$name %in% tf.list,
        label = V(g)$name
      ),
      repel = TRUE,
      hjust = "inward",
      color = "#ff7f00",
      size = 5,
      show.legend = FALSE,
      max.overlaps = Inf
    )
  }


  # highlight some genes
  if (!is.null(genes.highlight)) {
    p <-
      p + geom_node_label(
        aes(
          filter = V(g)$name %in% genes.highlight,
          label = V(g)$name
        ),
        repel = TRUE,
        hjust = "inward",
        size = 5,
        color = cols.highlight,
        show.legend = FALSE
      )

  }

  p <- p + theme_void()

  return(p)
}

########
########

CorrelationHeatmap_custom <- function(trajectory1,
                               trajectory2,
                               name1 = NULL,
                               name2 = NULL,
                               labelRows1 = TRUE,
                               labelRows2 = TRUE,
                               labelTop1 = 100,
                               labelTop2 =100,
                               limits1 = c(-2, 2),
                               limits2 = c(-2, 2)) {
  trajCombined <- trajectory1

  assay(trajCombined, withDimnames = FALSE) <-
    t(apply(assay(trajectory2), 1, scale)) +
    t(apply(assay(trajectory1), 1, scale))

  combinedMat <- TrajectoryHeatmap(trajCombined,
                                   returnMatrix = TRUE,
                                   varCutOff = 0)

  rowOrder <- match(rownames(combinedMat), rownames(trajectory1))

  ht1 <- TrajectoryHeatmap(
    trajectory1,
    pal = huntrix_pal,
#     paletteContinuous(set = "solarExtra"),
    varCutOff = 0,
    maxFeatures = nrow(trajectory1),
    rowOrder = rowOrder,
    limits = limits1,
    labelRows = labelRows1,
    labelTop = labelTop1,
    name = name1
  )

  ht2 <- TrajectoryHeatmap(
    trajectory2,
    pal = paletteContinuous(set = "horizonExtra"),
    varCutOff = 0,
    maxFeatures = nrow(trajectory2),
    rowOrder = rowOrder,
    limits = limits2,
    labelRows = labelRows2,
    labelTop = labelTop2,
    name = name2
  )

  ht <- ht1 + ht2

  return(ht)

}

########
########

SelectTFs_custom <- function(object,
                      tf.assay = "chromvar",
                      rna.assay = "RNA",
                      atac.assay= "ATAC",
                      trajectory.name = "Trajectory",
                      groupEvery = 1,
                      p.cutoff = 0.01,
                      cor.cutoff = 0.3,
                      return.heatmap = TRUE) {

    trajMM <- GetTrajectory(
    object,
    assay = tf.assay,
    trajectory.name=trajectory.name,
    groupEvery=groupEvery,
    slot = "data",
    smoothWindow = 7,
    log2Norm = FALSE
  )

  rownames(trajMM) <- object@assays[[atac.assay]]@motifs@motif.names

  trajRNA <- GetTrajectory(
    object,
    assay = rna.assay,
    trajectory.name=trajectory.name,
    groupEvery=groupEvery,
    slot = "data",
    smoothWindow = 7,
    log2Norm = TRUE
  )

rownames(trajRNA) = gene_name_ID$name[match(rownames(trajRNA), gene_name_ID$id)]


  df.cor <- GetCorrelation(trajMM, trajRNA)

  # we only select TFs that show significant correlation
  df.cor <- df.cor[df.cor$adj_p < p.cutoff & df.cor$correlation > cor.cutoff, ]

  matMM <- (TrajectoryHeatmap(
    trajMM,
    varCutOff = 0,
    pal = paletteContinuous(set = "solarExtra"),
    limits = c(-2, 2),
    name = "TF activity",
    returnMatrix = TRUE
  ))

  df_tf_time_point <- data.frame(tfs = rownames(matMM),
                                 time_point = seq(1, 100, length.out = nrow(matMM)))
  rownames(df_tf_time_point) <- df_tf_time_point$tfs

  df_tf_time_point <- df_tf_time_point[df.cor$tfs,]
  df.cor$time_point <- df_tf_time_point$time_point
  df.cor <- df.cor[order(df.cor$time_point),]

  trajMM <- trajMM[df.cor$tfs,]
  trajRNA <- trajRNA[df.cor$tfs,]


  if (return.heatmap) {
    ht <- (
      CorrelationHeatmap_custom(
        trajectory1 = trajMM,
        trajectory2 = trajRNA,
        name1 = "TF activity",
        name2 = "Gene expression"
      )
    )

    res <- list("tfs" = df.cor, "heatmap" = ht)

  } else{
    res <- list("tfs" = df.cor)

  }

  return(res)

}

########
########

GRNPlot_custom <- function(df.grn,
                    tfs.use = NULL,
                    show.tf.labels = TRUE,
                    tfs.timepoint = NULL,
                    genes.cluster = NULL,
                    genes.use = NULL,
                    genes.highlight = NULL,
                    cols.highlight = "#984ea3",
                    seed = 42,
                    plot.importance = TRUE,
                    min.importance = 2,
                    remove.isolated = FALSE) {
  if (is.null(tfs.timepoint)) {
    stop("Need time point for each TF!")
  }

  if (!is.null(tfs.use)){
    df.grn <- subset(df.grn, tf %in% tfs.use)

  }
  if (!is.null(genes.use)){
    df.grn <- subset(df.grn, gene %in% genes.use)
  }

  tf.list <- unique(df.grn$tf)
  gene.list <- setdiff(unique(df.grn$gene), tf.list)

  # create graph from data frame
  g <- igraph::graph_from_data_frame(df.grn, directed = TRUE)

  # remove the isolated if indicated
  if(remove.isolated){
    isolated <- which(degree(g)==0)
    g <- igraph::delete.vertices(g, isolated)
  }

  # compute pagerank and betweenness
  pagerank <- page_rank(g, weights = E(g)$weights)
  bet <-
    betweenness(g,
                weights = E(g)$weights,
                normalized = TRUE)

  df_measure <- data.frame(
    tf = V(g)$name,
    pagerank = pagerank$vector,
    betweenness = bet
  ) %>%
    subset(tf %in% df.grn$tf) %>%
    mutate(pagerank = scale(pagerank)[, 1]) %>%
    mutate(betweenness = scale(betweenness)[, 1])

  # compute importance only for TFs based on centrality and betweenness
  min.page <- min(df_measure$pagerank)
  min.bet <- min(df_measure$betweenness)
  df_measure$importance <-
    sqrt((df_measure$pagerank - min.page) ** 2 +
           (df_measure$betweenness - min.bet) ** 2)

  if (plot.importance) {
    p <- ggplot(data = df_measure) + aes(x = reorder(tf, -importance),
                                         y = importance) +
      geom_point() +
      xlab("TFs") + ylab("Importance") +
      cowplot::theme_cowplot() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))

    print(p)
  }

  df_measure_sub <- subset(df_measure, importance > 2)

  # assign size to each node
  # for TFs, the size is proportional to the importance
  tf_size <- df_measure$importance
  names(tf_size) <- df_measure$tf

  ## for genes, we use the minimum size of TFs
  gene_size <-
    rep(min(df_measure$importance), length(unique(df.grn$gene)))
  names(gene_size) <- gene.list
  v_size <- c(tf_size, gene_size)
  V(g)$size <- v_size[V(g)$name]

  # assign color to each node
  ## TFs are colored by pseudotime point
#   cols.tf <- ArchR::paletteContinuous(set = "blueYellow",
#                                       n = length(tfs.timepoint))
  cols.tf <- rumi_pal
  names(cols.tf) <- names(tfs.timepoint)

  ## genes are colored based on the clustering
  if (is.null(genes.cluster)) {
    cols.gene <- rep("gray", length(gene.list))
    names(cols.gene) <- gene.list
  } else{
    genes.cluster <- genes.cluster %>%
      subset(gene %in% gene.list)

    cols <-
      ArchR::paletteDiscrete(values = as.character(genes.cluster$cluster))

    df.gene <- lapply(1:length(cols), function(x) {
      df <- subset(genes.cluster, cluster == x)
      df$color <- rep(cols[[x]], nrow(df))
      return(df)

    }) %>% Reduce(rbind, .)

    cols.gene <- df.gene$color
    names(cols.gene) <- df.gene$gene
  }
  v_color <- c(cols.tf, cols.gene)
  v_color <- v_color[V(g)$name]

  ## assign alpha
  tf_alpha <- rep(1, length(tf.list))
  gene_alpha <- rep(0.5, length(gene.list))
  names(tf_alpha) <- tf.list
  names(gene_alpha) <- gene.list
  v_alpha <- c(tf.list, gene.list)
  V(g)$alpha <- v_alpha[V(g)$name]

  # compute layout
  set.seed(seed)
  layout <- layout_with_fr(
    g,
    weights = E(g)$weights,
    dim = 2,
    niter = 1000
  )

  p <- ggraph(g, layout = layout) +
    geom_edge_link(edge_colour = "gray", edge_alpha = 0.25) +
    geom_node_point(aes(
      size = V(g)$size,
      color = as.factor(name),
      alpha = V(g)$alpha
    ),
    show.legend = FALSE) +
    scale_size(range = c(1, 10)) +
    scale_color_manual(values = v_color)

  if(show.tf.labels){
    p <- p +geom_node_label(
      aes(
        #filter = V(g)$name %in% df_measure_sub$tf,
        filter = V(g)$name %in% tf.list,
        label = V(g)$name
      ),
      repel = TRUE,
      hjust = "inward",
      color = "#ff7f00",
      size = 5,
      show.legend = FALSE,
      max.overlaps = Inf
    )
  }


  # highlight some genes
  if (!is.null(genes.highlight)) {
    p <-
      p + geom_node_label(
        aes(
          filter = V(g)$name %in% genes.highlight,
          label = V(g)$name
        ),
        repel = TRUE,
        hjust = "inward",
        size = 5,
        color = cols.highlight,
        show.legend = FALSE
      )

  }

  p <- p + theme_void()

  return(p)
}



###### begin analysis


geneAnnotation = createGeneAnnotation(TxDb = TxDb.Athaliana.BioMart.plantsmart51, OrgDb = org.At.tair.db, annoStyle = "TAIR")

geneAnnotation$genes <- renameSeqlevels(geneAnnotation$genes, c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrM", "ChrC"))
geneAnnotation$exons <- renameSeqlevels(geneAnnotation$exons, c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrM", "ChrC"))
geneAnnotation$TSS <- renameSeqlevels(geneAnnotation$TSS, c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrM", "ChrC"))

geneAnno <- geneAnnotation$gene
exonAnno <- geneAnnotation$exon
TSSAnno <- geneAnnotation$TSS

geneAnno$symbol = gsub("NA_", "", geneAnno$symbol)
exonAnno$symbol = gsub("NA_", "", exonAnno$symbol)


geneAnnoAth <- createGeneAnnotation(
genome = NULL,
TxDb = NULL,
OrgDb = NULL,
genes = geneAnno,
exons = exonAnno,
TSS = TSSAnno,
annoStyle = "TAIR"
)

geneAnnoAth = GRangesList(genes = geneAnno, exons = exonAnno, TSS = TSSAnno, compress=FALSE)


gene_name_ID = data.frame(name = geneAnnoAth$genes$symbol, id = geneAnnoAth$genes$gene_id)


jaspar <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))

pfm <- TFBSTools::getMatrixSet(sq24, list(collection = "CORE", tax_group = 'plants', species = "3702", all_versions = FALSE))


#### manually add BBM motif from Horstman et al., 2017

ANT_matrix = Matrix(pfm[names(pfm) %in% "MA0571.2"])$MA0571.2

ANT_matrix[,1] = c(6,6,16,6)
ANT_matrix[,3] = c(25,3,3,3)
ANT_matrix[,5] = c(23,3,5,3)
ANT_matrix[,6] = c(9,8,9,8)
ANT_matrix[,7] = c(9,8,9,8)
ANT_matrix[,8] = c(9,8,9,8)
ANT_matrix[,9] = c(9,8,9,8)
ANT_matrix[,10] = c(4,16,4,10)
ANT_matrix[,11] = c(0,34,0,0)
ANT_matrix[,12] = c(9,8,9,8)
ANT_matrix[,13] = c(34,0,0,0)

BBM_matrix = ANT_matrix[,1:13]

BBM = PFMatrix(ID = "BBM", name = "BBM", matrixClass = "Unknown",
  strand = "+", bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
  tags = list(remap_tf_name = "BBM"), profileMatrix = BBM_matrix)

# pdf("BBM_motif.pdf")
# seqLogo::seqLogo(Matrix(toPWM(BBM, type="prob")), ic.scale=TRUE)
# seqLogo::seqLogo(Matrix(toPWM(BBM, type="prob")),ic.scale=FALSE)
# 
# seqLogo(toICM(BBM), ic.scale=TRUE)
# seqLogo(toICM(BBM), ic.scale=FALSE)
# dev.off()


### use GATA9 motif as a proxy for GATA2

GATA9_matrix = Matrix(pfm[names(pfm) %in% "MA1018.2"])$MA1018.2

GATA2 = PFMatrix(ID = "GATA2", name = "GATA2", matrixClass = "Unknown",
  strand = "+", bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
  tags = list(remap_tf_name = "GATA2"), profileMatrix = GATA9_matrix)


BBM_list = PFMatrixList("BBM" = BBM)
GATA2_list = PFMatrixList("GATA2" = GATA2)

pfm = c(pfm, BBM_list, GATA2_list)

### load c26 multiome data
c26 = readRDS("/anvil/projects/x-mcb130189/trlee/ethylene/GA_subclust26_multiome_imputed_250825.rds")

### subset only the multiome cels
c26=subset(c26, assay=="multiome")

c26 <- AddTrajectory(object = sub, 
                          trajectory = c(4, 0, 2, 1, 5),
                          group.by = "seurat_clusters", 
                          reduction = "umap",
                          dims = 1:2, 
                          use.all = FALSE)

obj <- c26[, !is.na(c26$Trajectory)]


rumi_list = list("rumi" = c("1" = "#ca6dca", "2"="#5826a3", "3"="#ffa22b", "4"="#b6a364", "5"="#f50501"))

p <- TrajectoryPlot_custom(object = obj, 
                    reduction = "umap",
                    continuousSet = "blueYellow",
                    size = 1,
                   addArrow = FALSE)

pdf("Fig4_UMAP_trajectory_recolored.pdf")
p
dev.off()



#### re-run chromvar on subset object
obj <- RegionStats(obj, genome = BSgenome.Athaliana.TAIR.TAIR9, assay="peaks_cluster")


obj <- AddMotifs(
  object = obj,
  genome = BSgenome.Athaliana.TAIR.TAIR9,
  pfm = pfm,
  assay="peaks_cluster"
)


obj <- RunChromVAR(
  object = obj,
  genome = BSgenome.Athaliana.TAIR.TAIR9,
    assay = "peaks_cluster"
)


DefaultAssay(obj) = "RNA"

obj = NormalizeData(obj)
obj = ScaleData(obj)
obj = JoinLayers(obj, layers=c("data1.2", "data.2"))




### identify TF regulators that positively correlate promoter chromatin accessibility and expression level

huntrix = c("#252157", "#6ea0b4", "#8bd282", "#91beaa", "#d29f41", "#fadc82",  "#aa82c8", "#cc1be3")
huntrix_pal = colorRampPalette(huntrix)(n = 256)

res <- SelectTFs_custom(object = obj, return.heatmap = TRUE, atac.assay="peaks_cluster", rna.assay="alra", p.cutoff = 0.05)

df.cor <- res$tfs
ht <- res$heatmap

pdf("Fig4_TFs_ATAC_RNA_correlation_heatmap_recolor.pdf", height=12)
draw(ht)
dev.off()


#### identify all putative target genes with positively correlated promoter chromatin accessibility and gene expression
res <- SelectGenes(object = obj,
                  labelTop1 = 50,
                  labelTop2 = 50,
                  atac.assay = "peaks_cluster",
                  rna.assay = "alra",
                  genome="Ath",
                  return.heatmap=TRUE,
                  distance.cutoff = 0)



df.p2g <- res$p2g
ht <- res$heatmap

write.table(df.p2g, "Fig4_ATAC_RNA_correlation_heatmap_table.xls", sep="\t", row.names=TRUE)


####### correlation heatmap


### identify TF - target gene interactions
tf.gene.cor <- GetTFGeneCorrelation(object = obj, 
                                    tf.use = df.cor$tfs, 
                                    gene.use = unique(df.p2g$gene),
                                    tf.assay = "chromvar", 
                                    gene.assay = "RNA",
                                    atac.assay = "peaks_cluster",
                                    trajectory.name = "Trajectory")



ht <- GRNHeatmap(tf.gene.cor, 
                 tf.timepoint = df.cor$time_point)

pdf("Fig4_GRN_ATAC_RNA_correlation_heatmap.pdf", width=12)
ht
dev.off()



motif.matching <- obj@assays$peaks_cluster@motifs@data
colnames(motif.matching) <- obj@assays$peaks_cluster@motifs@motif.names
motif.matching <- motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]


df.grn <- GetGRN(motif.matching = motif.matching, 
                 df.cor = tf.gene.cor, 
                 df.p2g = df.p2g)


# define colors for nodes representing TFs (i.e., regulators)
df.cor <- df.cor[order(df.cor$time_point), ]
tfs.timepoint <- df.cor$time_point
names(tfs.timepoint) <- df.cor$tfs

# plot the graph, here we can highlight some genes
df.grn2 <- df.grn %>%
    subset(correlation > 0.4) %>%
    dplyr::select(c(tf, gene, correlation)) %>%
    dplyr::rename(weights = correlation)


saveRDS(df.grn2,"final_grn.Rds")    




pdf("Fig4_GRN_network_full.pdf", width=12)

p <- GRNPlot(df.grn, 
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             seed = 42, 
             plot.importance = FALSE,
            min.importance = 2,
            remove.isolated = FALSE,
            genes.highlight=list("GATA2", "ARF6", "HDG12", "BBM", "HB-3", "ROT3"))

print(p)
dev.off()




GA_genes = c("GASA4")

selected_TFs = c("GATA10", "GATA11", "GATA6", "GATA12", "LBD13", "DOF3.4", "DF1", "AT5G47660", "AIL6", "BIM1", "ERF035", "AP3", "ERF104")

other_genes = c("ARF6", "ACO3", "GATA2")
highlight_genes = c(GA_genes, selected_TFs, other_genes)

highlight_colors = c(rep("magenta", length(GA_genes)), rep("purple", length(selected_TFs)), rep("orange", length(other_genes)))

pdf("Fig4_GRN_network_subset_labeled.pdf")

p <- GRNPlot(df.grn2, 
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = FALSE,
             seed = 42, 
             plot.importance = FALSE,
            min.importance = 2,
            remove.isolated = FALSE,
            genes.highlight=highlight_genes)

print(p)
dev.off()



write.table(df.grn, "Fig4_ATAC_RNA_GRN_full_table.xls", sep="\t", row.names=FALSE)
write.table(df.grn2, "Fig4_ATAC_RNA_GRN_subset_table.xls", sep="\t", row.names=FALSE)




##### GRN of GATA2 regulators and targets

GATA2_regulators = c("AIL6", "AP3", "AT5G05090", "AT5G47660", "ATHB-12", "ATHB-7", "ATMYB31", "AZF2", "BBM",  "CAMTA2", "CAMTA3", "DF1", "DOF3.4", "DOF5.6", "ERF018", "ERF035", "ERF104", "ERF115", "GATA11", "GATA12", "GATA2", "GATA6", "GBF3", "HAT1", "HAT2", "HYH", "JKD", "LBD13", "LBD18", "LHY", "MYB39", "MYC4", "MYR2", "NID1", "PHL2", "RVE1", "RVE5", "RVE6", "SPL10", "TCP9", "TCX6", "TGA1", "TGA3", "TGA7", "WRKY26", "WRKY3")

GATA2_targets = c("AGP1", "ANP3", "APOK3", "AR192", "AT1G01300", "AT1G06930", "AT1G12064", "AT1G12070", "AT1G14450", "AT1G49700", "AT1G54740", "AT1G56020", "AT1G67148", "AT1G74870", "AT1G75335", "AT1G75360", "AT2G21160", "AT2G30840", "AT2G34185", "AT2G42960", "AT3G06035", "AT3G07320", "AT3G16900", "AT3G55690", "AT3G56870", "AT4G08910", "AT4G22230", "AT4G38540", "AT5G07340", "AT5G08180", "AT5G15220", "AT5G28920", "AT5G56120", "ATFRO1", "ATGIP1", "ATGPAT4", "ATICS1", "ATMCM2", "ATMYB23", "ATWRKY48", "AtHB28", "AthA1-1", "CHR27", "CYCD4;1", "EDA14", "EFD", "GASA4", "GATA2", "GGL28", "ICK1", "IQD29", "IQD31", "KAN2", "KOM", "LAZY5", "MUB5", "PAI1", "PAP10", "SAUR17", "SAUR79", "SHM4", "cycp3;1", "eL29y", "eL34x", "eL41z", "mMDH1", "uS4z")

GATA2_regulators_grn = df.grn[df.grn$gene %in% GATA2_regulators,]
GATA2_regulators_grn = GATA2_regulators_grn[GATA2_regulators_grn$gene == "GATA2",]
GATA2_regulators_grn <- GATA2_regulators_grn %>%
    subset(correlation > 0.2) %>%
    dplyr::select(c(tf, gene, correlation)) %>%
    dplyr::rename(weights = correlation)



GATA2_targets_grn = df.grn[df.grn$tf %in% GATA2_targets,]
GATA2_targets_grn <- GATA2_targets_grn %>%
    subset(correlation > 0.2) %>%
    dplyr::select(c(tf, gene, correlation)) %>%
    dplyr::rename(weights = correlation)


rumi = c("#ca6dca", "#5826a3", "#ffa22b", "#b6a364", "#f50501")
rumi_pal = colorRampPalette(rumi)(n = length(df.cor$time_point))



GATA2_network = rbind(GATA2_regulators_grn, GATA2_targets_grn)

p <- GRNPlot_custom(GATA2_network, 
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             seed = 40, 
             plot.importance = FALSE,
            min.importance = 0,
            remove.isolated = FALSE,
            genes.highlight=GATA2_targets)


pdf("Fig4_GRN_GATA2_network.pdf", height=4, width=4)
print(p)
dev.off()




