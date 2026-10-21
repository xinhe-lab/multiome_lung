
source("code/utils_mapgen.R")

#' @title Make gene track plot using Gviz
#'
#' @param finemapstats A GRanges object of processed finemapping summary statistics
#' @param region The genomic region to visualize in the format of "chr:start-end".
#' @param gene.annots A GRanges object of gene annotations
#' @param bigSNP A `bigsnpr` object attached via bigsnpr::snp_attach()
#' @param txdb A `txdb` object of gene annotations
#' @param counts A list of counts data to display as quantitativ data tracks
#' @param peaks A list of peaks to display as binary data tracks
#' @param loops A list of chromatin loops, e.g. PC-HiC, ABC, etc.
#' @param genome Genome assembly version (default: "hg19").
#' @param filter_loop_genes A vector of gene names. Only show loops connected to the genes.
#' @param filter_loop_snps A vector of SNP IDs. Only show loops connected to the SNPs.
#' @param filter_protein_coding_genes Logical. If TRUE, only shows protein coding gene.
#' @param color_piptrack_by Color SNPs in the PIP track by
#' `locus`, `cs` (credible sets), or `none` (same color).
#' @param highlight_snps SNPs (rsIDs) to highlight. Or highlight top SNP ('topSNP').
#' @param counts.ylim ylim range (default: between 0 and 1) for the `counts` tracks
#' @param counts.color Colors for the `counts` tracks
#' @param peaks.color Colors for the `peak` tracks
#' @param highlight.color Colors for the highlighted SNPs
#' @param genelabel.side Side to put gene labels,
#' options are: 'above' (default), 'below', 'left', 'right'.
#' @param track.sizes Sizes of the tracks
#' @param rotation.title The rotation angle for the text in the track title
#' @param verbose If TRUE, print detail messages for plotting
#' @import GenomicRanges
#' @import tidyverse
#' @export
track_plot <- function(finemapstats,
                       region,
                       gene.annots,
                       rna.gr,
                       rna.counts,
                       genes = NULL,
                       bigSNP = NULL,
                       txdb  = NULL,
                       counts  = NULL,
                       peaks  = NULL,
                       loops  = NULL,
                       genome = 'hg19',
                       filter_loop_genes = NULL,
                       filter_loop_snps = NULL,
                       filter_protein_coding_genes = TRUE,
                       color_piptrack_by = c('locus', 'cs', 'none'),
                       highlight_snps = NULL,
                       counts.ylim = c(0,1),
                       counts.color,
                       peaks.color = 'navy',
                       rna.color = NULL,
                       loops.color = 'gray',
                       highlight.color = "orange",
                       genelabel.side = c('above', 'below', 'left', 'right'),
                       track.sizes = NULL,
                       rotation.title = 0,
                       verbose = FALSE,
                       ...) {

  genelabel.side <- match.arg(genelabel.side)
  color_piptrack_by <- match.arg(color_piptrack_by)

  seqlevelsStyle(finemapstats) <- 'UCSC'

  if( min(finemapstats$pval) >=0 && max(finemapstats$pval) <= 1 ){
    if(verbose){
      cat('Convert GWAS p-value to -log10(pvalue). \n')
    }
    finemapstats$pval <- -log10(finemapstats$pval)
  }

  # Limit to genomic region to visualize
  region <- as(region, 'GRanges')
  seqlevelsStyle(region) <- 'UCSC'

  rows.in <- which((as.character(seqnames(finemapstats)) == as.character(seqnames(region))) &
                     (finemapstats$pos >= start(region)) &
                     (finemapstats$pos <= end(region)))

  curr.finemapstats <- as.data.frame(finemapstats[rows.in,  ])

  # p-value track
  if(!is.null(bigSNP)){
    # Color SNPs by r2
    r2.breaks <- c(0, 0.1, 0.25, 0.75, 0.9, 1)
    r2.labels <- c('0-0.1','0.1-0.25','0.25-0.75','0.75-0.9','0.9-1')
    r2.colors <- c('black','blue','green','orange','red')
    names(r2.colors) <- r2.labels

    curr.finemapstats <- add_LD_bigSNP(curr.finemapstats, bigSNP, r2.breaks, r2.labels)
    if(verbose){ cat(nrow(curr.finemapstats), 'snps included.\n')}

    pval.df <- curr.finemapstats %>% dplyr::select(chr, pos, pval, r2) %>%
      dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
    pval.df <- pval.df %>% tidyr::pivot_wider(names_from = r2, names_sort = TRUE, values_from = 'pval')
    pval.gr <- makeGRangesFromDataFrame(pval.df, keep.extra.columns = T)
    seqlevelsStyle(pval.gr) <- 'UCSC'
    avail.r2.groups <- names(mcols(pval.gr))

    pval.track <- Gviz::DataTrack(range = pval.gr,
                                  genome = genome,
                                  groups = avail.r2.groups,
                                  col = r2.colors[avail.r2.groups],
                                  name = '-log10 P',
                                  legend = TRUE,
                                  box.legend = FALSE)
  }else{
    if(verbose){ cat(nrow(curr.finemapstats), 'snps included.\n')}
    pval.df <- curr.finemapstats %>% dplyr::select(chr, pos, pval) %>%
      dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
    pval.gr <- makeGRangesFromDataFrame(pval.df, keep.extra.columns = T)
    seqlevelsStyle(pval.gr) <- 'UCSC'

    pval.track <- Gviz::DataTrack(range = pval.gr,
                                  genome = genome,
                                  name = '-log10 P',
                                  col = 'black',
                                  legend = FALSE)
  }

  dpars.pval <- list(col.title = 'black',
                     col.axis = 'black',
                     background.title = "transparent",
                     #col.border.title = 'lightgray',
                     #col.frame = 'lightgray',
                     # rotation.title = rotation.title,
                     frame = TRUE,
                     cex.axis = 0.8,#1.2,
                     cex.legend = 0.8)
  Gviz::displayPars(pval.track) <- dpars.pval

  # PIP track
  if(color_piptrack_by == 'locus'){
    if(verbose){ cat('Color SNPs in PIP track by loci. \n')}
    pip.df <- curr.finemapstats %>% dplyr::select(chr, pos, pip, locus) %>%
      dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
    pip.df <- pip.df %>% tidyr::pivot_wider(names_from = locus, names_sort = TRUE, values_from = 'pip')
    pip.gr <- makeGRangesFromDataFrame(pip.df, keep.extra.columns = T)
    seqlevelsStyle(pip.gr) <- 'UCSC'
    avail.pip.groups <- names(mcols(pip.gr))
  }else if(color_piptrack_by == 'cs'){
    if(verbose){ cat('Color SNPs in PIP track by credible sets.\n')}
    pip.df <- curr.finemapstats %>% dplyr::select(chr, pos, pip, cs) %>%
      dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
    pip.df <- pip.df %>% tidyr::pivot_wider(names_from = cs, names_sort = TRUE, values_from = 'pip')
    pip.gr <- makeGRangesFromDataFrame(pip.df, keep.extra.columns = T)
    names(mcols(pip.gr)) <- paste0("CS", names(mcols(pip.gr)))
    seqlevelsStyle(pip.gr) <- 'UCSC'
    avail.pip.groups <- names(mcols(pip.gr))
    
    # Color SNPs by CS
    pip.breaks <- seq(0, 3, 1)
    pip.labels <- avail.pip.groups
    pip.colors <- c("#D3D3D3", "brown", "#0072B2", "#D55E00")
    names(pip.colors) <- pip.labels
    
  }else{
    pip.df <- curr.finemapstats %>% dplyr::select(chr, pos, pip) %>%
      dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
    pip.gr <- makeGRangesFromDataFrame(pip.df, keep.extra.columns = T)
    seqlevelsStyle(pip.gr) <- 'UCSC'
    avail.pip.groups <- NULL
  }

  pip.track <- Gviz::DataTrack(range = pip.gr,
                               genome = genome,
                               groups = avail.pip.groups,
                               col = pip.colors[avail.pip.groups],
                               name = 'PIP',
                               legend = TRUE,
                               box.legend = FALSE)

  dpars.pip <- list(col.title = 'black',
                    col.axis = 'black',
                    #col.border.title = 'lightgray',
                    #col.frame = 'lightgray',
                    # rotation.title = rotation.title,
                    background.title = "transparent",
                    frame = TRUE,
                    cex.axis = 0.8,#0.8,
                    cex.legend = 0.8)
  Gviz::displayPars(pip.track) <- dpars.pip

  # Data tracks
  if (length(counts) > 0) {
    dpars.data <- list(col.title = 'black',
                       col.axis = 'black',
                       background.title = "transparent",
                       #col.border.title = 'lightgray',
                       #col.frame = 'lightgray',
                       rotation.title = rotation.title,
                       frame = TRUE,
                       cex.axis = 0.2)

    if(missing(counts.color)){
      counts.color <- seq_along(counts)
    }

    counts.tracks <- lapply(names(counts), function(x){
      if(verbose){ cat('Adding', x, 'track...\n') }

      counts.gr <- counts[[x]]
      i <- which(names(counts) == x)
      seqlevelsStyle(counts.gr) <- 'UCSC'
      seqlevels(counts.gr, pruning.mode = 'coarse') <- paste0('chr',1:22)
      counts.track <- Gviz::DataTrack(range = counts.gr,
                                      type = 'h',
                                      genome = genome,
                                      col = counts.color[i],
                                      name = x,
                                      showAxis=FALSE,
                                      ylim = counts.ylim)

      Gviz::displayPars(counts.track) <- dpars.data
      counts.track
    })
  }else{
    counts.tracks <- NULL
  }

  # RNA-seq stacked tracks (added by Jing)
  if (length(genes) > 0) {
    dpars.rna <- list(col.title = 'black',
                        col.axis = 'black',
                        background.title = "transparent",
                        #col.border = "transparent",
                        #col.border.title = 'white',
                        #col.frame = 'white',
                        rotation.title = rotation.title,
                        frame = FALSE,
                        cex.axis = 1.0)
    
  # if(missing(rna.color)){
  #     rna.color <- seq_along(genes)
  #   }
  rna.track <-
      Gviz::DataTrack(
        range = rna.gr[which(rna.gr$symbol %in% genes)],
        data = matrix(rna.counts[rownames(rna.counts) %in% genes, ], 
                      byrow = TRUE,
                      nrow = ncol(rna.counts)),
        name = "RNA",
        genome = "hg38",
        type = "histogram",
        stackedBars = FALSE,
        col = rna.color,
        groups = factor(colnames(rna.counts),
                        levels = c("Naive_B","Memory_B","NK","CD8_T","CD4_T", "Th17","Treg")),
        legend = FALSE)
  }
  
  Gviz::displayPars(rna.track) <- dpars.rna

  
  # Peak annotation tracks
  if (length(peaks) > 0) {
    dpars.peaks <- list(col.title = 'black',
                        col.axis = 'black',
                        background.title = "transparent",
                        #col.border.title = 'lightgray',
                        #col.frame = 'lightgray',
                        rotation.title = rotation.title,
                        frame = FALSE,
                        cex.axis = 0.2)

    if(missing(peaks.color)){
      peaks.color <- seq_along(peaks)
    }

    peaks.tracks <- lapply(names(peaks), function(x){
      if(verbose){ cat('Adding', x, 'track...\n') }
      peaks.gr <- peaks[[x]]
      i <- which(names(peaks) == x)
      seqlevelsStyle(peaks.gr) <- 'UCSC'
      seqlevels(peaks.gr, pruning.mode = 'coarse') <- paste0('chr',1:22)
      peaks.gr$score <- 1
      peaks.gr <- peaks.gr[countOverlaps(peaks.gr, region) > 0]
      peaks.track <- Gviz::DataTrack(range = peaks.gr,
                                     type = 'h',
                                     genome = genome,
                                     name = x,
                                     col = peaks.color[i],
                                     showAxis = FALSE,
                                     ylim = c(0,1))
      Gviz::displayPars(peaks.track) <- dpars.peaks
      peaks.track
    })
  }else{
    peaks.tracks <- NULL
  }

  # Chromatin loop tracks
  if (length(loops) > 0) {
    # dpars.loops <- list(#col.interactions = loops.color,
    #                     col.anchors.fill = 'blue',
    #                     col.anchors.line = 'black',
    #                     interaction.dimension = 'width',
    #                     interaction.measure = 'score',
    #                     plot.trans = FALSE,
    #                     plot.outside = FALSE,
    #                     col.outside='lightblue',
    #                     anchor.height = 0.2,
    #                     rotation.title = rotation.title,
    #                     col.title = 'black',
    #                     col.axis = 'black',
    #                     background.title = "transparent",
    #                     #col.border.title = 'lightgray',
    #                     #col.frame = 'lightgray',
    #                     frame = FALSE,
    #                     cex.axis = 0.2)

    gene.annots$chr <- as.character(seqnames(gene.annots))
    gene.annots$tss <- start(resize(gene.annots, width = 1))

    loops.tracks <- vector("list", length(loops))
    names(loops.tracks) <- names(loops)
    for(x in names(loops)){
      if(verbose){ cat('Adding', x, 'track...\n') }
      loops.gr <- loops[[x]]
      if(filter_protein_coding_genes){
        loops.gr <- loops.gr[loops.gr$gene_name %in% gene.annots$gene_name]
      }

      loops_promoters.gr <- GRanges(seqnames = seqnames(loops.gr),
                                    ranges = IRanges(start = loops.gr$promoter_start,
                                                     end = loops.gr$promoter_end),
                                    score = loops.gr$score,
                                    gene = loops.gr$gene_name)
      loops_enhancers.gr <- GRanges(seqnames = seqnames(loops.gr),
                                    ranges = IRanges(start = start(loops.gr),
                                                     end = end(loops.gr)))
      loops.obj <- GenomicInteractions::GenomicInteractions(anchor1 = loops_promoters.gr,
                                                            anchor2 = loops_enhancers.gr)
      loops.obj$counts <- round(loops.obj$anchor1.score)

      if (length(filter_loop_genes) >0 ) {
        cat(sprintf('Only show %s loops linked to gene: %s \n',
                    x, paste(filter_loop_genes, collapse = ',')))
        loops.obj <- loops.obj[which(loops.obj$anchor1.gene %in% filter_loop_genes),]
      }

      if (length(filter_loop_snps) > 0){
        cat(sprintf('Only show %s loops linked to SNP: %s \n',
                    x, paste(filter_loop_snps, collapse = ',')))
        highlighted.snps.gr <- finemapstats[finemapstats$snp %in% filter_loop_snps]
        loops.obj <- subsetByOverlaps(loops.obj, highlighted.snps.gr)
      }
      
      loops.track <- GenomicInteractions::InteractionTrack(loops.obj, name = x)
      
      dpars.loops <- list(
          col.interactions = loops.color[[x]][findOverlaps(loops[[x]], region)@from],
          col.anchors.fill = 'blue',
          col.anchors.line = 'black',
          interaction.dimension = 'width',
          interaction.measure = 'score',
          plot.trans = FALSE,
          plot.outside = FALSE,
          col.outside='lightblue',
          anchor.height = 0.2,
          rotation.title = rotation.title,
          col.title = 'black',
          col.axis = 'black',
          background.title = "transparent",
          #col.border.title = 'lightgray',
          #col.frame = 'lightgray',
          frame = FALSE,
          cex.axis = 0.2)
      
      Gviz::displayPars(loops.track) <- dpars.loops
      loops.tracks[[x]] <- loops.track
    }
  }else{
    loops.tracks <- NULL
  }

  # gene track
  gene.track <- make_genetrack_obj(region,
                                   txdb,
                                   gene.annots,
                                   genome,
                                   name = 'Gene',
                                   filter_protein_coding_genes = filter_protein_coding_genes)

  dpars.genes <- list(col.title = 'black',
                      col.axis = 'black',
                      background.title = "transparent",
                      #col.border.title = 'lightgray',
                      #col.frame = 'lightgray',
                      rotation.title = rotation.title,
                      frame = FALSE,
                      cex.item = 1.2,
                      cex = 0.5)
  displayPars(gene.track) <- dpars.genes

  # Restrict to protein coding genes
  if (isTRUE(filter_protein_coding_genes)){
    gene.track <- gene.track[symbol(gene.track) %in% gene.annots$gene_name]
  }

  # Genome axis track
  axisTrack <- Gviz::GenomeAxisTrack()

  # List all tracks
  list.of.tracks <- c(pval.track,
                      pip.track,
                      counts.tracks,
                      peaks.tracks,
                      gene.track,
                      loops.tracks,
                      #rna.track,
                      axisTrack
                      )

  if (is.null(track.sizes)){
    track.sizes <- c(1,
                     0.6,
                     rep(0.3, length(counts.tracks)),
                     rep(0.2, length(peaks.tracks)),
                     rep(0.6, length(loops.tracks)),
                     #0.5, # rna.track
                     0.5,
                     0.4)
  }

  # Highlight SNPs
  if (length(highlight_snps) > 0){

    if('topSNP' %in% highlight_snps){
      highlight_snps[which(highlight_snps == "topSNP")] <- curr.finemapstats$snp[which.max(curr.finemapstats$pip)]
    }

    highlight_snps <- intersect(highlight_snps, curr.finemapstats$snp)
    highlight.pos <- curr.finemapstats$pos[match(highlight_snps, curr.finemapstats$snp)]
    cat('Highlight SNP:', highlight_snps[which(highlight_snps %in% curr.finemapstats$snp)], '\n')

    if(length(highlight.pos) == 1){
      list.of.tracks <- Gviz::HighlightTrack(trackList = list.of.tracks,
                                             start = highlight.pos, #c(highlight.pos-500),
                                             width = 1,#1000,
                                             chromosome = as.character(seqnames(region)),
                                             col = highlight.color)
    }else if(length(highlight.pos) > 1){
      cat('Highlight position:', highlight.pos, '\n')

      if(length(highlight.color) > 1){
        highlight.color <- highlight.color[which(highlight_snps %in% curr.finemapstats$snp)]
      }
      list.of.tracks <- Gviz::HighlightTrack(trackList = list.of.tracks,
                                             start = highlight.pos,
                                             width = 1,
                                             chromosome = as.character(seqnames(region)),
                                             col = highlight.color)
    }
  }

  if(verbose){ cat('Making track plot ...\n') }

  Gviz::plotTracks(list.of.tracks,
                   chromosome = as.character(seqnames(region)),
                   transcriptAnnotation = 'symbol',
                   collapseTranscripts = 'longest',
                   background.title = "transparent", # remove grey background
                   from = start(region),
                   to = end(region),
                   sizes = track.sizes,
                   fontsize = 10,
                   just.group = genelabel.side,
                   ...)

}


# Make gene track object for Gviz
make_genetrack_obj <- function(region,
                               txdb = NULL,
                               gene.annots = NULL,
                               genome = 'hg19',
                               keytype = 'ENSEMBL',
                               name = 'Gene',
                               filter_protein_coding_genes = TRUE){

  if(!is.null(txdb)){
    # use supplied txdb if available.
    cat('Making gene track object using txdb database ...\n')
    grtrack <- Gviz::GeneRegionTrack(range = txdb,
                                     genome = genome,
                                     chromosome = as.character(seqnames(region)),
                                     start = start(region),
                                     end = end(region),
                                     name = name)

    symbol(grtrack) <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                             keys=sub('\\.\\d.*$', '', gene(grtrack)),
                                             keytype=keytype,
                                             column='SYMBOL')

    symbol(grtrack) <- ifelse(is.na(symbol(grtrack)), '', symbol(grtrack))

  }else if(!is.null(gene.annots)){
    # use supplied gene.annots if available.
    cat('txdb not available. Making gene track object using gene.annots ...\n')
    grtrack <- Gviz::GeneRegionTrack(range = gene.annots,
                                     genome = genome,
                                     chromosome = as.character(seqnames(region)),
                                     start = start(region),
                                     end = end(region),
                                     name = name,
                                     symbol = gene.annots$gene_name)
  }

  # restrict to protein coding genes
  if (isTRUE(filter_protein_coding_genes)){
    grtrack <- grtrack[symbol(grtrack) %in% gene.annots$gene_name]
  }

  return(grtrack)
}

