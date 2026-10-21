

EnrichmentTestE2G <- function(eqtl_gr, multiome_gr) {
  eqtl_gr <- eqtl_gr[eqtl_gr$distance <= 100000 & eqtl_gr$distance > 10000]
  multiome_gr <- multiome_gr[multiome_gr$distance <= 100000 & multiome_gr$distance > 10000]
  positive_gr <- multiome_gr[multiome_gr$E2G.Score.qnorm >= 0.2]
  positive_links <- data.frame(peak = as.character(positive_gr),
                               gene = positive_gr$TargetGene)
  # matching SNP-gene pairs
  overlaps <-
    findOverlaps(eqtl_gr, multiome_gr+250) 
  overlapped_links <- data.frame(
    peak = as.character(multiome_gr[overlaps@to]),
    gene = as.character(eqtl_gr$targetGene[overlaps@from]))
  
  out <-
    EnrichmentTest(unique(paste(positive_links$peak, positive_links$gene, sep = ":")),
                   unique(paste(overlapped_links$peak, overlapped_links$gene, sep = ":")),
                   unique(paste(as.character(multiome_gr), multiome_gr$TargetGene, sep = ":")))
  return(out)
}

# Make a list of fine-mapped eQTL-gene pairs (PIP>=0.3) from a public dataset
ref <- fread("tabix_ftp_paths.tsv")
ftp_path <- 
  ref %>% filter(study_label == "OneK1K" & quant_method == "ge") %>%
  filter(sample_group %in% c("B_memory", "B_naive", "CD4_Naive", 
                             "CD4_TCM", "CD4_TEM", "CD8_Naive",
                             "CD8_TCM", "CD8_TEM", "NK", "NK_Proliferating", "Treg")) 

eqtl_list <- 
  lapply(c("B_memory", "B_naive", "CD4", "CD8", "NK", "Treg"), function(i) {
    idx = grep(i, ftp_path$sample_group)
    print(idx)
    eqtl <- c()
    for(i in idx){
      f <- fread(paste0("/your path/", basename(ftp_path$ftp_cs_path[i])))
      eqtl <- rbind(eqtl, f[f$pip >= 0.3, ])}
    eqtl <- eqtl[eqtl$variant %in% all_variants$variant] #keep common variants
    # make a GRange object
    eqtl$position <-
      unlist(lapply(eqtl$variant, function(i) {
        paste(strsplit(i, "_")[[1]][1:2], collapse = ":")
      }))
    eqtl_gr <- GRanges(eqtl$position)
    mcols(eqtl_gr) <- eqtl
    # convert to gene symbols
    mapped_genes <- mapGeneIDs(
      query = eqtl_gr$gene_id,
      query_type = "ENSEMBL",
      output_types = "SYMBOL"
    )
    eqtl_gr$targetGene <-
      mapped_genes$SYMBOL[match(eqtl_gr$gene_id, mapped_genes$ENSEMBL)]
    eqtl_gr <-
      eqtl_gr[!is.na(eqtl_gr$targetGene) &
                eqtl_gr$targetGene %in% proj@geneAnnotation$genes$symbol]
    # compute the distance of variants to gene body
    eqtl_gr$distance <-
      unlist(lapply(seq_along(eqtl_gr), function(i) {
        gene = proj@geneAnnotation$genes[which(proj@geneAnnotation$genes$symbol == eqtl_gr$targetGene[i])]
        return(distance(eqtl_gr[i], gene))
      }))
    return(eqtl_gr) 
  })
names(eqtl_list) <- c("B_memory", "B_naive", "CD4", "CD8", "NK", "Treg")




df <- data.frame(eqtl = c("CD4", "CD8", "B_naive", "B_memory", "NK"),
                 scE2G = names(p2g_list)[c(1:4,8)])
enrich_out <-
  apply(df, 1, function(i) {
    return(EnrichmentTestE2G(eqtl_list[[i[1]]], p2g_list[[i[[2]]]])[[2]])
  })

gginput <- data.frame(CellType = df$scE2G, 
                      do.call(rbind, enrich_out))

ggplot(gginput %>% filter(CellType != "Memory_B"),
       aes(
         y = odds_ratio,
         x = CellType,
         fill = CellType)) + 
  geom_col(width = 0.4) + 
  geom_errorbar(aes(ymin = conf_low , 
                    ymax = conf_high),
                width = 0.1) + 
  theme_pub() + xlab("") + ylab("Enrichment of eQTL-gene pairs \nin E2G links vs. non-E2G links \n(within 10-100Kb of genes)") + 
  theme(legend.position = "none")
