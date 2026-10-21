# Validate GRNs with marker genes of CD4 T cells and B cell subsets
ct <- "CD4_T"
marker_expr <- readRDS(sprintf("%s/lung_CD4_T_vs_Mem_B_DEG.RDS", data_dir))
DEGenes <- 
  marker_expr %>% 
  mutate(gene = rownames(marker_expr)) %>%
  filter(p_val_adj <= 0.05 & avg_log2FC >= 0.58) %>%
  pull(gene)


enrich_out <-
  lapply(names(target_list), function(i) {
    out_tbl <- c()
    for(j in unique(unlist(select_TFs))) {
      if (j %in% names(target_list[[as.character(i)]][[ct]])) {
        target <- unique(target_list[[as.character(i)]][[ct]][[j]])
        out_tbl = rbind(out_tbl, EnrichmentTest(target, DEGenes, all_genes[[ct]]))
      }
    } 
    return(data.frame(out_tbl, TF = unique(unlist(select_TFs)), cutoff = i))
  })

enrich_out <-
  do.call(rbind, enrich_out) %>% 
  na.replace(0) 

count_tbl <- 
  enrich_out %>% 
  group_by(cutoff) %>% 
  dplyr::summarise(numEnriched= sum(p.adjust(p_value, method = "fdr") <= 0.05 & odds_ratio > 1),
                   numNotEnriched = sum(p.adjust(p_value, method = "fdr") >0.05 | odds_ratio <= 1))

count_tbl  

# Identify master TF regulators for asthma 
sparsity = "0.015"
idx <- which(sapply(target_list[[sparsity]][["CD4_T"]], function(i){length(i) >= 100}) == TRUE)
regTFs <- unique(names(target_list[[sparsity]][["CD4_T"]])[idx])

tf_tbl <- c()
for(ct in ct_order) {
  #regTFs <- novel_TFs
  # test enrichment for asthma risk genes across marker TFs
  idx <- which(sapply(target_list$`0.015`[[ct]], function(i){length(i) >= 100}) == TRUE)
  regTFs <- unique(names(target_list$`0.015`[[ct]])[idx])
  gene_set <- gene_scores$Gene
  out_tbl <-
    lapply(regTFs, function(i){
      target = unique(target_list[[sparsity]][[ct]][[i]])
      return(EnrichmentTest(var1 = target, var2 = gene_set, bkg = all_list[[sparsity]][[ct]]))
    })
  out_tbl <- do.call(rbind, out_tbl)
  #rownames(out_tbl) <- sprintf("%s(%s)", regTFs$gene, regTFs$CellType)
  out_tbl <- data.frame(TF = regTFs, out_tbl)
  tf_tbl <- rbind(
    tf_tbl,
    data.frame(
      out_tbl,
      CellType = ct,
      FDR = p.adjust(out_tbl$p_value, method = "fdr"),
      row.names = NULL
    ))
}

asthma_TFs <- 
  unique(unlist(lapply(names(target_list$`0.015`), function(i){
    intersect(c(unique(gene_scores$Gene), novel_TFs), names(target_list$`0.015`[[i]]))}
  )))

asthma_out <- c()
for(i in asthma_TFs) {
  enrich_out <-
    lapply(ct_order[ct_order!= "Other"], function(j) {
      EnrichmentTest(var1 = target_list$`0.015`[[j]][[i]],
                     var2 = gene_scores$Gene,
                     bkg = all_genes[[j]])
    })
  enrich_out <- do.call(rbind, enrich_out)
  enrich_out$CellType <- ct_order[ct_order!= "Other"]
  asthma_out <- rbind(asthma_out, data.frame(TF = i, enrich_out))
}

asthma_out$FDR <- p.adjust(asthma_out$p_value, method = "fdr")
asthma_out %>% filter(FDR <= 0.1)

selectedTFs <-
  unique(c(asthma_out %>% filter(TF %in% asthma_out$TF[asthma_out$FDR <= 0.1]) %>% pull(TF),
           tf_tbl %>% filter(TF %in% tf_tbl$TF[tf_tbl$FDR <= 0.1]) %>% pull(TF)))
  
