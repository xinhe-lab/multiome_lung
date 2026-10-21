# Identify asthma-associated cis-regulatory elements

## Two sources
## 1) Overlapping with confident variants that contribute to asthma risk enhancers 
## confident enhancer sequences: ePIPs $\ge$ 0.5
## peaks need to be within 250 bp from confidently fine-mapped variants (PIP $\ge$ 0.3)

## 2) Identify peaks near confidently fine-mapped variants 
## peaks need to be within 250 bp from confidently fine-mapped variants (PIP $\ge$ 0.3)

f <- fread("data/CRE_0_summary.txt")

## process GWAS summary statistics with fine-mapped results

f_aoa <- data.frame(f) %>% filter(aoa_epip > 0) %>%
  select(-c(peak_cell_type, lineage, coa_epip, coa_snps_pip, coa_snps, snps)) %>%
  separate_longer_delim(c(aoa_snps, aoa_snps_pip), delim = ",") %>%
  mutate(epip = aoa_epip, susie_pip = aoa_snps_pip, 
         snp = aoa_snps, trait = "aoa")

f_coa <- f %>% filter(coa_epip > 0) %>%
  select(-c(peak_cell_type, lineage, aoa_epip, aoa_snps_pip, aoa_snps, snps)) %>%
  separate_longer_delim(c(coa_snps, coa_snps_pip), delim = ",") %>%
  mutate(epip = coa_epip, susie_pip = coa_snps_pip, 
         snp = coa_snps, trait = "coa")

f_aoa$pos <- ukbb_snps$pos[match(f_aoa$aoa_snps, ukbb_snps$id)]
f_coa$pos <- ukbb_snps$pos[match(f_coa$coa_snps, ukbb_snps$id)]
f_coa$pos[is.na(f_coa$pos)] <- 137605991

f_comb <- rbind(f_coa %>% select(-c(coa_snps, coa_snps_pip, coa_epip)),
                f_aoa %>% select(-c(aoa_snps, aoa_snps_pip, aoa_epip)))

## Identify asthma CREs using the first approach

colnames(f_comb)[7:8] <- c("susie_pip", "snp")
f_comb$susie_pip <- as.numeric(f_comb$susie_pip)
f_comb <- rbind(left_join(f_comb[f_comb$trait == "coa", ], coa[, c("snp", "pval")]),
                left_join(f_comb[f_comb$trait == "aoa", ], aoa[, c("snp", "pval")]))
# Overlap asthma fine-mapped regions with OCRs in our dataset
gr <- GRanges(seqnames=paste0("chr", f_comb$chr),
              IRanges(start = f_comb$start,
                      end = f_comb$end),
              epip = f_comb$epip,
              snps_pip = f_comb$susie_pip,
              snps = f_comb$snp,
              pos = f_comb$pos,
              #pval = f_comb$pval,
              trait = f_comb$trait)

# liftover to hg38
ch=import.chain("~/resources/genomes/chain_files/hg19ToHg38.over.chain")

seqlevelsStyle(gr) <- "UCSC"
gr_hg38 <- unlist(liftOver(gr, ch))
gr_hg38$labels <- paste(as.character(gr_hg38))

gr_hg38$type <- ifelse(gr_hg38$epip >= 0.5, "high", "low")
gr_hg38$snp_type <- ifelse(gr_hg38$snps_pip >= 0.3, "high", "low")

union_peaks <- GRanges(seqnames(proj@peakSet),
                       IRanges(start = start(proj@peakSet),
                               end = end(proj@peakSet)),
                       nearestGene = proj@peakSet$nearestGene)

# overlap with high ePIP enhancers
paired <- findOverlaps(union_peaks, gr_hg38)
gr_peaks <-union_peaks[paired@from]

# calculate the distances btw fine-mapped SNPs and the overlapped peak sets
gr_peaks$distToSNP <- unlist(
  lapply(1:length(paired), function(i){
    gr_snps <- 
      GRanges(seqnames = as.character(gr_hg38[paired@to[i]]@seqnames), 
              IRanges(start = gr_hg38[paired@to[i]]$pos))
    dist_gr = distanceToNearest(union_peaks[paired@from[i]], gr_snps)
    return(dist_gr@elementMetadata$distance)
  }))

# add metadata
mcols(gr_peaks) <- 
  cbind(mcols(gr_peaks),
        mcols(gr_hg38[paired@to]))
gr_peaks$peakID <- as.character(gr_peaks)
gr_hg38_paired <- gr_hg38[paired@to]
#gr_peaks_paired <- gr_peaks

# filter out peaks that are not near any SNPs
paired <- paired[gr_peaks$distToSNP <= 250]
gr_peaks <- gr_peaks[gr_peaks$distToSNP <= 250]
ref_peaks <- gr_peaks

## Identify asthma CREs using second approach
coa <- readRDS("data/coa_gp3_finemapping_gwas_L5.rds")
aoa <- readRDS("data/aoa_gp3_finemapping_gwas_L5.rds")

aoa_gr <- GRanges(seqnames = paste0("chr", aoa$chr),
                  IRanges(start = aoa$pos,
                          end = aoa$pos+1))
mcols(aoa_gr) <- data.frame(aoa[, -4:-1], trait = "aoa")
coa_gr <- GRanges(seqnames = paste0("chr", coa$chr),
                  IRanges(start = coa$pos,
                          end = coa$pos+1))
mcols(coa_gr) <- data.frame(coa[, -4:-1], trait = "coa")
merge_gr <- c(coa_gr, aoa_gr)
#liftOver to hg38
ch=import.chain("data/hg19ToHg38.over.chain")
seqlevelsStyle(merge_gr) <- "UCSC"
snp_gr <- unlist(liftOver(merge_gr, ch))
paired <- findOverlaps(proj@peakSet + 250, snp_gr)
gr_peaks <- proj@peakSet[paired@from[which(paired@from %in% peak_idx)]]
gr_peaks$labels <- paste(as.character(gr_peaks),
                         gr_peaks$nearestGene,
                         sep = "_")
gr_peaks$peakID <- peak_idx
mcols(gr_peaks) <- data.frame(mcols(gr_peaks), 
                              mcols(snp_gr)[paired@to[which(paired@from %in% peak_idx)], c("snp", "susie_pip", "pval", "trait")])
gr_peaks$type <- ifelse(gr_peaks$susie_pip >=0.3, "high", "low")


# combine two approaches
ref_peaks <- ref_peaks[ref_peaks$type == "high" & ref_peaks$snp_type == "high"]
paired <- findOverlaps(ref_peaks, gr_peaks)
diff <- ref_peaks[setdiff(1:length(ref_peaks), paired@from)]
colnames(mcols(diff))[4:5] <- c("susie_pip", "snp")
diff_gr <- proj@peakSet[findOverlaps(proj@peakSet, diff)@from]
mcols(diff_gr) <- 
  data.frame(mcols(diff_gr), 
             mcols(diff)[findOverlaps(proj@peakSet, diff)@to, c("snp", "susie_pip", "pval", "trait", "type")])
union_peaks <- c(gr_peaks, diff_gr)


# filter CREs that contain SNPs not in a credible set
gr_peaks$CS <- 0
gr_peaks$pval <- 0
coa_index <- which(gr_peaks$trait == "coa")
aoa_index <- which(gr_peaks$trait == "aoa")

gr_peaks$CS[coa_index] <- coa$CS[match(gr_peaks$snp[coa_index], coa$snp)]
gr_peaks$CS[aoa_index] <- aoa$CS[match(gr_peaks$snp[aoa_index], aoa$snp)]
gr_peaks$pval[coa_index] <- coa$pval[match(gr_peaks$snp[coa_index], coa$snp)]
gr_peaks$pval[aoa_index] <- aoa$pval[match(gr_peaks$snp[aoa_index], aoa$snp)]

gr_peaks_sub <- gr_peaks[gr_peaks$CS > 0]

# Identify cell types where asthma CREs were called as peaks
peak_dir = "output/u19_peaks_tissue/PeakCalls/"
peaks<-list()
for(i in list.files(peak_dir, pattern = ".rds")[1:9]){
  term <- strsplit(i, "-")[[1]][1]
  peaks[[term]] <- readRDS(paste0(peak_dir, i))
}

celltype_list <- 
  lapply(1:length(gr_peaks_sub), function(i){
    selected_ct = c()
    for(j in names(peaks)){
      res = findOverlaps(gr_peaks_sub[i], peaks[[j]])
      if(length(res) > 0){
        selected_ct = gsub("Lung_", "", c(selected_ct, j))
      }
    }
    return(selected_ct)
  })
gr_peaks_sub$Peak <- as.character(gr_peaks_sub)
gr_peaks_sub$CellTypes <-
  unlist(lapply(celltype_list, function(i){
    paste(i[i != "Other"], collapse = ",")
  }))

celltype_list <- 
  lapply(1:length(gr_peaks_sub), function(i){
    selected_ct = c()
    for(j in names(peaks)){
      res = findOverlaps(gr_peaks_sub[i], peaks[[j]])
      if(length(res) > 0){
        if(grepl("CD4_T", j) | grepl("Treg", j)){
          j = "CD4_T"
        }else if(grepl("B", j)){
          j = "B"
        }
        selected_ct = gsub("Lung_", "", c(selected_ct, j))
      }
    }
    return(unique(selected_ct))
  })

gr_peaks_sub <- gr_peaks_sub[gr_peaks_sub$CellTypes!= ""]
