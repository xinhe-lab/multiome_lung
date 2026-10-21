

make_summary_table <- 
  #make a summary table for LDSC outputs from individual test
  function(traits, 
           referenceTab, # a table of two columns: input annotations and displayed annotations
           baseline="baseline_v1.2"
  ){
    
    output_df<-c()
    for (trait in traits){
      for (annot in referenceTab$annot_id){
        f<-read.table(
          sprintf("~/cluster/projects/ldsc_enrichment/%s/%s/%s.results", 
                  trait, baseline, annot), header = T)
        if(TRUE %in% grepl("^L2_1", f$Category)){
          f$Category[!grepl("L2_0", f$Category)] <- paste0(annot, "L2_1")
        }
        output_df<-rbind(output_df, cbind(trait, f[!grepl("L2_0", f$Category), ]))
      }
    }
    
    # trim off suffix of labels
    output_df$Category <- unlist(
      lapply(output_df$Category, function(i){strsplit(i, "L2")[[1]][1]})
    )
    
    return(output_df) 
  }


estimate_sd_annot <- function(cell_path){
  ## A function to estimate standard deviation for a given annotation
  ldsc_dir <- "/home/jinggu/cluster/projects/ldsc_enrichment/ldscore"
  Mref = 5961159
  ll <- list.files(
    paste(ldsc_dir, cell_path, sep="/"),
    pattern = "annot.gz")
  num = 0
  sd_annot = rep(0, 22)
  for(m in 1:length(ll)){
    f = read.table(paste(ldsc_dir, cell_path, ll[m], sep="/"), skip=1)
    num = num + nrow(f)
    sd_annot[m] = (nrow(f)-1) * var(f$V1)
  }
  est_sd_annot = sqrt(sum(sd_annot)/(num-1))
  return(est_sd_annot)
}



compute_tau_star <- function(ldsc_tbl, annotations){
  # compute tau* enrichment and p-value given a summary table and a vector of annotations
  sd_annot <-
    data.frame(Category = annotations,
               annot_std_error = unlist(lapply(annotations, function(i){estimate_sd_annot(i)})))
  ldsc_tbl <- left_join(ldsc_tbl, trait_snp_h2g, by = "trait")
  ldsc_tbl <- left_join(ldsc_tbl, sd_annot, by = "Category")
  ldsc_tbl$tau_star <-  ldsc_tbl$Coefficient/ldsc_tbl$perSNPh2g * ldsc_tbl$annot_std_error
  corr_ldsc_tbl <-
    ldsc_tbl %>%
    mutate(
      tar_star_se = Coefficient_std_error / perSNPh2g *
        annot_std_error,
      
      tau_star_p = pnorm(-abs(Coefficient_z.score)) * 2)
      #tau_star_p = format(pnorm(-abs(Coefficient_z.score)) * 2,
                          #scientific = TRUE, digit = 1),
      #tau_star = round(tau_star, 1))
  return(corr_ldsc_tbl)
}


trait_dict <- function(trait_list){
  # retrieve trait ID from a dictionary
  dict = list(
    height = "Height",
    bmi = "BMI",
    LDL = "LDL",
    HDL = "HDL",
    rbc = "rbc",
    wbc = "wbc",
    monocyte = "monocyte",
    lymph = "lymphocyte",
    scz = "SCZ",
    aFib = "AFib",
    AllergicRhinitis = "AR",
    AtopicDermatitis = "AtopicDerm",
    aoa_ukb = "AOA",
    coa_ukb = "COA",
    IPF = "IPF",
    COPD = "COPD",
    COVID19 = "COVID19"
  )
  return(unlist(lapply(trait_list, function(i){dict[[i]]})))
}

# save a list of tables into excel sheets
saveToExcelSheets <-
  function(input_list,
           file_name,
           rowname_include = TRUE) {
    library(openxlsx)
    wb <- createWorkbook()
    for (name in names(input_list)) {
      addWorksheet(wb, sheetName = name)
      writeData(wb,
                sheet = name,
                x = input_list[[name]],
                rowNames = TRUE)
    }
    saveWorkbook(wb, sprintf("output/%s.xlsx", file_name), overwrite = TRUE)
  }

EnrichmentTest <- function(var1, var2, bkg){
  # fisher exact test
  var2 <- intersect(var2, bkg) # this step makes sure var2 were all included in bkg
  mat <- matrix(
    c(length(intersect(var1, var2)),
      length(intersect(setdiff(bkg, var2), var1)),
      length(intersect(setdiff(bkg, var1), var2)),
      length(intersect(setdiff(bkg, var1), setdiff(bkg, var2)))
    ), nrow = 2, byrow = TRUE)
  #print(mat)
  res <- fisher.test(mat)
  out <- data.frame(
    p_value = format(res$p.value, scientific = TRUE, digits = 2),
    odds_ratio = round(res$estimate, 2),
    conf_low = round(res$conf.int[1], 2),
    conf_high = round(res$conf.int[2], 2)
    #method = res$method
  )
  return(list(mat = mat,
              out = out))
}

# map gene id from one type to another
mapGeneIDs <- function(query, query_type,
                       output_types) {
  
  library(org.Hs.eg.db)
  mapped = AnnotationDbi::select(org.Hs.eg.db,
                                 keys = query,
                                 columns = output_types,
                                 keytype = query_type)
  return(mapped)
}


# Liftover a GRange from one genome build to the other
library(rtracklayer)
makeLiftOver <- function(input_gr, chain) {
  #Given an input of GRanges object, liftover it to a different genome build.
  
  # This outputs a GRangesList object
  lifted_list <- liftOver(input_gr, chain)
  
  # Step 5: Flatten the list back into a standard GRanges object
  lifted_ranges <- unlist(lifted_list)
  return(lifted_ranges)
  
  # Step 6: (Optional) Convert the results back to a standard R data frame
  #lifted_df <- as.data.frame(lifted_ranges)
}

