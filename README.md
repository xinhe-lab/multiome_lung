# Single-cell lung multiome paper

The manuscript is currently under review by Genome Biology. Please check out the preprint at https://www.medrxiv.org/content/10.64898/2026.02.05.26345013v1.

## Analyses pipelines and codes

* LDSC workflow: used to assess trait heritability partitioned by different functional annotations 
* scE2G workflow: used to identify enhancer-to-gene links from single-cell multiome data
* Processing scRNA-seq datasets:
  - scRNA_dublets_removal.ipynb
  - scRNA_seq_integration_Seurat.ipynb
  - validate_cell_type_annotations.ipynb
  - cross-tissue_differential_gene_expression.ipynb
* Processing scATAC-seq datasets:
  - analyze_scATAC_seq_ArchR.ipynb
  - validate_cell_type_annotations_with_ATAC-seq.Rmd
* Interpreting hertiability enrichment results
  - heritability_enrichment_for_lung_OCRs.Rmd
* Interpreting scE2G results and link them to asthma genetics
  - link_CREs_to_genes.Rmd
  - summarize_candidate_genes_43CREs.Rmd
  - make_locus_plots.Rmd
* Interpreting GRN results and link them to asthma genetics
  - GRN_plotting.ipynb
  - TF_GRN_analysis.Rmd

