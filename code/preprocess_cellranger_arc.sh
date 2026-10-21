#!/bin/bash

# Run Cell Ranger ARC's pipelines

## Generate single cell feature counts for a single library
arc_dir=""  # path to cellranger-arc software
sample_id=""  # unique ID for each sample
reference_dir=""  # path to a reference package downloade from "https://www.10xgenomics.com/support/software/cell-ranger-arc/downloads"
library_file=""  # full path to a an aggregation CSV with the required columns: library_id, atac_fragments, per_barcode_metrics, gex_molecule_info

${arc_dir}/cellranger-arc count --id=${sample_id} \
                       --reference=${reference_dir} \
                       --libraries=${library_file} \
                       --localcores=16 \
                       --localmem=64

## Get an aggregated matrix by pooling multiple libraries from each tissue
group_id=""  # unique ID for each group
AGG_library=""  # full path to a 3-column libraries CSV file that specifies the location of fastq files, sample_ID, and library type.

 ${arc_dir}/cellranger-arc aggr --id=${group_id} \
	    --csv=${AGG_library} \
      --reference=${reference_dir} \	  
      --normalize=depth \ 
		  ---nosecondary 
  
                     
                       
                