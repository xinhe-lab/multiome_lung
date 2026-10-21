
## Load packages
load_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}
load_package("WebGestaltR")
load_package("dplyr")
load_package("optparse")
## Process command-line arguments 
#args <- commandArgs(trailingOnly = TRUE)
parser <- OptionParser()
parser <- add_option(parser, c("--geneset"), type = "character", default = "")  # a list of genesets saved in RDS
parser <- add_option(parser, c("--bkg"), type = "character", default = "")  # background genes saved in RDS
parser <- add_option(parser, c("--output"), type = "character", default = "")  # prefix of the file name (XX_GO_enrichment.RDS)
out <- parse_args(parser)

## Check the list of database
listGeneSet(
  organism = "hsapiens",
  hostName = "https://www.webgestalt.org/",
  cache = NULL
)

genesets <- readRDS(out$geneset)
bkgGenes <- readRDS(out$bkg)
database <- c(#"pathway_KEGG",
               "geneontology_Biological_Process_noRedundant",
               "geneontology_Cellular_Component_noRedundant",
               "geneontology_Molecular_Function_noRedundant")
               #"community-contributed_Hallmark50")

## Run GO enrichment test for each gene set in the input list
outputs<-list()

for(k in 1:length(genesets)){
    go_out<-WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                     #interestGene=top_genes[[k]][!grepl("MT-", top_genes[[k]])],
                     interestGene = genesets[[k]],
                     referenceGene = bkgGenes[[k]],
                     enrichDatabase= database,
                     interestGeneType="genesymbol",
                     referenceGeneType="genesymbol",
                     projectName = names(genesets)[k],
                     fdrThr= 1,
                     isOutput= F,
                     nThreads = 10)
    #print(colnames(out))
    if(!is.null(go_out)){
      outputs[[k]]<-data.frame(k = names(genesets)[k], 
                               go_out[, c("geneSet", "description", "database", "enrichmentRatio", "FDR", "userId")]
      )
    }
}

saveRDS(outputs, paste0(out$output, "_GO_enrichment.RDS"))