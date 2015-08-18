\name{sleuthData}
\alias{sleuthData}
\docType{data}
\title{RNAseq data for the `sleuth` example}
\description{
This package contains hdf5 output from running Kallisto 0.42.2.1 against the 
samples in GSE37704, accompanying Trapnell et al, Nat Biotech 2013.  
}
\examples{

library(dplyr)
library(sleuth)
library(sleuthData)

## find the data
base_dir <- system.file("extdata", "", package="sleuthData")
sample_id <- dir(file.path(base_dir,"results"))
kal_dirs <- sapply(sample_id, 
                   function(id) file.path(base_dir, "results", id, "kallisto"))

## map samples to covariates
s2c <- read.table(file.path(base_dir, "hiseq_info.txt"), 
                  header=TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample=run_accession, condition)

## run sleuth with no annotation
so <- sleuth_prep(kal_dirs, s2c, ~ condition)
so <- sleuth_fit(so)
so <- sleuth_test(so, which_beta = 'conditionscramble')
models(so)

## annotate transcripts
library(biomaRt)
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", 
                                     "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
so <- sleuth_prep(kal_dirs, s2c, ~ condition, target_mapping = t2g)
so <- sleuth_fit(so)
so <- sleuth_test(so, which_beta = 'conditionscramble')

## visualize the results
sleuth_live(so)

## a table for downstream analysis
results_table <- sleuth_results(so, 'conditionscramble') 

}
\keyword{datasets}
