requireNamespace("biomaRt", quietly = TRUE)

mart <- biomaRt::useMart("ensembl", dataset = "drerio_gene_ensembl")
zeb_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id", 
                                           "transcript_length"), 
                            mart = mart)
rm(mart)