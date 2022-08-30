requireNamespace("biomaRt", quietly = TRUE)

mart <- biomaRt::useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
droso_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id", 
                                             "ensembl_transcript_id",
                                             "external_gene_name",
                                             "transcript_end", "transcript_start"),
                              mart = mart)
droso_genes$transcript_length <- droso_genes$transcript_end - droso_genes$transcript_start
droso_genes <- droso_genes[,c(1:3,6)]
colnames(droso_genes)[1:3] <- c("fb_id", "transcript_id", "gene_name")

rm(mart)