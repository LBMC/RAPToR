p_url_dswhite2017 <- paste0("http://europepmc.org/articles/PMC5690287/",
                            "bin/elife-30860-supp1.tsv")
g_url_dswhite2017 <- apste0("http://europepmc.org/articles/PMC5690287/",
                            "bin/elife-30860-supp2.tsv")
g_file_dswhite2017 <- paste0(data_folder, "dswhite2017.tsv")
utils::download.file(g_url_dswhite2017, destfile = g_file_dswhite2017)

X_dswhite2017 <- read.table(g_file_dswhite2017, h = T, sep  ="\t", 
                            as.is = T, quote = "\"")
rownames(X_dswhite2017) <- X_dswhite2017$Gene.ID
X_dswhite2017 <- X_dswhite2017[,-(1:8)]

# convert to tpm & ensembl_id
X_dswhite2017 <- X_dswhite2017[
  rownames(X_dswhite2017)%in%zeb_genes$ensembl_gene_id,]
X_dswhite2017 <- raw2tpm(
  rawcounts = X_dswhite2017, 
  genelengths = zeb_genes$transcript_length[
    match(rownames(X_dswhite2017), zeb_genes$ensembl_gene_id)])

# pheno data
P_dswhite2017 <- read.table(p_url_dswhite2017, h = T, sep = "\t", as.is = T)
P_dswhite2017 <- P_dswhite2017[P_dswhite2017$sequencing == "RNASeq", 
                               c("sample", "accession_number", "stage", 
                                 "stageName", "sampleName")]

# timings of stages (from White et al. eLife (2017)).
# in hours post-fertilization
timepoints <- data.frame(stage = unique(P_dswhite2017$stageName), 
                         hours_pf = c(0, .75, 2.25, 3, 4.3, 5.25, 6, 8, 10.3, 
                                      16, 19, 24, 30, 36, 48, 72, 96, 120),
                         stringsAsFactors = F, row.names = "stage")
P_dswhite2017$age <- timepoints[P_dswhite2017$stageName, "hours_pf"]

# formatting
P_dswhite2017$batch <- factor(gsub(".*-(\\d)$", "\\1", 
                                   P_dswhite2017$sampleName))
X_dswhite2017 <- X_dswhite2017[, P_dswhite2017$sample]

# save data
dswhite2017 <- list(g = X_dswhite2017, p = P_dswhite2017)
save(dswhite2017, file = paste0(data_folder, "dswhite2017.RData"), 
     compress = "xz")

# cleanup
file.remove(g_file_dswhite2017)
rm(p_url_dswhite2017, g_url_dswhite2017, g_file_dswhite2017, 
   X_dswhite2017, P_dswhite2017, timepoints)