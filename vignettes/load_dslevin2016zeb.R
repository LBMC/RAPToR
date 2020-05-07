geo_dslevin2016zeb <- "GSE60619"

g_url_dslevin2016zeb <- GEOquery::getGEOSuppFiles(geo_dslevin2016zeb, makeDirectory = FALSE, fetch_files = FALSE)
g_file_dslevin2016zeb <- paste0(data_folder, "dslevin2016zeb.txt.gz")
utils::download.file(url = as.character(g_url_dslevin2016zeb$url[2]), destfile = g_file_dslevin2016zeb)

X_dslevin2016zeb <- read.table(gzfile(g_file_dslevin2016zeb), h = T, sep = '\t', as.is = T, row.names = 1, comment.char = "")

# convert to tpm & ensembl_id
X_dslevin2016zeb <- X_dslevin2016zeb[rownames(X_dslevin2016zeb)%in%zeb_genes$ensembl_gene_id,]
X_dslevin2016zeb <- raw2tpm(rawcounts = X_dslevin2016zeb, 
                            genelengths = zeb_genes$transcript_length[match(rownames(X_dslevin2016zeb),
                                                                            zeb_genes$ensembl_gene_id)])


# pheno data
P_dslevin2016zeb <- Biobase::pData(GEOquery::getGEO(geo_dslevin2016zeb, getGPL = F)[[1]])

# filter relevant fields/samples
P_dslevin2016zeb <- P_dslevin2016zeb[, c("title", "geo_accession", "time (min after fertilization):ch1")]
colnames(P_dslevin2016zeb)[3] <- "time"
P_dslevin2016zeb$title <- as.character(P_dslevin2016zeb$title)

P_dslevin2016zeb <- P_dslevin2016zeb[P_dslevin2016zeb$title %in% colnames(X_dslevin2016zeb),]
X_dslevin2016zeb <- X_dslevin2016zeb[, P_dslevin2016zeb$title]

# formatting
P_dslevin2016zeb$title <- gsub('Metazome_ZF_timecourse_', '', P_dslevin2016zeb$title)
colnames(X_dslevin2016zeb) <- P_dslevin2016zeb$title

P_dslevin2016zeb$age <- as.numeric(P_dslevin2016zeb$time) / 60

dslevin2016zeb <- list(g = X_dslevin2016zeb, p = P_dslevin2016zeb)
save(dslevin2016zeb, file = paste0(data_folder, "dslevin2016zeb.RData"), compress = "xz")

# cleanup
file.remove(g_file_dslevin2016zeb)
rm(geo_dslevin2016zeb, g_url_dslevin2016zeb, g_file_dslevin2016zeb, X_dslevin2016zeb, P_dslevin2016zeb)
