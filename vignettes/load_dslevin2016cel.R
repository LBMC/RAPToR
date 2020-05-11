geo_dslevin2016cel <- "GSE60755"

g_url_dslevin2016cel <- GEOquery::getGEOSuppFiles(geo_dslevin2016cel, makeDirectory = FALSE, fetch_files = FALSE)
g_file_dslevin2016cel <- paste0(data_folder, "dslevin2016cel.txt.gz")
utils::download.file(url = as.character(g_url_dslevin2016cel$url[1]), destfile = g_file_dslevin2016cel)

X_dslevin2016cel <- read.table(gzfile(g_file_dslevin2016cel), h = T, sep = '\t', as.is = T, row.names = 1, comment.char = "")

# filter poor quality samples
cm_dslevin2016cel <- RAPToR::cor.gene_expr(X_dslevin2016cel, X_dslevin2016cel)
f_dslevin2016cel <- which(0.67 > apply(cm_dslevin2016cel, 1, quantile, probs = .99))
X_dslevin2016cel <- X_dslevin2016cel[, -f_dslevin2016cel]


# convert to tpm & FBgn

X_dslevin2016cel <- X_dslevin2016cel[rownames(X_dslevin2016cel)%in%wormRef::Cel_genes$sequence_name,]
X_dslevin2016cel <- raw2tpm(rawcounts = X_dslevin2016cel, 
                             genelengths = wormRef::Cel_genes$transcript_length[match(rownames(X_dslevin2016cel),
                                                                                      wormRef::Cel_genes$sequence_name)])
X_dslevin2016cel <- RAPToR::format_ids(X_dslevin2016cel, wormRef::Cel_genes, from = "sequence_name", to = "wb_id")


# pheno data
P_dslevin2016cel <- Biobase::pData(GEOquery::getGEO(geo_dslevin2016cel, getGPL = F)[[1]])

# filter relevant fields/samples
P_dslevin2016cel <- P_dslevin2016cel[, c("title", "geo_accession", "time point (minutes after 4-cell):ch1")]
colnames(P_dslevin2016cel)[3] <- "time"
P_dslevin2016cel$title <- as.character(P_dslevin2016cel$title)

P_dslevin2016cel <- P_dslevin2016cel[P_dslevin2016cel$title %in% colnames(X_dslevin2016cel),]

# formatting
P_dslevin2016cel$age <- as.numeric(P_dslevin2016cel$time) / 60
P_dslevin2016cel <- P_dslevin2016cel[order(P_dslevin2016cel$age),]
X_dslevin2016cel <- X_dslevin2016cel[, P_dslevin2016cel$title]

P_dslevin2016cel$title <- gsub('Metazome_CE_timecourse_', '', P_dslevin2016cel$title)
colnames(X_dslevin2016cel) <- P_dslevin2016cel$title

X_dslevin2016cel <- X_dslevin2016cel[,-127] # remove extra outlier
P_dslevin2016cel <- P_dslevin2016cel[-127,]


dslevin2016cel <- list(g = X_dslevin2016cel, p = P_dslevin2016cel)
save(dslevin2016cel, file = paste0(data_folder, "dslevin2016cel.RData"), compress = "xz")

# cleanup
file.remove(g_file_dslevin2016cel)
rm(geo_dslevin2016cel, g_url_dslevin2016cel, g_file_dslevin2016cel, 
   f_dslevin2016cel, cm_dslevin2016cel, X_dslevin2016cel, P_dslevin2016cel)
