geo_dslevin2016dmel <- "GSE60471"
g_url_dslevin2016dmel <- GEOquery::getGEOSuppFiles(geo_dslevin2016dmel, 
                                                   makeDirectory = FALSE, 
                                                   fetch_files = FALSE)
g_file_dslevin2016dmel <- paste0(data_folder, "dslevin2016dmel.txt.gz")
utils::download.file(url = as.character(g_url_dslevin2016dmel$url[3]), 
                     destfile = g_file_dslevin2016dmel)

X_dslevin2016dmel <- read.table(gzfile(g_file_dslevin2016dmel), h = T, 
                                sep = '\t', as.is = T, row.names = 1, 
                                comment.char = "")

# filter poor quality samples
cm_dslevin2016dmel <- RAPToR::cor.gene_expr(X_dslevin2016dmel, 
                                            X_dslevin2016dmel)
f_dslevin2016dmel <- which(0.6 > apply(cm_dslevin2016dmel, 1, 
                                       quantile, probs = .99))
X_dslevin2016dmel <- X_dslevin2016dmel[, -f_dslevin2016dmel]

# convert to tpm & FBgn
X_dslevin2016dmel <- X_dslevin2016dmel[
  rownames(X_dslevin2016dmel)%in%droso_genes$fb_id,]
X_dslevin2016dmel <- raw2tpm(
  rawcounts = X_dslevin2016dmel, 
  genelengths = droso_genes$transcript_length[
    match(rownames(X_dslevin2016dmel), droso_genes$fb_id)])

# pheno data
P_dslevin2016dmel <- Biobase::pData(GEOquery::getGEO(geo_dslevin2016dmel,
                                                     getGPL = F)[[1]])

# filter relevant fields/samples
P_dslevin2016dmel <- P_dslevin2016dmel[, 
  c("title", "geo_accession", "time (minutes cellularization stage):ch1")]
colnames(P_dslevin2016dmel)[3] <- "time"
P_dslevin2016dmel$title <- as.character(P_dslevin2016dmel$title)

P_dslevin2016dmel <- P_dslevin2016dmel[
  P_dslevin2016dmel$title %in% colnames(X_dslevin2016dmel),]
X_dslevin2016dmel <- X_dslevin2016dmel[, P_dslevin2016dmel$title]

# formatting
P_dslevin2016dmel$title <- gsub('Metazome_Drosophila_timecourse_', '', 
                                P_dslevin2016dmel$title)
colnames(X_dslevin2016dmel) <- P_dslevin2016dmel$title
P_dslevin2016dmel$age <- as.numeric(P_dslevin2016dmel$time) / 60

# save data
dslevin2016dmel <- list(g = X_dslevin2016dmel, p = P_dslevin2016dmel)
save(dslevin2016dmel, file = paste0(data_folder, "dslevin2016dmel.RData"), compress = "xz")

# cleanup
file.remove(g_file_dslevin2016dmel)
rm(geo_dslevin2016dmel, g_url_dslevin2016dmel, g_file_dslevin2016dmel, 
   X_dslevin2016dmel, P_dslevin2016dmel, 
   cm_dslevin2016dmel, f_dslevin2016dmel)
