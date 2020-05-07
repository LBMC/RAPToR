g_url_dsgraveley2011 <- "ftp://ftp.fruitfly.org/pub/download/modencode_expression_scores/Celniker_Drosophila_Annotation_20120616_1428_allsamps_MEAN_gene_expression.csv.gz"
g_file_dsgraveley2011 <- paste0(data_folder, "dsgraveley2011.csv.gz")
utils::download.file(g_url_dsgraveley2011, destfile = g_file_dsgraveley2011)


X_dsgraveley2011 <- read.table(gzfile(g_file_dsgraveley2011), sep = ',', row.names = 1, h = T)

# convert gene ids to FBgn
X_dsgraveley2011 <- RAPToR::format_ids(X_dsgraveley2011, droso_genes, from = "gene_name", to = "fb_id")

# select embryo time series samples
X_dsgraveley2011 <- X_dsgraveley2011[,1:12]

P_dsgraveley2011 <- data.frame(sname = colnames(X_dsgraveley2011),
                               age = as.numeric(gsub("em(\\d+)\\.\\d+hr", "\\1", colnames(X_dsgraveley2011))),
                               stringsAsFactors = FALSE)

dsgraveley2011 <- list(g = X_dsgraveley2011, p = P_dsgraveley2011)

save(dsgraveley2011, file = paste0(data_folder, "dsgraveley2011.RData"), compress = "xz")

# cleanup
file.remove(g_file_dsgraveley2011)
rm(g_url_dsgraveley2011, g_file_dsgraveley2011, X_dsgraveley2011, P_dsgraveley2011)