geo_dsmiki2017 <- "GSE97775"
url_dsmiki2017 <- as.character(GEOquery::getGEOSuppFiles(geo_dsmiki2017, makeDirectory = FALSE, fetch_files = FALSE)[9,"url"])
tmpf <- file.path(data_folder, "miki2017.txt.gz")
utils::download.file(url_dsmiki2017, destfile = tmpf)
g <- read.table(gzfile(tmpf), h=T, as.is = T, row.names = 1)
file.remove(tmpf)


# format genes to WB ids
g <- RAPToR::format_ids(g, wormRef::Cel_genes, from = "wb_id", to = "wb_id", aggr.fun = sum)


# store raw counts and convert to log(TPM+1)
g.raw <- g
g <- raw2tpm(g.raw, wormRef::Cel_genes$transcript_length[match(rownames(g.raw), wormRef::Cel_genes$wb_id)])
g <- log1p(g)


# sample metadata 
p <- Biobase::pData(GEOquery::getGEO(geo_dsmiki2017, getGPL = F)[[1]])[14:35, c(1,2, 45)]
colnames(p)[3] <- "strain_long"
p[,3] <- factor(p[, 3], levels = unique(p[,3]), labels = c("N2 (wild-type)", "HW1660 (xrn-2(xe31))"))
p$strain <- p[,3]
levels(p$strain) <- c("wt", "xrn2")

# get chronological age from sample name
p$age <- as.numeric(gsub(".*_(\\d+)h", "\\1", p$title))
p$title <- colnames(g)


dsmiki2017 <- list(g = g, p = p, g.raw = g.raw)

save(dsmiki2017, file = file.path(data_folder, "dsmiki2017.RData"), compress = "xz")
rm(url_dsmiki2017, geo_dsmiki2017, tmpf, g, g.raw, p, raw2tpm)
