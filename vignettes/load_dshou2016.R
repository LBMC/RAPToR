geo_id <- "GSE77110"
geo_obj <- GEOquery::getGEO(geo_id)[[1]]

# get pheno data
p <- Biobase::pData(geo_obj)
p <- p[, c("title", "geo_accession", "age:ch1", 
           "diet:ch1", "source_name_ch1")]
p$age <- as.numeric(gsub("adult day ", "", as.character(p$`age:ch1`)))
p <- p[, -3]
colnames(p) <- c("title", "geo_accession", "diet", "source_name", "age")
p$diet <- factor(p$diet, levels = c("ad libitum", 
                                    "calorie restriction", 
                                    "intermittent fasting"))

# get microarray probe IDs
mart <- biomaRt::useMart(biomart = "ensembl", 
                         dataset = "celegans_gene_ensembl")
probe_ids <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                           "affy_c_elegans"),
                            mart = mart)

# download data
sfile <- GEOquery::getGEOSuppFiles(geo_id, makeDirectory = F, fetch_files = F)
# You may need to increase the allowed connection size to download the dataset
# Sys.setenv("VROOM_CONNECTION_SIZE") <- 131072*4

tarfolder <- paste0(data_folder,"raw_hou")
dir.create(tarfolder, showWarnings = F)
tarfile <- paste0(tarfolder,'/', as.character(sfile$fname[1]))
utils::download.file(url = as.character(sfile$url[1]), destfile = tarfile)
untar(tarfile = tarfile, exdir = tarfolder)

flist <- paste0(tarfolder,'/',list.files(tarfolder))
for (f in flist[grepl(".gz", flist)]){
  GEOquery::gunzip(filename = f, destname = gsub("(.*).gz", "\\1", f), 
         overwrite = T, remove = T)
}

# load and format gdata
flist <- list.files(tarfolder)
flist <- flist[sapply(p$geo_accession, function(gg) which(grepl(gg, flist)))]
g <- affy::ReadAffy(filenames = flist, celfile.path = tarfolder, phenoData = p)
g <- affy::expresso(g, bg.correct = F, normalize = F,
                    pmcorrect.method = "pmonly", summary.method = "median")
g <- 2^Biobase::exprs(g) # expresso log2s the data
g <- RAPToR::format_ids(g, probe_ids, from = 2, to = 1)

# save data
dshou2016 <- list(g=g, p=p)
save(dshou2016, file = file.path(data_folder, "dshou2016.RData"), 
     compress = "xz")

# cleanup
file.remove(file.path(tarfolder, list.files(tarfolder)), tarfolder)
file.remove(tarfile)

rm(geo_id, geo_obj, g, p, probe_ids, mart, 
   tarfolder, tarfile, flist, sfile, f)
