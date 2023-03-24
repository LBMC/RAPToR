geo_id <- "GSE93826"
geo_obj <- GEOquery::getGEO(geo_id)[[1]]

# get expression data
sfs <- GEOquery::getGEOSuppFiles(geo_id, fetch_files = F, makeDirectory = F)

tmpfolder <- file.path(data_folder, "byrne2020tmp")
dir.create(tmpfolder)
destfile <- file.path(tmpfolder, sfs$fname[1])
download.file(sfs$url[1], destfile = destfile)

untar(destfile, exdir = tmpfolder)

flist <- list.files(tmpfolder)
g <- lapply(seq_along(flist)[-1], function(i){
  gi <- read.table(gzfile(file.path(tmpfolder,flist[i])), 
                   h=F, row.names = 1)
  colnames(gi) <- gsub("(GSM246\\d+)_.*", "\\1", flist[i])
  return(gi)
})

g <- do.call(cbind, g)

# format and convert to tpm
g <- RAPToR::format_ids(g, wormRef::Cel_genes, from = "wb_id", 
                        to = "wb_id", aggr.fun = sum)

g <- raw2tpm(g, wormRef::Cel_genes$transcript_length[
  match(rownames(g), wormRef::Cel_genes$wb_id)])
g <- g[apply(g, 1, sum) > 0, ] # keep expressed genes only
g <- g[-3123, ] # remove 0-variance artefact gene

# get pheno data
p <-  Biobase::pData(geo_obj)
p <- p[, c("title", "geo_accession", "age:ch1")]
colnames(p)[3] <- "age"
p$age <- as.numeric(as.character(
  factor(p$age, levels = paste0('day ', c(3,6,9,12,15,18)),
         labels = c(3,6,9,12,15,18))
  ))

p <- p[order(p$age), ]
g <- g[, p$geo_accession]

# save data
dsbyrne2020 <- list(g=g, p=p)
save(dsbyrne2020, file = file.path(data_folder, "dsbyrne2020.RData"), 
     compress = "xz")

# cleanup
file.remove(file.path(tmpfolder, flist))
file.remove(tmpfolder, destfile)
rm(g, p, tmpfolder, flist, geo_id, geo_obj, sfs, destfile)

