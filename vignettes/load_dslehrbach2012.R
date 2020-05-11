# get probe ids from the arrayexpress
probe_ids <- read.table("https://www.ebi.ac.uk/arrayexpress/files/A-AFFY-60/A-AFFY-60.adf.txt", h=T,
                        comment.char = "", sep = "\t", skip = 17, as.is = T)


# pheno data
p_url <- "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1333/E-MTAB-1333.sdrf.txt"
p <- read.table(p_url, h = T, sep = "\t", comment.char = "", as.is = T)

# geno data
g_url <- unique(p$Comment..ArrayExpress.FTP.file.)
g_dir <- paste0(data_folder, "dslehrbach2012_raw/")
g_file <- paste0(g_dir, "raw.zip")

if(!file.exists(g_file)){
  dir.create(g_dir)
  utils::download.file(g_url, destfile = g_file)
}
utils::unzip(g_file, exdir = g_dir)

g <- affy::ReadAffy(filenames = p$Array.Data.File, celfile.path = g_dir, phenoData = p)
g <- affy::expresso(g, bg.correct = F, normalize = F,
                    pmcorrect.method = "pmonly", summary.method = "median")
g <- 2^exprs(g) # expresso log2s the data
colnames(g) <- p$Source.Name

# format ids
g <- RAPToR::format_ids(g, probe_ids, from = 1, to = 7)
g <- RAPToR::format_ids(g, wormRef::Cel_genes, from = 3, to = 1)


# filter relevant fields
p <- p[, c(1,22,24)]
colnames(p) <- c("title", "tpastL4", "strain")
p$strain <- factor(p$strain, levels = c('wild type', 'pash-1(mj100)'), labels = c('wt', 'strain'))
p$rep <- factor(gsub("\\w+_h\\d+\\.(\\d)","\\1", p$title))

dslehrbach2012 <- list(g = g, p = p)
save(dslehrbach2012, file = paste0(data_folder, "dslerhbach2012.RData"), compress = "xz")

# cleanup
file.remove(g_file)
unlink(g_dir, recursive = T)

rm(g, g_url, g_dir, g_file, p, p_url, probe_ids)