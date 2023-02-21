geo_dsrockman2010 <- "GSE23857"
geo_dsrockman2010 <- GEOquery::getGEO(geo_dsrockman2010, GSEMatrix = F)

# get pdata
P_dsrockman2010 <- do.call(
  rbind, 
  lapply(GEOquery::GSMList(geo_dsrockman2010), function(go){
    unlist(GEOquery::Meta(go)[
      c("geo_accession", "characteristics_ch1", "characteristics_ch2",
        "label_ch1", "label_ch2")]
    )
  })
  )

P_dsrockman2010 <- as.data.frame(P_dsrockman2010, stringsAsFactors = F)

# get RG data from microarray
RG_dsrockman2010 <- list(
  R = do.call(cbind, lapply(GEOquery::GSMList(geo_dsrockman2010), function(go){
    GEOquery::Table(go)[, "R_MEAN_SIGNAL"]
  })),
  G = do.call(cbind, lapply(GEOquery::GSMList(geo_dsrockman2010), function(go){
    GEOquery::Table(go)[, "G_MEAN_SIGNAL"]
  }))
)

# normalize within channels
RG_dsrockman2010 <- limma::normalizeWithinArrays(RG_dsrockman2010, 
                                                 method = "loess")
RG_dsrockman2010 <- limma::RG.MA(RG_dsrockman2010) # convert back from MA to RG

# get only the RIALs (not the mixed stage controls)
X_dsrockman2010 <- RG_dsrockman2010$R
X_dsrockman2010[, P_dsrockman2010$label_ch1 == "Cy5"] <- 
  RG_dsrockman2010$G[, P_dsrockman2010$label_ch1 == "Cy5"]

# format probe/gene ids
gpl <- GEOquery::getGEO(GEOquery::Meta(
  GEOquery::GSMList(geo_dsrockman2010)[[1]]
  )$platform_id)
gpl <- GEOquery::Table(gpl)
gpl <- gpl[as.character(gpl$ID) %in% as.character(GEOquery::Table(
  GEOquery::GSMList(geo_dsrockman2010)[[1]]
  )$ID_REF), ]

sel <- gpl$ORF%in%wormRef::Cel_genes$sequence_name
gpl <- gpl[sel,]
X_dsrockman2010 <- X_dsrockman2010[sel,]


# filter bad quality samples
cm_dsrockman2010 <- cor(log1p(X_dsrockman2010), method = 'spearman')
f_dsrockman2010 <- which(0.95 > apply(cm_dsrockman2010, 1, 
                                      quantile, probs = .95))

X_dsrockman2010 <- X_dsrockman2010[,-f_dsrockman2010]
P_dsrockman2010 <- P_dsrockman2010[-f_dsrockman2010,]

# format ids
rownames(X_dsrockman2010) <- gpl$ORF
X_dsrockman2010 <- RAPToR::format_ids(X_dsrockman2010, wormRef::Cel_genes, 
                                      from = "sequence_name", to = "wb_id")


dsrockman2010 <- list(g = X_dsrockman2010, p = P_dsrockman2010)
save(dsrockman2010, 
     file = paste0(data_folder, "dsrockman2010.RData"), compress = "xz")

rm(geo_dsrockman2010, gpl, sel, 
   RG_dsrockman2010, X_dsrockman2010, P_dsrockman2010, 
   cm_dsrockman2010, f_dsrockman2010)