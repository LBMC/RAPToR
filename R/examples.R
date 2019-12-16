ae_example <- function(){
  ex <- "
@examples
\\donttest{
requireNamespace('wormRef', quietly = TRUE)

# get some samples to stage
samp <- wormRef::Cel_larval$g[,13:15]

# load interpolated reference
r_larv <- prepare_refdata('Cel_larval', n.inter = 200)

# perform age estimate
ae_test <- ae(samp = samp, 
              refdata = r_larv$interpGE, 
              ref.time_series = r_larv$time.series)

# check output
summary(ae_test)
plot(ae_test, show.boot_estimates = TRUE) # plot all sample estimates
plot_cor.ae(ae_test) # plot individual correlation profiles of samples

# get results
ae_test$age.estimates
}
"
  
  return(strsplit(ex, split = "\n")[[1]])
}


interpol_example <- function(){
  ex <- "
@examples
\\donttest{
requireNamespace('wormRef', quietly = TRUE)

# gene expression data
G <- wormRef::Cel_larval$g

# pheno data (e.g time, batch)
P <- wormRef::Cel_larval$p



# find optimal spline degree of freedom for PLSR interpolation
dfCV_G <- df_CV(X = G, 
                time.series = P$age,
                covar = P$cov,
                dfs = 10:20)

# check minimum of CV error
summary(dfCV_G)
plot(dfCV_G) # min for df = 17

# perform PLSR interpolation using optimal parameter
r_G <- plsr_interpol(X = G, 
                     time.series = P$age,
                     covar = P$cov,
                     df = 17) 

}
"

return(strsplit(ex, split = "\n")[[1]])
}