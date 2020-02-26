ae_example <- function(){
  ex <- "
@examples
\\donttest{
requireNamespace('wormRef', quietly = TRUE)

# get some samples to stage
samp <- wormRef::Cel_larval$g[,13:15]

# load interpolated reference
r_larv <- prepare_refdata(ref = 'Cel_larval', datapkg = 'wormRef' , n.inter = 200)

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
requireNamespace('stats', quietly = TRUE)

# gene expression data
X <- wormRef::Cel_larval$g

# pheno data (e.g time, batch)
p <- wormRef::Cel_larval$p

# do a pca & select nb of components to use for interpol
pca <- stats::prcomp(X, rank = 20)
nc <- sum(summary(pca)$importance[3, ] < .999) + 1


# find optimal spline type
# setup formulas
smooths <- c('bs', 'tp', 'cr', 'ds')
flist <- as.list(paste0('X ~ s(age, bs = \\'', smooths, '\\') + cov'))
# do CV
cvres <- ge_imCV(X = scale(X), p = p, formula_list = flist,
                 cv.n = 20, nc = nc)
# check results
plot(cvres, names.arrange = 4) # lowest pred error with 'ds' spline

# build model & interpolation data
m <- ge_im(X = X, p = p, formula = 'X ~ s(age, bs = \\'ds\\') + cov', nc = nc)
n.inter = 100
ndat <- data.frame(age = seq(min(p$age), max(p$age), l = n.inter),
                   cov = rep(p$cov[1], n.inter))

# check interpolation on pca components
pred_pca <- predict(m, ndat, as.c = T)

par(mfrow = c(2,2))
invisible(sapply(seq_len(nc), function(i){
  plot(p$age, pca$rotation[, i], xlab = 'age', ylab = 'PC', main = paste0('PC',i),
       lwd = (p$cov == p$cov[1]) + 1)
  points(ndat$age, pred_pca[, i], type = 'l', lwd = 2)
}))

# get interpolated GE matrix, as a reference
r_X <- list(interpGE = predict(m, ndat), time.series = ndat$age)

# test
ae_X <- ae(X, r_X$interpGE, r_X$time.series)
par(mfrow = c(1,1))
plot(p$age, ae_X$age.estimates[,1])


}
"

return(strsplit(ex, split = "\n")[[1]])
}