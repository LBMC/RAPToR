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
ae_test <- ae(samp, r_larv)

# check output
summary(ae_test)
plot(ae_test) # plot all sample estimates
plot_cor(ae_test) # plot individual correlation profiles of samples

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

# build model & make reference
m <- ge_im(X = X, p = p, formula = 'X ~ s(age, bs = \\'ds\\') + cov', nc = nc)

ref <- make_ref(m, cov.levels = list('cov'='O.20'), n.inter = 100, 
                t.unit='h past egg-laying (20C)')

# check model interpolation on pca components
par(mfrow = c(2,2))
plot(m, ref, ncs=1:4) # showing first 4 PCs


# test
ae_X <- ae(X, ref)
par(mfrow = c(1,2))
plot(p$age, ae_X$age.estimates[,1])
plot(ae_X, groups = p$cov)


}
"

return(strsplit(ex, split = "\n")[[1]])
}


rc_example <- function(){
  ex <- "
@examples
\\donttest{
requireNamespace('wormRef', quietly = TRUE)

# get sample gene expression data
X <- wormRef::Cel_larval$g[,1:9]

# get reference
ref <- prepare_refdata(ref = 'Cel_larval', datapkg = 'wormRef' , n.inter = 200)

# define groups
fac <- factor(c('a', 'a', 'b', 
                'a', 'b', 'b',
                'c', 'c', 'c'))

# estimate sample age
ae_X <- ae(X, ref)

# compare group diff. expr. with matching reference
rc <- ref_compare(X, ref, fac, ae_X)
print(rc)

# get sample and reference (ie. development) logFCs between groups
lfc_a_vs_b <- get_logFC(rc)
lfc_a_vs_c <- get_logFC(rc, l = 'c')
lfc_b_vs_c <- get_logFC(rc, l0='b', l = 'c')

# plot sample vs. reference logFCs
par(mfrow = c(2,2))
plot(ae_X, groups = fac)
plot(lfc_a_vs_b, main = 'a vs. b logFC')
plot(lfc_a_vs_c, main = 'a vs. c logFC')
plot(lfc_b_vs_c, main = 'b vs. c logFC')

}
" 

return(strsplit(ex, split = "\n")[[1]])
}