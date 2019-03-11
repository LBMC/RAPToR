library(wormAge)
par(mfrow=c(1,1))


# load oudenaarden ref data
data(oud_ref)

# compute interpolation
interp.oud <- interpol_refdata(oud_ref$X, 150,
                               time.series = oud_ref$time.series,
                               t.min = min(oud_ref$time.series),
                               t.max = max(oud_ref$time.series))

# get sample dataset (here, Gorrepati et al. 2015)
gorrepati.expr_data <- GPL200_to_WB("GSE52747")

overlap.set <- format_to_ref(gorrepati.expr_data, interp.oud$interpol.gene_expr)

age.est <- estimate.worm_age(samp = overlap.set$samp, 
                             refdata = overlap.set$refdata, 
                             ref.time_series = interp.oud$time.series,
                             time.sd = 7, bootstrap.n = 50,
                             est.time = 27, bootstrap.time_window = 5)

age.est$age.estimates
plot(age.est)

plot(age.est$init.est.times, age.est$age.estimates[,1])



