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
                             est.time = 27)

age.est$age.estimates


par(mfrow=c(4,2), mar=c(4,3,3,1))
pb <- sapply(1:ncol(age.est$age.estimates), function(i){
  
  plot(interp.oud$time.series, age.est$cors[,i], type = 'l', lwd=2,
       main=colnames(age.est$cors)[i], xlab = 'time', ylab='corr.score')
  
  ae <- t(age.est$age.estimates[c(1,2),i])
  points(ae, pch='|', cex=2, col='firebrick')
  
  text(ae, pos=1, 
       labels = paste(round(ae[1], 2), 'h', sep=''),
       offset = 1)
})



