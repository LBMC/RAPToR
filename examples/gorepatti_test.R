library(wormAge)


# get sample dataset (here, Gorrepati et al. 2015)
gorrepati.expr_data <- GPL200_to_WB("GSE52747")

gorrepati.pheno <- read.table(textConnection(
  "id 	batch 	strain
  GSM2092609 	1 	N2
  GSM2092610 	1 	N2
  GSM2092611 	1 	N2
  GSM2092612 	2 	Huis1
  GSM2092613 	2 	DeltaNTPOP1
  GSM2092614 	2 	Huis1
  GSM2092615 	2 	DeltaNTPOP1
  GSM2092616 	2 	Huis1
  GSM2092617 	2 	DeltaNTPOP1
  "), sep='\t', h=T, row.names=1)
gorrepati.pheno <- gorrepati.pheno[colnames(gorrepati.expr_data),]

guess <- 26 # the paper states the worms were grown for 26 hrs at 20C



# load and prepare oudenaarden ref data
interp.oud <- prepare_refdata("oudenaarden")


# Do age estimation
age.est <- estimate.worm_age(samp = gorrepati.expr_data, 
                             refdata = interp.oud$interpol.gene_expr, 
                             ref.time_series = interp.oud$time.series,
                             time.sd = 5, bootstrap.n = 50,
                             est.time = guess, bootstrap.time_window = 1)

age.est$age.estimates

# plot correlation profiles
par(mfrow=c(3,2))
plot(age.est, show.init_estimate = T)




# plot age vs pheno data
p <- par(mfrow=c(1,2), mar=c(7,4,3,1))
boxplot(age.est$age.estimates[,1]~gorrepati.pheno$strain,
        main="Estimated age per strain", las=2)
require(beeswarm)
beeswarm(age.est$age.estimates[,1]~gorrepati.pheno$strain, 
         col = 2, pch=16, cex=1.5, lwd=2, add=T)


boxplot(age.est$age.estimates[,1]~gorrepati.pheno$batch, 
        main="Estimated age per batch")
beeswarm(age.est$age.estimates[,1]~gorrepati.pheno$batch, 
         col = 2, pch=16, cex=1.5, lwd=2, add=T)
par(p)


# print time window
diff(range(age.est$age.estimates[,1]))

