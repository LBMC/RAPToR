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
plot_cor.ae(age.est, show.init_estimate = T)


# plot age estimates
par(mfrow=c(1,1))
plot(age.est, col.i=4)


# plot age vs pheno data
par(mfrow=c(2,1))

plot(age.est, groups = gorrepati.pheno$strain, main="By strain",
     col.i=4, col= as.numeric(gorrepati.pheno$strain))

plot(age.est, groups = gorrepati.pheno$batch, main="By batch",
     col.i=4, col=gorrepati.pheno$batch)


# print time window
diff(range(age.est$age.estimates[,1]))


