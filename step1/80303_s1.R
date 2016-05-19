rm(list=ls())

library(GenABEL)
library(survival)


################ input #####################

tfamfile="/lustre/scr/j/i/jinjin/GRA/meta_analysis/80303_geno.tfam"
tools::md5sum(tfamfile)

tpedfile="/lustre/scr/j/i/jinjin/GRA/meta_analysis/80303_geno.tped"
tools::md5sum(tpedfile)

phenofile="/lustre/scr/j/i/jinjin/GRA/meta_analysis/80303_pheno.txt"
tools::md5sum(phenofile)

convert.snp.tped(tpedfile = tpedfile, tfamfile = tfamfile, outfile = "genotype.raw")
df <- load.gwaa.data(phe = phenofile, gen = "genotype.raw", force = T)

mc = check.marker(df, callrate = 0.95, extr.call = 0.95, maf = 0.05,
                  p.level = 1e-08, het.fdr = 0, ibs.exclude = "lower")
df = df[mc$idok, mc$snpok]

class(df)
names(phdata(df))
descriptives.trait(df)
descriptives.marker(df)
nsnps <- nsnps(gtdata(df))

############ MAF ##########################################
#mc = check.marker(df, callrate = 0.95, extr.call = 0.95, 
#                  p.level = 1e-08, het.fdr = 0, maf = 0.01)
#df2 = df[, !is.element(snpnames(gtdata(df)), mc$nofreq)]
#nsnps(gtdata(df2))
#
#mc = check.marker(df, callrate = 0.95, extr.call = 0.95, 
#                  p.level = 1e-08, het.fdr = 0, maf = 0.01)
#df2 = df[, !is.element(snpnames(gtdata(df)), mc$nofreq)]
#nsnps(gtdata(df2))
#############################################################

################ get statistics #############################

snps <- as.numeric(gtdata(df))
chr <- c(as.numeric(chromosome(gtdata(df))))
snp <- c(colnames(snps))
bp <- c(as.numeric(map(gtdata(df))))

time <- phdata(df)$ostime
status <- phdata(df)$osevent

coef <- rep(0, nsnps)
v <- rep(0, nsnps)
p <- rep(0, nsnps)
ss <- rep(0, nsnps)
maf <- rep(0, nsnps) 
maincoef <- rep(0, nsnps)
mainv <- rep(0, nsnps)

save("df" = df, file = "/lustre/scr/j/i/jinjin/GRA/meta_analysis/80303.Rdata")

arm <- phdata(df)$arm - 1
fit1 <- survfit(Surv(time, status) ~ arm)
png("/lustre/scr/j/i/jinjin/GRA/meta_analysis/result/80303.png", width=1000, height=1000, pointsize=18)
plot(fit1, lty=as.numeric(as.factor(fit1$strata)), xlab = "Time since Randomization(Months)", ylab = "Probability", cex = 1)
legend(30, 0.9, substring(names(fit1$strata), 6, nchar(names(fit1$strata))), lty=as.numeric(as.factor(fit1$strata)))
dev.off()

for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()))
  int <- snps[,i]*arm
  fit <- summary(coxph(Surv(time, status) ~ int + arm 
                       + snps[,i] + strata(phdata(df)$priorrt, phdata(df)$dzext,
                                           phdata(df)$ps) ))
  coef[i] <- fit$coef[1, 1]
  v[i] <- fit$coef[1,3]
  p[i] <- fit$coef[1,5]
  ss[i] <- fit$n
  maincoef[i] <- fit$coef[2, 1]
  mainv[i] <- fit$coef[2, 3]
}

refallele <- refallele(gtdata(df))
effallele <- effallele(gtdata(df))
sumdf <- summary(gtdata(df))
afr <- sumdf[, "Q.2"]
maf <- pmin(afr, (1. - afr))

result <- data.frame(snp = snp, chr = chr, refallele = refallele, effallele = effallele,
                     bp = bp, coef = coef, sd = v, maincoef = maincoef, mainsd = mainv,
                     pvalue = p, n = ss, maf = maf)

################ output ######################################

write.table(result, "/lustre/scr/j/i/jinjin/GRA/meta_analysis/result/80303.txt", 
            sep="\t")
