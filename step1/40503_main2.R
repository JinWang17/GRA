######### test main effect of snp with interaction in the model ##############

rm(list=ls())

library(GenABEL)
library(survival)

################ input #####################

load("/lustre/scr/j/i/jinjin/GRA/meta_analysis/40503.Rdata")

################ analysis ##################

nsnps <- nsnps(gtdata(df))

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

for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()));
  int <- snps[,i]*phdata(df)$arm
  fit <- summary(coxph(Surv(time, status) ~ snps[,i] 
                       + phdata(df)$arm + strata(phdata(df)$stra2, phdata(df)$stra3) ))
  coef[i] <- fit$coef[1, 1]
  v[i] <- fit$coef[1, 3]
  p[i] <- fit$coef[1, 5]
  ss[i] <- fit$n
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

write.table(result, "/lustre/scr/j/i/jinjin/GRA/meta_analysis/result_main2/40503.txt", 
            sep="\t")
