rm(list=ls())

library(GenABEL)
library(survival)

################ input #####################

load("/lustre/scr/j/i/jinjin/GRA/meta_analysis/80405.Rdata")
load("/lustre/scr/j/i/jinjin/GRA/meta_analysis/80405_1.Rdata")
load("/lustre/scr/j/i/jinjin/GRA/meta_analysis/80405_2.Rdata")


################ analysis ###################

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

loc1 <- as.numeric(phdata(df)$tumorside == "Left")
loc2 <- as.numeric(phdata(df)$tumorside == "Right")
loc3 <- as.numeric(phdata(df)$tumorside == "Multiple")
loc4 <- as.numeric(phdata(df)$tumorside == "Transverse")



for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()))
  int <- snps[,i]*phdata(df)$Bev
  fit <- summary(coxph(Surv(time, status) ~ snps[,i] + int + phdata(df)$Bev
                       + loc1 + loc2 + loc3 + loc4 
                       + strata(phdata(df)$prot_chemo, phdata(df)$pr_adj, phdata(df)$pr_rad)))
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
                     bp = bp, coef = coef, sd = v, intcoef = maincoef, intsd = mainv,
                     pvalue = p, n = ss, maf = maf)

################ output ######################################

write.table(result, "/lustre/scr/j/i/jinjin/GRA/meta_analysis/result_main/80405.txt", 
            sep="\t")

################ get statistics #############################

class(df1)
names(phdata(df1))
descriptives.trait(df1)
descriptives.marker(df1)
nsnps <- nsnps(gtdata(df1))


snps <- as.numeric(gtdata(df1))
chr <- c(as.numeric(chromosome(gtdata(df1))))
snp <- c(colnames(snps))
bp <- c(as.numeric(map(gtdata(df1))))

time <- phdata(df1)$ostime
status <- phdata(df1)$osevent

coef <- rep(0, nsnps)
v <- rep(0, nsnps)
p <- rep(0, nsnps)
ss <- rep(0, nsnps)
maf <- rep(0, nsnps)
maincoef <- rep(0, nsnps)
mainv <- rep(0, nsnps)

loc1 <- as.numeric(phdata(df1)$tumorside == "Left")
loc2 <- as.numeric(phdata(df1)$tumorside == "Right")
loc3 <- as.numeric(phdata(df1)$tumorside == "Multiple")
loc4 <- as.numeric(phdata(df1)$tumorside == "Transverse")

for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()))
  int <- snps[,i]*phdata(df1)$Bev
  fit <- summary(coxph(Surv(time, status) ~ snps[,i] + int + phdata(df1)$Bev
                       + loc1 + loc2 + loc3 + loc4 
                       + strata(phdata(df1)$prot_chemo, phdata(df1)$pr_adj, phdata(df1)$pr_rad)))
  coef[i] <- fit$coef[1, 1]
  v[i] <- fit$coef[1,3]
  p[i] <- fit$coef[1,5]
  ss[i] <- fit$n
  maincoef[i] <- fit$coef[2, 1]
  mainv[i] <- fit$coef[2, 3]
}

refallele <- refallele(gtdata(df1))
effallele <- effallele(gtdata(df1))

sumdf1 <- summary(gtdata(df1))
afr <- sumdf1[, "Q.2"]
maf <- pmin(afr, (1. - afr))

result <- data.frame(snp = snp, chr = chr, refallele = refallele, effallele = effallele,
                     bp = bp, coef = coef, sd = v, intcoef = maincoef, intsd = mainv,
                     pvalue = p, n = ss, maf = maf)

################ output ######################################

write.table(result, "/lustre/scr/j/i/jinjin/GRA/meta_analysis/result_main/80405_1.txt", 
            sep="\t")

################ get statistics #############################

class(df2)
names(phdata(df2))
descriptives.trait(df2)
descriptives.marker(df2)
nsnps <- nsnps(gtdata(df2))


snps <- as.numeric(gtdata(df2))
chr <- c(as.numeric(chromosome(gtdata(df2))))
snp <- c(colnames(snps))
bp <- c(as.numeric(map(gtdata(df2))))

time <- phdata(df2)$ostime
status <- phdata(df2)$osevent

coef <- rep(0, nsnps)
v <- rep(0, nsnps)
p <- rep(0, nsnps)
ss <- rep(0, nsnps)
maf <- rep(0, nsnps)
maincoef <- rep(0, nsnps)
mainv <- rep(0, nsnps)

loc1 <- as.numeric(phdata(df2)$tumorside == "Left")
loc2 <- as.numeric(phdata(df2)$tumorside == "Right")
loc3 <- as.numeric(phdata(df2)$tumorside == "Multiple")
loc4 <- as.numeric(phdata(df2)$tumorside == "Transverse")

for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()))
  int <- snps[,i]*phdata(df2)$Bev
  fit <- summary(coxph(Surv(time, status) ~ snps[,i] + int + phdata(df2)$Bev
                       + loc1 + loc2 + loc3 + loc4 
                       + strata(phdata(df2)$prot_chemo, phdata(df2)$pr_adj, phdata(df2)$pr_rad)))
  coef[i] <- fit$coef[1, 1]
  v[i] <- fit$coef[1,3]
  p[i] <- fit$coef[1,5]
  ss[i] <- fit$n
  maincoef[i] <- fit$coef[2, 1]
  mainv[i] <- fit$coef[2, 3]
}

refallele <- refallele(gtdata(df2))
effallele <- effallele(gtdata(df2))

sumdf2 <- summary(gtdata(df2))
afr <- sumdf2[, "Q.2"]
maf <- pmin(afr, (1. - afr))

result <- data.frame(snp = snp, chr = chr, refallele = refallele, effallele = effallele,
                     bp = bp, coef = coef, sd = v, intcoef = maincoef, intsd = mainv,
                     pvalue = p, n = ss, maf = maf)

################ output ######################################

write.table(result, "/lustre/scr/j/i/jinjin/GRA/meta_analysis/result_main/80405_2.txt", 
            sep="\t")
