rm(list=ls())

library(GenABEL)
library(survival)

################ input #####################

tfamfile1="/lustre/scr/j/i/jinjin/GRA/meta_analysis/80405_geno_1.tfam"
tools::md5sum(tfamfile1)

tpedfile1="/lustre/scr/j/i/jinjin/GRA/meta_analysis/80405_geno_1.tped"
tools::md5sum(tpedfile1)

tfamfile2="/lustre/scr/j/i/jinjin/GRA/meta_analysis/80405_geno_2.tfam"
tools::md5sum(tfamfile2)

tpedfile2="/lustre/scr/j/i/jinjin/GRA/meta_analysis/80405_geno_2.tped"
tools::md5sum(tpedfile2)

phenofile="/lustre/scr/j/i/jinjin/GRA/meta_analysis/80405_pheno.txt"
tools::md5sum(phenofile)

convert.snp.tped(tpedfile = tpedfile1, tfamfile = tfamfile1, outfile = "genotype1.raw")
convert.snp.tped(tpedfile = tpedfile2, tfamfile = tfamfile2, outfile = "genotype2.raw")
#merge.snp.data("genotype1.raw", "genotype2.raw")
df1 <- load.gwaa.data(phe = phenofile, gen = "genotype1.raw", force = T)
df2 <- load.gwaa.data(phe = phenofile, gen = "genotype2.raw", force = T)
df1 <- df1[phdata(df1)$cau == 1, ] 
df2 <- df2[phdata(df2)$cau == 1, ] 
df1 <- df1[phdata(df1)$primary == 1, ]
df2 <- df2[phdata(df2)$primary == 1, ] 
mc1 = check.marker(df1, callrate = 0.95, extr.call = 0.95, maf = 0.05, 
                  p.level = 1e-08, het.fdr = 0, ibs.exclude = "lower")
df1 = df1[mc1$idok, mc1$snpok]
mc2 = check.marker(df2, callrate = 0.95, extr.call = 0.95, maf = 0.05,
                  p.level = 1e-08, het.fdr = 0, ibs.exclude = "lower")
df2 = df2[mc2$idok, mc2$snpok]
nsnps(gtdata(df1))
nsnps(gtdata(df2))

df0 = merge(df1, df2)

df1_m = df1[, is.element(snpnames(gtdata(df1)), snpnames(gtdata(df2)))]
df2_m = df2[, is.element(snpnames(gtdata(df2)), snpnames(gtdata(df1)))]
df = merge(df1_m, df2_m)
mc = check.marker(df, callrate = 0.95, extr.call = 0.95, maf = 0.05,
                  p.level = 1e-08, het.fdr = 0, ibs.exclude = "lower")
df = df[mc$idok, mc$snpok]
nsnps(gtdata(df))

df1 = df1[, !is.element(snpnames(gtdata(df1)), snpnames(gtdata(df1_m)))]
df2 = df2[, !is.element(snpnames(gtdata(df2)), snpnames(gtdata(df2_m)))]
nsnps(gtdata(df1))
nsnps(gtdata(df2))
mc1 = check.marker(df1, callrate = 0.95, extr.call = 0.95, maf = 0.05,
                   p.level = 1e-08, het.fdr = 0, ibs.exclude = "lower")
df1 = df1[mc1$idok, mc1$snpok]
mc2 = check.marker(df2, callrate = 0.95, extr.call = 0.95, maf = 0.05,
                   p.level = 1e-08, het.fdr = 0, ibs.exclude = "lower")
df2 = df2[mc2$idok, mc2$snpok]
nsnps(gtdata(df1))
nsnps(gtdata(df2))

rm(df1_m)
rm(df2_m)


#############################################################

################ plot #######################################

time <- phdata(df0)$ostime
status <- phdata(df0)$osevent

fit1 <- survfit(Surv(time, status) ~ phdata(df0)$Bev)
png("/lustre/scr/j/i/jinjin/GRA/meta_analysis/result/80405.png", width=1000, height=1000, pointsize=18)
plot(fit1, lty=as.numeric(as.factor(fit1$strata)), xlab = "Time since Randomization(Months)", ylab = "Probability", cex = 1)
legend(30, 0.9, substring(names(fit1$strata), 6, nchar(names(fit1$strata))), lty=as.numeric(as.factor(fit1$strata)))
dev.off()

################ get statistics #############################

class(df)
names(phdata(df))
descriptives.trait(df)
descriptives.marker(df)
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

save("df" = df, "df0" = df0, 
     file = "/lustre/scr/j/i/jinjin/GRA/meta_analysis/80405.Rdata")


for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()))
  int <- snps[,i]*phdata(df)$Bev
  fit <- summary(coxph(Surv(time, status) ~ int + phdata(df)$Bev
                       + snps[,i] + loc1 + loc2 + loc3 + loc4 
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
                     bp = bp, coef = coef, sd = v, maincoef = maincoef, mainsd = mainv,
                     pvalue = p, n = ss, maf = maf)

################ output ######################################

write.table(result, "/lustre/scr/j/i/jinjin/GRA/meta_analysis/result/80405.txt", 
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

save("df1" = df1, file = "/lustre/scr/j/i/jinjin/GRA/meta_analysis/80405.Rdata")
for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()))
  int <- snps[,i]*phdata(df1)$Bev
  fit <- summary(coxph(Surv(time, status) ~ int + phdata(df1)$Bev
                       + snps[,i] + loc1 + loc2 + loc3 + loc4 
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
                     bp = bp, coef = coef, sd = v, maincoef = maincoef, mainsd = mainv,
                     pvalue = p, n = ss, maf = maf)

################ output ######################################

write.table(result, "/lustre/scr/j/i/jinjin/GRA/meta_analysis/result/80405_1.txt", 
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

save("df2" = df2, file = "/lustre/scr/j/i/jinjin/GRA/meta_analysis/80405.Rdata")
for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()))
  int <- snps[,i]*phdata(df2)$Bev
  fit <- summary(coxph(Surv(time, status) ~ int + phdata(df2)$Bev
                       + snps[,i] + loc1 + loc2 + loc3 + loc4 
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
                     bp = bp, coef = coef, sd = v, maincoef = maincoef, mainsd = mainv,
                     pvalue = p, n = ss, maf = maf)

################ output ######################################

write.table(result, "/lustre/scr/j/i/jinjin/GRA/meta_analysis/result/80405_2.txt", 
            sep="\t")
