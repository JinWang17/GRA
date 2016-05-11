rm(list=ls())

library(GenABEL)
library(survival)

################ input #####################

tfamfile="/lustre/scr/j/i/jinjin/GRA/meta_analysis/40503_geno.tfam"
tools::md5sum(tfamfile)

tpedfile="/lustre/scr/j/i/jinjin/GRA/meta_analysis/40503_geno.tped"
tools::md5sum(tpedfile)

phenofile="/lustre/scr/j/i/jinjin/GRA/meta_analysis/40503_pheno.txt"
tools::md5sum(phenofile)

convert.snp.tped(tpedfile = tpedfile, tfamfile = tfamfile, outfile = "genotype.raw")
df <- load.gwaa.data(phe = phenofile, gen = "genotype.raw", force = T)
df <- df[phdata(df)$cau == 1, ]
df <- df[phdata(df)$treated == 1, ]

mc = check.marker(df, callrate = 0.95, extr.call = 0.95, 
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

for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()));
  fit <- summary(coxph(Surv(time, status) ~ snps[,i]*phdata(df)$arm + snps[,i] 
                       + phdata(df)$arm + strata(phdata(df)$stra2, phdata(df)$stra3) ))
  coef[i] <- fit$coef[1, 1]
  v[i] <- fit$coef[1,3]^2
  p[i] <- fit$coef[1,5]
  ss[i] <- fit$n
}

refallele <- refallele(gtdata(df))
effallele <- effallele(gtdata(df))

result <- data.frame(snp = snp, chr = chr, refallele = refallele, effallele = effallele,
                     bp = bp, coef = coef, sd = sqrt(v), pvalue = p, n = ss)

################ output ######################################

write.table(result, "/lustre/scr/j/i/jinjin/GRA/meta_analysis/result/40503.txt", 
            sep="\t")
