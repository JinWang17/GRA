rm(list=ls())

library(GenABEL)
library(cmprsk)
library(qqman)
require(graphics)

######### For subjects and snps alreday filtered and QCed ##########################

tfamfile="/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/plink.tfam"
tools::md5sum(tfamfile)

tpedfile="/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/step5.tped"
tools::md5sum(tpedfile)

phenofile="/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/neut3_slim.txt"
tools::md5sum(phenofile)

expidfile="/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/expidlist_294.txt"
tools::md5sum(expidfile)

convert.snp.tped(tpedfile = tpedfile, tfamfile = tfamfile, 
                 outfile = "genotype.raw")
df <- load.gwaa.data(phe = phenofile, gen = "genotype.raw", force = T)


class(df)
names(phdata(df))
attach(phdata(df))
descriptives.trait(df)
descriptives.marker(df)

nsnps <- nsnps(gtdata(df))


############## MAF > 0.015 #########################################################

mc = check.marker(df, callrate = 0.95, extr.call = 0.95, 
                  p.level = 1e-08, het.fdr = 0, maf = 0.015)
df2 = df[, !is.element(snpnames(gtdata(df)), mc$nofreq)]
nsnps(gtdata(df2))


############## preparing data ######################################################

nsnps <- nsnps(gtdata(df2))
nids <- nids(gtdata(df2))

arm <- phdata(df2)$Arm
str(arm)

snps <- as.numeric(gtdata(df2))
time <- phdata(df2)$event_time
status_1 <- phdata(df2)$event4
status_2 <- phdata(df2)$event3
status_0 <- phdata(df2)$event1

pvalue1 <- rep(0, nsnps)
converge1 <- rep(0, nsnps)
pvalue2 <- pvalue1
converge2 <- converge1
chr <- c(as.numeric(chromosome(gtdata(df2))))
snp <- c(colnames(snps))
bp <- c(as.numeric(map(gtdata(df2))))

pvalue1_en <- rep(0, nsnps)
converge1_en <- rep(0, nsnps)
pvalue2_en <- pvalue1_en
converge2_en <- converge1_en

before3m <- (time <= 91.3125)
time_en <- time
time_en[before3m] <- time[before3m]
time_en[!before3m] <- 91.3125

status_1_en <- status_1
status_2_en <- status_2
status_1_en[!before3m] <- 0
status_2_en[!before3m] <- 0

p_mar <- pvalue1
p_mar_en <- pvalue1_en

############## Survival plot ######################################################

# Figure 11
xx <- cuminc(time, status_2, , cencode=0)
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/neutropenia/fig11_1.png", 
    width=2000, height=1000, pointsize=22)
plot(xx, color = 1:3, xlab = "Time since Start of Treatment (Days)", 
     ylab = "Cumulative Incidence",
     curvlab = c("Neutropenia", "Progression/Death", "TTAE"))
dev.off()

xx <- cuminc(time, status_2, arm, cencode=0) 
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/neutropenia/fig11_2.png", 
    width=1000, height=1000, pointsize=22)
plot(xx, color = c(1, 1, 2, 2, 3, 3), 
     xlab = "Time since Start of Treatment (Days)", 
     ylab = "Cumulative Incidence",
     curvlab = c("Neutropenia", "Progression/Death", "TTAE"))
dev.off()

# Figure 13
snps13 <- c("rs10516812", "rs7372391", "rs13112901", "rs3748396")
for (s in snps13){
  xx <- cuminc(time_en, status_2, snps[,match(s, snp)], cencode=0)
  png(paste("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/neutropenia/fig13_", s, 
            "_en.png"), 
      width=1000, height=1000, pointsize=22)
  plot(xx, color = c(2, 2, 2, 3, 3, 3, 4, 4, 4),
       xlab = "Time since Start of Treatment (Days)", 
       ylab = "Cumulative Incidence")
  dev.off()
}

# Figure 14
snps14 <- c("rs6599244", "rs7431144", "rs2799083", "rs10509655")
for (s in snps14){
  xx <- cuminc(time_en, status_2, snps[,match(s, snp)], cencode=0)
  png(paste("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/neutropenia/fig14_", s, 
            "_en.png"), 
      width=1000, height=1000, pointsize=22)
  plot(xx, color = c(2, 2, 2, 3, 3, 3, 4, 4, 4),
       xlab = "Time since Start of Treatment (Days)", 
       ylab = "Cumulative Incidence")
  dev.off()
}

# Figure 15
snps15 <- c("rs10516812", "rs3111779", "rs9614417", "rs7372391")
for (s in snps15){
  xx <- cuminc(time, status_2, snps[,match(s, snp)], cencode=0)
  png(paste("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/neutropenia/fig15_", s, 
            ".png"), 
      width=1000, height=1000, pointsize=22)
  plot(xx, color = c(2, 2, 2, 3, 3, 3, 4, 4, 4),
       xlab = "Time since Start of Treatment (Days)", 
       ylab = "Cumulative Incidence")
  dev.off()
}

# Figure 16
snps16 <- c("rs7431144", "rs4847428", "rs2580817", "rs7707427")
for (s in snps16){
  xx <- cuminc(time, status_2, snps[,match(s, snp)], cencode=0)
  png(paste("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/neutropenia/fig16_", s, 
            ".png"), 
      width=1000, height=1000, pointsize=22)
  plot(xx, color = c(2, 2, 2, 3, 3, 3, 4, 4, 4),
       xlab = "Time since Start of Treatment (Days)", 
       ylab = "Cumulative Incidence")
  dev.off()
}


############## Marginal association for neutropina ################################

#Figure 12
for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()));
  fit <- summary(coxph(Surv(time, status_0) ~ snps[,i]))
  p_mar[i] <- fit$sctest[3]
}

score <- data.frame(CHR=chr, SNP=snp, BP=bp, P=p_mar)

index <- order(p_mar)[1:60]
print("scrore test top 60 hits with marginal association for neutropina")
score[index,]
snpsOfInterest1 <- snp[index]
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/neutropenia/mar.png", width=2000,
    height=1000, pointsize=22)
manhattan(score, main = "Manhattan Plot - cox model with score test for 
          neutropenia", ylim = c(0, 8), highlight = snpsOfInterest)
dev.off()

png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/neutropenia/qqplot.png", 
    width=1000, height=1000, pointsize=22)
qq(score$P, xlim=c(0,8), ylim=c(0,8))
dev.off()

############## Marginal association for early neutropina ##########################

for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()));
  fit <- summary(coxph(Surv(time_en, status_0) ~ snps[,i]))
  p_mar_en[i] <- fit$sctest[3]
}

score_en <- data.frame(CHR=chr, SNP=snp, BP=bp, P=p_mar_en)

index_en <- order(p_mar_en)[1:60]
print("scrore test top 60 hits with marginal association for early neutropina")
score_en[index_en,]
snpsOfInterest <- snp[index_en]
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/neutropenia/mar_en.png", 
    width=2000, height=1000, pointsize=22)
manhattan(score_en, main = "Manhattan Plot - cox model with score test for 
          early neutropenia", ylim = c(0, 8), highlight = snpsOfInterest)
dev.off()

png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/neutropenia/qqplot_en.png", 
    width=1000, height=1000, pointsize=22)
qq(score_en$P, xlim=c(0,8), ylim=c(0,8))
dev.off()

############## Fine and Gray model for neutropina #################################

for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()));
  #sink('/dev/null')
  sink(tempfile())
  fit1 <- summary(crr(time, status_1, snps[,i], failcode=1, cencode=0))
  fit2 <- summary(crr(time, status_2, snps[,i], failcode=1, cencode=0))
  pvalue1[i] <- fit1$coef[5]
  converge1[i] <- fit1$converged
  pvalue2[i] <- fit2$coef[5]
  converge2[i] <- fit2$converged
  sink() 
}

cnv1 <- (converge1 == 1)
chr1 <- chr[cnv1]
snp1 <- snp[cnv1]
bp1 <- bp[cnv1]
pvalue1 <- pvalue1[cnv1]
test1 <- data.frame(CHR=chr1, SNP=snp1, BP=bp1, P=pvalue1)
index1 <- order(pvalue1)[1:60]
print(c("test top 60 for neutropenia with 1st coding scheme after removing ",
        length(converge1) - length(cnv1), " non-converging cases"))
test1[index1,]
snpsOfInterest1 <- snp1[index1]

cnv2 <- (converge2 == 1)
chr2 <- chr[cnv2]
snp2 <- snp[cnv2]
bp2 <- bp[cnv2]
pvalue2 <- pvalue2[cnv2]
test2 <- data.frame(CHR=chr2, SNP=snp2, BP=bp2, P=pvalue2)
index2 <- order(pvalue2)[1:60]
print(c("test top 60 for neutropenia with 2nd coding scheme after removing ",
        length(converge2) - length(cnv1), " non-converging cases"))
test2[index2,]
snpsOfInterest2 <- snp2[index2]

############## Fine and Gray model for early neutropina ############################

for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()));
  #sink('/dev/null')
  sink(tempfile())
  fit1_en <- summary(crr(time_en, status_1_en, snps[,i], failcode=1, cencode=0))
  fit2_en <- summary(crr(time_en, status_2_en, snps[,i], failcode=1, cencode=0))
  pvalue1_en[i] <- fit1_en$coef[5]
  converge1_en[i] <- fit1_en$converged
  pvalue2_en[i] <- fit2_en$coef[5]
  converge2_en[i] <- fit2_en$converged
  sink() 
}

cnv1_en <- (converge1_en == 1)
chr1_en <- chr[cnv1_en]
snp1_en <- snp[cnv1_en]
bp1_en <- bp[cnv1_en]
pvalue1_en <- pvalue1_en[cnv1_en]
test1_en <- data.frame(CHR=chr1_en, SNP=snp1_en, BP=bp1_en, P=pvalue1_en)
index1_en <- order(pvalue1_en)[1:60]
print(c("test top 60 for early neutropenia with 1st coding scheme after removing ",
        length(converge1_en) - length(cnv1_en), " non-converging cases"))
test1_en[index1_en,]
snpsOfInterest1 <- snp1_en[index1_en]

cnv2_en <- (converge2_en == 1)
chr2_en <- chr[cnv2_en]
snp2_en <- snp[cnv2_en]
bp2_en <- bp[cnv2_en]
pvalue2_en <- pvalue2_en[cnv2_en]
test2_en <- data.frame(CHR=chr2_en, SNP=snp2_en, BP=bp2_en, P=pvalue2_en)
index2_en <- order(pvalue2_en)[1:60]
print(c("test top 60 for early neutropenia with 2nd coding scheme after removing ",
        length(converge2_en) - length(cnv1_en), " non-converging cases"))
test2_en[index2_en,]
snpsOfInterest2 <- snp2[index2_en]

save(test1, test2, test1_en, test2_en, score, score_en, df2, 
     file = "/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/neutropenia/data.Rdata")

