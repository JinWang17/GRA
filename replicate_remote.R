rm(list=ls())


library(GenABEL)
library(survival)
library(SNPassoc)
library(xtable)
library(qqman)
library(gap)


######### For subjects and snps alreday filtered and QCed ###############################

tfamfile="/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/plink.tfam"
tools::md5sum(tfamfile)

tpedfile="/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/step5.tped"
tools::md5sum(tpedfile)

phenofile="/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/phenowithPC_new.txt"
tools::md5sum(phenofile)

tfamtopfile="/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/tophits.tfam"
tools::md5sum(tfamfile)

tpedtopfile="/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/tophits.tped"
tools::md5sum(tpedfile)

expidfile="/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/expidlist_294.txt"
tools::md5sum(expidfile)

convert.snp.tped(tpedfile = tpedfile, tfamfile = tfamfile, outfile = "genotype.raw")
df <- load.gwaa.data(phe = phenofile, gen = "genotype.raw", force = T)

convert.snp.tped(tpedfile = tpedtopfile, tfamfile = tfamtopfile, outfile = "tophits.raw")
tophits <- load.gwaa.data(phe = phenofile, gen = "tophits.raw", force = T)


load("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/Pat803722Data.RData", verbose=TRUE)
pheno_top <- read.csv("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/CALGB80303_338_pheno.csv")
pheno_top <- pheno_top[pheno_top$id == 803722, ]

class(df)
names(phdata(df))
attach(phdata(df))
descriptives.trait(df)
descriptives.marker(df)

nsnps <- nsnps(gtdata(df))

############## MGC > 9 #################################################################

mgcfun<-function(x){ 
  x0=x[1]
  x1=x[2]
  x2=x[3]
  if(x0*x1*x2>0)
    MGC=min(x0,x1,x2)
  else if(x0==0)
    MGC=min(x1,x2)
  else if(x1==0)
    MGC=min(x0,x2)
  else if(x2==0)
    MGC=min(x1,x0)
  MGC 
}

temp <- summary(df[, 1:nsnps])
MGC <- apply(temp[,c("P.11","P.12","P.22")],1,mgcfun)
sum(MGC>9)

df2 <-df[,(MGC>9)]


################# check result ##############################################################

if(1 < 0){

select <-c(365014, 365005, 160998, 216329, 41180, 203715, 308532, 40897, 272293,
           216046, 170654, 190554, 303203, 416787, 32893, 216043, 292594, 116434,
           32890, 332935)
dta <- df[, select]
snpss <- as.numeric(gtdata(dta))
times <- phdata(dta)$ostime
statuss <- phdata(dta)$osevent
p <- rep(0, nsnps(gtdata(dta)))
chr <- c(as.numeric(chromosome(gtdata(dta))))
snp <- c(colnames(snpss))
bp <- c(as.numeric(map(gtdata(dta))))
for(i in 1:nsnps(gtdata(dta))){
  fit <- survdiff(Surv(times, statuss) ~ snpss[,i])
  print(fit)
  p[i] <- pchisq(fit$chisq,length(fit$n)-1, lower.tail = FALSE)
}
snp
p



select2 <-c(182427, 182428, 334345, 429939, 86710, 332254, 11374, 231507,
            416784, 269305, 328676, 483284, 262276, 408773, 402912, 157560,
            431274, 63831, 88585, 271608)
dta2 <- df[, select2]
snpss <- as.numeric(gtdata(dta2))
times <- phdata(dta2)$ostime
statuss <- phdata(dta2)$osevent
p2 <- rep(0, nsnps(gtdata(dta2)))
chr <- c(as.numeric(chromosome(gtdata(dta2))))
snp <- c(colnames(snpss))
bp <- c(as.numeric(map(gtdata(dta2))))
for(i in 1:nsnps(gtdata(dta2))){
  fit <- survdiff(Surv(times, statuss) ~ snpss[,i])
  print(fit)
  p2[i] <- pchisq(fit$chisq,length(fit$n)-1, lower.tail = FALSE)
}
snp
p2

for(i in 1:nsnps(gtdata(dta2))){
  fit <- survdiff(Surv(times, statuss) ~ snpss[,i])
  print(fit)
  p2[i] <- 1 - pchisq(fit$chisq,length(fit$n)-1)
}
snp
p2
}

############## Before MGC Filtering ##############################################################

if(1 < 0){
nsnps <- nsnps(gtdata(df))
nids <- nids(gtdata(df))


# primary analysis - log rank test 
snps <- as.numeric(gtdata(df))
time <- phdata(df)$ostime
status <- phdata(df)$osevent

p <- rep(0, nsnps)
chr <- c(as.numeric(chromosome(gtdata(df))))
snp <- c(colnames(snps))
bp <- c(as.numeric(map(gtdata(df))))
for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()));
  fit <- survdiff(Surv(time, status) ~ snps[,i])
  p[i] <- pchisq(fit$chisq,length(fit$n)-1, lower.tail = FALSE)
}

logrank <- data.frame(CHR=chr, SNP=snp, BP=bp, P=p)
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/logrank.png", width=2000, height=1000, pointsize=18)
manhattan(logrank, main = "Manhattan Plot - Primary Analysis", ylim = c(0, 14), cex = 0.6, 
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F)
dev.off()
index <- order(p)[1:20]
logrank[index,]

# validation - cox proportional hazard regression model
ph <- mlreg(GASurv(time,status)~1,df, trait.type="survival")
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/ph.png", width=2000, height=1000, pointsize=18)
plot(ph, ylim = c(0, 14), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"))
dev.off()
descriptives.scan(data = ph, top = 20, sortby = "P1df")

}

############ After MGC filtering #################################################################


nsnps <- nsnps(gtdata(df2))
nids <- nids(gtdata(df2))

# primary analysis - log rank test 
snps <- as.numeric(gtdata(df2))
time <- phdata(df2)$ostime
status <- phdata(df2)$osevent


p.wald <- rep(0, nsnps)
p.ll <- p.wald
p.score <- p.wald
chr <- c(as.numeric(chromosome(gtdata(df2))))
snp <- c(colnames(snps))
bp <- c(as.numeric(map(gtdata(df2))))
for(i in 1:nsnps){
  if (i%%10000==0) print(paste(i, date()));
  fit <- summary(coxph(Surv(time, status) ~ snps[,i]))
  #p.wald[i] <- fit$waldtest[3]
  #p.ll[i] <- fit$logtest[3]
  p.score[i] <- fit$sctest[3]
}
warnings()

#waldtest <- data.frame(CHR=chr, SNP=snp, BP=bp, P=p.wald)
#lltest <- data.frame(CHR=chr, SNP=snp, BP=bp, P=p.ll)
scoretest <- data.frame(CHR=chr, SNP=snp, BP=bp, P=p.score)

if(1<0){
index <- order(p.wald)[1:20]
print("Wald test top 20")
waldtest[index,]
snpsOfInterest <- snp[index]
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/replicate/primary_ph_wald.png", width=2000, height=1000, pointsize=18)
manhattan(waldtest, main = "Manhattan Plot - Primary Analysis - cox model with Wald test", ylim = c(0, 8), highlight = snpsOfInterest) 
dev.off()

index <- order(p.ll)[1:20]
print("Log likelihood test top 20")
lltest[index,]
snpsOfInterest <- snp[index]
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/replicate/primary_ph_ll.png", width=2000, height=1000, pointsize=18)
manhattan(lltest, main = "Manhattan Plot - Primary Analysis - cox model with log likelihood test", ylim = c(0, 8), highlight = snpsOfInterest)
dev.off()
}

# figure 2

index <- order(p.score)[1:20]
print("scrore test top 20")
scoretest[index,]
snpsOfInterest <- snp[index]
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/replicate/primary_ph_score.png", width=2000, height=1000, pointsize=22)
manhattan(scoretest, main = "Manhattan Plot - Primary Analysis - cox model with score test", ylim = c(0, 8), 
          highlight = snpsOfInterest)
dev.off()

# Chen
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/replicate/primary_ph_score_Chen.png", width=2000, height=1000, pointsize=22)
mhtplot(scoretest[, c("CHR","BP","P")], labels=1:22, base=10, pch=19, ylim=c(0,9), 
        control=mht.control(type="p", colors=rep(gray(2:12/15),2), usepos=TRUE, logscale=TRUE, cutoffs=c(4,5,6,7,8)))
dev.off()

# produce table 2 needs using this: http://www.scandb.org/newinterface/index_v1.html

scoretest$strand <- strand(df2)
scoretest$effa <- effallele(df2)
scoretest$refa <- refallele(df2)

index <- order(p.score)[1:20]
top20 <- scoretest[index, ]

write.table(top20, "/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/replicate/top20.txt", sep= "\t")
top20$SNP

# supplementary material table 2
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/replicate/qqplot.png", width=1000, height=1000, pointsize=22)
qq(scoretest$P, xlim=c(0,8), ylim=c(0,8))
dev.off()

# figure 3

gene1 <- "rs763780"
pos1 <- match(gene1, snp)
snp1 <- as.character(gtdata(df2[, gene1]))
fit1 <- survfit(Surv(time, status) ~ snp1)
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/replicate/km_rs763780.png", width=1000, height=1000, pointsize=18)
plot(fit1, lty=as.numeric(as.factor(fit1$strata)), xlab = "Time since Randomization(Months)", ylab = "Probability", cex = 1)
legend(30, 0.9, substring(names(fit1$strata), 6, nchar(names(fit1$strata))), lty=as.numeric(as.factor(fit1$strata)))
dev.off()

gene2 <- "rs11644322"
pos2 <- match(gene2, snp)
snp2 <- as.character(gtdata(df2[, gene2]))
fit2 <- survfit(Surv(time, status) ~ snp2)
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/replicate/km_rs11644322.png", width=1000, height=1000, pointsize=18)
plot(fit2, lty=as.numeric(as.factor(fit2$strata)), xlab = "Time since Randomization(Months)", ylab = "Probability", cex = 1)
legend(30, 0.9, substring(names(fit2$strata), 6, nchar(names(fit2$strata))), lty=as.numeric(as.factor(fit2$strata)))
dev.off()

gene3 <- "rs10883617"
pos3 <- match(gene3, snp)
snp3 <- as.character(gtdata(df2[, gene3]))
fit3 <- survfit(Surv(time, status) ~ snp3)
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/replicate/km_rs10883617.png", width=1000, height=1000, pointsize=18)
plot(fit3, lty=as.numeric(as.factor(fit3$strata)), xlab = "Time since Randomization(Months)", ylab = "Probability", cex = 1)
legend(30, 0.9, substring(names(fit3$strata), 6, nchar(names(fit3$strata))), lty=as.numeric(as.factor(fit3$strata)))
dev.off()


# OS result on page 5
ph <- phdata(df2)
fit_arm1 <- survfit(Surv(ph$ostime[ph$arm==1], ph$osevent[ph$arm==1]) ~ 1)
fit_arm2 <- survfit(Surv(ph$ostime[ph$arm==2], ph$osevent[ph$arm==2]) ~ 1)

gene4 <- "rs7771466"
pos4 <- match(gene4, snp)
snp4 <- as.character(gtdata(df2[, gene4]))
fit4 <- survfit(Surv(time, status) ~ snp4)

summary(fit_arm1)
summary(fit_arm2)
summary(fit1)
summary(fit4)

# supplementary material table 3
bb <- summary(df2[, snpsOfInterest])
####print(xtable(bb, caption = "snpsummary", label="tab:snpsumary"), size = "\\tiny")
print(bb)


# supplementary material table 1
eff1 <- rep(0, 20)
cn_u1 <- rep(0, 20)
cn_l1 <- rep(0, 20)
pval1 <- rep(0, 20)

eff2 <- eff1
eff3 <- eff1
eff4 <- eff1
eff5 <- eff1

cn_u2 <- cn_u1
cn_u3 <- cn_u1
cn_u4 <- cn_u1
cn_u5 <- cn_u1

cn_l2 <- cn_l1
cn_l3 <- cn_l1
cn_l4 <- cn_l1
cn_l5 <- cn_l1

pval2 <- pval1
pval3 <- pval1
pval4 <- pval1
pval5 <- pval1

snp_top <- snpnames(tophits)
pos <- match(top20$SNP, snp_top)
snps_top <- as.numeric(gtdata(tophits))

snptocheck <- snpnames(gtdata(kk))
summary(is.element(snpsOfInterest, snptocheck))
pos_803722 <- match(snpsOfInterest, snptocheck)
snps_803722 <- as.numeric(gtdata(kk))

snps_338 <- rbind(snps_top[, pos], snps_803722[, pos_803722])

pheno_338 <- data.frame(ostime = phdata(tophits)$ostime, osevent = phdata(tophits)$osevent, 
                        PC1 = phdata(tophits)$PC1_338, PC2 = phdata(tophits)$PC2_338, 
                        PC3 = phdata(tophits)$PC3_338)
pheno_338 <- rbind(pheno_338, c(pheno_top$ostime, pheno_top$osevent, pheno_top$PC1,
                                pheno_top$PC2, pheno_top$PC3))

j <- 0

for (i in index){
  j <- j + 1
  fit1 <- summary(coxph(Surv(time, status) ~ snps[,i]))
  eff1[j] <- fit1$coef[1, 2]
  pval1[j] <- fit1$sctest[3]
  cn_l1[j] <- fit1$conf.int[1, 3]
  cn_u1[j] <- fit1$conf.int[1, 4]
  if (eff1[j] < 1){
    eff1[j] <- 1/fit1$coef[1, 2]
    cn_l1[j] <- 1/fit1$conf.int[1, 4]
    cn_u1[j] <- 1/fit1$conf.int[1, 3]
  }
  
  #  fit2 <- summary(coxph(Surv(time, status) ~ snps[,i] + phdata(df2)$priorrt_new + phdata(df2)$dzext_new 
  #                        + phdata(df2)$ps_new))  
  fit2 <- summary(coxph(Surv(time, status) ~ snps[,i] + phdata(df2)$priorrt + phdata(df2)$dzext 
                        + phdata(df2)$ps))
  eff2[j] <- fit2$coef[1, 2]
  pval2[j] <- fit2$coef[1, 5]
  cn_l2[j] <- fit2$conf.int[1, 3]
  cn_u2[j] <- fit2$conf.int[1, 4]
  if (eff2[j] < 1){
    eff2[j] <- 1/fit2$coef[1, 2]
    cn_l2[j] <- 1/fit2$conf.int[1, 4]
    cn_u2[j] <- 1/fit2$conf.int[1, 3]
  }
  
  #  fit3 <- summary(coxph(Surv(time, status) ~ snps[,i] + phdata(df2)$priorrt_new + phdata(df2)$dzext_new 
  #                        + phdata(df2)$ps_new + phdata(df2)$arm + phdata(df2)$PC1 + phdata(df2)$PC2 + phdata(df2)$PC3))
  fit3 <- summary(coxph(Surv(time, status) ~ snps[,i] + phdata(df2)$priorrt 
                        + phdata(df2)$dzext + phdata(df2)$ps + phdata(df2)$arm 
                        + phdata(df2)$PC1 + phdata(df2)$PC2 + phdata(df2)$PC3))
  eff3[j] <- fit3$coef[1, 2]
  pval3[j] <- fit3$coef[1, 5]
  cn_l3[j] <- fit3$conf.int[1, 3]
  cn_u3[j] <- fit3$conf.int[1, 4]
  if (eff3[j] < 1){
    eff3[j] <- 1/fit3$coef[1, 2]
    cn_l3[j] <- 1/fit3$conf.int[1, 4]
    cn_u3[j] <- 1/fit3$conf.int[1, 3]
  }
  
  fit4 <- summary(coxph(Surv(pheno_338$ostime, pheno_338$osevent) ~ 
                          snps_338[,j] + pheno_338$PC1 + pheno_338$PC2 
                        + pheno_338$PC3))
  eff4[j] <- fit4$coef[1, 2]
  pval4[j] <- fit4$coef[1, 5]
  cn_l4[j] <- fit4$conf.int[1, 3]
  cn_u4[j] <- fit4$conf.int[1, 4]
  if (eff4[j] < 1){
    eff4[j] <- 1/fit4$coef[1, 2]
    cn_l4[j] <- 1/fit4$conf.int[1, 4]
    cn_u4[j] <- 1/fit4$conf.int[1, 3]
  }
  
  fit5 <- summary(coxph(Surv(time, status) ~ snps[,i] + phdata(df2)$PC1 
                        + phdata(df2)$PC2 + phdata(df2)$PC3))
  eff5[j] <- fit5$coef[1, 2]
  pval5[j] <- fit5$coef[1, 5]
  cn_l5[j] <- fit5$conf.int[1, 3]
  cn_u5[j] <- fit5$conf.int[1, 4]
  if (eff5[j] < 1){
    eff5[j] <- 1/fit5$coef[1, 2]
    cn_l5[j] <- 1/fit5$conf.int[1, 4]
    cn_u5[j] <- 1/fit5$conf.int[1, 3]
  } 
}


aa <- data.frame(snps=snpsOfInterest,
                 pC1=pval1, HR1=eff1, HR1_L=cn_l1, HR1_U=cn_u1, 
                 pC2=pval2, HR2=eff2, HR2_L=cn_l2, HR2_U=cn_u2, 
                 pC3=pval3, HR3=eff3, HR3_L=cn_l3, HR3_U=cn_u3, 
                 pC4=pval4, HR4=eff4, HR4_L=cn_l4, HR4_U=cn_u4, 
                 pC5=pval5, HR5=eff5, HR5_L=cn_l5, HR5_U=cn_u5)

##print(xtable(aa, caption = "supplementary table 1", label="tabtophits"), size = "\\tiny")
print(aa)


# check using priorrt_new, dzzext_new, ps_new
if(1<0)
{
j <- 0

for (i in index){
  fit2 <- summary(coxph(Surv(time, status) ~ snps[,i] + phdata(df2)$priorrt_new 
                        + phdata(df2)$dzext_new + phdata(df2)$ps_new))
  eff2[j] <- fit2$coef[1, 2]
  pval2[j] <- fit2$coef[1, 5]
  cn_l2[j] <- fit2$conf.int[1, 3]
  cn_u2[j] <- fit2$conf.int[1, 4]
  if (eff2[j] < 1){
    eff2[j] <- 1/fit2$coef[1, 2]
    cn_l2[j] <- 1/fit2$conf.int[1, 4]
    cn_u2[j] <- 1/fit2$conf.int[1, 3]
  }
  
  fit3 <- summary(coxph(Surv(time, status) ~ snps[,i] + phdata(df2)$priorrt_new 
                        + phdata(df2)$dzext_new + phdata(df2)$ps + phdata(df2)$arm + phdata(df2)$PC1 + phdata(df2)$PC2 + phdata(df2)$PC3))
  eff3[j] <- fit3$coef[1, 2]
  pval3[j] <- fit3$coef[1, 5]
  cn_l3[j] <- fit3$conf.int[1, 3]
  cn_u3[j] <- fit3$conf.int[1, 4]
  if (eff3[j] < 1){
    eff3[j] <- 1/fit3$coef[1, 2]
    cn_l3[j] <- 1/fit3$conf.int[1, 4]
    cn_u3[j] <- 1/fit3$conf.int[1, 3]
  }
  
}

aa_check <- data.frame(snps=snpsOfInterest, 
                 pC2=pval2, HR2=eff2, HR2_L=cn_l2, HR2_U=cn_u2, 
                 pC3=pval3, HR3=eff3, HR3_L=cn_l3, HR3_U=cn_u3)

##print(xtable(aa, caption = "supplementary table 1", label="tabtophits"), size = "\\tiny")
print(aa_check)
}

# table 1 

ph_338 <- data.frame(sex = phdata(tophits)$sex, age = phdata(tophits)$age, 
                     race = phdata(tophits)$race, dzext = phdata(tophits)$dzext, 
                     priorrt = phdata(tophits)$priorrt, ps = phdata(tophits)$ps, 
                     arm = phdata(tophits)$arm)
ph_338 <- rbind(ph_338, c(pheno_top$sex, pheno_top$age, pheno_top$race, pheno_top$dzext, 
                          pheno_top$priorrt, pheno_top$ps, pheno_top$arm))
ph_338$race[338] <- pheno_top$race
table(ph_338$sex)
table(ph_338$race)
table(ph_338$dzext)
table(ph_338$priorrt)
table(ph_338$ps)
table(ph_338$arm)
mean(ph_338$age)
sqrt(var(ph_338$age))
quantile(ph_338$age, c(0.05, 0.5, 0.95))

ph_294 <- data.frame(sex = phdata(df2)$sex, age = phdata(df2)$age, 
                  race = phdata(df2)$race, dzext = phdata(df2)$dzext, 
                  priorrt = phdata(df2)$priorrt, ps = phdata(df2)$ps, 
                  arm = phdata(df2)$arm)
table(ph_294$sex)
table(ph_294$race)
table(ph_294$dzext)
table(ph_294$priorrt)
table(ph_294$ps)
table(ph_294$arm)
mean(ph_294$age)
sqrt(var(ph_294$age))
quantile(ph_294$age, c(0.05, 0.5, 0.95))

########## check the ph result by genAbel ##########################################################  
if(1<0){
ph <- mlreg(GASurv(time,status)~1,df2, trait.type="survival")
png("/lustre/scr/j/i/jinjin/GRA/CALGB80303/data/replicate/ph_check.png", width=2000, height=1000, pointsize=18)
plot(ph, ylim = c(0, 8))
dev.off()
descriptives.scan(data = ph, top = 20, sortby = "P1df")

#cox1 <- mlreg(GASurv(time,status)~1,df, trait.type="survival")
#cox2 <- mlreg(GASurv(time,status)~1,df, trait.type="survival")
#cox3 <- mlreg(GASurv(time,status)~1,df, trait.type="survival")
#cox4 <- mlreg(GASurv(time,status)~1,df, trait.type="survival")
}


rm(list=ls())






