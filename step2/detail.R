library(qqman)

wd <- "/lustre/scr/j/i/jinjin/GRA/meta_analysis/result_main/"
top100 <- read.table(paste(wd, "top100.txt", sep = ""), header = TRUE)
f80303 <- read.table(paste(wd, "80303.txt", sep = ""), header = TRUE)
f80405 <- read.table(paste(wd, "80405.txt", sep = ""), header = TRUE)
f40503 <- read.table(paste(wd, "40503.txt", sep = ""), header = TRUE)
full <- read.table(paste(wd, "result.txt", sep = ""), header = TRUE)

snp <- top100$MarkerName
eff.80303 <- rep(".", 100)
eff.80405 <- rep(".", 100)
eff.40503 <- rep(".", 100)
sd.80303 <- rep(".", 100)
sd.80405 <- rep(".", 100)
sd.40503 <- rep(".", 100)
maf.80303 <- rep(".", 100)
maf.80405 <- rep(".", 100)
maf.40503 <- rep(".", 100)
maineff.80303 <- rep(0, 100)
maineff.80405 <- rep(0, 100)
maineff.40503 <- rep(0, 100)
mainsd.80303 <- rep(0, 100)
mainsd.80405 <- rep(0, 100)
mainsd.40503 <- rep(0, 100)
main <- rep(0, 100)
mainsd <- rep(0, 100)

for (i in 1:100){
  tempbeta <- 0
  tempsd <- 0
  if(is.element(snp[i], f80303$snp)){
    index <- match(snp[i], f80303$snp)
    if(top100$Allele1[i] == tolower(f80303$effallele[index])){
      eff.80303[i] = -f80303$coef[index] 
    } else{
      eff.80303[i] = f80303$coef[index]
    }
    sd.80303[i] = f80303$sd[index]
    maf.80303[i] = f80303$maf[index]
  }
  if(is.element(snp[i], f80405$snp)){
    index <- match(snp[i], f80405$snp)
    if(top100$Allele1[i] == tolower(f80405$effallele[index])){
      eff.80405[i] = -f80405$coef[index] 
    } else{
      eff.80405[i] = f80405$coef[index]
    }
    sd.80405[i] = f80405$sd[index]
    maf.80405[i] = f80405$maf[index]
  }
  if(is.element(snp[i], f40503$snp)){
    index <- match(snp[i], f40503$snp)
    if(top100$Allele1[i] == tolower(f40503$effallele[index])){
      eff.40503[i] = -f40503$coef[index] 
    } else{
      eff.40503[i] = f40503$coef[index]
    }
    sd.40503[i] = f40503$sd[index]
    maf.40503[i] = f40503$maf[index]
  }
}

detail <- data.frame(top100, eff.80303, sd.80303, maf.80303, 
eff.40503, 
                     sd.40503, maf.40503, eff.80405, sd.80405, maf.80405)

write.csv(detail, file = paste(wd, "top100_detail.csv"))

#### QQ plot ############################################################

png(paste(wd, "qq80303.png", sep = ""), width=1000, height=1000, pointsize=22)
qq(f80303$pvalue, xlim=c(0,8), ylim=c(0,8))
dev.off()

png(paste(wd, "qq80405.png", sep = ""), width=1000, height=1000, pointsize=22)
qq(f80405$pvalue, xlim=c(0,8), ylim=c(0,8))
dev.off()

png(paste(wd, "qq40503.png", sep = ""), width=1000, height=1000, pointsize=22)
qq(f40503$pvalue, xlim=c(0,8), ylim=c(0,8))
dev.off()

png(paste(wd, "qq_meta.png", sep = ""), width=1000, height=1000, pointsize=22)
qq(full$P.value, xlim=c(0,8), ylim=c(0,8))
dev.off()

