a <- read.table("METAANALYSIS1.TBL", header = TRUE)
str(a)
newa <- a[order(a$P.value), ] 
head(newa)
write.table(newa, "meta_result.txt", quote = FALSE)