#seqinr

a <- read.table("dataset.01.23.txt", skip = 13, sep= "\t", header = T)

length(unique(a$PtID))
colnames(a)

summary(a$RegimenName)
summary(a$Subtype)
#housekeeping on regimen names

#more mutations in monotherapy/old
#Are we seeing DRMs in PR when not given PIs
