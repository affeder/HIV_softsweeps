library(seqinr)

list.files()
a <- read.table("~/Dropbox/HIV_softsweeps/dataset.01.24.txt", skip = 13, sep= "\t", header = T)

length(unique(a$PtID))
colnames(a)

summary(a$RegimenName)
summary(a$Subtype)
#housekeeping on regimen names

#more mutations in monotherapy/old
#Are we seeing DRMs in PR when not given PIs

conv <- function(x){
    return(strsplit(as.vector(x), split = "")[[1]])
   }

translate(conv(a$RTSeq[2]))

colnames(a)
table(a$RegimenName)

standardizeRegimenName <- function(x){
    sort(strsplit(as.vector(x), split = '\\+')[[1]])
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

tmp <- strsplit(as.vector(a$RegimenName), split = '\\+')

newRegimen <- c()
numMeds <- c()
for(i in 1:length(tmp)){
    newRegimen[i] <- paste(sort(trim(tmp[[i]])), collapse = "+")
    numMeds[i] <- length(tmp[[i]])
}

table(newRegimen)
a <- cbind(a, newRegimen, numMeds)


colnames(a)
sort(unique(a$newRegimen))

#people with 8 meds?
a[a$numMeds == 8, ]

#what therapies are monotherapies?
monodrugs_factor <- sort(unique(subset(a, a$numMeds == 1)$newRegimen))
monodrugs <- c()
for(i in 1:length(monodrugs_factor)){
    monodrugs[i] <- paste(monodrugs_factor[i])
}
monodrugs

#How many of each monotherapy do we have?
monodrugcount <- c()
for(i in 1:length(monodrugs)){
    monodrugcount[i] <- sum(a$newRegimen == monodrugs[i])
}
monos<- cbind(monodrugcount, monodrugs)
monos_high <- subset(monos, as.numeric(monos[,1]) > 100)
monos

#let's list every combination of those two
doubdrugs <- c()
for(i in 1:(length(monodrugs)-1)){
    for(j in (i+1):length(monodrugs)){
        doubdrugs[length(doubdrugs) + 1] <- paste(monodrugs[i], monodrugs[j], sep = "+")
    }
}
doubdrugs

#How many of the doubles do we actually have?
doubdrugcount <- c()
for(i in 1:length(doubdrugs)){
    doubdrugcount[i] <- sum(a$newRegimen == doubdrugs[i])
}

doubles <- cbind(doubdrugcount, doubdrugs)
doubs <- subset(doubles, as.numeric(doubles[,1]) > 100)

doubs
monos

#The most promising looks like AZT+DDI
AZT <- subset(a, a$newRegimen == "AZT")
DDI <- subset(a, a$newRegimen == "DDI")
AZTDDI <- subset(a, a$newRegimen == "AZT+DDI")

dim(AZT)
dim(DDI)
dim(AZTDDI)

colnames(a)
AZT$RTFirstAA

#we'll need to do the alignment
x <- AZT[1,]

alignRT<- function(x){
    
    first <- x$RTFirstAA
    last <- x$RTLastAA
    aligned <- x$RTSeq

    longest <- max(a$RTLastAA)
    
    while(first > 1){
        aligned <- paste(".", aligned, sep = "")
        first <- first - 1
    }

    while(last < longest){
        aligned <- paste(aligned, ".", sep = "")
        last <- last + 1
    }

    return(aligned)
}


align <- c()
for(i in 1:nrow(AZT)){
    align[i] <- alignRT(AZT[i,])
}

