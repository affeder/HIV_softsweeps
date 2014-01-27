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

# ACT TGC ATT TAG
# (1) (4)  (4) -> 
#         ATT TAG AGG GAT
# (3) (6)
#
AZT[400,]$RTFirstAA
AZT[100,]$RTFirstAA

paste(AZT[100,]$RTSeq)

alignRT<- function(x){

    first <- (x$RTFirstAA - 1) *3 + 1
    last <- (x$RTLastAA) * 3 

    aligned <- paste(x$RTSeq)

    longest <- max(a$RTLastAA)*3


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

nchar(alignRT(x))

#AZT alignment
align <- c()
for(i in 1:nrow(AZT)){
    align[i] <- alignRT(AZT[i,])
}
dim(AZT)
length(align)
AZT <- cbind(AZT, align)

#DDI alignment
align <- c()
for(i in 1:nrow(DDI)){
    align[i] <- alignRT(DDI[i,])
}
dim(DDI)
length(align)
DDI <- cbind(DDI, align)

#AZTDDI alignment
align <- c()
for(i in 1:nrow(AZTDDI)){
    align[i] <- alignRT(AZTDDI[i,])
}
dim(AZTDDI)
length(align)
AZTDDI <- cbind(AZTDDI, align)


#Now, let's put all of this info into an array that we can query for polymorphisms


#First, let's just check that our alignments look ok

tmp <- matrix(data = NA, nrow = 2, ncol = 1680)
for(i in 1:nrow(AZT)){
    i <- 400
    tmp[1,] <- c((strsplit(as.vector(AZT$align)[1], split = "")[[1]]))
    tmp[2,] <- (strsplit(as.vector(AZT$align)[400], split = "")[[1]])
}

t(tmp[1:2, 1:100])


#Let's make one array that has all alignments
allaligns <- matrix(data = NA, nrow = nrow(AZT) + nrow(DDI) + nrow(AZTDDI), ncol =  max(a$RTLastAA)*3)
dim(allaligns)

idents <- c()

for(i in 1:nrow(AZT)){
    allaligns[i,] <- c((strsplit(as.vector(AZT$align)[i], split = "")[[1]]))
    idents <- c(idents, "AZT")
}

startind <- which(is.na(allaligns[,1]))[1] -1
for(i in 1:nrow(DDI)){
    allaligns[startind + i,] <- c((strsplit(as.vector(DDI$align)[i], split = "")[[1]]))
    idents <- c(idents, "DDI")
}

startind <- which(is.na(allaligns[,1]))[1] -1
for(i in 1:nrow(AZTDDI)){
    allaligns[startind + i,] <- c((strsplit(as.vector(AZTDDI$align)[i], split = "")[[1]]))
    idents <- c(idents, "AZTDDI")
}

#what does the 100th position look like?
allaligns[,100]

#Ok, let's find our polymorphic reads
#Let's make a list of all positions that have at least two reads (not .s)

x <- allaligns[,100]
polycheck <- function(x){

    if(length(unique(x[x != "."]) > 1)){ return(1) };
    return(false);
}

polypos <- apply(allaligns, 2, polycheck)
table(polypos)                                     

#EVERY position is polymorphic?

#Ok, so, let's try looking at positions in the treatment divided populations

AZTseq <- allaligns[idents == "AZT", ]
DDIseq <- allaligns[idents == "DDI",]
AZTDDIseq <- allaligns[idents == "AZTDDI",]
#First question, more polymorphism in single drug regimes?

#let's start by just counting the ambiguous reads by sequence
x <- allaligns[100,]
numAmbig <- function(x){

    nondot <- (x[x != "."])
#    size <- length(nondot)
    return(length(grep('[ATCG]', nondot, perl = T, invert = T)))#/size)
    
}

hist(apply(AZTseq, 1, numAmbig))
hist(apply(DDIseq, 1, numAmbig))
hist(apply(AZTDDIseq, 1, numAmbig))
