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

#translate(conv(a$RTSeq[2]))

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

#nchar(alignRT(x))

#AZT alignment


#Ok, let's see if we can find any mutations specific to subclasses

nucOpts <- c("A","C","T","G","B","D","H","K","M","N","R","S","V","W","Y","a","c","g","k","m","r","s","t","w","y","~") 
mutCount <- matrix(data = NA, nrow = max(a$RTLastAA)*3, ncol = length(nucOpts))
dim(mutCount)

listmuts <- function(x){

    toReturn <- c()
    for(i in nucOpts){
        toReturn[length(toReturn) +1 ] <- sum(x == i)
    }
    return(toReturn)
}

#OK, I am not sure this is yielding anything. Let's try looking at those specifically marked mutations

b <- read.table("~/Dropbox/HIV_softsweeps/DRM_file.txt", header = T)
unique(b[,7])
colnames(b)

#AZT = Zidovudine
#DDI = Didanosine

#number of ambiguous reads per bp
numAmbigperbp <- function(x){

    tmp <- listmuts(strsplit(x, split = "")[[1]])
    c(sum(tmp[c(-1:-4, -11)])/sum(tmp))
    
}

