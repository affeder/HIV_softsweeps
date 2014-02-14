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
    return(length(grep('[atcgATCGN]', nondot, perl = T, invert = T)))#/size)
    
}
layout(matrix(1:3, ncol = 3))
hist(apply(AZTseq, 1, numAmbig))
hist(apply(DDIseq, 1, numAmbig))
hist(apply(AZTDDIseq, 1, numAmbig))

#There are some really large numbers of ambiguous reads in AZT
#What are these? Ns and atcgs

AZTnumAmbig <- apply(AZTseq, 1, numAmbig)
AZTnumAmbig
AZTseq[201,]

ks.test(apply(DDIseq,1,numAmbig), apply(AZTDDIseq,1,numAmbig))

AZTseq[110:120,100:110]



AZTdist <- apply(AZTseq, 2, listmuts)
DDIdist <- apply(DDIseq, 2, listmuts)
AZTDDIdist <- apply(AZTDDIseq, 2, listmuts)

#par(ask = T)

for(i in 100:200){

    tab <- cbind(AZTdist[,i]/sum(AZTdist[,i]),DDIdist[,i]/sum(DDIdist[,i]),AZTDDIdist[,i]/sum(AZTDDIdist[,i]))
    barplot(tab, col = rainbow(18))

}



subset(b, b$Drug == "Didanosine")
subset(b, b$Drug == "Zidovudine")


AAAZTseq <- apply(AZTseq, 1, translate)
AAAZTseq[2,]

#how about 210? 210 -> AA 70
AAAZTseq[70,]
#Ks and Rs. Is that what we'd expect?
#Yup, it's a K -> R transformation
#Do we see anything else?

unique(AAAZTseq[70,])
#Not at the amino acid level
#How about the nucleotide level?

trips <- apply((cbind(AZTseq[,210], AZTseq[,211], AZTseq[,212])), 1, paste, collapse = "")
unique(trips)
#[1] "ATG" "RTG" "GTG" "rTG" "ATg" "NNN"
#    M      M/V    V    M/V    M     -                                        

#Let's try this 



#AZ#Analyze out
#For the moment, let's make this AZT specific, and then we can generalize
analyzeAZT<- function(x){

    pos1 <- AZTseq[,x]
    pos2 <- AZTseq[,x+1]
    pos3 <- AZTseq[,x+2]

    
    x <- 67
    pos1 <- AZTseq[,3*x]
    pos2 <- AZTseq[,3*x+1]
    pos3 <- AZTseq[,3*x+2]

    
    i <- 382
    print(c(pos1[i], pos2[i], pos3[i]))
    print(translate(c(pos1[i], pos2[i], pos3[i])))
    
    print(unique(pos1))
    print(unique(pos2))
    print(unique(pos3))

    translate(rbind(pos1, pos2, pos3))
    for(i in 1:nrow(AZTseq)){
        print(c(pos1[i], pos2[i], pos3[i]))
        print(translate(c(pos1[i], pos2[i], pos3[i])))
    }
}




#diversity in monotherapies

drug1 <- (a[a$numMeds == 1, ])
drug2 <- (a[a$numMeds == 2, ])
drug3 <- (a[a$numMeds == 3, ])

numAmbig(drug1$RTSeq[3])
table(strsplit(as.vector(drug1$RTSeq[3]), split = "")[[1]])
numAmbig(strsplit(as.vector(drug1$RTSeq[3]), split = "")[[1]])

d1align <- c()
for(i in 1:nrow(drug1)){
    d1align[i] <- alignRT(drug1[i,])
}
d2align <- c()
for(i in 1:nrow(drug2)){
    d2align[i] <- alignRT(drug2[i,])
}
d3align <- c()
for(i in 1:nrow(drug3)){
    d3align[i] <- alignRT(drug3[i,])
}




numAmbigperbp(d1align[5])

d1amb <- c()
for(i in 1:nrow(drug1)){
    d1amb[i] <- numAmbigperbp(d1align[i])
}
d2amb <- c()
for(i in 1:nrow(drug2)){
    d2amb[i] <- numAmbigperbp(d2align[i])
}
d3amb <- c()
for(i in 1:nrow(drug3)){
    d3amb[i] <- numAmbigperbp(d3align[i])
}

d1ambtmp <- d1amb
d2ambtmp <- d2amb
d3ambtmp <- d3amb

d1amb[d1amb == 0] <- 10^-5
d2amb[d2amb == 0] <- 10^-5
d3amb[d3amb == 0] <- 10^-5
par(ask = F)
par(mar = c(4,4, 3.5,0))
layout(matrix(1:3, ncol = 3))
breaks <- seq(-13, -.5, by = .5)
hist(log(d1amb), main = "One drug", xlab = "log(Ambiguous/bp)", breaks = breaks)
hist(log(d2amb), main = "Two drugs", xlab = "log(Ambiguous/bp)", breaks = breaks)
hist(log(d3amb), main = "Three drugs", xlab = "log(Ambiguous/bp)", breaks = breaks)

plot(density(log(d3ambtmp)), col = "red", xlab = "log(ambiguous/bp)", main = "")
lines(density(log(d1ambtmp)), col = "blue")
lines(density(log(d2ambtmp)), col = "green")
legend("topright", c("1 drug", "2 drugs", "3 drugs"), col = c("blue", "green", "red"), lty = ("solid"))


#drug1$Weeks

length(d2amb)
length(d3amb)





#Let's look at the years with treatment now

colnames(a)
a$IsolateYear

