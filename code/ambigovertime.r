numAmbigperbp <- function(x){

    tmp <- listmuts(strsplit(x, split = "")[[1]])
    c(sum(tmp[c(-1:-4, -11)])/sum(tmp))
    
}

drug3 <- (a[a$numMeds == 3, ])
drugs <- names(summary(drug3$newRegimen))

drugreg <- drugs[6]
ther <- subset(a, newRegimen == drugreg)
years <- rownames(table(ther$IsolateYear))
avambig <- matrix(data = NA, nrow = length(years), ncol = 1000)

for(i in 1:length(years)){

    subther <- subset(ther, ther$IsolateYear == years[i])
    tmp <- c()

    for(j in 1:nrow(subther)){
       
       tmp[j] <- numAmbigperbp(as.vector(subther[j,]$RTSeq)[[1]])
       avambig[i,j] <- tmp[j]
   }

}

k <- 6
subset(avambig[k,], !is.na(avambig[k,]))

layout(matrix(1:4, nrow = 2))
par(mar = c(4,4,2,1))
plot(0, type = "n", xlim = range(as.numeric(years)), ylim = c(0, .15), xlab = "Year", ylab = "Ambig/bp",  main = paste(drugreg))
for(i in 1:length(years)){
    tmp <- subset(avambig[i,], !is.na(avambig[i,]))
    points(rep(as.numeric(years[i]), length(tmp)) + runif( length(tmp), -.1, .1), tmp , pch = 16, col = rgb(0,0,0,.1))
}
plot(years, apply(avambig, 1, function(x){ c(sum(x >= 0, na.rm = T)) }), xlab = "Year", pch = 16, ylab = "Number of Sequences",  main = paste(drugreg))

plot(years, apply(avambig, 1, function(x){ c(mean(x, na.rm = T))}), xlab = "Year", pch = 16, ylab = "Average ambig/bp",  main = paste(drugreg))
plot(years, apply(avambig, 1, function(x){ c(sum(x > 0, na.rm = T)/sum(x >= 0, na.rm = T)) }), xlab = "Year", pch = 16, ylab = "Prop. seqs with >0 ambig/bp",  main = paste(drugreg))



