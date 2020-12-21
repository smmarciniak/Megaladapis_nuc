#!/usr/bin/env Rscript

#Generate individual histograms of read frequency based on read length (bp) using quality filtered (min 20bp) non-redundant trimmed/merged reads

reads11N <- read.table("bar11_not_unique_min20.txt", sep="", header=FALSE)
reads12N <- read.table("bar12_not_unique_min20.txt", sep="", header=FALSE)
reads13Nd <- read.table("bar13_dec_not_unique_min20.txt", sep="", header=FALSE)
reads13Nj <- read.table("bar13_jun2_not_unique_min20.txt", sep="", header=FALSE)
reads14Nd <- read.table("bar14_dec_not_unique_min20.txt", sep="", header=FALSE)
reads14Nj <- read.table("bar14_jun2_not_unique_min20.txt", sep="", header=FALSE)
reads41Nj <- read.table("bar41_jul2_not_unique_min20.txt", sep="", header=FALSE)
reads41Nn <- read.table("bar41_nov2_not_unique_min20.txt", sep="", header=FALSE)
reads42Nj <- read.table("bar42_jul2_not_unique_min20.txt", sep="", header=FALSE)
reads42Nn <- read.table("bar42_nov2_not_unique_min20.txt", sep="", header=FALSE)
readspjm <- read.table("PJX039_may2_not_unique_min20.txt", sep="", header=FALSE)
readspjj <- read.table("PJX039_jun2_not_unique_min20.txt", sep="", header=FALSE)
readssN <- read.table("schuster_CORR_not_unique_min20.txt", sep="", header=FALSE)

hist(reads11N$V1, col="#E58606", main="UA5180_11", yaxt="n", xlab="Read length (bp)", ylim=c(0,1000000)) 
axis(2, axTicks(2), format(axTicks(2), scientific = F)) 

hist(reads12N$V1, col="#5D69B1", main="UA5180_12", yaxt="n", xlab="Read length (bp)", ylim=c(0,1000000)) 
axis(2, axTicks(2), format(axTicks(2), scientific = F))

hist(reads13Nd$V1, col="#99C945", main="UA5180_13B", yaxt="n", xlab="Read length (bp)", ylim=c(0,10000000)) 
axis(2, axTicks(2), format(axTicks(2), scientific = F)) 

hist(reads13Nj$V1, col="#52BCA3", main="UA5180_13A", yaxt="n", xlab="Read length (bp)", ylim=c(0,10000000)) 
axis(2, axTicks(2), format(axTicks(2), scientific = F)) 

hist(reads14Nd$V1, col="#24796C", main="UA5180_14B", yaxt="n", xlab="Read length (bp)", ylim=c(0,10000000)) 
axis(2, axTicks(2), format(axTicks(2), scientific = F)) 

hist(reads14Nj$V1, col="#CC61B0", main="UA5180_14A", yaxt="n", xlab="Read length (bp)", ylim=c(0,10000000)) 
axis(2, axTicks(2), format(axTicks(2), scientific = F))

hist(reads41Nj$V1, col="#DAA51B", main="UA5180_41A", yaxt="n", xlab="Read length (bp)") 
axis(2, axTicks(2), format(axTicks(2), scientific = F)) 

hist(reads41Nn$V1, col="#2F8AC4", main="UA5180_41B", yaxt="n", xlab="Read length (bp)") 
axis(2, axTicks(2), format(axTicks(2), scientific = F)) 

hist(reads42Nj$V1, col="#764E9F", main="UA5180_42A", yaxt="n", xlab="Read length (bp)") 
axis(2, axTicks(2), format(axTicks(2), scientific = F)) 

hist(reads42Nn$V1, col="#ED645A", main="UA5180_42B", yaxt="n", xlab="Read length (bp)") 
axis(2, axTicks(2), format(axTicks(2), scientific = F)) 

hist(readspjm$V1, col="#CC3A8E", main="UA5180_39A", yaxt="n", xlab="Read length (bp)") 
axis(2, axTicks(2), format(axTicks(2), scientific = F)) 

hist(readspjj$V1, col="#94346E", main="UA5180_39B", yaxt="n", xlab="Read length (bp)", ylim=c(0,3000000)) 
axis(2, axTicks(2), format(axTicks(2), scientific = F)) 

hist(readssN$V1, col="#332288", main="UA5180_3", yaxt="n", xlab="Read length (bp)", xlim=c(20,140)) 
axis(2, axTicks(2), format(axTicks(2), scientific = F))

