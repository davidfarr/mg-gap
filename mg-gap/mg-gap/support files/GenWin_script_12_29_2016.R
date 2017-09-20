#Script to take B* results from JK's python script & window size of 1 to call window breakpoints
#A. Scoville 12_29_2016

#download and install packages the first time
install.packages("tidyr")
install.packages("GenWin")
install.packages("ggplot2")

#load packages
library(tidyr)  #for pre-processing; contains function separate
library(GenWin) #for spline-based window analysis
#library(ggplot2) #not necessary yet, but we may want to add graphing later 

#Set appropriate working directory and read in B1.txt file containing B* values based on a window size of 1
#Note: B1.txt is output from JK's python script
#setwd("N:\\app dev\\scoville research\\program files\\github repo\\mg-gap\\mg-gap\\mg-gap\\mg-gap\\bin\\Debug")
setwd("C:\\Users\\David\\Documents\\GitHub\\mg-gap\\mg-gap\\mg-gap\\bin\\Debug")

file <- read.table("B1_new.txt", header=FALSE)
  
#give the file appropriate column names
# colnames(file) <- c("Scaffoldraw","B","Bs","P") #this was original
colnames(file) <- c("Scaffoldraw","B")
  
#split the location into separate columns for SNP, chromosome, and base pair and give appropriate column headers
file <- separate(data=file, col = Scaffoldraw, into = c("SNP","CHR","BP"),sep = "\\_") 

#trim out the first row, since the first BP is repeated and will cause GenWin to choke
file <- file[-1,]

#set CHR and BP to numeric mode, and then trim all data not associated with scaffolds 1-14 (note to Ali: any useful data in CHR 15-18? ask JK)
file$CHR <- as.numeric(file$CHR)
file$BP <- as.numeric(file$BP)
file <- subset(file, CHR<15)

#set CHR back to factor mode in order to identify the chromosomes in the .txt file
file$CHR <- as.factor(file$CHR)
CHRnames <- levels(file$CHR)

#set CHRnames to numeric mode so that a logical statement in for loop will work
CHRnames <- as.numeric(CHRnames)
file$CHR <- as.numeric(file$CHR)

#Break B1.txt into different chromosomes
#for each chromosome, use GenWin to identify window breaks

#create an empty list
windowsList <- list()

for(input in CHRnames){
  
  #take data from a single chromosome and make sure it is ordered by base pair
  window1 <- subset(file, CHR==input)
  window1 <- window1[order(window1$BP),]
  
  #use GenWin program to identify windows.  We may want to try different smoothness levels
  spline <- splineAnalyze(Y=window1$B,map=window1$BP,smoothness=100, plotRaw=TRUE,plotWindows=TRUE,method=4)
  
  #add a column to the spline$windowData output to identify the chromosome under consideration
  CHRcol <- rep(input,dim(spline$windowData)[1])
  spline$windowData <- cbind(CHRcol,spline$windowData)
  
  #add the output to a list of output from all chromosomes
  windowsList[[input]] <- spline$windowData
}

#take the output from all chromosomes and bind together into a single matrix called splinewindows
splinewindows = do.call(rbind,windowsList)
  
#write the splinewindows matrix to a txt file
write.table(splinewindows, file = "splinewindows.txt", quote=FALSE, sep="\t", row.names = FALSE)