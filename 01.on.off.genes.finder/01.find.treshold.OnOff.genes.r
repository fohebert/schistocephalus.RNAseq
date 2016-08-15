### REQUIRED LIBRARIES ###
library(edgeR)
library(limma) # limma is supposed to be part of edgeR though

#########################
### 1. IMPORTING DATA ###
#########################

# Setting working directory
setwd("~/YOUR/WORKING/DIRECTORY")

# The data that is imported here is a read-count matrix, i.e. rows = genes & columns = samples
# Read-count matrices can be obtained using various methods, among others we find
# HTseq or custom scripts based on counting the number of RNA-seq reads that align onto
# each reference sequence (by parsing alignement SAM file(s)). In this case, I obtained
# the read-count matrix by running the program "Corset" (https://github.com/Oshlack/Corset/wiki)
# on my raw data. You can find the pipeline I used on my github page
# (https://github.com/fohebert/corset_pipeline).
data <- read.delim("ssol.OliWorms.read-count.2016-04-29.txt", row.names=1)

# Telling R which sampl belongs to which group
group <- factor(c("NonInfect","NonInfect","NonInfect","NonInfect",
	"NonInfect","NonInfect","NonInfect","Infect","Infect","Infect",
	"Infect","Adult","Adult","Adult"))
dge <- DGEList(counts=data,group=group)

###################################
### 2. PLOTTING RAW READ COUNTS ###
###################################

# Storing all of the log2(CPM) values for the read counts into a vector.
# This vector will later be used to produce a histogram.
all.counts.log2CPM <- as.vector(log2(cpm(dge)))
all.counts.CPM <- as.vector(cpm(dge))
all.counts.raw <- as.vector(dge$counts)

# Building the histogram
hist(all.counts.log2CPM,
	main = "Read count distribution (all samples)",
	xlab = "log2(CPM)",
	ylab = "Frequency")

# Getting the 5th percentile of the distribution generated with the histogram
# This will indicate which log2(CPM) value corresponds to the bottom 5% of the
# distribution. This log2(CPM) would be the threshold to use in order to assess
# if a transcript can be considered not expressed at all in a given sample. This
# would be the very basics towards identifying on/off genes.
quantile(all.counts.log2CPM, 0.05)
#       5% 
# -0.857688

# So the result we get is log2(CPM) = -0.857688. This means that below this value
# we can consider that the transcript is not expressed and its expression is 0.
# By calculating the inverse of log2: http://www.rapidtables.com/calc/math/anti-log-calculator.htm
# I found that log2(CPM) = -0.857688 corresponds to a CPM value of ~0.55. So if
# we want to use a CPM threshold to label "zero expression transcripts", we could
# use the value of 0.55; everything below that is considered like not expressed.

######################################################
### 3. PREPARING FOR THE FILTERING OF ON/OFF GENES ###
######################################################

# Converting the dge object into a data frame in order to be able to
# adequately subset the object, i.e. take a subsample based on the
# threshold conditions previously established.
df <- as.data.frame(cpm(dge))
colnames(df) <- c("NI.3","NI.4","NI.1","NI.5","NI.6","NI.7",
    "NI.2","I.1","I.2","I.4","I.3","A.1","A.2","A.3")
# Writing data frame to output file
write.table(df, file = "ssol.CPM.table.all-sp.2016-06-06.txt")