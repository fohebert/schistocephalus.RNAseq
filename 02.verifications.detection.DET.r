###################################
# PACKAGES THAT NEED TO BE UPLOADED
###################################
library(limma)
library(edgeR)
library(VennDiagram)
library(calibrate)
library(Tmisc)

#############################
# -      DESCRIPTION      - #
#############################
# This script looks at the impact of various detection thresholds for
# differentially expressed transcripts. Parameters that are targeted by
# this scripts are: minimum CPM value, minimum number of individuals in
# which the CPM threshold is respected. 

####################################
# -      PREPARING THE DATA      - #
####################################

# Setting working directory on my computer
setwd("~/YOUR/WORKING/DIRECTORY/")
# Importing the data
data <- read.delim("Schs.all-sp.counts.customSPnames.txt", row.names=1)
# First look at the dataset
head(data, n = 10)
# Dataset structure
str(data)

# Creating a vector with ALL of the read counts, i.e. every single
# cell in the original count table called "data". We want to see
# what the distribution of these read counts looks like.
all.read.counts <- c(
	data$NI.1, data$NI.2, data$NI.3, data$NI.4, 
	data$NI.5, data$NI.6, data$NI.7, data$I.1, data$I.2, data$I.3, 
	data$I.4, data$A.1, data$A.2, data$A.3
	)
# Check the vector
str(all.read.counts)
# Quantiles, mean, median, min/max values for that vectoro
summary(all.read.counts)
# Transforming all of the read counts in the vector into
# log values because there are extreme values in the
# dataset, i.e. values ranging from 0 to 454,200. It is
# super hard to visualize if we do not convert into log
# values.
all.read.counts.log <- log(all.read.counts)

#########################################
### -       PLOTTING THE DATA       - ###
###     (read count distribution)     ###
#########################################
# These graphs are made with 
# RAW DATA (i.e. raw counts)

# Boxplot
boxplot(
	all.read.counts.log, 
	main = "Read count distribution (log-transformed)", ylab = "log(read.count)", 
	ylim = c(0,14), 
	xlab = "Unfiltered unigenes from Corset's output\n(115 702 transcripts in 14 samples)"
	)

# Histogram
hist(
	all.read.counts.log, 
	main = "Read count distribution (log-transformed)\n(115 702 transcripts in 14 samples)", 
	ylim = c(0,0.45), 
	prob = T, 
	xlab = "log(read.count)"
	)

####################################################
### COMPARISON FILTER BEFORE/AFTER NORMALIZATION ###
####################################################

# The idea is to go through the whole DGE analysis
# by creating 2 different DGE objects: one with
# filtration process BEFORE normalizing for library
# size and one with filtration process AFTER 
# normalizing for library size.

######################################################
# (1) Preparing the data and producing the initial DGE
# ----------------------------------------------------

# Vector of the groups i.e. non-infective, infective, adults)
group <- factor(c("NonInfect","NonInfect","NonInfect",
	"NonInfect","NonInfect","NonInfect","NonInfect",
	"Infect","Infect","Infect","Infect","Adult","Adult","Adult"))
# Creating edgeR's DGE object
dge <- DGEList(counts=data,group=group)

###########################################
# (2) NORMALIZATION AND FILTERING PROCESSES
# -----------------------------------------
###########################################

######################################
# (A) TRANSITION = ADULT vs. INFECTIVE
# ------------------------------------

# -----------------------------
# Filtering BEFORE normalizing:
# -----------------------------

# All possible combinations for thresholds
# on CPM and no. of samples:

# ------ AvsI | filter FIRST | CPM >= 0 | in at least 1-3 samples ------

# CPM >= 0 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge)>=0) >= 1
AvsI.filtFirst.cpm0.sp1.dge <- dge[keep.sp1,]
AvsI.filtFirst.cpm0.sp1.dge$samples$lib.size <- colSums(AvsI.filtFirst.cpm0.sp1.dge$counts)
# Normalization
AvsI.filtFirst.cpm0.sp1.dge <- calcNormFactors(AvsI.filtFirst.cpm0.sp1.dge)

# CPM >= 0 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge)>=0) >= 2
AvsI.filtFirst.cpm0.sp2.dge <- dge[keep.sp2,]
AvsI.filtFirst.cpm0.sp2.dge$samples$lib.size <- colSums(AvsI.filtFirst.cpm0.sp2.dge$counts)
# Normalization
AvsI.filtFirst.cpm0.sp2.dge <- calcNormFactors(AvsI.filtFirst.cpm0.sp2.dge)

# CPM >= 0 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge)>=0) >= 3
AvsI.filtFirst.cpm0.sp3.dge <- dge[keep.sp3,]
AvsI.filtFirst.cpm0.sp3.dge$samples$lib.size <- colSums(AvsI.filtFirst.cpm0.sp3.dge$counts)
# Normalization
AvsI.filtFirst.cpm0.sp3.dge <- calcNormFactors(AvsI.filtFirst.cpm0.sp3.dge)

# ------ AvsI | filter FIRST | CPM >= 5 | in at least 1-3 samples ------

# CPM >= 5 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge)>=5) >= 1
AvsI.filtFirst.cpm05.sp1.dge <- dge[keep.sp1,]
AvsI.filtFirst.cpm05.sp1.dge$samples$lib.size <- colSums(AvsI.filtFirst.cpm05.sp1.dge$counts)
# Normalization
AvsI.filtFirst.cpm05.sp1.dge <- calcNormFactors(AvsI.filtFirst.cpm05.sp1.dge)

# CPM >= 5 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge)>=5) >= 2
AvsI.filtFirst.cpm05.sp2.dge <- dge[keep.sp2,]
AvsI.filtFirst.cpm05.sp2.dge$samples$lib.size <- colSums(AvsI.filtFirst.cpm05.sp2.dge$counts)
# Normalization
AvsI.filtFirst.cpm05.sp2.dge <- calcNormFactors(AvsI.filtFirst.cpm05.sp2.dge)

# CPM >= 5 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge)>=5) >= 3
AvsI.filtFirst.cpm05.sp3.dge <- dge[keep.sp3,]
AvsI.filtFirst.cpm05.sp3.dge$samples$lib.size <- colSums(AvsI.filtFirst.cpm05.sp3.dge$counts)
# Normalization
AvsI.filtFirst.cpm05.sp3.dge <- calcNormFactors(AvsI.filtFirst.cpm05.sp3.dge)

# ------ AvsI | filter FIRST | CPM >= 10 | in at least 1-3 samples ------

# CPM >= 10 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge)>=10) >= 1
AvsI.filtFirst.cpm10.sp1.dge <- dge[keep.sp1,]
AvsI.filtFirst.cpm10.sp1.dge$samples$lib.size <- colSums(AvsI.filtFirst.cpm10.sp1.dge$counts)
# Normalization
AvsI.filtFirst.cpm10.sp1.dge <- calcNormFactors(AvsI.filtFirst.cpm10.sp1.dge)

# CPM >= 10 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge)>=10) >= 2
AvsI.filtFirst.cpm10.sp2.dge <- dge[keep.sp2,]
AvsI.filtFirst.cpm10.sp2.dge$samples$lib.size <- colSums(AvsI.filtFirst.cpm10.sp2.dge$counts)
# Normalization
AvsI.filtFirst.cpm10.sp2.dge <- calcNormFactors(AvsI.filtFirst.cpm10.sp2.dge)

# CPM >= 10 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge)>=10) >= 3
AvsI.filtFirst.cpm10.sp3.dge <- dge[keep.sp3,]
AvsI.filtFirst.cpm10.sp3.dge$samples$lib.size <- colSums(AvsI.filtFirst.cpm10.sp3.dge$counts)
# Normalization
AvsI.filtFirst.cpm10.sp3.dge <- calcNormFactors(AvsI.filtFirst.cpm10.sp3.dge)

# ------ AvsI | filter FIRST | CPM >= 15 | in at least 1-3 samples ------

# CPM >= 15 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge)>=15) >= 1
AvsI.filtFirst.cpm15.sp1.dge <- dge[keep.sp1,]
AvsI.filtFirst.cpm15.sp1.dge$samples$lib.size <- colSums(AvsI.filtFirst.cpm15.sp1.dge$counts)
# Normalization
AvsI.filtFirst.cpm15.sp1.dge <- calcNormFactors(AvsI.filtFirst.cpm15.sp1.dge)

# CPM >= 15 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge)>=15) >= 2
AvsI.filtFirst.cpm15.sp2.dge <- dge[keep.sp2,]
AvsI.filtFirst.cpm15.sp2.dge$samples$lib.size <- colSums(AvsI.filtFirst.cpm15.sp2.dge$counts)
# Normalization
AvsI.filtFirst.cpm15.sp2.dge <- calcNormFactors(AvsI.filtFirst.cpm15.sp2.dge)

# CPM >= 15 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge)>=15) >= 3
AvsI.filtFirst.cpm15.sp3.dge <- dge[keep.sp3,]
AvsI.filtFirst.cpm15.sp3.dge$samples$lib.size <- colSums(AvsI.filtFirst.cpm15.sp3.dge$counts)
# Normalization
AvsI.filtFirst.cpm15.sp3.dge <- calcNormFactors(AvsI.filtFirst.cpm15.sp3.dge)

# ----------------------------
# Filtering AFTER normalizing:
# ----------------------------

# All possible combinations for thresholds 
# on CPM and no. of samples:

# ------ AvsI | filter AFTER | CPM >= 0 | in at least 1-3 samples ------

# Normalization
dge.norm <- calcNormFactors(dge)

# CPM >= 0 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge.norm)>=0) >= 1
AvsI.filtAfter.cpm0.sp1.dge <- dge.norm[keep.sp1,]
AvsI.filtAfter.cpm0.sp1.dge$samples$lib.size <- colSums(AvsI.filtAfter.cpm0.sp1.dge$counts)

# CPM >= 0 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge.norm)>=0) >= 2
AvsI.filtAfter.cpm0.sp2.dge <- dge.norm[keep.sp2,]
AvsI.filtAfter.cpm0.sp2.dge$samples$lib.size <- colSums(AvsI.filtAfter.cpm0.sp2.dge$counts)

# CPM >= 0 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge.norm)>=0) >= 3
AvsI.filtAfter.cpm0.sp3.dge <- dge.norm[keep.sp3,]
AvsI.filtAfter.cpm0.sp3.dge$samples$lib.size <- colSums(AvsI.filtAfter.cpm0.sp3.dge$counts)

# ------ AvsI | filter AFTER | CPM >= 5 | in at least 1-3 samples ------

# CPM >= 5 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge.norm)>=5) >= 1
AvsI.filtAfter.cpm05.sp1.dge <- dge.norm[keep.sp1,]
AvsI.filtAfter.cpm05.sp1.dge$samples$lib.size <- colSums(AvsI.filtAfter.cpm05.sp1.dge$counts)

# CPM >= 5 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge.norm)>=5) >= 2
AvsI.filtAfter.cpm05.sp2.dge <- dge.norm[keep.sp2,]
AvsI.filtAfter.cpm05.sp2.dge$samples$lib.size <- colSums(AvsI.filtAfter.cpm05.sp2.dge$counts)

# CPM >= 5 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge.norm)>=5) >= 3
AvsI.filtAfter.cpm05.sp3.dge <- dge.norm[keep.sp3,]
AvsI.filtAfter.cpm05.sp3.dge$samples$lib.size <- colSums(AvsI.filtAfter.cpm05.sp3.dge$counts)

# ------ AvsI | filter AFTER | CPM >= 10 | in at least 1-3 samples ------

# CPM >= 10 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge.norm)>=10) >= 1
AvsI.filtAfter.cpm10.sp1.dge <- dge.norm[keep.sp1,]
AvsI.filtAfter.cpm10.sp1.dge$samples$lib.size <- colSums(AvsI.filtAfter.cpm10.sp1.dge$counts)

# CPM >= 10 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge.norm)>=10) >= 2
AvsI.filtAfter.cpm10.sp2.dge <- dge.norm[keep.sp2,]
AvsI.filtAfter.cpm10.sp2.dge$samples$lib.size <- colSums(AvsI.filtAfter.cpm10.sp2.dge$counts)

# CPM >= 10 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge.norm)>=10) >= 3
AvsI.filtAfter.cpm10.sp3.dge <- dge.norm[keep.sp3,]
AvsI.filtAfter.cpm10.sp3.dge$samples$lib.size <- colSums(AvsI.filtAfter.cpm10.sp3.dge$counts)

# ------ AvsI | filter AFTER | CPM >= 15 | in at least 1-3 samples ------

# CPM >= 15 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge.norm)>=15) >= 1
AvsI.filtAfter.cpm15.sp1.dge <- dge.norm[keep.sp1,]
AvsI.filtAfter.cpm15.sp1.dge$samples$lib.size <- colSums(AvsI.filtAfter.cpm15.sp1.dge$counts)

# CPM >= 15 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge.norm)>=15) >= 2
AvsI.filtAfter.cpm15.sp2.dge <- dge.norm[keep.sp2,]
AvsI.filtAfter.cpm15.sp2.dge$samples$lib.size <- colSums(AvsI.filtAfter.cpm15.sp2.dge$counts)

# CPM >= 15 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge.norm)>=15) >= 3
AvsI.filtAfter.cpm15.sp3.dge <- dge.norm[keep.sp3,]
AvsI.filtAfter.cpm15.sp3.dge$samples$lib.size <- colSums(AvsI.filtAfter.cpm15.sp3.dge$counts)


##############################################
# (B) TRANSITION = INFECTIVE vs. NON-INFECTIVE
# --------------------------------------------

# -----------------------------
# Filtering BEFORE normalizing:
# -----------------------------

# All possible combinations for thresholds 
# on CPM and no. of samples:

# ------ IvsNI | filter FIRST | CPM >= 0 | in at least 1-3 samples ------

# CPM >= 0 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge)>=0) >= 1
IvsNI.filtFirst.cpm0.sp1.dge <- dge[keep.sp1,]
IvsNI.filtFirst.cpm0.sp1.dge$samples$lib.size <- colSums(IvsNI.filtFirst.cpm0.sp1.dge$counts)
# Normalization
IvsNI.filtFirst.cpm0.sp1.dge <- calcNormFactors(IvsNI.filtFirst.cpm0.sp1.dge)

# CPM >= 0 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge)>=0) >= 2
IvsNI.filtFirst.cpm0.sp2.dge <- dge[keep.sp2,]
IvsNI.filtFirst.cpm0.sp2.dge$samples$lib.size <- colSums(IvsNI.filtFirst.cpm0.sp2.dge$counts)
# Normalization
IvsNI.filtFirst.cpm0.sp2.dge <- calcNormFactors(IvsNI.filtFirst.cpm0.sp2.dge)

# CPM >= 0 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge)>=0) >= 3
IvsNI.filtFirst.cpm0.sp3.dge <- dge[keep.sp3,]
IvsNI.filtFirst.cpm0.sp3.dge$samples$lib.size <- colSums(IvsNI.filtFirst.cpm0.sp3.dge$counts)
# Normalization
IvsNI.filtFirst.cpm0.sp3.dge <- calcNormFactors(IvsNI.filtFirst.cpm0.sp3.dge)

# ------ IvsNI | filter FIRST | CPM >= 5 | in at least 1-3 samples ------

# CPM >= 5 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge)>=5) >= 1
IvsNI.filtFirst.cpm05.sp1.dge <- dge[keep.sp1,]
IvsNI.filtFirst.cpm05.sp1.dge$samples$lib.size <- colSums(IvsNI.filtFirst.cpm05.sp1.dge$counts)
# Normalization
IvsNI.filtFirst.cpm05.sp1.dge <- calcNormFactors(IvsNI.filtFirst.cpm05.sp1.dge)

# CPM >= 5 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge)>=5) >= 2
IvsNI.filtFirst.cpm05.sp2.dge <- dge[keep.sp2,]
IvsNI.filtFirst.cpm05.sp2.dge$samples$lib.size <- colSums(IvsNI.filtFirst.cpm05.sp2.dge$counts)
# Normalization
IvsNI.filtFirst.cpm05.sp2.dge <- calcNormFactors(IvsNI.filtFirst.cpm05.sp2.dge)

# CPM >= 5 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge)>=5) >= 3
IvsNI.filtFirst.cpm05.sp3.dge <- dge[keep.sp3,]
IvsNI.filtFirst.cpm05.sp3.dge$samples$lib.size <- colSums(IvsNI.filtFirst.cpm05.sp3.dge$counts)
# Normalization
IvsNI.filtFirst.cpm05.sp3.dge <- calcNormFactors(IvsNI.filtFirst.cpm05.sp3.dge)

# ------ IvsNI | filter FIRST | CPM >= 10 | in at least 1-3 samples ------

# CPM >= 10 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge)>=10) >= 1
IvsNI.filtFirst.cpm10.sp1.dge <- dge[keep.sp1,]
IvsNI.filtFirst.cpm10.sp1.dge$samples$lib.size <- colSums(IvsNI.filtFirst.cpm10.sp1.dge$counts)
# Normalization
IvsNI.filtFirst.cpm10.sp1.dge <- calcNormFactors(IvsNI.filtFirst.cpm10.sp1.dge)

# CPM >= 10 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge)>=10) >= 2
IvsNI.filtFirst.cpm10.sp2.dge <- dge[keep.sp2,]
IvsNI.filtFirst.cpm10.sp2.dge$samples$lib.size <- colSums(IvsNI.filtFirst.cpm10.sp2.dge$counts)
# Normalization
IvsNI.filtFirst.cpm10.sp2.dge <- calcNormFactors(IvsNI.filtFirst.cpm10.sp2.dge)

# CPM >= 10 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge)>=10) >= 3
IvsNI.filtFirst.cpm10.sp3.dge <- dge[keep.sp3,]
IvsNI.filtFirst.cpm10.sp3.dge$samples$lib.size <- colSums(IvsNI.filtFirst.cpm10.sp3.dge$counts)
# Normalization
IvsNI.filtFirst.cpm10.sp3.dge <- calcNormFactors(IvsNI.filtFirst.cpm10.sp3.dge)

# ------ IvsNI | filter FIRST | CPM >= 15 | in at least 1-3 samples ------

# CPM >= 15 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge)>=15) >= 1
IvsNI.filtFirst.cpm15.sp1.dge <- dge[keep.sp1,]
IvsNI.filtFirst.cpm15.sp1.dge$samples$lib.size <- colSums(IvsNI.filtFirst.cpm15.sp1.dge$counts)
# Normalization
IvsNI.filtFirst.cpm15.sp1.dge <- calcNormFactors(IvsNI.filtFirst.cpm15.sp1.dge)

# CPM >= 15 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge)>=15) >= 2
IvsNI.filtFirst.cpm15.sp2.dge <- dge[keep.sp2,]
IvsNI.filtFirst.cpm15.sp2.dge$samples$lib.size <- colSums(IvsNI.filtFirst.cpm15.sp2.dge$counts)
# Normalization
IvsNI.filtFirst.cpm15.sp2.dge <- calcNormFactors(IvsNI.filtFirst.cpm15.sp2.dge)

# CPM >= 15 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge)>=15) >= 3
IvsNI.filtFirst.cpm15.sp3.dge <- dge[keep.sp3,]
IvsNI.filtFirst.cpm15.sp3.dge$samples$lib.size <- colSums(IvsNI.filtFirst.cpm15.sp3.dge$counts)
# Normalization
IvsNI.filtFirst.cpm15.sp3.dge <- calcNormFactors(IvsNI.filtFirst.cpm15.sp3.dge)

# ----------------------------
# Filtering AFTER normalizing:
# ----------------------------

# All possible combinations for thresholds 
# on CPM and no. of samples:

# ------ IvsNI | filter AFTER | CPM >= 0 | in at least 1-3 samples ------

# CPM >= 0 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge.norm)>=0) >= 1
IvsNI.filtAfter.cpm0.sp1.dge <- dge.norm[keep.sp1,]
IvsNI.filtAfter.cpm0.sp1.dge$samples$lib.size <- colSums(IvsNI.filtAfter.cpm0.sp1.dge$counts)

# CPM >= 0 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge.norm)>=0) >= 2
IvsNI.filtAfter.cpm0.sp2.dge <- dge.norm[keep.sp2,]
IvsNI.filtAfter.cpm0.sp2.dge$samples$lib.size <- colSums(IvsNI.filtAfter.cpm0.sp2.dge$counts)

# CPM >= 0 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge.norm)>=0) >= 3
IvsNI.filtAfter.cpm0.sp3.dge <- dge.norm[keep.sp3,]
IvsNI.filtAfter.cpm0.sp3.dge$samples$lib.size <- colSums(IvsNI.filtAfter.cpm0.sp3.dge$counts)

# ------ IvsNI | filter AFTER | CPM >= 5 | in at least 1-3 samples ------

# CPM >= 5 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge.norm)>=5) >= 1
IvsNI.filtAfter.cpm05.sp1.dge <- dge.norm[keep.sp1,]
IvsNI.filtAfter.cpm05.sp1.dge$samples$lib.size <- colSums(IvsNI.filtAfter.cpm05.sp1.dge$counts)

# CPM >= 5 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge.norm)>=5) >= 2
IvsNI.filtAfter.cpm05.sp2.dge <- dge.norm[keep.sp2,]
IvsNI.filtAfter.cpm05.sp2.dge$samples$lib.size <- colSums(IvsNI.filtAfter.cpm05.sp2.dge$counts)

# CPM >= 5 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge.norm)>=5) >= 3
IvsNI.filtAfter.cpm05.sp3.dge <- dge.norm[keep.sp3,]
IvsNI.filtAfter.cpm05.sp3.dge$samples$lib.size <- colSums(IvsNI.filtAfter.cpm05.sp3.dge$counts)

# ------ IvsNI | filter AFTER | CPM >= 10 | in at least 1-3 samples ------

# CPM >= 10 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge.norm)>=10) >= 1
IvsNI.filtAfter.cpm10.sp1.dge <- dge.norm[keep.sp1,]
IvsNI.filtAfter.cpm10.sp1.dge$samples$lib.size <- colSums(IvsNI.filtAfter.cpm10.sp1.dge$counts)

# CPM >= 10 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge.norm)>=10) >= 2
IvsNI.filtAfter.cpm10.sp2.dge <- dge.norm[keep.sp2,]
IvsNI.filtAfter.cpm10.sp2.dge$samples$lib.size <- colSums(IvsNI.filtAfter.cpm10.sp2.dge$counts)

# CPM >= 10 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge.norm)>=10) >= 3
IvsNI.filtAfter.cpm10.sp3.dge <- dge.norm[keep.sp3,]
IvsNI.filtAfter.cpm10.sp3.dge$samples$lib.size <- colSums(IvsNI.filtAfter.cpm10.sp3.dge$counts)

# ------ IvsNI | filter AFTER | CPM >= 15 | in at least 1-3 samples ------

# CPM >= 15 in at least 1 sample
keep.sp1 <- rowSums(cpm(dge.norm)>=15) >= 1
IvsNI.filtAfter.cpm15.sp1.dge <- dge.norm[keep.sp1,]
IvsNI.filtAfter.cpm15.sp1.dge$samples$lib.size <- colSums(IvsNI.filtAfter.cpm15.sp1.dge$counts)

# CPM >= 15 in at least 2 samples
keep.sp2 <- rowSums(cpm(dge.norm)>=15) >= 2
IvsNI.filtAfter.cpm15.sp2.dge <- dge.norm[keep.sp2,]
IvsNI.filtAfter.cpm15.sp2.dge$samples$lib.size <- colSums(IvsNI.filtAfter.cpm15.sp2.dge$counts)

# CPM >= 15 in at least 3 samples
keep.sp3 <- rowSums(cpm(dge.norm)>=15) >= 3
IvsNI.filtAfter.cpm15.sp3.dge <- dge.norm[keep.sp3,]
IvsNI.filtAfter.cpm15.sp3.dge$samples$lib.size <- colSums(IvsNI.filtAfter.cpm15.sp3.dge$counts)

###################
# (3) DESIGN MATRIX
# -----------------
###################

# Producing a matrix design needed for the GLM. 
# NB: this design will be the same for all
# DGE objects.
d <- model.matrix(~0+group, data=dge$samples)

#########################
# (4) VOOM TRANSFORMATION
# -----------------------
#########################

# -----------------------------------
# (A) TRANSITION: ADULT vs. INFECTIVE
# -----------------------------------

# DGE with filtering BEFORE normalization:
# ----------------------------------------

# ---- CPM >= 0 | in >= 1 sample ----
v.AvsI.filtFirst.cpm0.sp1.dge <- voom(AvsI.filtFirst.cpm0.sp1.dge, d, plot=F)
# ---- CPM >= 0 | in >= 2 samples ----
v.AvsI.filtFirst.cpm0.sp2.dge <- voom(AvsI.filtFirst.cpm0.sp2.dge, d, plot=F)
# ---- CPM >= 0 | in >= 3 samples ----
v.AvsI.filtFirst.cpm0.sp3.dge <- voom(AvsI.filtFirst.cpm0.sp3.dge, d, plot=F)

# ---- CPM >= 5 | in >= 1 sample ----
v.AvsI.filtFirst.cpm05.sp1.dge <- voom(AvsI.filtFirst.cpm05.sp1.dge, d, plot=F)
# ---- CPM >= 5 | in >= 2 samples ----
v.AvsI.filtFirst.cpm05.sp2.dge <- voom(AvsI.filtFirst.cpm05.sp2.dge, d, plot=F)
# ---- CPM >= 5 | in >= 3 samples ----
v.AvsI.filtFirst.cpm05.sp3.dge <- voom(AvsI.filtFirst.cpm05.sp3.dge, d, plot=F)

# ---- CPM >= 10 | in >= 1 sample ----
v.AvsI.filtFirst.cpm10.sp1.dge <- voom(AvsI.filtFirst.cpm10.sp1.dge, d, plot=F)
# ---- CPM >= 10 | in >= 2 samples ----
v.AvsI.filtFirst.cpm10.sp2.dge <- voom(AvsI.filtFirst.cpm10.sp2.dge, d, plot=F)
# ---- CPM >= 10 | in >= 3 samples ----
v.AvsI.filtFirst.cpm10.sp3.dge <- voom(AvsI.filtFirst.cpm10.sp3.dge, d, plot=F)

# ---- CPM >= 15 | in >= 1 sample ----
v.AvsI.filtFirst.cpm15.sp1.dge <- voom(AvsI.filtFirst.cpm15.sp1.dge, d, plot=F)
# ---- CPM >= 15 | in >= 2 samples ----
v.AvsI.filtFirst.cpm15.sp2.dge <- voom(AvsI.filtFirst.cpm15.sp2.dge, d, plot=F)
# ---- CPM >= 15 | in >= 3 samples ----
v.AvsI.filtFirst.cpm15.sp3.dge <- voom(AvsI.filtFirst.cpm15.sp3.dge, d, plot=F)

# DGE with filtering AFTER normalization:
# ----------------------------------------

# ---- CPM >= 0 | in >= 1 sample ----
v.AvsI.filtAfter.cpm0.sp1.dge <- voom(AvsI.filtAfter.cpm0.sp1.dge, d, plot=F)
# ---- CPM >= 0 | in >= 2 samples ----
v.AvsI.filtAfter.cpm0.sp2.dge <- voom(AvsI.filtAfter.cpm0.sp2.dge, d, plot=F)
# ---- CPM >= 0 | in >= 3 samples ----
v.AvsI.filtAfter.cpm0.sp3.dge <- voom(AvsI.filtAfter.cpm0.sp3.dge, d, plot=F)

# ---- CPM >= 5 | in >= 1 sample ----
v.AvsI.filtAfter.cpm05.sp1.dge <- voom(AvsI.filtAfter.cpm05.sp1.dge, d, plot=F)
# ---- CPM >= 5 | in >= 2 samples ----
v.AvsI.filtAfter.cpm05.sp2.dge <- voom(AvsI.filtAfter.cpm05.sp2.dge, d, plot=F)
# ---- CPM >= 5 | in >= 3 samples ----
v.AvsI.filtAfter.cpm05.sp3.dge <- voom(AvsI.filtAfter.cpm05.sp3.dge, d, plot=F)

# ---- CPM >= 10 | in >= 1 sample ----
v.AvsI.filtAfter.cpm10.sp1.dge <- voom(AvsI.filtAfter.cpm10.sp1.dge, d, plot=F)
# ---- CPM >= 10 | in >= 2 samples ----
v.AvsI.filtAfter.cpm10.sp2.dge <- voom(AvsI.filtAfter.cpm10.sp2.dge, d, plot=F)
# ---- CPM >= 10 | in >= 3 samples ----
v.AvsI.filtAfter.cpm10.sp3.dge <- voom(AvsI.filtAfter.cpm10.sp3.dge, d, plot=F)

# ---- CPM >= 15 | in >= 1 sample ----
v.AvsI.filtAfter.cpm15.sp1.dge <- voom(AvsI.filtAfter.cpm15.sp1.dge, d, plot=F)
# ---- CPM >= 15 | in >= 2 samples ----
v.AvsI.filtAfter.cpm15.sp2.dge <- voom(AvsI.filtAfter.cpm15.sp2.dge, d, plot=F)
# ---- CPM >= 15 | in >= 3 samples ----
v.AvsI.filtAfter.cpm15.sp3.dge <- voom(AvsI.filtAfter.cpm15.sp3.dge, d, plot=F)

# -------------------------------------------
# (B) TRANSITION: INFECTIVE vs. NON-INFECTIVE
# -------------------------------------------

# DGE with filtering BEFORE normalization:
# ----------------------------------------

# ---- CPM >= 0 | in >= 1 sample ----
v.IvsNI.filtFirst.cpm0.sp1.dge <- voom(IvsNI.filtFirst.cpm0.sp1.dge, d, plot=F)
# ---- CPM >= 0 | in >= 2 samples ----
v.IvsNI.filtFirst.cpm0.sp2.dge <- voom(IvsNI.filtFirst.cpm0.sp2.dge, d, plot=F)
# ---- CPM >= 0 | in >= 3 samples ----
v.IvsNI.filtFirst.cpm0.sp3.dge <- voom(IvsNI.filtFirst.cpm0.sp3.dge, d, plot=F)

# ---- CPM >= 5 | in >= 1 sample ----
v.IvsNI.filtFirst.cpm05.sp1.dge <- voom(IvsNI.filtFirst.cpm05.sp1.dge, d, plot=F)
# ---- CPM >= 5 | in >= 2 samples ----
v.IvsNI.filtFirst.cpm05.sp2.dge <- voom(IvsNI.filtFirst.cpm05.sp2.dge, d, plot=F)
# ---- CPM >= 5 | in >= 3 samples ----
v.IvsNI.filtFirst.cpm05.sp3.dge <- voom(IvsNI.filtFirst.cpm05.sp3.dge, d, plot=F)

# ---- CPM >= 10 | in >= 1 sample ----
v.IvsNI.filtFirst.cpm10.sp1.dge <- voom(IvsNI.filtFirst.cpm10.sp1.dge, d, plot=F)
# ---- CPM >= 10 | in >= 2 samples ----
v.IvsNI.filtFirst.cpm10.sp2.dge <- voom(IvsNI.filtFirst.cpm10.sp2.dge, d, plot=F)
# ---- CPM >= 10 | in >= 3 samples ----
v.IvsNI.filtFirst.cpm10.sp3.dge <- voom(IvsNI.filtFirst.cpm10.sp3.dge, d, plot=F)

# ---- CPM >= 15 | in >= 1 sample ----
v.IvsNI.filtFirst.cpm15.sp1.dge <- voom(IvsNI.filtFirst.cpm15.sp1.dge, d, plot=F)
# ---- CPM >= 15 | in >= 2 samples ----
v.IvsNI.filtFirst.cpm15.sp2.dge <- voom(IvsNI.filtFirst.cpm15.sp2.dge, d, plot=F)
# ---- CPM >= 15 | in >= 3 samples ----
v.IvsNI.filtFirst.cpm15.sp3.dge <- voom(IvsNI.filtFirst.cpm15.sp3.dge, d, plot=F)

# DGE with filtering AFTER normalization:
# ----------------------------------------

# ---- CPM >= 0 | in >= 1 sample ----
v.IvsNI.filtAfter.cpm0.sp1.dge <- voom(IvsNI.filtAfter.cpm0.sp1.dge, d, plot=F)
# ---- CPM >= 0 | in >= 2 samples ----
v.IvsNI.filtAfter.cpm0.sp2.dge <- voom(IvsNI.filtAfter.cpm0.sp2.dge, d, plot=F)
# ---- CPM >= 0 | in >= 3 samples ----
v.IvsNI.filtAfter.cpm0.sp3.dge <- voom(IvsNI.filtAfter.cpm0.sp3.dge, d, plot=F)

# ---- CPM >= 5 | in >= 1 sample ----
v.IvsNI.filtAfter.cpm05.sp1.dge <- voom(IvsNI.filtAfter.cpm05.sp1.dge, d, plot=F)
# ---- CPM >= 5 | in >= 2 samples ----
v.IvsNI.filtAfter.cpm05.sp2.dge <- voom(IvsNI.filtAfter.cpm05.sp2.dge, d, plot=F)
# ---- CPM >= 5 | in >= 3 samples ----
v.IvsNI.filtAfter.cpm05.sp3.dge <- voom(IvsNI.filtAfter.cpm05.sp3.dge, d, plot=F)

# ---- CPM >= 10 | in >= 1 sample ----
v.IvsNI.filtAfter.cpm10.sp1.dge <- voom(IvsNI.filtAfter.cpm10.sp1.dge, d, plot=F)
# ---- CPM >= 10 | in >= 2 samples ----
v.IvsNI.filtAfter.cpm10.sp2.dge <- voom(IvsNI.filtAfter.cpm10.sp2.dge, d, plot=F)
# ---- CPM >= 10 | in >= 3 samples ----
v.IvsNI.filtAfter.cpm10.sp3.dge <- voom(IvsNI.filtAfter.cpm10.sp3.dge, d, plot=F)

# ---- CPM >= 15 | in >= 1 sample ----
v.IvsNI.filtAfter.cpm15.sp1.dge <- voom(IvsNI.filtAfter.cpm15.sp1.dge, d, plot=F)
# ---- CPM >= 15 | in >= 2 samples ----
v.IvsNI.filtAfter.cpm15.sp2.dge <- voom(IvsNI.filtAfter.cpm15.sp2.dge, d, plot=F)
# ---- CPM >= 15 | in >= 3 samples ----
v.IvsNI.filtAfter.cpm15.sp3.dge <- voom(IvsNI.filtAfter.cpm15.sp3.dge, d, plot=F)

# ###################################################
# (5) FIT TO GLM MODEL & PRODUCE TOP TABLE WITH DET
#   (DET = Differentially Expressed Transcripts)  #
# -------------------------------------------------
###################################################

# Building contrast matrix that will be used
# (same matrix for all of the groups)
contrast.matrix <- makeContrasts(groupAdult-groupNonInfect, 
	groupAdult-groupInfect, groupInfect-groupNonInfect, levels=d)

# -----------------------------------
# (A) TRANSITION: ADULT vs. INFECTIVE
# -----------------------------------

# DGE with filtering BEFORE normalization:
# ----------------------------------------

# ---- CPM >= 0 | in >= 1 sample ----
fitLim.AvsI.filtFirst.cpm0.sp1.dge <- lmFit(v.AvsI.filtFirst.cpm0.sp1.dge, d)
fit.AvsI.filtFirst.cpm0.sp1.dge <- contrasts.fit(fitLim.AvsI.filtFirst.cpm0.sp1.dge, contrast.matrix)
fit.AvsI.filtFirst.cpm0.sp1.dge <- eBayes(fitLim.AvsI.filtFirst.cpm0.sp1.dge)
AvsI.filtFirst.cpm0.sp1.DET <- topTable(fit.AvsI.filtFirst.cpm0.sp1.dge, 
	coef=2, n=dim(fit.AvsI.filtFirst.cpm0.sp1.dge)[1]) # DET table
# ---- CPM >= 0 | in >= 2 samples ----
fitLim.AvsI.filtFirst.cpm0.sp2.dge <- lmFit(v.AvsI.filtFirst.cpm0.sp2.dge, d)
fit.AvsI.filtFirst.cpm0.sp2.dge <- contrasts.fit(fitLim.AvsI.filtFirst.cpm0.sp2.dge, contrast.matrix)
fit.AvsI.filtFirst.cpm0.sp2.dge <- eBayes(fitLim.AvsI.filtFirst.cpm0.sp2.dge)
AvsI.filtFirst.cpm0.sp2.DET <- topTable(fit.AvsI.filtFirst.cpm0.sp2.dge, 
	coef=2, n=dim(fit.AvsI.filtFirst.cpm0.sp2.dge)[1]) # DET table
# ---- CPM >= 0 | in >= 3 samples ----
fitLim.AvsI.filtFirst.cpm0.sp3.dge <- lmFit(v.AvsI.filtFirst.cpm0.sp3.dge, d)
fit.AvsI.filtFirst.cpm0.sp3.dge <- contrasts.fit(fitLim.AvsI.filtFirst.cpm0.sp3.dge, contrast.matrix)
fit.AvsI.filtFirst.cpm0.sp3.dge <- eBayes(fitLim.AvsI.filtFirst.cpm0.sp3.dge)
AvsI.filtFirst.cpm0.sp3.DET <- topTable(fit.AvsI.filtFirst.cpm0.sp3.dge, 
	coef=2, n=dim(fit.AvsI.filtFirst.cpm0.sp3.dge)[1]) # DET table

# ---- CPM >= 5 | in >= 1 sample ----
fitLim.AvsI.filtFirst.cpm05.sp1.dge <- lmFit(v.AvsI.filtFirst.cpm05.sp1.dge, d)
fit.AvsI.filtFirst.cpm05.sp1.dge <- contrasts.fit(fitLim.AvsI.filtFirst.cpm05.sp1.dge, contrast.matrix)
fit.AvsI.filtFirst.cpm05.sp1.dge <- eBayes(fitLim.AvsI.filtFirst.cpm05.sp1.dge)
AvsI.filtFirst.cpm05.sp1.DET <- topTable(fit.AvsI.filtFirst.cpm05.sp1.dge, 
	coef=2, n=dim(fit.AvsI.filtFirst.cpm05.sp1.dge)[1]) # DET table
# ---- CPM >= 5 | in >= 2 samples ----
fitLim.AvsI.filtFirst.cpm05.sp2.dge <- lmFit(v.AvsI.filtFirst.cpm05.sp2.dge, d)
fit.AvsI.filtFirst.cpm05.sp2.dge <- contrasts.fit(fitLim.AvsI.filtFirst.cpm05.sp2.dge, contrast.matrix)
fit.AvsI.filtFirst.cpm05.sp2.dge <- eBayes(fitLim.AvsI.filtFirst.cpm05.sp2.dge)
AvsI.filtFirst.cpm05.sp2.DET <- topTable(fit.AvsI.filtFirst.cpm05.sp2.dge, 
	coef=2, n=dim(fit.AvsI.filtFirst.cpm05.sp2.dge)[1]) # DET table
# ---- CPM >= 5 | in >= 3 samples ----
fitLim.AvsI.filtFirst.cpm05.sp3.dge <- lmFit(v.AvsI.filtFirst.cpm05.sp3.dge, d)
fit.AvsI.filtFirst.cpm05.sp3.dge <- contrasts.fit(fitLim.AvsI.filtFirst.cpm05.sp3.dge, contrast.matrix)
fit.AvsI.filtFirst.cpm05.sp3.dge <- eBayes(fitLim.AvsI.filtFirst.cpm05.sp3.dge)
AvsI.filtFirst.cpm05.sp3.DET <- topTable(fit.AvsI.filtFirst.cpm05.sp3.dge, 
	coef=2, n=dim(fit.AvsI.filtFirst.cpm05.sp3.dge)[1]) # DET table

# ---- CPM >= 10 | in >= 1 sample ----
fitLim.AvsI.filtFirst.cpm10.sp1.dge <- lmFit(v.AvsI.filtFirst.cpm10.sp1.dge, d)
fit.AvsI.filtFirst.cpm10.sp1.dge <- contrasts.fit(fitLim.AvsI.filtFirst.cpm10.sp1.dge, contrast.matrix)
fit.AvsI.filtFirst.cpm10.sp1.dge <- eBayes(fitLim.AvsI.filtFirst.cpm10.sp1.dge)
AvsI.filtFirst.cpm10.sp1.DET <- topTable(fit.AvsI.filtFirst.cpm10.sp1.dge, 
	coef=2, n=dim(fit.AvsI.filtFirst.cpm10.sp1.dge)[1]) # DET table
# ---- CPM >= 10 | in >= 2 samples ----
fitLim.AvsI.filtFirst.cpm10.sp2.dge <- lmFit(v.AvsI.filtFirst.cpm10.sp2.dge, d)
fit.AvsI.filtFirst.cpm10.sp2.dge <- contrasts.fit(fitLim.AvsI.filtFirst.cpm10.sp2.dge, contrast.matrix)
fit.AvsI.filtFirst.cpm10.sp2.dge <- eBayes(fitLim.AvsI.filtFirst.cpm10.sp2.dge)
AvsI.filtFirst.cpm10.sp2.DET <- topTable(fit.AvsI.filtFirst.cpm10.sp2.dge, 
	coef=2, n=dim(fit.AvsI.filtFirst.cpm10.sp2.dge)[1]) # DET table
# ---- CPM >= 10 | in >= 3 samples ----
fitLim.AvsI.filtFirst.cpm10.sp3.dge <- lmFit(v.AvsI.filtFirst.cpm10.sp3.dge, d)
fit.AvsI.filtFirst.cpm10.sp3.dge <- contrasts.fit(fitLim.AvsI.filtFirst.cpm10.sp3.dge, contrast.matrix)
fit.AvsI.filtFirst.cpm10.sp3.dge <- eBayes(fitLim.AvsI.filtFirst.cpm10.sp3.dge)
AvsI.filtFirst.cpm10.sp3.DET <- topTable(fit.AvsI.filtFirst.cpm10.sp3.dge, 
	coef=2, n=dim(fit.AvsI.filtFirst.cpm10.sp3.dge)[1]) # DET table

# ---- CPM >= 15 | in >= 1 sample ----
fitLim.AvsI.filtFirst.cpm15.sp1.dge <- lmFit(v.AvsI.filtFirst.cpm15.sp1.dge, d)
fit.AvsI.filtFirst.cpm15.sp1.dge <- contrasts.fit(fitLim.AvsI.filtFirst.cpm15.sp1.dge, contrast.matrix)
fit.AvsI.filtFirst.cpm15.sp1.dge <- eBayes(fitLim.AvsI.filtFirst.cpm15.sp1.dge)
AvsI.filtFirst.cpm15.sp1.DET <- topTable(fit.AvsI.filtFirst.cpm15.sp1.dge, 
	coef=2, n=dim(fit.AvsI.filtFirst.cpm15.sp1.dge)[1]) # DET table
# ---- CPM >= 15 | in >= 2 samples ----
fitLim.AvsI.filtFirst.cpm15.sp2.dge <- lmFit(v.AvsI.filtFirst.cpm15.sp2.dge, d)
fit.AvsI.filtFirst.cpm15.sp2.dge <- contrasts.fit(fitLim.AvsI.filtFirst.cpm15.sp2.dge, contrast.matrix)
fit.AvsI.filtFirst.cpm15.sp2.dge <- eBayes(fitLim.AvsI.filtFirst.cpm15.sp2.dge)
AvsI.filtFirst.cpm15.sp2.DET <- topTable(fit.AvsI.filtFirst.cpm15.sp2.dge, 
	coef=2, n=dim(fit.AvsI.filtFirst.cpm15.sp2.dge)[1]) # DET table
# ---- CPM >= 15 | in >= 3 samples ----
fitLim.AvsI.filtFirst.cpm15.sp3.dge <- lmFit(v.AvsI.filtFirst.cpm15.sp3.dge, d)
fit.AvsI.filtFirst.cpm15.sp3.dge <- contrasts.fit(fitLim.AvsI.filtFirst.cpm15.sp3.dge, contrast.matrix)
fit.AvsI.filtFirst.cpm15.sp3.dge <- eBayes(fitLim.AvsI.filtFirst.cpm15.sp3.dge)
AvsI.filtFirst.cpm15.sp3.DET <- topTable(fit.AvsI.filtFirst.cpm15.sp3.dge, 
	coef=2, n=dim(fit.AvsI.filtFirst.cpm15.sp3.dge)[1]) # DET table

# DGE with filtering AFTER normalization:
# ----------------------------------------

# ---- CPM >= 0 | in >= 1 sample ----
fitLim.AvsI.filtAfter.cpm0.sp1.dge <- lmFit(v.AvsI.filtAfter.cpm0.sp1.dge, d)
fit.AvsI.filtAfter.cpm0.sp1.dge <- contrasts.fit(fitLim.AvsI.filtAfter.cpm0.sp1.dge, contrast.matrix)
fit.AvsI.filtAfter.cpm0.sp1.dge <- eBayes(fitLim.AvsI.filtAfter.cpm0.sp1.dge)
AvsI.filtAfter.cpm0.sp1.DET <- topTable(fit.AvsI.filtAfter.cpm0.sp1.dge, 
	coef=2, n=dim(fit.AvsI.filtAfter.cpm0.sp1.dge)[1]) # DET table
# ---- CPM >= 0 | in >= 2 samples ----
fitLim.AvsI.filtAfter.cpm0.sp2.dge <- lmFit(v.AvsI.filtAfter.cpm0.sp2.dge, d)
fit.AvsI.filtAfter.cpm0.sp2.dge <- contrasts.fit(fitLim.AvsI.filtAfter.cpm0.sp2.dge, contrast.matrix)
fit.AvsI.filtAfter.cpm0.sp2.dge <- eBayes(fitLim.AvsI.filtAfter.cpm0.sp2.dge)
AvsI.filtAfter.cpm0.sp2.DET <- topTable(fit.AvsI.filtAfter.cpm0.sp2.dge, 
	coef=2, n=dim(fit.AvsI.filtAfter.cpm0.sp2.dge)[1]) # DET table
# ---- CPM >= 0 | in >= 3 samples ----
fitLim.AvsI.filtAfter.cpm0.sp3.dge <- lmFit(v.AvsI.filtAfter.cpm0.sp3.dge, d)
fit.AvsI.filtAfter.cpm0.sp3.dge <- contrasts.fit(fitLim.AvsI.filtAfter.cpm0.sp3.dge, contrast.matrix)
fit.AvsI.filtAfter.cpm0.sp3.dge <- eBayes(fitLim.AvsI.filtAfter.cpm0.sp3.dge)
AvsI.filtAfter.cpm0.sp3.DET <- topTable(fit.AvsI.filtAfter.cpm0.sp3.dge, 
	coef=2, n=dim(fit.AvsI.filtAfter.cpm0.sp3.dge)[1]) # DET table

# ---- CPM >= 5 | in >= 1 sample ----
fitLim.AvsI.filtAfter.cpm05.sp1.dge <- lmFit(v.AvsI.filtAfter.cpm05.sp1.dge, d)
fit.AvsI.filtAfter.cpm05.sp1.dge <- contrasts.fit(fitLim.AvsI.filtAfter.cpm05.sp1.dge, contrast.matrix)
fit.AvsI.filtAfter.cpm05.sp1.dge <- eBayes(fitLim.AvsI.filtAfter.cpm05.sp1.dge)
AvsI.filtAfter.cpm05.sp1.DET <- topTable(fit.AvsI.filtAfter.cpm05.sp1.dge, 
	coef=2, n=dim(fit.AvsI.filtAfter.cpm05.sp1.dge)[1]) # DET table
# ---- CPM >= 5 | in >= 2 samples ----
fitLim.AvsI.filtAfter.cpm05.sp2.dge <- lmFit(v.AvsI.filtAfter.cpm05.sp2.dge, d)
fit.AvsI.filtAfter.cpm05.sp2.dge <- contrasts.fit(fitLim.AvsI.filtAfter.cpm05.sp2.dge, contrast.matrix)
fit.AvsI.filtAfter.cpm05.sp2.dge <- eBayes(fitLim.AvsI.filtAfter.cpm05.sp2.dge)
AvsI.filtAfter.cpm05.sp2.DET <- topTable(fit.AvsI.filtAfter.cpm05.sp2.dge, 
	coef=2, n=dim(fit.AvsI.filtAfter.cpm05.sp2.dge)[1]) # DET table
# ---- CPM >= 5 | in >= 3 samples ----
fitLim.AvsI.filtAfter.cpm05.sp3.dge <- lmFit(v.AvsI.filtAfter.cpm05.sp3.dge, d)
fit.AvsI.filtAfter.cpm05.sp3.dge <- contrasts.fit(fitLim.AvsI.filtAfter.cpm05.sp3.dge, contrast.matrix)
fit.AvsI.filtAfter.cpm05.sp3.dge <- eBayes(fitLim.AvsI.filtAfter.cpm05.sp3.dge)
AvsI.filtAfter.cpm05.sp3.DET <- topTable(fit.AvsI.filtAfter.cpm05.sp3.dge, 
	coef=2, n=dim(fit.AvsI.filtAfter.cpm05.sp3.dge)[1]) # DET table

# ---- CPM >= 10 | in >= 1 sample ----
fitLim.AvsI.filtAfter.cpm10.sp1.dge <- lmFit(v.AvsI.filtAfter.cpm10.sp1.dge, d)
fit.AvsI.filtAfter.cpm10.sp1.dge <- contrasts.fit(fitLim.AvsI.filtAfter.cpm10.sp1.dge, contrast.matrix)
fit.AvsI.filtAfter.cpm10.sp1.dge <- eBayes(fitLim.AvsI.filtAfter.cpm10.sp1.dge)
AvsI.filtAfter.cpm10.sp1.DET <- topTable(fit.AvsI.filtAfter.cpm10.sp1.dge, 
	coef=2, n=dim(fit.AvsI.filtAfter.cpm10.sp1.dge)[1]) # DET table
# ---- CPM >= 10 | in >= 2 samples ----
fitLim.AvsI.filtAfter.cpm10.sp2.dge <- lmFit(v.AvsI.filtAfter.cpm10.sp2.dge, d)
fit.AvsI.filtAfter.cpm10.sp2.dge <- contrasts.fit(fitLim.AvsI.filtAfter.cpm10.sp2.dge, contrast.matrix)
fit.AvsI.filtAfter.cpm10.sp2.dge <- eBayes(fitLim.AvsI.filtAfter.cpm10.sp2.dge)
AvsI.filtAfter.cpm10.sp2.DET <- topTable(fit.AvsI.filtAfter.cpm10.sp2.dge, 
	coef=2, n=dim(fit.AvsI.filtAfter.cpm10.sp2.dge)[1]) # DET table
# ---- CPM >= 10 | in >= 3 samples ----
fitLim.AvsI.filtAfter.cpm10.sp3.dge <- lmFit(v.AvsI.filtAfter.cpm10.sp3.dge, d)
fit.AvsI.filtAfter.cpm10.sp3.dge <- contrasts.fit(fitLim.AvsI.filtAfter.cpm10.sp3.dge, contrast.matrix)
fit.AvsI.filtAfter.cpm10.sp3.dge <- eBayes(fitLim.AvsI.filtAfter.cpm10.sp3.dge)
AvsI.filtAfter.cpm10.sp3.DET <- topTable(fit.AvsI.filtAfter.cpm10.sp3.dge, 
	coef=2, n=dim(fit.AvsI.filtAfter.cpm10.sp3.dge)[1]) # DET table

# ---- CPM >= 15 | in >= 1 sample ----
fitLim.AvsI.filtAfter.cpm15.sp1.dge <- lmFit(v.AvsI.filtAfter.cpm15.sp1.dge, d)
fit.AvsI.filtAfter.cpm15.sp1.dge <- contrasts.fit(fitLim.AvsI.filtAfter.cpm15.sp1.dge, contrast.matrix)
fit.AvsI.filtAfter.cpm15.sp1.dge <- eBayes(fitLim.AvsI.filtAfter.cpm15.sp1.dge)
AvsI.filtAfter.cpm15.sp1.DET <- topTable(fit.AvsI.filtAfter.cpm15.sp1.dge, 
	coef=2, n=dim(fit.AvsI.filtAfter.cpm15.sp1.dge)[1]) # DET table
# ---- CPM >= 15 | in >= 2 samples ----
fitLim.AvsI.filtAfter.cpm15.sp2.dge <- lmFit(v.AvsI.filtAfter.cpm15.sp2.dge, d)
fit.AvsI.filtAfter.cpm15.sp2.dge <- contrasts.fit(fitLim.AvsI.filtAfter.cpm15.sp2.dge, contrast.matrix)
fit.AvsI.filtAfter.cpm15.sp2.dge <- eBayes(fitLim.AvsI.filtAfter.cpm15.sp2.dge)
AvsI.filtAfter.cpm15.sp2.DET <- topTable(fit.AvsI.filtAfter.cpm15.sp2.dge, 
	coef=2, n=dim(fit.AvsI.filtAfter.cpm15.sp2.dge)[1]) # DET table
# ---- CPM >= 15 | in >= 3 samples ----
fitLim.AvsI.filtAfter.cpm15.sp3.dge <- lmFit(v.AvsI.filtAfter.cpm15.sp3.dge, d)
fit.AvsI.filtAfter.cpm15.sp3.dge <- contrasts.fit(fitLim.AvsI.filtAfter.cpm15.sp3.dge, contrast.matrix)
fit.AvsI.filtAfter.cpm15.sp3.dge <- eBayes(fitLim.AvsI.filtAfter.cpm15.sp3.dge)
AvsI.filtAfter.cpm15.sp3.DET <- topTable(fit.AvsI.filtAfter.cpm15.sp3.dge, 
	coef=2, n=dim(fit.AvsI.filtAfter.cpm15.sp3.dge)[1]) # DET table

# -------------------------------------------
# (B) TRANSITION: INFECTIVE vs. NON-INFECTIVE
# -------------------------------------------

# DGE with filtering BEFORE normalization:
# ----------------------------------------

# ---- CPM >= 0 | in >= 1 sample ----
fitLim.IvsNI.filtFirst.cpm0.sp1.dge <- lmFit(v.IvsNI.filtFirst.cpm0.sp1.dge, d)
fit.IvsNI.filtFirst.cpm0.sp1.dge <- contrasts.fit(fitLim.IvsNI.filtFirst.cpm0.sp1.dge, contrast.matrix)
fit.IvsNI.filtFirst.cpm0.sp1.dge <- eBayes(fitLim.IvsNI.filtFirst.cpm0.sp1.dge)
IvsNI.filtFirst.cpm0.sp1.DET <- topTable(fit.IvsNI.filtFirst.cpm0.sp1.dge, 
	coef=3, n=dim(fit.IvsNI.filtFirst.cpm0.sp1.dge)[1]) # DET table
# ---- CPM >= 0 | in >= 2 samples ----
fitLim.IvsNI.filtFirst.cpm0.sp2.dge <- lmFit(v.IvsNI.filtFirst.cpm0.sp2.dge, d)
fit.IvsNI.filtFirst.cpm0.sp2.dge <- contrasts.fit(fitLim.IvsNI.filtFirst.cpm0.sp2.dge, contrast.matrix)
fit.IvsNI.filtFirst.cpm0.sp2.dge <- eBayes(fitLim.IvsNI.filtFirst.cpm0.sp2.dge)
IvsNI.filtFirst.cpm0.sp2.DET <- topTable(fit.IvsNI.filtFirst.cpm0.sp2.dge, 
	coef=3, n=dim(fit.IvsNI.filtFirst.cpm0.sp2.dge)[1]) # DET table
# ---- CPM >= 0 | in >= 3 samples ----
fitLim.IvsNI.filtFirst.cpm0.sp3.dge <- lmFit(v.IvsNI.filtFirst.cpm0.sp3.dge, d)
fit.IvsNI.filtFirst.cpm0.sp3.dge <- contrasts.fit(fitLim.IvsNI.filtFirst.cpm0.sp3.dge, contrast.matrix)
fit.IvsNI.filtFirst.cpm0.sp3.dge <- eBayes(fitLim.IvsNI.filtFirst.cpm0.sp3.dge)
IvsNI.filtFirst.cpm0.sp3.DET <- topTable(fit.IvsNI.filtFirst.cpm0.sp3.dge, 
	coef=3, n=dim(fit.IvsNI.filtFirst.cpm0.sp3.dge)[1]) # DET table

# ---- CPM >= 5 | in >= 1 sample ----
fitLim.IvsNI.filtFirst.cpm05.sp1.dge <- lmFit(v.IvsNI.filtFirst.cpm05.sp1.dge, d)
fit.IvsNI.filtFirst.cpm05.sp1.dge <- contrasts.fit(fitLim.IvsNI.filtFirst.cpm05.sp1.dge, contrast.matrix)
fit.IvsNI.filtFirst.cpm05.sp1.dge <- eBayes(fitLim.IvsNI.filtFirst.cpm05.sp1.dge)
IvsNI.filtFirst.cpm05.sp1.DET <- topTable(fit.IvsNI.filtFirst.cpm05.sp1.dge, 
	coef=3, n=dim(fit.IvsNI.filtFirst.cpm05.sp1.dge)[1]) # DET table
# ---- CPM >= 5 | in >= 2 samples ----
fitLim.IvsNI.filtFirst.cpm05.sp2.dge <- lmFit(v.IvsNI.filtFirst.cpm05.sp2.dge, d)
fit.IvsNI.filtFirst.cpm05.sp2.dge <- contrasts.fit(fitLim.IvsNI.filtFirst.cpm05.sp2.dge, contrast.matrix)
fit.IvsNI.filtFirst.cpm05.sp2.dge <- eBayes(fitLim.IvsNI.filtFirst.cpm05.sp2.dge)
IvsNI.filtFirst.cpm05.sp2.DET <- topTable(fit.IvsNI.filtFirst.cpm05.sp2.dge, 
	coef=3, n=dim(fit.IvsNI.filtFirst.cpm05.sp2.dge)[1]) # DET table
# ---- CPM >= 5 | in >= 3 samples ----
fitLim.IvsNI.filtFirst.cpm05.sp3.dge <- lmFit(v.IvsNI.filtFirst.cpm05.sp3.dge, d)
fit.IvsNI.filtFirst.cpm05.sp3.dge <- contrasts.fit(fitLim.IvsNI.filtFirst.cpm05.sp3.dge, contrast.matrix)
fit.IvsNI.filtFirst.cpm05.sp3.dge <- eBayes(fitLim.IvsNI.filtFirst.cpm05.sp3.dge)
IvsNI.filtFirst.cpm05.sp3.DET <- topTable(fit.IvsNI.filtFirst.cpm05.sp3.dge, 
	coef=3, n=dim(fit.IvsNI.filtFirst.cpm05.sp3.dge)[1]) # DET table

# ---- CPM >= 10 | in >= 1 sample ----
fitLim.IvsNI.filtFirst.cpm10.sp1.dge <- lmFit(v.IvsNI.filtFirst.cpm10.sp1.dge, d)
fit.IvsNI.filtFirst.cpm10.sp1.dge <- contrasts.fit(fitLim.IvsNI.filtFirst.cpm10.sp1.dge, contrast.matrix)
fit.IvsNI.filtFirst.cpm10.sp1.dge <- eBayes(fitLim.IvsNI.filtFirst.cpm10.sp1.dge)
IvsNI.filtFirst.cpm10.sp1.DET <- topTable(fit.IvsNI.filtFirst.cpm10.sp1.dge, 
	coef=3, n=dim(fit.IvsNI.filtFirst.cpm10.sp1.dge)[1]) # DET table
# ---- CPM >= 10 | in >= 2 samples ----
fitLim.IvsNI.filtFirst.cpm10.sp2.dge <- lmFit(v.IvsNI.filtFirst.cpm10.sp2.dge, d)
fit.IvsNI.filtFirst.cpm10.sp2.dge <- contrasts.fit(fitLim.IvsNI.filtFirst.cpm10.sp2.dge, contrast.matrix)
fit.IvsNI.filtFirst.cpm10.sp2.dge <- eBayes(fitLim.IvsNI.filtFirst.cpm10.sp2.dge)
IvsNI.filtFirst.cpm10.sp2.DET <- topTable(fit.IvsNI.filtFirst.cpm10.sp2.dge, 
	coef=3, n=dim(fit.IvsNI.filtFirst.cpm10.sp2.dge)[1]) # DET table
# ---- CPM >= 10 | in >= 3 samples ----
fitLim.IvsNI.filtFirst.cpm10.sp3.dge <- lmFit(v.IvsNI.filtFirst.cpm10.sp3.dge, d)
fit.IvsNI.filtFirst.cpm10.sp3.dge <- contrasts.fit(fitLim.IvsNI.filtFirst.cpm10.sp3.dge, contrast.matrix)
fit.IvsNI.filtFirst.cpm10.sp3.dge <- eBayes(fitLim.IvsNI.filtFirst.cpm10.sp3.dge)
IvsNI.filtFirst.cpm10.sp3.DET <- topTable(fit.IvsNI.filtFirst.cpm10.sp3.dge, 
	coef=3, n=dim(fit.IvsNI.filtFirst.cpm10.sp3.dge)[1]) # DET table

# ---- CPM >= 15 | in >= 1 sample ----
fitLim.IvsNI.filtFirst.cpm15.sp1.dge <- lmFit(v.IvsNI.filtFirst.cpm15.sp1.dge, d)
fit.IvsNI.filtFirst.cpm15.sp1.dge <- contrasts.fit(fitLim.IvsNI.filtFirst.cpm15.sp1.dge, contrast.matrix)
fit.IvsNI.filtFirst.cpm15.sp1.dge <- eBayes(fitLim.IvsNI.filtFirst.cpm15.sp1.dge)
IvsNI.filtFirst.cpm15.sp1.DET <- topTable(fit.IvsNI.filtFirst.cpm15.sp1.dge, 
	coef=3, n=dim(fit.IvsNI.filtFirst.cpm15.sp1.dge)[1]) # DET table
# ---- CPM >= 15 | in >= 2 samples ----
fitLim.IvsNI.filtFirst.cpm15.sp2.dge <- lmFit(v.IvsNI.filtFirst.cpm15.sp2.dge, d)
fit.IvsNI.filtFirst.cpm15.sp2.dge <- contrasts.fit(fitLim.IvsNI.filtFirst.cpm15.sp2.dge, contrast.matrix)
fit.IvsNI.filtFirst.cpm15.sp2.dge <- eBayes(fitLim.IvsNI.filtFirst.cpm15.sp2.dge)
IvsNI.filtFirst.cpm15.sp2.DET <- topTable(fit.IvsNI.filtFirst.cpm15.sp2.dge, 
	coef=3, n=dim(fit.IvsNI.filtFirst.cpm15.sp2.dge)[1]) # DET table
# ---- CPM >= 15 | in >= 3 samples ----
fitLim.IvsNI.filtFirst.cpm15.sp3.dge <- lmFit(v.IvsNI.filtFirst.cpm15.sp3.dge, d)
fit.IvsNI.filtFirst.cpm15.sp3.dge <- contrasts.fit(fitLim.IvsNI.filtFirst.cpm15.sp3.dge, contrast.matrix)
fit.IvsNI.filtFirst.cpm15.sp3.dge <- eBayes(fitLim.IvsNI.filtFirst.cpm15.sp3.dge)
IvsNI.filtFirst.cpm15.sp3.DET <- topTable(fit.IvsNI.filtFirst.cpm15.sp3.dge, 
	coef=3, n=dim(fit.IvsNI.filtFirst.cpm15.sp3.dge)[1]) # DET table

# DGE with filtering AFTER normalization:
# ----------------------------------------

# ---- CPM >= 0 | in >= 1 sample ----
fitLim.IvsNI.filtAfter.cpm0.sp1.dge <- lmFit(v.IvsNI.filtAfter.cpm0.sp1.dge, d)
fit.IvsNI.filtAfter.cpm0.sp1.dge <- contrasts.fit(fitLim.IvsNI.filtAfter.cpm0.sp1.dge, contrast.matrix)
fit.IvsNI.filtAfter.cpm0.sp1.dge <- eBayes(fitLim.IvsNI.filtAfter.cpm0.sp1.dge)
IvsNI.filtAfter.cpm0.sp1.DET <- topTable(fit.IvsNI.filtAfter.cpm0.sp1.dge, 
	coef=3, n=dim(fit.IvsNI.filtAfter.cpm0.sp1.dge)[1]) # DET table
# ---- CPM >= 0 | in >= 2 samples ----
fitLim.IvsNI.filtAfter.cpm0.sp2.dge <- lmFit(v.IvsNI.filtAfter.cpm0.sp2.dge, d)
fit.IvsNI.filtAfter.cpm0.sp2.dge <- contrasts.fit(fitLim.IvsNI.filtAfter.cpm0.sp2.dge, contrast.matrix)
fit.IvsNI.filtAfter.cpm0.sp2.dge <- eBayes(fitLim.IvsNI.filtAfter.cpm0.sp2.dge)
IvsNI.filtAfter.cpm0.sp2.DET <- topTable(fit.IvsNI.filtAfter.cpm0.sp2.dge, 
	coef=3, n=dim(fit.IvsNI.filtAfter.cpm0.sp2.dge)[1]) # DET table
# ---- CPM >= 0 | in >= 3 samples ----
fitLim.IvsNI.filtAfter.cpm0.sp3.dge <- lmFit(v.IvsNI.filtAfter.cpm0.sp3.dge, d)
fit.IvsNI.filtAfter.cpm0.sp3.dge <- contrasts.fit(fitLim.IvsNI.filtAfter.cpm0.sp3.dge, contrast.matrix)
fit.IvsNI.filtAfter.cpm0.sp3.dge <- eBayes(fitLim.IvsNI.filtAfter.cpm0.sp3.dge)
IvsNI.filtAfter.cpm0.sp3.DET <- topTable(fit.IvsNI.filtAfter.cpm0.sp3.dge, 
	coef=3, n=dim(fit.IvsNI.filtAfter.cpm0.sp3.dge)[1]) # DET table

# ---- CPM >= 5 | in >= 1 sample ----
fitLim.IvsNI.filtAfter.cpm05.sp1.dge <- lmFit(v.IvsNI.filtAfter.cpm05.sp1.dge, d)
fit.IvsNI.filtAfter.cpm05.sp1.dge <- contrasts.fit(fitLim.IvsNI.filtAfter.cpm05.sp1.dge, contrast.matrix)
fit.IvsNI.filtAfter.cpm05.sp1.dge <- eBayes(fitLim.IvsNI.filtAfter.cpm05.sp1.dge)
IvsNI.filtAfter.cpm05.sp1.DET <- topTable(fit.IvsNI.filtAfter.cpm05.sp1.dge, 
	coef=3, n=dim(fit.IvsNI.filtAfter.cpm05.sp1.dge)[1]) # DET table
# ---- CPM >= 5 | in >= 2 samples ----
fitLim.IvsNI.filtAfter.cpm05.sp2.dge <- lmFit(v.IvsNI.filtAfter.cpm05.sp2.dge, d)
fit.IvsNI.filtAfter.cpm05.sp2.dge <- contrasts.fit(fitLim.IvsNI.filtAfter.cpm05.sp2.dge, contrast.matrix)
fit.IvsNI.filtAfter.cpm05.sp2.dge <- eBayes(fitLim.IvsNI.filtAfter.cpm05.sp2.dge)
IvsNI.filtAfter.cpm05.sp2.DET <- topTable(fit.IvsNI.filtAfter.cpm05.sp2.dge, 
	coef=3, n=dim(fit.IvsNI.filtAfter.cpm05.sp2.dge)[1]) # DET table
# ---- CPM >= 5 | in >= 3 samples ----
fitLim.IvsNI.filtAfter.cpm05.sp3.dge <- lmFit(v.IvsNI.filtAfter.cpm05.sp3.dge, d)
fit.IvsNI.filtAfter.cpm05.sp3.dge <- contrasts.fit(fitLim.IvsNI.filtAfter.cpm05.sp3.dge, contrast.matrix)
fit.IvsNI.filtAfter.cpm05.sp3.dge <- eBayes(fitLim.IvsNI.filtAfter.cpm05.sp3.dge)
IvsNI.filtAfter.cpm05.sp3.DET <- topTable(fit.IvsNI.filtAfter.cpm05.sp3.dge, 
	coef=3, n=dim(fit.IvsNI.filtAfter.cpm05.sp3.dge)[1]) # DET table

# ---- CPM >= 10 | in >= 1 sample ----
fitLim.IvsNI.filtAfter.cpm10.sp1.dge <- lmFit(v.IvsNI.filtAfter.cpm10.sp1.dge, d)
fit.IvsNI.filtAfter.cpm10.sp1.dge <- contrasts.fit(fitLim.IvsNI.filtAfter.cpm10.sp1.dge, contrast.matrix)
fit.IvsNI.filtAfter.cpm10.sp1.dge <- eBayes(fitLim.IvsNI.filtAfter.cpm10.sp1.dge)
IvsNI.filtAfter.cpm10.sp1.DET <- topTable(fit.IvsNI.filtAfter.cpm10.sp1.dge, 
	coef=3, n=dim(fit.IvsNI.filtAfter.cpm10.sp1.dge)[1]) # DET table
# ---- CPM >= 10 | in >= 2 samples ----
fitLim.IvsNI.filtAfter.cpm10.sp2.dge <- lmFit(v.IvsNI.filtAfter.cpm10.sp2.dge, d)
fit.IvsNI.filtAfter.cpm10.sp2.dge <- contrasts.fit(fitLim.IvsNI.filtAfter.cpm10.sp2.dge, contrast.matrix)
fit.IvsNI.filtAfter.cpm10.sp2.dge <- eBayes(fitLim.IvsNI.filtAfter.cpm10.sp2.dge)
IvsNI.filtAfter.cpm10.sp2.DET <- topTable(fit.IvsNI.filtAfter.cpm10.sp2.dge, 
	coef=3, n=dim(fit.IvsNI.filtAfter.cpm10.sp2.dge)[1]) # DET table
# ---- CPM >= 10 | in >= 3 samples ----
fitLim.IvsNI.filtAfter.cpm10.sp3.dge <- lmFit(v.IvsNI.filtAfter.cpm10.sp3.dge, d)
fit.IvsNI.filtAfter.cpm10.sp3.dge <- contrasts.fit(fitLim.IvsNI.filtAfter.cpm10.sp3.dge, contrast.matrix)
fit.IvsNI.filtAfter.cpm10.sp3.dge <- eBayes(fitLim.IvsNI.filtAfter.cpm10.sp3.dge)
IvsNI.filtAfter.cpm10.sp3.DET <- topTable(fit.IvsNI.filtAfter.cpm10.sp3.dge, 
	coef=3, n=dim(fit.IvsNI.filtAfter.cpm10.sp3.dge)[1]) # DET table

# ---- CPM >= 15 | in >= 1 sample ----
fitLim.IvsNI.filtAfter.cpm15.sp1.dge <- lmFit(v.IvsNI.filtAfter.cpm15.sp1.dge, d)
fit.IvsNI.filtAfter.cpm15.sp1.dge <- contrasts.fit(fitLim.IvsNI.filtAfter.cpm15.sp1.dge, contrast.matrix)
fit.IvsNI.filtAfter.cpm15.sp1.dge <- eBayes(fitLim.IvsNI.filtAfter.cpm15.sp1.dge)
IvsNI.filtAfter.cpm15.sp1.DET <- topTable(fit.IvsNI.filtAfter.cpm15.sp1.dge, 
	coef=3, n=dim(fit.IvsNI.filtAfter.cpm15.sp1.dge)[1]) # DET table
# ---- CPM >= 15 | in >= 2 samples ----
fitLim.IvsNI.filtAfter.cpm15.sp2.dge <- lmFit(v.IvsNI.filtAfter.cpm15.sp2.dge, d)
fit.IvsNI.filtAfter.cpm15.sp2.dge <- contrasts.fit(fitLim.IvsNI.filtAfter.cpm15.sp2.dge, contrast.matrix)
fit.IvsNI.filtAfter.cpm15.sp2.dge <- eBayes(fitLim.IvsNI.filtAfter.cpm15.sp2.dge)
IvsNI.filtAfter.cpm15.sp2.DET <- topTable(fit.IvsNI.filtAfter.cpm15.sp2.dge, 
	coef=3, n=dim(fit.IvsNI.filtAfter.cpm15.sp2.dge)[1]) # DET table
# ---- CPM >= 15 | in >= 3 samples ----
fitLim.IvsNI.filtAfter.cpm15.sp3.dge <- lmFit(v.IvsNI.filtAfter.cpm15.sp3.dge, d)
fit.IvsNI.filtAfter.cpm15.sp3.dge <- contrasts.fit(fitLim.IvsNI.filtAfter.cpm15.sp3.dge, contrast.matrix)
fit.IvsNI.filtAfter.cpm15.sp3.dge <- eBayes(fitLim.IvsNI.filtAfter.cpm15.sp3.dge)
IvsNI.filtAfter.cpm15.sp3.DET <- topTable(fit.IvsNI.filtAfter.cpm15.sp3.dge, 
	coef=3, n=dim(fit.IvsNI.filtAfter.cpm15.sp3.dge)[1]) # DET table

################################################################
# (6) EXTRACTING DIFF. EXPR. TRANSCRIPTS WITH ADJ. P-VAL < 0.001
# --------------------------------------------------------------
################################################################

# TRANSITION: ADULT vs. INFECTIVE
# -------------------------------

# NB: each object created here contains a subset of the original
# 'topTable' object from the edgeR/limma packages, but only for
# for the transcripts with adj.P.Val < 0.001. So if the names of
# the transcripts is the required information, then one must
# use the 'row.names()' function to get them and do something
# with them (e.g. finding the intercept, i.e. common names).

# ---- AvsI, filtFirst | CPM >= 0 | in >= 1 sample ---- #
AvsI.filtFirst.subset.cpm0.sp1 <- subset(x = AvsI.filtFirst.cpm0.sp1.DET,
	subset = AvsI.filtFirst.cpm0.sp1.DET$adj.P.Val < 0.001)
# ---- AvsI, filtFirst | CPM >= 0 | in >= 2 samples ---- #
AvsI.filtFirst.subset.cpm0.sp2 <- subset(x = AvsI.filtFirst.cpm0.sp2.DET,
	subset = AvsI.filtFirst.cpm0.sp2.DET$adj.P.Val < 0.001)
# ---- AvsI, filtFirst | CPM >= 0 | in >= 3 samples ---- #
AvsI.filtFirst.subset.cpm0.sp3 <- subset(x = AvsI.filtFirst.cpm0.sp3.DET,
	subset = AvsI.filtFirst.cpm0.sp3.DET$adj.P.Val < 0.001)

# ---- AvsI, filtFirst | CPM >= 5 | in >= 1 sample ---- #
AvsI.filtFirst.subset.cpm05.sp1 <- subset(x = AvsI.filtFirst.cpm05.sp1.DET,
	subset = AvsI.filtFirst.cpm05.sp1.DET$adj.P.Val < 0.001)
# ---- AvsI, filtFirst | CPM >= 5 | in >= 2 samples ---- #
AvsI.filtFirst.subset.cpm05.sp2 <- subset(x = AvsI.filtFirst.cpm05.sp2.DET,
	subset = AvsI.filtFirst.cpm05.sp2.DET$adj.P.Val < 0.001)
# ---- AvsI, filtFirst | CPM >= 5 | in >= 3 samples ---- #
AvsI.filtFirst.subset.cpm05.sp3 <- subset(x = AvsI.filtFirst.cpm05.sp3.DET,
	subset = AvsI.filtFirst.cpm05.sp3.DET$adj.P.Val < 0.001)

# ---- AvsI, filtFirst | CPM >= 10 | in >= 1 sample ---- #
AvsI.filtFirst.subset.cpm10.sp1 <- subset(x = AvsI.filtFirst.cpm10.sp1.DET,
	subset = AvsI.filtFirst.cpm10.sp1.DET$adj.P.Val < 0.001)
# ---- AvsI, filtFirst | CPM >= 10 | in >= 2 samples ---- #
AvsI.filtFirst.subset.cpm10.sp2 <- subset(x = AvsI.filtFirst.cpm10.sp2.DET,
	subset = AvsI.filtFirst.cpm10.sp2.DET$adj.P.Val < 0.001)
# ---- AvsI, filtFirst | CPM >= 10 | in >= 3 samples ---- #
AvsI.filtFirst.subset.cpm10.sp3 <- subset(x = AvsI.filtFirst.cpm10.sp3.DET,
	subset = AvsI.filtFirst.cpm10.sp3.DET$adj.P.Val < 0.001)

# ---- AvsI, filtFirst | CPM >= 15 | in >= 1 sample ---- #
AvsI.filtFirst.subset.cpm15.sp1 <- subset(x = AvsI.filtFirst.cpm15.sp1.DET,
	subset = AvsI.filtFirst.cpm15.sp1.DET$adj.P.Val < 0.001)
# ---- AvsI, filtFirst | CPM >= 15 | in >= 2 samples ---- #
AvsI.filtFirst.subset.cpm15.sp2 <- subset(x = AvsI.filtFirst.cpm15.sp2.DET,
	subset = AvsI.filtFirst.cpm15.sp2.DET$adj.P.Val < 0.001)
# ---- AvsI, filtFirst | CPM >= 15 | in >= 3 samples ---- #
AvsI.filtFirst.subset.cpm15.sp3 <- subset(x = AvsI.filtFirst.cpm15.sp3.DET,
	subset = AvsI.filtFirst.cpm15.sp3.DET$adj.P.Val < 0.001)

# ----------------------------------------------------------

# ---- AvsI, filtAfter | CPM >= 0 | in >= 1 sample ---- #
AvsI.filtAfter.subset.cpm0.sp1 <- subset(x = AvsI.filtAfter.cpm0.sp1.DET,
	subset = AvsI.filtAfter.cpm0.sp1.DET$adj.P.Val < 0.001)
# ---- AvsI, filtAfter | CPM >= 0 | in >= 2 samples ---- #
AvsI.filtAfter.subset.cpm0.sp2 <- subset(x = AvsI.filtAfter.cpm0.sp2.DET,
	subset = AvsI.filtAfter.cpm0.sp2.DET$adj.P.Val < 0.001)
# ---- AvsI, filtAfter | CPM >= 0 | in >= 3 samples ---- #
AvsI.filtAfter.subset.cpm0.sp3 <- subset(x = AvsI.filtAfter.cpm0.sp3.DET,
	subset = AvsI.filtAfter.cpm0.sp3.DET$adj.P.Val < 0.001)

# ---- AvsI, filtAfter | CPM >= 5 | in >= 1 sample ---- #
AvsI.filtAfter.subset.cpm05.sp1 <- subset(x = AvsI.filtAfter.cpm05.sp1.DET,
	subset = AvsI.filtAfter.cpm05.sp1.DET$adj.P.Val < 0.001)
# ---- AvsI, filtAfter | CPM >= 5 | in >= 2 samples ---- #
AvsI.filtAfter.subset.cpm05.sp2 <- subset(x = AvsI.filtAfter.cpm05.sp2.DET,
	subset = AvsI.filtAfter.cpm05.sp2.DET$adj.P.Val < 0.001)
# ---- AvsI, filtAfter | CPM >= 5 | in >= 3 samples ---- #
AvsI.filtAfter.subset.cpm05.sp3 <- subset(x = AvsI.filtAfter.cpm05.sp3.DET,
	subset = AvsI.filtAfter.cpm05.sp3.DET$adj.P.Val < 0.001)

# ---- AvsI, filtAfter | CPM >= 10 | in >= 1 sample ---- #
AvsI.filtAfter.subset.cpm10.sp1 <- subset(x = AvsI.filtAfter.cpm10.sp1.DET,
	subset = AvsI.filtAfter.cpm10.sp1.DET$adj.P.Val < 0.001)
# ---- AvsI, filtAfter | CPM >= 10 | in >= 2 samples ---- #
AvsI.filtAfter.subset.cpm10.sp2 <- subset(x = AvsI.filtAfter.cpm10.sp2.DET,
	subset = AvsI.filtAfter.cpm10.sp2.DET$adj.P.Val < 0.001)
# ---- AvsI, filtAfter | CPM >= 10 | in >= 3 samples ---- #
AvsI.filtAfter.subset.cpm10.sp3 <- subset(x = AvsI.filtAfter.cpm10.sp3.DET,
	subset = AvsI.filtAfter.cpm10.sp3.DET$adj.P.Val < 0.001)

# ---- AvsI, filtAfter | CPM >= 15 | in >= 1 sample ---- #
AvsI.filtAfter.subset.cpm15.sp1 <- subset(x = AvsI.filtAfter.cpm15.sp1.DET,
	subset = AvsI.filtAfter.cpm15.sp1.DET$adj.P.Val < 0.001)
# ---- AvsI, filtAfter | CPM >= 15 | in >= 2 samples ---- #
AvsI.filtAfter.subset.cpm15.sp2 <- subset(x = AvsI.filtAfter.cpm15.sp2.DET,
	subset = AvsI.filtAfter.cpm15.sp2.DET$adj.P.Val < 0.001)
# ---- AvsI, filtAfter | CPM >= 15 | in >= 3 samples ---- #
AvsI.filtAfter.subset.cpm15.sp3 <- subset(x = AvsI.filtAfter.cpm15.sp3.DET,
	subset = AvsI.filtAfter.cpm15.sp3.DET$adj.P.Val < 0.001)

# TRANSITION: INFECTIVE vs. NON-INFECTIVE
# ---------------------------------------

# ---- IvsNI, filtFirst | CPM >= 0 | in >= 1 sample ---- #
IvsNI.filtFirst.subset.cpm0.sp1 <- subset(x = IvsNI.filtFirst.cpm0.sp1.DET,
	subset = IvsNI.filtFirst.cpm0.sp1.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtFirst | CPM >= 0 | in >= 2 samples ---- #
IvsNI.filtFirst.subset.cpm0.sp2 <- subset(x = IvsNI.filtFirst.cpm0.sp2.DET,
	subset = IvsNI.filtFirst.cpm0.sp2.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtFirst | CPM >= 0 | in >= 3 samples ---- #
IvsNI.filtFirst.subset.cpm0.sp3 <- subset(x = IvsNI.filtFirst.cpm0.sp3.DET,
	subset = IvsNI.filtFirst.cpm0.sp3.DET$adj.P.Val < 0.001)

# ---- IvsNI, filtFirst | CPM >= 5 | in >= 1 sample ---- #
IvsNI.filtFirst.subset.cpm05.sp1 <- subset(x = IvsNI.filtFirst.cpm05.sp1.DET,
	subset = IvsNI.filtFirst.cpm05.sp1.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtFirst | CPM >= 5 | in >= 2 samples ---- #
IvsNI.filtFirst.subset.cpm05.sp2 <- subset(x = IvsNI.filtFirst.cpm05.sp2.DET,
	subset = IvsNI.filtFirst.cpm05.sp2.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtFirst | CPM >= 5 | in >= 3 samples ---- #
IvsNI.filtFirst.subset.cpm05.sp3 <- subset(x = IvsNI.filtFirst.cpm05.sp3.DET,
	subset = IvsNI.filtFirst.cpm05.sp3.DET$adj.P.Val < 0.001)

# ---- IvsNI, filtFirst | CPM >= 10 | in >= 1 sample ---- #
IvsNI.filtFirst.subset.cpm10.sp1 <- subset(x = IvsNI.filtFirst.cpm10.sp1.DET,
	subset = IvsNI.filtFirst.cpm10.sp1.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtFirst | CPM >= 10 | in >= 2 samples ---- #
IvsNI.filtFirst.subset.cpm10.sp2 <- subset(x = IvsNI.filtFirst.cpm10.sp2.DET,
	subset = IvsNI.filtFirst.cpm10.sp2.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtFirst | CPM >= 10 | in >= 3 samples ---- #
IvsNI.filtFirst.subset.cpm10.sp3 <- subset(x = IvsNI.filtFirst.cpm10.sp3.DET,
	subset = IvsNI.filtFirst.cpm10.sp3.DET$adj.P.Val < 0.001)

# ---- IvsNI, filtFirst | CPM >= 15 | in >= 1 sample ---- #
IvsNI.filtFirst.subset.cpm15.sp1 <- subset(x = IvsNI.filtFirst.cpm15.sp1.DET,
	subset = IvsNI.filtFirst.cpm15.sp1.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtFirst | CPM >= 15 | in >= 2 samples ---- #
IvsNI.filtFirst.subset.cpm15.sp2 <- subset(x = IvsNI.filtFirst.cpm15.sp2.DET,
	subset = IvsNI.filtFirst.cpm15.sp2.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtFirst | CPM >= 15 | in >= 3 samples ---- #
IvsNI.filtFirst.subset.cpm15.sp3 <- subset(x = IvsNI.filtFirst.cpm15.sp3.DET,
	subset = IvsNI.filtFirst.cpm15.sp3.DET$adj.P.Val < 0.001)

# -----------------------------------------------------------

# ---- IvsNI, filtAfter | CPM >= 0 | in >= 1 sample ---- #
IvsNI.filtAfter.subset.cpm0.sp1 <- subset(x = IvsNI.filtAfter.cpm0.sp1.DET,
	subset = IvsNI.filtAfter.cpm0.sp1.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtAfter | CPM >= 0 | in >= 2 samples ---- #
IvsNI.filtAfter.subset.cpm0.sp2 <- subset(x = IvsNI.filtAfter.cpm0.sp2.DET,
	subset = IvsNI.filtAfter.cpm0.sp2.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtAfter | CPM >= 0 | in >= 3 samples ---- #
IvsNI.filtAfter.subset.cpm0.sp3 <- subset(x = IvsNI.filtAfter.cpm0.sp3.DET,
	subset = IvsNI.filtAfter.cpm0.sp3.DET$adj.P.Val < 0.001)

# ---- IvsNI, filtAfter | CPM >= 5 | in >= 1 sample ---- #
IvsNI.filtAfter.subset.cpm05.sp1 <- subset(x = IvsNI.filtAfter.cpm05.sp1.DET,
	subset = IvsNI.filtAfter.cpm05.sp1.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtAfter | CPM >= 5 | in >= 2 samples ---- #
IvsNI.filtAfter.subset.cpm05.sp2 <- subset(x = IvsNI.filtAfter.cpm05.sp2.DET,
	subset = IvsNI.filtAfter.cpm05.sp2.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtAfter | CPM >= 5 | in >= 3 samples ---- #
IvsNI.filtAfter.subset.cpm05.sp3 <- subset(x = IvsNI.filtAfter.cpm05.sp3.DET,
	subset = IvsNI.filtAfter.cpm05.sp3.DET$adj.P.Val < 0.001)

# ---- IvsNI, filtAfter | CPM >= 10 | in >= 1 sample ---- #
IvsNI.filtAfter.subset.cpm10.sp1 <- subset(x = IvsNI.filtAfter.cpm10.sp1.DET,
	subset = IvsNI.filtAfter.cpm10.sp1.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtAfter | CPM >= 10 | in >= 2 samples ---- #
IvsNI.filtAfter.subset.cpm10.sp2 <- subset(x = IvsNI.filtAfter.cpm10.sp2.DET,
	subset = IvsNI.filtAfter.cpm10.sp2.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtAfter | CPM >= 10 | in >= 3 samples ---- #
IvsNI.filtAfter.subset.cpm10.sp3 <- subset(x = IvsNI.filtAfter.cpm10.sp3.DET,
	subset = IvsNI.filtAfter.cpm10.sp3.DET$adj.P.Val < 0.001)

# ---- IvsNI, filtAfter | CPM >= 15 | in >= 1 sample ---- #
IvsNI.filtAfter.subset.cpm15.sp1 <- subset(x = IvsNI.filtAfter.cpm15.sp1.DET,
	subset = IvsNI.filtAfter.cpm15.sp1.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtAfter | CPM >= 15 | in >= 2 samples ---- #
IvsNI.filtAfter.subset.cpm15.sp2 <- subset(x = IvsNI.filtAfter.cpm15.sp2.DET,
	subset = IvsNI.filtAfter.cpm15.sp2.DET$adj.P.Val < 0.001)
# ---- IvsNI, filtAfter | CPM >= 15 | in >= 3 samples ---- #
IvsNI.filtAfter.subset.cpm15.sp3 <- subset(x = IvsNI.filtAfter.cpm15.sp3.DET,
	subset = IvsNI.filtAfter.cpm15.sp3.DET$adj.P.Val < 0.001)

################################################################
# (7) PRODUCING COMPARISON MATRIX BETWEEN 'DET' AMONG CONDITIONS
# --------------------------------------------------------------
################################################################

# Creating a list that contains, for each filtering condition, 
# all of the transcript names that are significantly diff.
# expressed. This represents the "reference object" from which
# intersect() will be computed in an "all vs. all" design.
DET.names.all.conditions.ref <- list(
	"AvsI.filtFirst.subset.cpm0.sp1" = row.names(AvsI.filtFirst.subset.cpm0.sp1),
	"AvsI.filtFirst.subset.cpm0.sp2" = row.names(AvsI.filtFirst.subset.cpm0.sp2),
	"AvsI.filtFirst.subset.cpm0.sp3" = row.names(AvsI.filtFirst.subset.cpm0.sp3),
	"AvsI.filtFirst.subset.cpm05.sp1" = row.names(AvsI.filtFirst.subset.cpm05.sp1),
	"AvsI.filtFirst.subset.cpm05.sp2" = row.names(AvsI.filtFirst.subset.cpm05.sp2),
	"AvsI.filtFirst.subset.cpm05.sp3" = row.names(AvsI.filtFirst.subset.cpm05.sp3),
	"AvsI.filtFirst.subset.cpm10.sp1" = row.names(AvsI.filtFirst.subset.cpm10.sp1),
	"AvsI.filtFirst.subset.cpm10.sp2" = row.names(AvsI.filtFirst.subset.cpm10.sp2),
	"AvsI.filtFirst.subset.cpm10.sp3" = row.names(AvsI.filtFirst.subset.cpm10.sp3),
	"AvsI.filtFirst.subset.cpm15.sp1" = row.names(AvsI.filtFirst.subset.cpm15.sp1),
	"AvsI.filtFirst.subset.cpm15.sp2" = row.names(AvsI.filtFirst.subset.cpm15.sp2),
	"AvsI.filtFirst.subset.cpm15.sp3" = row.names(AvsI.filtFirst.subset.cpm15.sp3),
	"AvsI.filtAfter.subset.cpm0.sp1" = row.names(AvsI.filtAfter.subset.cpm0.sp1),
	"AvsI.filtAfter.subset.cpm0.sp2" = row.names(AvsI.filtAfter.subset.cpm0.sp2),
	"AvsI.filtAfter.subset.cpm0.sp3" = row.names(AvsI.filtAfter.subset.cpm0.sp3),
	"AvsI.filtAfter.subset.cpm05.sp1" = row.names(AvsI.filtAfter.subset.cpm05.sp1),
	"AvsI.filtAfter.subset.cpm05.sp2" = row.names(AvsI.filtAfter.subset.cpm05.sp2),
	"AvsI.filtAfter.subset.cpm05.sp3" = row.names(AvsI.filtAfter.subset.cpm05.sp3),
	"AvsI.filtAfter.subset.cpm10.sp1" = row.names(AvsI.filtAfter.subset.cpm10.sp1),
	"AvsI.filtAfter.subset.cpm10.sp2" = row.names(AvsI.filtAfter.subset.cpm10.sp2),
	"AvsI.filtAfter.subset.cpm10.sp3" = row.names(AvsI.filtAfter.subset.cpm10.sp3),
	"AvsI.filtAfter.subset.cpm15.sp1" = row.names(AvsI.filtAfter.subset.cpm15.sp1),
	"AvsI.filtAfter.subset.cpm15.sp2" = row.names(AvsI.filtAfter.subset.cpm15.sp2),
	"AvsI.filtAfter.subset.cpm15.sp3" = row.names(AvsI.filtAfter.subset.cpm15.sp3),
	"IvsNI.filtFirst.subset.cpm0.sp1" = row.names(IvsNI.filtFirst.subset.cpm0.sp1),
	"IvsNI.filtFirst.subset.cpm0.sp2" = row.names(IvsNI.filtFirst.subset.cpm0.sp2),
	"IvsNI.filtFirst.subset.cpm0.sp3" = row.names(IvsNI.filtFirst.subset.cpm0.sp3),
	"IvsNI.filtFirst.subset.cpm05.sp1" = row.names(IvsNI.filtFirst.subset.cpm05.sp1),
	"IvsNI.filtFirst.subset.cpm05.sp2" = row.names(IvsNI.filtFirst.subset.cpm05.sp2),
	"IvsNI.filtFirst.subset.cpm05.sp3" = row.names(IvsNI.filtFirst.subset.cpm05.sp3),
	"IvsNI.filtFirst.subset.cpm10.sp1" = row.names(IvsNI.filtFirst.subset.cpm10.sp1),
	"IvsNI.filtFirst.subset.cpm10.sp2" = row.names(IvsNI.filtFirst.subset.cpm10.sp2),
	"IvsNI.filtFirst.subset.cpm10.sp3" = row.names(IvsNI.filtFirst.subset.cpm10.sp3),
	"IvsNI.filtFirst.subset.cpm15.sp1" = row.names(IvsNI.filtFirst.subset.cpm15.sp1),
	"IvsNI.filtFirst.subset.cpm15.sp2" = row.names(IvsNI.filtFirst.subset.cpm15.sp2),
	"IvsNI.filtFirst.subset.cpm15.sp3" = row.names(IvsNI.filtFirst.subset.cpm15.sp3),
	"IvsNI.filtAfter.subset.cpm0.sp1" = row.names(IvsNI.filtAfter.subset.cpm0.sp1),
	"IvsNI.filtAfter.subset.cpm0.sp2" = row.names(IvsNI.filtAfter.subset.cpm0.sp2),
	"IvsNI.filtAfter.subset.cpm0.sp3" = row.names(IvsNI.filtAfter.subset.cpm0.sp3),
	"IvsNI.filtAfter.subset.cpm05.sp1" = row.names(IvsNI.filtAfter.subset.cpm05.sp1),
	"IvsNI.filtAfter.subset.cpm05.sp2" = row.names(IvsNI.filtAfter.subset.cpm05.sp2),
	"IvsNI.filtAfter.subset.cpm05.sp3" = row.names(IvsNI.filtAfter.subset.cpm05.sp3),
	"IvsNI.filtAfter.subset.cpm10.sp1" = row.names(IvsNI.filtAfter.subset.cpm10.sp1),
	"IvsNI.filtAfter.subset.cpm10.sp2" = row.names(IvsNI.filtAfter.subset.cpm10.sp2),
	"IvsNI.filtAfter.subset.cpm10.sp3" = row.names(IvsNI.filtAfter.subset.cpm10.sp3),
	"IvsNI.filtAfter.subset.cpm15.sp1" = row.names(IvsNI.filtAfter.subset.cpm15.sp1),
	"IvsNI.filtAfter.subset.cpm15.sp2" = row.names(IvsNI.filtAfter.subset.cpm15.sp2),
	"IvsNI.filtAfter.subset.cpm15.sp3" = row.names(IvsNI.filtAfter.subset.cpm15.sp3)
)

####################
# INTERSECT MATRIX #

# Building a matrix with the INTERSECT of each pairwise comparison between 'treatments'
# (1 threshold combination = 1 treatment):

# Getting all the possible pairwise comparisons in the list
nms <- combn(names(DET.names.all.conditions.ref),
	2, # all combinations when comparing list elements 2 by 2
	FUN = paste0,
	collapse = "",
	simplify = FALSE)

# Putting together the information necessary to perform each pairwise comparison
all.combn <- combn(DET.names.all.conditions.ref,
	2,
	simplify = FALSE)
# Setting names for each comparison
all.combn <- setNames(all.combn, nms)

# Computing the intersect of each pairwise comparison
out.int <- lapply(all.combn, function(x) length(intersect(x[[1]], x[[2]])))
# Creating a vector with the results. This vector corresponds to the lower off
# diagonal in the final matrix.
out.int.v <- unlist(out.int, use.names = F)

# Creating the matrix, with respective values off diagonal. The diagonal will be = 100
# ------------------------------------------------------------------------------------
# Diagonal = 48 times the number 100. 48 because there are 48 comparisons
diag <- rep(100, 48)
# Empty matrix, 48x48
m.int <- matrix(NA, ncol = length(diag), nrow = length(diag))
# Filling in the lower part of the matrix with intersect values
m.int[lower.tri(m.int)] <- out.int.v
# Filling in the upper part of the matrix with the same intersect values, but
# transposed in order to output the numbers in the proper order off diagonal
m.int[upper.tri(m.int)] <- t(m.int)[upper.tri(t(m.int))]
# Defining the diagonal
diag(m.int) <- diag

################
# UNION MATRIX #

# Building a matrix with the UNION of each pairwise comparison between 'treatments'
# (1 threshold combination = 1 treatment):

# Computing the union of each pairwise comparison
out.un <- lapply(all.combn, function(x) length(union(x[[1]], x[[2]])))
# Creating a vector with the results. This vector corresponds to the lower off
# diagonal in the final matrix.
out.un.v <- unlist(out.un, use.names = F)

# Creating the matrix, with respective values off diagonal
# --------------------------------------------------------
# Empty matrix (NA), 48x48
m.un <- matrix(NA, ncol = length(diag), nrow = length(diag))
# Filling in the lower part of the matrix with intersect values
m.un[lower.tri(m.un)] <- out.un.v
# Filling in the upper part of the matrix with the same intersect values, but
# transposed in order to output the numbers in the proper order off diagonal
m.un[upper.tri(m.un)] <- t(m.un)[upper.tri(t(m.un))]
# Defining the diagonal
diag(m.un) <- diag

#########################################
# FINAL MATRIX - DIVIDE INTERSECT/UNION #

# Divide the matrix with intersect by the matrix with union to get
# the proportion of overlapping DE transcripts. The resulting 
# matrix = proportion of overlapping DE transcripts between 2  
# threshold combinations ('treatments').
final.mat <- (m.int/m.un)
# Setting column + row names
colnames(final.mat) <- names(DET.names.all.conditions.ref)
row.names(final.mat) <- names(DET.names.all.conditions.ref)
final.mat[upper.tri(final.mat)] <- 0 # Placing 0s above diagonal for better readability
diag(final.mat) <- 1

# Writing the final matrix to output file for visual 
# exploration/inspection in Excel for instance.
write.table(final.mat, file = "~/Desktop/RNAseq.filterThreshold.comparison.041116.txt")

########################################################
# 8.PRODUCING A CLUSTER GRAPH WITH SIMILARITY MATRIX - #
# ------------------------------------------------------
########################################################

# Packages required for the production of clustering graphs
library(cluster)

# Producing one matrix for each transition. The final matrix
# produced above contains all the information, but can be
# separated to simplify the graphical analysis.
AvsI.sim.mat <- final.mat[1:24, 1:24] # Adult vs Infective
IvsNI.sim.mat <- final.mat[25:48, 25:48] # Infective vs Non-Infective

# Converting similarities into dissimilarities (1-sim)
# and expressing these dissimilarities into distances.
# Distance method = euclidean.
AvsI.dist.mat <- as.dist(1 - AvsI.sim.mat) # Adult vs Infective
IvsNI.dist.mat <- as.dist(1 - IvsNI.sim.mat) # Infective vs Non-Infective

# Using hierarchical clustering to create groups of treatments with
# similar distance values. The 'Ward' clustering method is used here.
clust.AvsI <- hclust(AvsI.dist.mat, method="ward.D") # Adult vs Infective
clust.IvsNI <- hclust(IvsNI.dist.mat, method="ward.D") # Infective vs Non-Infective

# Clusters generated based on distance measures are plotted. This will
# define how many clusters there are among the groups.
plot(clust.AvsI, main = "Adult vs Infective") # Adult vs Infective
plot(clust.IvsNI, main = "Infective vs Non-infective") # Infective vs Non-Infective

# Cutting tree into 4 clusters, according to the results shown in 
# the plots generated above (cluster plots). Both dendrograms return
# the same pattern, i.e. 4 predominant groups.
groups.AvsI <- cutree(clust.AvsI, k = 4) # AvsI transition
groups.IvsNI <- cutree(clust.IvsNI, k = 4) # IvsNI transition

# Red rectangles around each cluster in each transition for 
# visualization purposes.
rect.hclust(clust.AvsI, k = 4, border = "red") # AvsI transition
rect.hclust(clust.IvsNI, k = 4, border = "red") # IvsNI transition

# Ploting cluster solutions on a PCA-like graph:
# ------------------
# Adult vs Infective
# ------------------
clusplot(as.matrix(AvsI.dist.mat),
	groups.AvsI, 
	color = T, 
	shade = T, 
	labels = 4, 
	lines = 0,
	cex.axis = 1.2,
	cex.lab = 1.2,
	cex = 1.2,
	main = "Clusters of filtering conditions\nAdult vs Infective worms")

# Adding labels for each point on the graph:
# Cluster no.2 (blue)
text(0, 1.98, label = c("CPM > 05 in > 1 sp"))
text(-0.42, 1.78, label = c("CPM > 05 in > 2 sp"))
text(-0.8, 1.45, label = c("CPM > 05 in > 3 sp"))
# Cluster no.3 (red)
text(-1.4, 0.4261749, label = c("CPM > 10 in > 1 sp"))
text(-1.4, -0.4316064, label = c("CPM > 10 in > 2 sp"))
text(-1.4, -0.8781314, label = c("CPM > 10 in > 3 sp"))
text(-1.4, -1.0298834, label = c("CPM > 15 in > 1 sp"))
# Cluster no.4 (green)
text(-1, -1.50506406, label = c("CPM > 15 in > 2 sp"))
text(-0.5, -1.74830822, label = c("CPM > 15 in > 3 sp"))
# Cluster no.1 (pink)
text(6.5, -0.06, label = c("CPM > 0 in > 1 sp"))
text(6.5, -0.2, label = c("CPM > 0 in > 2 sp"))
text(6.5, -0.34, label = c("CPM > 0 in > 3 sp"))

# --------------------------
# Infective vs Non-Infective
# --------------------------
clusplot(as.matrix(IvsNI.dist.mat), 
	groups.IvsNI,
	color = T,
	shade = T,
	labels = 4,
	lines = 0,
	cex.axis = 1.2,
	cex.lab = 1.2,
	cex = 1.2,
	main = "Clusters of filtering conditions\nInfective vs Non-Infective worms")

# Adding labels for each point on the graph:
# Cluster no.2 (red)
text(0, 1.98, label = c("CPM > 05 in > 1 sp"))
text(-0.46, 1.72, label = c("CPM > 05 in > 2 sp"))
text(-0.88, 1.35, label = c("CPM > 05 in > 3 sp"))
# Cluster no.3 (blue)
text(-1.4, 0.4261749, label = c("CPM > 10 in > 1 sp"))
text(-1.4, -0.4316064, label = c("CPM > 10 in > 2 sp"))
text(-1.4, -0.8781314, label = c("CPM > 15 in > 3 sp"))
text(-1.4, -1.0298834, label = c("CPM > 10 in > 1 sp"))
# Cluster no.4 (green)
text(-1, -1.45, label = c("CPM > 15 in > 2 sp"))
text(-0.55, -1.65, label = c("CPM > 15 in > 3 sp"))
# Cluster no.1 (pink)
text(6.5, -0.06, label = c("CPM > 0 in > 1 sp"))
text(6.5, -0.2, label = c("CPM > 0 in > 2 sp"))
text(6.5, -0.34, label = c("CPM > 0 in > 3 sp"))