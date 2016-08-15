#######################
### 1 GENERAL SETUP ###
#######################
library(edgeR) # Importing edgeR package
# # Changing the working directory to the directory 
# where the count/cluster files are stored
setwd("~/YOUR/WORKING/DIRECTORY/")

#############################
### 2 CREATING A DGE-LIST ###
#############################

# The data that is imported here is a read-count matrix, i.e. rows = genes & columns = samples
# Read-count matrices can be obtained using various methods, among others we find
# HTseq or custom scripts based on counting the number of RNA-seq reads that align onto
# each reference sequence (by parsing alignement SAM file(s)). In this case, I obtained
# the read-count matrix by running the program "Corset" (https://github.com/Oshlack/Corset/wiki)
# on my raw data. You can find the pipeline I used on my github page
# (https://github.com/fohebert/corset_pipeline).
data <- read.delim("Schs.all-sp.counts.customSPnames.txt", row.names=1)
# Telling R which sampl belongs to which group
group <- factor(c("NonInfect","NonInfect","NonInfect","NonInfect","NonInfect","NonInfect","NonInfect","Infect","Infect","Infect","Infect","Adult","Adult","Adult"))
dge <- DGEList(counts=data,group=group)

#########################################
### 3 FILTERING LOWLY EXPRESSED GENES ###
#########################################
# Keeping genes with CPM of at least 15 in at least 3 
# samples of a given treatment. This can vary depending 
# on the dataset used. In my case, one of the treatments 
# has only 3 samples, so a CPM of 15 in at least 3 
# individuals ensures me that each transcript kept in the
# analysis will have a fair amount of reads. This threshold 
# was chosen upon an empirical test conducted on multiple
# threshold values ("verif.detection.DE.transcripts" directory)
keep <- rowSums(cpm(dge[,12:14])>15) >= 3 | # Adult
	rowSums(cpm(dge[,8:11])>15) >= 3 | # Infective
	rowSums(cpm(dge[,1:7])>15) >= 3 # Non-Infective
dge <- dge[keep,]
dge$samples$lib.size <- colSums(dge$counts)

###############################################
### 4 NORMALIZING ACCORDING TO LIBRARY SIZE ###
###############################################
# Performing normalization based on library size. This is the TMM normalization, which
# is the default one. There exists other methods, such as the one implemented in 
# DESeq, or RPKM, Total Count & Quantile.
dge <- calcNormFactors(dge)

##########################################
#### 5 PRODUCING MDS PLOT USING edgeR ####
##########################################
# This produces a PCA-type of graph, use parameter 'top=' to change number 
# of top diff. genes used (default = 500)
plotMDS(dge, 
	top=1000, 
	cex.axis=1.2, 
	col = c("darkgreen","orangered","blue")[unclass(dge$samples[,1])], 
	xlim=c(-2.5,6.5), 
	ylim=c(-2.5,6.5), 
	cex.lab=1.2, 
	cex = 1.3
	)
# Adding a legend to the figure with color code for the points on the graph
legend(3.5, 6.5, 
	c("NI = Non-Infective", "I = Infective", "A = Adult"), 
	col = c("blue", "orangered", "darkgreen"), 
	fill = c("blue", "orangered", "darkgreen"), 
	border = "black"
	)

#########################
#### 6 DESIGN MATRIX ####
#########################
# Producing a matrix design needed for the GLM. NB: the design contains an intercept
d <- model.matrix(~0+group, data=dge$samples)
# d2 <- model.matrix(~group, data=dge$samples) # Alternatively, this line would put the adult group as the intercept (i.e. "control gr.")

#################################
#### 7.A VOOM transformation ####
#################################
# Normal 'voom' function
v <- voom(dge, d, plot=F)
# 

# Converting the 'v' object into a matrix.
# Will be used to produce the heatmaps.
all.sp.mat <- as.matrix(v)

# Outputting the matrix in a file and this file
# will be used to extract, for both transitions,
# the DET. This matrix will be hand-modified in
# Excel in order to keep the appropriate columns
# in each comparison being made (either AvsI or
# IvsNI) and to keep the appropriate transcripts
# that show an adj.P.Val > 0.01. In the end, 2
# matrices will be produced: i) all DET for AvsI
# and ii) all DET for IvsNI. These log(CPM) values
# will serve as the input for the heatmap R code.
write.table(all.sp.mat,
	file= "/OUTPUT/DIRECTORY/file.name.txt", 
	quote=FALSE, 
	sep="\t", 
	col.names=NA)

# NB: now, this general matrix, as in the output file, can be
# used to produce 2 files, i.e. one matrix per transition that
# contains the columns that correspond to the samples being
# compared and the rows that correspond to the transcripts that
# have an adj.P.Val (FDR) < 0.01. In order to do that, the
# Python script < extract.lines.from.names.py > is required.
# This script will extract the lines that correspond to the
# transcripts that have an adj.P.Val < 0.01 (names of these
# transcripts must be fed to the script in the form of a
# list, i.e. one transcript name per line) for a given transition
# and then the resulting output file can be hand-modified to
# delete the columns that correspond to the samples not being
# compared (e.g if I want to produce the table for AvsI, I will
# manually remove the columns that correspond to NI samples).
# The resulting 2 files (2 log(CPM) tables for DET) will be
# used in the R script that produces the heatmaps: 
# < hm.AvsI.euclidean.r > and < hm.IvsNI.euclidean.r >.

######################################################################
### 8 CLASSIC GLM APPROACH, SIMILAR TO A ONE-WAY ANOVA USING limma ###
######################################################################
# Building a contrast matrix to make all the possible comparisons (A vs. NI, A vs. I, NI vs. I)
fit_limma <- lmFit(v, d)
contrast.matrix <- makeContrasts(groupAdult-groupNonInfect, groupAdult-groupInfect, groupInfect-groupNonInfect, levels=d)
fit <- contrasts.fit(fit_limma, contrast.matrix)
fit <- eBayes(fit)

############################################################################################################
# ADULT vs. NON-INFECTIVE TopTable
# --------------------------------

# Producing a full table with all the differentially expressed unigenes
DE_genes_AvsNI <- topTable(fit, coef=1, n=dim(fit)[1])
# Writing to output file
write.table(DE_genes_AvsNI, file='DGE.AvsNI.CPM15.sp3.limmaV.all-trans.tsv', quote=FALSE, sep='\t', col.names=NA)

############################################################################################################
# ADULT vs. INFECTIVE TopTable
# ----------------------------

# Producing a full table with all the differentially expressed unigenes
DE_genes_AvsI <- topTable(fit,coef=2, n=dim(fit)[1])
# Writing to output file
write.table(DE_genes_AvsI, file='DGE.AvsI.CPM15.sp3.limmaV.all-trans.tsv', quote=FALSE, sep='\t', col.names=NA)

############################################################################################################
# INFECTIVE vs. NON-INFECTIVE TopTable
# ------------------------------------

# Producing a full table with all the differentially expressed unigenes
DE_genes_IvsNI <- topTable(fit,coef=3, n=dim(fit)[1])
# Writing to output file
write.table(DE_genes_IvsNI, file='DGE.IvsNI.CPM15.sp3.limmaV.all-trans.tsv', quote=FALSE, sep='\t', col.names=NA)



#######################################################################
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # -
#      DOWN HERE YOU HAVE SOME CODE FOR THE PRODUCTION OF GRAPHS      #
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # -
#######################################################################

# NB: this code can also be found in an independent file
# in the directory where the graphs can be found 
# (i.e. descriptive_graphs > volcano.plots.DET)


#########################
#### VOLCANO PLOTS ####
#########################

# Libraries needed for plotting the data and to add text on
# the graph to highlight significant points.
library(Tmisc)
library(calibrate)

# Volcano plot to visualize the number and extent of differentially expressed transcripts

####################
# Adult vs Infective
# ------------------

# Transcript names
Gene <- row.names(DE_genes_AvsI)
# logFC
log2FoldChange <- DE_genes_AvsI$logFC
# P-Value
P.Val <- DE_genes_AvsI$P.Value
# Adjusted P-Value
adj.P.Val <- DE_genes_AvsI$adj.P.Val

# Data frame with unigene names, log2(FC), P-value & adjusted P-value
AvsI.df <- data.frame(
	Gene,
	log2FoldChange,
	P.Val,
	adj.P.Val
	)

# Building Volcano plot:
# ----------------------

# Basic plot
par(cex = 1.2)
with(AvsI.df, 
	plot(log2FoldChange, 
		-log10(P.Val), 
		pch=20,
		xlim = c(-10,14),
		ylim = c(0,17)
		)
	)

# Adding colored points: red if adj.P.Val<0.001
with(subset(AvsI.df, adj.P.Val<0.001), 
	points(log2FoldChange, -log10(P.Val), pch=20, col="red")
	)
# Orange if adj.P.Val < 1.0e-05
with(subset(AvsI.df, adj.P.Val<1.0e-05), 
	points(log2FoldChange, -log10(P.Val), pch=20, col="orange")
	)
# Green if adj.P.Val < 1.0e-08
with(subset(AvsI.df, adj.P.Val<1.0e-08), 
	points(log2FoldChange, -log10(P.Val), pch=20, col="green")
	)
# Adding Legend
legend(
	x = "bottomright",
	legend = c("FDR < 0.001", "FDR < 1.0e-05", "FDR < 1.0e-08"),
	fill = c("red", "orange", "green"),
	bty = "n",
	cex = 0.8
	)

# NEXT CODE NOT NECESSARY IF ONE WANTS TO MANUALLY LABEL THE POINT ON THE GRAPH #
# (can be super messy with long names, the graph becomes ugly)

# Label points with the textxy function from the calibrate plot
# Transcripts with adj.P.val < 1.0e-08 AND log(FC) < 0
with(subset(AvsI.df, adj.P.Val<1.0e-08 & log2FoldChange<0), 
	textxy(log2FoldChange, -log10(P.Val), labs=Gene, cex=0.8)
	)
# Transcripts with adj.P.val < 1.0e-10 AND log(FC) > 11
with(subset(AvsI.df, adj.P.Val<1.0e-10 & log2FoldChange>11), 
	textxy(log2FoldChange, -log10(P.Val), labs=Gene, cex=0.8)
	)

############################
# INFECTIVE vs NON-INFECTIVE
# --------------------------

# Transcript names
Gene <- row.names(DE_genes_IvsNI)
# logFC
log2FoldChange <- DE_genes_IvsNI$logFC
# P-Value
P.Val <- DE_genes_IvsNI$P.Value
# Adjusted P-Value
adj.P.Val <- DE_genes_IvsNI$adj.P.Val

# Data frame with unigene names, log2(FC), P-value & adjusted P-value
IvsNI.df <- data.frame(
	Gene,
	log2FoldChange,
	P.Val,
	adj.P.Val
	)

# Building Volcano plot:
# ----------------------

# Basic plot
par(cex = 1.2)
with(IvsNI.df, 
	plot(log2FoldChange, 
		-log10(P.Val), 
		pch=20,
		xlim = c(-10,12.5),
		ylim = c(0,13)
		)
	)

# Adding colored points: red if adj.P.Val<0.001
with(subset(IvsNI.df, adj.P.Val<0.01), 
	points(log2FoldChange, -log10(P.Val), pch=20, col="red")
	)
# Orange if adj.P.Val < 1.0e-05
with(subset(IvsNI.df, adj.P.Val<1.0e-03), 
	points(log2FoldChange, -log10(P.Val), pch=20, col="orange")
	)
# Green if adj.P.Val < 1.0e-08
with(subset(IvsNI.df, adj.P.Val<1.0e-04), 
	points(log2FoldChange, -log10(P.Val), pch=20, col="green")
	)
# Adding Legend
legend(
	x = "bottomright",
	legend = c("FDR < 0.01", "FDR < 1.0e-03", "FDR < 1.0e-04"),
	fill = c("red", "orange", "green"),
	bty = "n",
	cex = 0.8
	)

# NEXT CODE IS NOT NECESSARY IF ONE WANTS TO MANUALLY LABEL THE POINT ON THE GRAPH #
# (can be super messy with long names, the graph becomes ugly)

# Label points with the textxy function from the calibrate plot
# Transcripts with adj.P.val < 1.0e-08 AND log(FC) < 0
with(subset(IvsNI.df, adj.P.Val<1.0e-05 & log2FoldChange<0), 
	textxy(log2FoldChange, -log10(P.Val), labs=Gene, cex=0.8)
	)
# Transcripts with adj.P.val < 1.0e-09 AND log(FC) > 10
with(subset(IvsNI.df, adj.P.Val<1.0e-05 & log2FoldChange>11), 
	textxy(log2FoldChange, -log10(P.Val), labs=Gene, cex=0.8)
	)