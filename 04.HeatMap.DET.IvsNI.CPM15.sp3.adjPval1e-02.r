#######################################################################
#-#-# Hierarchical clustering from log(CPM) voom-transformed data #-#-#
#######################################################################

# Required packages
library(gplots)
library(RColorBrewer)

# Setting color palette for the heatmap
myheatcol <- rev(redgreen(75))

# Setting up working directory
setwd("~/Dropbox/travail/data/chapitre_2_rna-seq_ver/Hi_Seq/descriptive_graphs/heatmaps/IvsNI/cpm15.sp3")

# Importing dataset. Table that contains the log(CPM) of differentially
# expressed transcripts (DET) for this specific transition. This matrix
# was obtained by selecting the transcripts that have an adj.P.Val < 0.01
# when running limma-voom. I simply extracted the significant transcripts
# from the original matrix of voom-transformed data ('v' object in limma-
# voom script).
m <- as.matrix(
	read.table("ssol.IvsNI.DET.logCPM.CPM15.sp3.adjPval1e-02.txt", 
		header=T, 
		row.names=1
		)
	)

# Building a distance matrix among samples and transcripts, then clustering the 
# data by using the hclust( ) function. Euclidean distance is computed and ward.D2
# method is used for the hierarchical clustering.
hr_eucl <- hclust(dist(m, method="euclidean"), method="ward.D2")

# Visually looking at the clustering dendrogram
plot(hr_eucl)

# According to the plot produced above, there seems to be ~ 10 clusters
rect.hclust(hr_eucl, k = 10) # k = number of clusters

# Cutting the dendrogram tree according to the threshold chosen above
mycl_eucl <- cutree(hr_eucl, k = 10)

# Getting a color palette equal to the number of clusters
clusterCols_eucl <- rainbow(length(unique(mycl_eucl))) 

# creating a vector of colors for side bar
myClusterSideBar_eucl <- clusterCols_eucl[mycl_eucl]

# Plotting heatmap with all samples
heatmap.2(m, 
	Rowv=as.dendrogram(hr_eucl), 
	dendrogram="both", # clusters (dendrograms) by transcripts AND samples
	scale="row", 
	col=myheatcol, 
	density.info="none", 
	trace="none", 
	RowSideColors=myClusterSideBar_eucl)

##########################################################################
# 						DOWNSTREAM VERIFICATIONS                         #
##########################################################################

# Add cluster ID to data (i.e. to the initial matrix)
m_clust_eucl <- cbind(m, clusterID = mycl_eucl) # All samples

# Visualize the number of sequences in each cluster (according to cluster #)
# The numerical order is not the order in which the clusters appear on the heatmap
# This might be confusing.
table(mycl_eucl)

##################
#    RESULTS:    #
# ----------------
mycl_eucl

  1   2   3   4   5   6   7   8   9  10 
240 478 278 338  20  69 602 335 214 388

# This command will give the matrix data with cluster IDs as the last column
# and the table is ordered like the heat map, reading from bottom to top.
# That means, the first cluster that appears in the table that you obtain
# with the code below is the cluster located at the bottom of the heatmap.
# The second cluster in the table is the second from the bottom on the
# heatmap and so on.
m_clust_eucl[hr_eucl$order,]

# In that particular case with the data here, the order of the cluster no.
# on the heatmap is, from bottom to top: 9 (purple), 3 (yellow), 2 (orange), 8 (darkBlue),
# 10 (pink), 7 (lightBlue), 1 (red), 4 (green), 5 (?), 6 (bluePastelle).

# Here is a way of saving the resulting table with cluster IDs as the last column
write.table(m_clust_eucl[hr_eucl$order,], file = "hm.IvsNI.withClustersID.CPM15.sp3.txt", sep="\t")
# NOTE: it is important to add a tab in the first row so that the header is well formated.

#################################################################
# BOOTSTRAP EVALUATION OF HIERARCHICAL CLUSTERING - 'fpc' PACKAGE
# ---------------------------------------------------------------

# Loading 'fpc' package
library(fpc)

# Set desired number of clusters. Here, I chose 6 clusters because that's
# the number that came out of the heatmap after visualization with the
# rect.hclust() and cutree() functions.
kbest.p <- 10

# Run clusterboot() with hclust ('clustermethod=hclustCBI')
# using the 'complete' method ('method="complete"') and kbest.p
# clusters ('k=kbest.p'). It returns the results in an object
# called 'cboot.hclust'. The original matrix called 'm' is used
# as the first argument. The 'B=' argument stands for the number
# of iterations performed by clustboot(), but this can be changed
# to any desired number (default value = 100).
cboot.hclust <- clusterboot(m, 
	B = 1000, # 1000 bootstraps
	clustermethod = hclustCBI, 
	method = "ward.D2", 
	k = kbest.p)

###############
### RESULTS ###
###############

# Just calling the resulting object will output on the scree
# the results.
cboot.hclust

# Here are the results:

* Cluster stability assessment *
Cluster method:  hclust/cutree 
Full clustering results are given as parameter result
of the clusterboot object, which also provides further statistics
of the resampling results.
Number of resampling runs:  1000 

Number of clusters found in data:  10 

 Clusterwise Jaccard bootstrap (omitting multiple points) mean:
 [1] 0.5936722 0.6319901 0.6300652 0.6767492 0.7412475 0.4087933 0.6520781 0.5326715
 [9] 0.6333642 0.4939861
dissolved:
 [1] 345 225 230 159 209 721 190 453 233 591
recovered:
 [1] 191 240 240 380 525  92 286  56 227 338

# Results of the clustering are in cboot.hclust$result.
# The output of hclust() is in cboot.hclust$result$result.
# cboot.hclust$result$partition returns a vector of clusterlabels.
# This will give you, for each sequence, the cluster to which it
# belongs according to the hclust() function used.
groups <- cboot.hclust$result$partition