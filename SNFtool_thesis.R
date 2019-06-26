### Master thesis 
### R script for subclassification
### Similarity network fusion algorithm
### Author: Daphne Leunissen
### Date: 12-04-2019

#-----------------------------------------------------------------------------#
# Install and load required packages
#-----------------------------------------------------------------------------#
# first install required packges if needed:
# install.packages("SNFtool")
# install.packages("mlr")
# install.packages("tidyverse")

library(SNFtool)
library(mlr)
library(tidyverse)

options(stringsAsFactors = FALSE)
#-----------------------------------------------------------------------------#
# Set working directory and import data files 
#-----------------------------------------------------------------------------#

setwd("C:\\Users\\daphn\\OneDrive\\Documenten\\Master Systems Biology\\Master thesis\\R")

# import smMIPs result
data_clinchar <- read.delim("SNF_clinchar.txt", stringsAsFactor = F, header=TRUE , sep="\t")

# import annovar results 
data_mutations <- read.delim("SNF_mutation.txt", stringsAsFactor = F, header=TRUE , sep="\t")

#-----------------------------------------------------------------------------#
# Data preparation 
#-----------------------------------------------------------------------------#

## create dummy variables for PDL1 
pdl1 <- createDummyFeatures(data_clinchar$PDL1)

## delete pdl.nb column, the pdl1 expression wasn't determined for these data subject
pdl1$nb <- NULL

# replace PDL1 column in clinchar data with pdl data columns
data_clinchar[["PDL1"]] <- pdl1

# remove rows of which the sum equals zero (1 in pdl.nb column) --> pdl1 expression not determined
data_clinchar1 <- data_clinchar[rowSums(data_clinchar$PDL1)!=0,]

# set MBT-number as rownames
rownames(data_clinchar1) <- data_clinchar1[,1]
data_clinchar1$MBT.number <- NULL

# change sex column content to F = 2, M = 1
data_sex <- str_replace_all(data_clinchar1$Sex,"F", "2")
data_sex1 <- str_replace_all(data_sex,"M", "1")
data_sex1 <- as.numeric(data_sex1)
data_clinchar1[["Sex"]] <- data_sex1

# remove PDL1 columns, these are used as the labels so they cant be used for clustering
data_clinchar1$PDL1 <- NULL

# set MBT-number as rownames
rownames(data_mutations) <- data_mutations[,1]
data_mutations$MBT.number <- NULL

# remove samples from data_mutations, with no PDL1 expression determined
remove <- c("1095", "1115", "1107", "1108", "1137", "1140", "1147","1217" ,"1231", "414", "420", "528", "578", "585", "638")
data_mutations1 <- data_mutations[ !(row.names(data_mutations) %in% remove), ]

#-----------------------------------------------------------------------------#
# Similarity Network Fusion 
#-----------------------------------------------------------------------------#
# parameter settings:
# number of neighbors
K <- 20         
# hyperparameter
alpha <- 0.5    
# number of iterations
T <- 10         

# PDL1 label: 1 = no expression, 2 = low expression, 3 = high expression. 
label = c(2,3,2,1,1,3,2,1,1,1,1,2,3,1,2,1,2,1,1,2,2,1,3,2,1,1,2,1,1,2,3,2,3,3,1)

# # normalize
data_clinchar_norm <- as.data.frame(standardNormalization(data_clinchar1))
data_mutations_norm <- as.data.frame(standardNormalization(data_mutations1))

# calculate pair-wise distance
dist_clinchar <- (dist2(as.matrix(data_clinchar1),as.matrix(data_clinchar1)))^(1/2)
dist_mutations <- (dist2(as.matrix(data_mutations1),as.matrix(data_mutations1)))^(1/2)

# construct similarity graph
sim_clinchar <- affinityMatrix(dist_clinchar, K, alpha)
sim_mutations <- affinityMatrix(dist_mutations, K, alpha)

# display similarity graphs
displayClusters(sim_clinchar, label)
displayClusters(sim_mutations, label)

## next, we fuse the two graphs
## then the overall matrix can be computed by similarity network fusion(SNF):
W = SNF(list(sim_clinchar,sim_mutations), K, T)

# the number of clusters:
C = 2 
# perform spectral clustering 
group = spectralClustering(W,C);

# Get a matrix containing the group information, including spectralclustering results and true label
M_label=cbind(group,label)
colnames(M_label)=c("spectralClustering","Label")

# assign a color to each group
M_label_colors=t(apply(M_label,1,getColorsForGroups))

## Visualize the clusters present in the given similarity matrix, together with some sample information
displayClustersWithHeatmap(W, group, M_label_colors)

#-----------------------------------------------------------------------------#
# Prediction model
#-----------------------------------------------------------------------------#
# make a list out of the two data frames
data <- list(data_clinchar1, data_mutations1)

# data labels 
label = c(2,3,2,1,1,3,2,1,1,1,1,2,3,1,2,1,2,1,1,2,2,1,3,2,1,1,2,1,1,2,3,2,3,3,1)
# number of training cases
n = floor(0.8*length(label))
# training samples 
trainSample = sample.int(length(label), n)
# divide data in training and test data set
train = lapply(data, function(x) x[trainSample, ]) 
test = lapply(data, function(x) x[-trainSample, ])
# labels of training samples 
groups <- label[trainSample]

# set number of neighbors 
K = 20
# hyperparameter 
alpha = 0.5
# number of iterations 
t = 20
method = TRUE

## Apply the prediction function to the data
newLabel = groupPredict(train,test,groups,K,alpha,t,method)
## Compare the prediction accuracy
accuracy = sum(label[-trainSample] == newLabel[-c(1:n)])/(length(label) - n)

