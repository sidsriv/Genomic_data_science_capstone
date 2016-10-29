
# set working directory
setwd("C:/Users/Felix/Dropbox/Coursera/Genomic-Data-Science/8-Genomic-datascience-capstone/github-files")

# update bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite()


# load datasets 
## gene expression dataset
gene_exp <- read.delim("gene_approx_counts_in_FPKM.txt")

## dataset with phenotypic information
phenotypic <- read.delim("gdc-sample-phenotypic-summalign.txt")

#-----------------------------------------------
# Data Cleaning
## explore colums names on datasets
names(gene_exp)

names(phenotypic)

## create tidy dataset
library (tidyr)
gene_exp_long <- gather (gene_exp, sample, fpkm, FPKM535:FPKM566)

## re-order the sample column

gene_exp_long$sample <- factor (gene_exp_long$sample, order=TRUE, levels=c("FPKM534", "FPKM535",
	"FPKM538", "FPKM541", "FPKM561", "FPKM566"))

#------------------------------------------------
# Exploratory Analysis

## descriptive statistics fpkms per sample
library(pander)
pander(summary (gene_exp_long))

## boxplot
boxplot (fpkm ~ sample, gene_exp_long, main= "Boxplots of FPKMs per Sample", ylab="FPKM", xlab="Samples", las=1)

## log2 transform FPLM
library (dplyr)
gene_exp_long <- mutate (gene_exp_long, logstrnf = log2(gene_exp_long$fpkm))

## boxplot of log2 transformed FPKM
boxplot(logstrnf ~ sample, gene_exp_long, main="Boxplots of Log2 Transformed FPKM per Sample", ylab="Log2 FPKM", xlab="Samples", las=1)

#----------------------------------------------
# Principal Commponent Analysis
## calculate the mean, including all fpkm
mean_all <- mean (gene_exp_long$fpkm)

## substract the mean of fpkms to each fpkm observation
gene_exp_long <- mutate (gene_exp_long, dfallmean = mean_all - gene_exp_long$fpkm)

## calculat the mean difference by run and create a data frame
pca <- gene_exp_long %>% group_by(sample) %>% summarise(dfmeanall_fpkm = mean(dfallmean))
pca <- as.data.frame(pca)

## calculate the mean fpkm for each sample
samples <- select(gene_exp, FPKM534, FPKM535, FPKM538, FPKM541, FPKM561, FPKM566)
mean_all_within <- apply(samples, 2, mean)

## create columns in which the mean fpkm for each sample is substracted to each sample
gene_exp <- mutate (gene_exp, FPKM534mindiff = FPKM534 - 28.5, FPKM535mindiff = FPKM535 - 28.5, FPKM538mindiff = FPKM538 - 56.4, FPKM541mindiff = FPKM541 - 56.4, FPKM561mindiff = FPKM561 - 14.4, FPKM566mindiff = FPKM566 - 37.5)

## calculate the mean for each sample mean difference columns 
gene_exp_mindiff_within <- select(gene_exp, FPKM534mindiff, FPKM535mindiff, FPKM538mindiff, FPKM541mindiff, FPKM561mindiff, FPKM566mindiff)
meandiff_within <- apply (gene_exp_mindiff_within, 2, mean)

## add the phenotypic variables
pca <- data.frame (pca, meandiff_within)
age <- c("adult", "adult", "fetus", "fetus", "adult", "fetus")
pca <- data.frame (pca, meandiff_within, age)
sample_id <- c("R2857", "R2869", "R3462", "R3485", "R4166", "R4706")
pca <- data.frame (pca, meandiff_within, age, sample_id)

## create a pca scatter plot
library (ggplot2)
pcaplot <- ggplot(pca, aes(meandiff_within, dfmeanall_fpkm, color=age, shape=sample)) + geom_point(size=3) + xlab(Variance FPKM withn Samples) + ylab(Overall FPKM Variance)
pcaplot