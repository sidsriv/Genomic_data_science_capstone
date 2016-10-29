# Phenotypic Information --------------------------------------
## load phenotypic dataset
library (xlsx)
pheno_align <- read.xlsx("gdc-sample-phenotypic-summalign.xlsx", 1)

#-------------------------------------------------------------

#
# load Raw Data for each of the runs selected. Each run is loaded into R
# in the same order as they appear in the phenotypic dataset

# Run SRR1554538-----------------------------------------------
## load dataset
run4538 <- read.xlsx("FastQC_on_SRR1554538__RawData.xlsx", 1)

## subset by per sequence quality
pbsq4538 <- run4538[c(13:67),] # subset by per sequence quality

### rename variables
names(pbsq4538) <- c("base", "mean", "median", "lowerq", "upperq", "10pctl", "90pctl")

### base and mean as numeric 
pbsq4538$base <- as.numeric(pbsq4538$base)
pbsq4538$mean <- as.numeric(pbsq4538$mean)

## subset by per sequence quality scores
psqc4538 <- run4538[c(71:109), c(1:2)] 
names(psqc4538) <- c("quality", "count")

### quality and count as numeric 
psqc4538$quality <- as.numeric(psqc4538$quality)
psqc4538$count <- as.numeric(psqc4538$count)

# Run SRR1554541-----------------------------------------------
## load dataset
run4541 <- read.xlsx("FastQC_on__SRR1554541__RawData.xlsx", 1)

## subset the dataset
pbsq4541 <- run4541[c(13:67),] # subset by per sequence quality
names(pbsq4541) <- c("base", "mean", "median", "lowerq", "upperq", "10pctl", "90pctl")

### base and mean as numeric 
pbsq4541$base <- as.numeric(pbsq4541$base)
pbsq4541$mean <- as.numeric(pbsq4541$mean)

psqc4541 <- run4541[c(71:109), c(1:2)] # subset by per sequence quality scores
names(psqc4541) <- c("quality", "count")

### quality and count as numeric 
psqc4541$quality <- as.numeric(psqc4541$quality)
psqc4541$count <- as.numeric(psqc4541$count)


# Run SRR1554566-----------------------------------------------
## load dataset
run4566 <- read.xlsx("FastQC_on__SRR1554566__RawData.xlsx", 1)

## subset the dataset
pbsq4566 <- run4566[c(13:67),] # subset by per sequence quality

### rename variables
names(pbsq4566) <- c("base", "mean", "median", "lowerq", "upperq", "10pctl", "90pctl")

### base and mean as numeric 
pbsq4566$base <- as.numeric(pbsq4566$base)
pbsq4566$mean <- as.numeric(pbsq4566$mean)

## subset by per sequence quality scores
psqc4566 <- run4566[c(71:109), c(1:2)] 
names(psqc4566) <- c("quality", "count")

### quality and count as numeric 
psqc4566$quality <- as.numeric(psqc4566$quality)
psqc4566$count <- as.numeric(psqc4566$count)


# Run SRR1554535-----------------------------------------------
## load dataset
run4535 <- read.xlsx("FastQC_on__SRR1554535__RawData.xlsx", 1)

## subset the dataset
pbsq4535 <- run4535[c(13:67),] # subset by per sequence quality

### rename variables
names(pbsq4535) <- c("base", "mean", "median", "lowerq", "upperq", "10pctl", "90pctl")

### base and mean as numeric 
pbsq4535$base <- as.numeric(pbsq4535$base)
pbsq4535$mean <- as.numeric(pbsq4535$mean)

## subset by per sequence quality scores
psqc4535 <- run4535[c(71:109), c(1:2)] 
names(psqc4535) <- c("quality", "count")

### quality and count as numeric 
psqc4535$quality <- as.numeric(psqc4535$quality)
psqc4535$count <- as.numeric(psqc4535$count)


# Run SRR1554534-----------------------------------------------
## load dataset
run4534 <- read.xlsx("FastQC_on__SRR1554534__RawData.xlsx", 1)

## subset the dataset
pbsq4534 <- run4534[c(13:67),] # subset by per sequence quality

names(pbsq4534) <- c("base", "mean", "median", "lowerq", "upperq", "10pctl", "90pctl")

### base and mean as numeric 
pbsq4534$base <- as.numeric(pbsq4534$base)
pbsq4534$mean <- as.numeric(pbsq4534$mean)

## subset by per sequence quality scores
psqc4534 <- run4534[c(71:109), c(1:2)] 
names(psqc4534) <- c("quality", "count")

### quality and count as numeric 
psqc4534$quality <- as.numeric(psqc4534$quality)
psqc4534$count <- as.numeric(psqc4534$count)

# Run SRR1554561-----------------------------------------------
## load dataset
run4561 <- read.xlsx("FastQC_on__SRR1554561__RawData.xlsx", 1)

## subset the dataset
pbsq4561 <- run4561[c(13:67),] # subset by per sequence quality

### rename variables
names(pbsq4561) <- c("base", "mean", "median", "lowerq", "upperq", "10pctl", "90pctl")

### base and mean as numeric 
pbsq4561$base <- as.numeric(pbsq4561$base)
pbsq4561$mean <- as.numeric(pbsq4561$mean)

## subset by per sequence quality scores
psqc4561 <- run4561[c(71:109), c(1:2)] 
names(psqc4561) <- c("quality", "count")

### quality and count as numeric 
psqc4561$quality <- as.numeric(psqc4561$quality)
psqc4561$count <- as.numeric(psqc4561$count)


#------Mean per sequante quality and per sequence quality score------
## extract per sequence quality for each run
psqall <- cbind (pbsq4538[2], pbsq4541[2], pbsq4566[2], pbsq4535[2], pbsq4534[2], pbsq4561[2])

### calculate the mean for each column (each column is a sample)
psqallmean <- apply(psqall, 2, mean, na.rm.=TRUE)

## extract per sequence quality scores for each run
psqscall <- cbind (psqc4538[2], psqc4541[2], psqc4566[2], psqc4535[2], psqc4534[2], psqc4561[2])

### calculate the mean for each column (each column is a sample).
psqscallmean <- apply (psqscall, 2, mean, na.rm=TRUE)


#---- GC % per sample -----------------------------------------------

## information extracted from row 10, column 2 of each dataset in the order
## in which each run has been loaded in this script 
gc_content <- c(0.47, 0.46, 0.49, 0.47, 0.51, 0.52)


#----Create dataframe with run, sample, agegroup, pct_mapped, psqallmean, psqscallmean, gc_content

## calculate the percentge of mapped reads per sample
### the average between the left and right reads was calculated, and this average divied by the input
library (dplyr)
pheno_align <- mutate (pheno_align, pct_mapped = ((lrmapped + rrmapped)/2) / input)

## extract the variables from pheno_align, and mean per sequence quality, mean quality scores
## and gc_content
alldata <- data.frame (pheno_align[,c(1:6, 8, 19)], psqallmean, psqscallmean, gc_content)


#---Subset and summary stats by age group----------------------------------------------
## subset by age group
allf <- subset (alldata, ageg == "fetus")
alld <- subset (alldata, ageg == "adult")

## extract last four columns
allf_last4 <- select (allf, pct_mapped:gc_content)

alla_last4 <- select (alla, pct_mapped:gc_content)

## summary stats
summary (allf_last4)

summary (alla_last4)


