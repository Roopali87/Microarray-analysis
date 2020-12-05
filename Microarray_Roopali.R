if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

#Package for human ref and annotation 
#BiocManager::install("pd.hugene.1.0.st.v1")
#library("pd.hugene.1.0.st.v1")
#BioacManager::install("hugene10stv1cdf")
#library("hugene10stv1cdf")
###########################################

BiocManager::install("limma")
library("limma")


## For wheat reference genome
BiocManager::install("pd.wheat")
library("pd.wheat")

## For wheat affymetrix genome array annotation
BiocManager::install("wheatcdf")
library("wheatcdf")

##
BiocManager::install("affyPLM")
library("affyPLM")

##
BiocManager::install("affy")
library("affy")

BiocManager::install("IRanges")
library("IRanges")
BiocManager::install("RColorBrewer")
library("RColorBrewer")
BiocManager::install("methods")
library("methods")
BiocManager::install("S4Vectors")
library("S4Vectors")
BiocManager::install("Hmisc")
library("Hmisc")

#To get Boxplot for Pre-Normalized Expression
targets <- readTargets("Target.txt")
#Read CEL Files
dat <- ReadAffy(filenames = targets$FileName) ###FileName is the first column name in txt
dat
eset<-rma(dat)
eset
normset<-eset
pData(normset)


###Oligo package, if required.
#BiocManager::install("oligo", version = "3.8")
#library("oligo")
#Oligo Read in the CEL files in directory
#celFiles<- list.celfiles()
#affyRaw<-read.celfiles(celFiles)
#affyRaw
#exprset <-affyRaw
#exprset
#pData(exprset)
#RMA Normalization
#BiocManager::install("gcrma")
#library("gcrma")
#exprset <- gcrma(exprset)

#Finally, save the data to an output file to be used by other programs, etc (Data will be log2 transformed and normalized)
write.exprs(eset,file="PostNormalisedData.txt")


#Boxplot Before Normalization ##############################################################
tiff(file="Control-treatment Pre-Normalization [BoxPlot].tiff", bg="transparent", width=600, height=600)
par(mar = c(7, 5, 3, 2) + 0.1); # This sets the plot margins #bottom,left,top,right
boxplot(dat,col="red", main="Pre-Normalization", las=2, cex.axis=0.74, ylab="Intensities")#, ylim=c(2,14))
title(xlab = "Sample Array", line = 6); # Add x axis title
dev.off()

#Boxplot After Normalization
tiff(file="Control-treatment Post-Normalization [BoxPlot].tiff", bg="transparent", width=600, height=600)
par(eset,mar = c(7, 5, 3, 2) + 0.1); # This sets the plot margins #bottom,left,top,right
boxplot(eset, col="blue",main="Post-Normalization", las=2, cex.axis=0.74, ylab="Intensities")#, ylim=c(2,14))
title(xlab = "Sample Array", line = 6); # Add x axis title
dev.off()
###################################################################################

#https://rpubs.com/ge600/limma
data<-read.table(file = "PostNormalisedData.txt", header = T, row.names=1)
groups = gsub("_.*", "", colnames(data)) #clip sample name after underscore(_)
groups <- factor(groups, levels = c("Control","Treatment") )
design <- model.matrix( ~ 0 + groups )
colnames(design) <- c("Control","Treatment") 
design


library(limma)
# Fits a linear model for each gene based on the given series of arrays.
fit <- lmFit(data, design) 
#write.csv(fit, "lmFit.csv", quote = F)


#Matrix
#design<-model.matrix(~factor(c("Control", "Control", "Treatment", "Treatment")))
#colnames(design)<-c("Control","Treatment")


head(data)



#Contrast Matrix Design
cont.matrix <- makeContrasts(contrasts = "Treatment-Control", levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
# Computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
fit2 <- eBayes(fit2, trend=FALSE)

# calls differential gene expression 1 for up, -1 for down
#results <- decideTests(fit2, p.value = 0.05, lfc= log2(2) )
topGenes =topTable(fit2, number = 1e12,sort.by="M" )
head(topGenes)
write.csv(topGenes, "Result Top Table Final.csv", quote = F)
############################################################################