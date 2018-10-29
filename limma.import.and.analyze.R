## this script will work with limma to analyze affymetrix data

################################################################################################################
################################### INSTALLATION AND ACTIVATION ################################################
################################################################################################################


#this installs bioconductor, assume it will upgrade as well. run sparingly
#source("https://bioconductor.org/biocLite.R")
#biocLite()


#install/activate(?) these packages; again, only needed first time or if attempting an upgrade
#biocLite("limma")
#biocLite("statmod")
#biocLite("affy")
#install.packages("RMySQL")
#biocLite("GenomicFeatures")
#biocLite("VariantAnnotation")
#biocLite("OrganismDbi")
#biocLite("ensembldb")
#biocLite("biovizBase")
#biocLite("ggbio")
#biocLite("ReportingTools")
#biocLite("affycoretools")
#biocLite("oligo")
#biocLite("hta20transcriptcluster.db")

######################################################################################################
######################  ONE TIME ACTIONS  ##############################################
######################################################################################################


library(limma)
library(statmod)
library(affy)
library(affycoretools)
library(oligo)
library(hta20transcriptcluster.db)

Annotation <- read.csv("HTA-2_0.na35.2.hg19.transcript.csv/HTA-2_0.na35.2.hg19.transcript.csv", header=TRUE, skip=19)

######################################################################################################
######################  READ IN ALL TEST MATRICES ##############################################
######################################################################################################

targets <- read.csv("targets.txt", as.is=T)
targets.tGF.cGF <- read.csv("targets.GFonly.txt", as.is=T)
targets.tGFR.cGFR <- read.csv("targets.GFRonly.txt", as.is=T)
targets.cGFR.cGF <- read.csv("targets.noDrug.Strain.txt", as.is=T)
targets.tGFR.tGF <- read.csv("targets.drugonly.txt", as.is=T)

######################################################################################################
######################  GET FULL DIFF EXPRESSION MATRIX ##############################################
######################################################################################################

#this imports the CEL files in the working directory
dat <- read.celfiles(targets$FileName)

#Background correct, normalize, and calculate gene expression
#these objects are large "eset"s
#Summarise over genes
project.bgcorrect.norm.avg <- rma(dat, background=TRUE, normalize=TRUE)
#Calculate expression over individual exons
#I don't have any next steps for the exons yet, so this is useless now
#project.bgcorrect.norm.avg.Exons <- rma(dat, background=TRUE, normalize=TRUE, target="probeset")

#this adds annotations
#this doesn't work, not sure why... have to assume its for an older version of affymetrix
#eset.main <- annotTreatmentateEset(project.bgcorrect.norm.avg, hta20transcriptcluster.db)

#I will likely need to mess with the database myself to refine 
#check with https://www.biostars.org/p/274570/
#Annotation.Exons <- read.csv("Annotation/HuGene-2_1-st-v1.na36.hg19.probeset.csv", header=TRUE, skip=22)

#this should remove control probes (I don't think it is working?)
project.bgcorrect.norm.avg.filt <- exprs(project.bgcorrect.norm.avg)[which(rownames(project.bgcorrect.norm.avg) %in% Annotation[,1]),]

#Replace the transcript IDs with rownames (you can change to anything in annotation file, change the last numbers)
#This will replace your unique probe/gene IDs, use with caution, probably best to do this later in the process since
#you will also likely want to map to pathways and get annotations
#rownames(project.bgcorrect.norm.avg.filt) <- Annotation[match(rownames(project.bgcorrect.norm.avg.filt), Annotation[,1]), 8]

#this will spit out your normalized background corrected data matrix which you can use in other ways if you want
#no stat testing has been done yet

write.table(project.bgcorrect.norm.avg.filt, "full.expression.results.tsv", sep="\t", row.names = TRUE, quote=FALSE)

######################################################################################################
#########################################################
######################################################################################################



######################################################################################################
######################  TESTING GROUP PAIRS #########################################################
######################################################################################################

######################  cGFR vs. cGF #########################################################

dat <- read.celfiles(targets.cGFR.cGF$FileName)
project.bgcorrect.norm.avg <- rma(dat, background=TRUE, normalize=TRUE)
project.bgcorrect.norm.avg.filt <- exprs(project.bgcorrect.norm.avg)[which(rownames(project.bgcorrect.norm.avg) %in% Annotation[,1]),]

######################  LOADING IN METADATA 
Exp <- factor(targets.cGFR.cGF$Experiment)
Strain <- factor(targets.cGFR.cGF$Strain, levels=c("GF","GFR"))
Treatment <- factor(targets.cGFR.cGF$Treatment, levels=c("control","drug"))

#this compares all drug vs no drug samples
design <- model.matrix(~Exp+Strain)

#I'm just redefining the eset name so taht Ch.17 is easier to copy
eset <- project.bgcorrect.norm.avg.filt

fit <- lmFit(eset, design)
fit <- eBayes(fit, trend=T, robust=T)
results <- decideTests(fit)
summary(results)
write.table(summary(results), "controls.summary", quote=F)

#this generates a simple table of all differentially expressed genes
diff.cGFR.cGF <- topTable(fit, coef="StrainGFR", n=73)
write.table(diff.cGFR.cGF, "diff.cGFR.cGF", quote = F)


#this makes a scatter plot of all relative expression with the significant ones colored
#plotMD(fit,coef="StrainGFR",status=results[,3],values=c(1,-1),hl.col=c("red","blue"))

#VennDiagram?  how to plot this?
#vennDiagram(summary(results)$Treatmentdrug)

##now the question is how to test small subsets, for example untreated strains?
##The simple answer is to change your targets.txt file, make sure only target data files are in the...
##working directory, change the analysis matrix, and run all again

######################  tGFR vs. tGF #########################################################

dat <- read.celfiles(targets.tGFR.tGF$FileName)
project.bgcorrect.norm.avg <- rma(dat, background=TRUE, normalize=TRUE)
project.bgcorrect.norm.avg.filt <- exprs(project.bgcorrect.norm.avg)[which(rownames(project.bgcorrect.norm.avg) %in% Annotation[,1]),]

######################  LOADING IN METADATA 
Exp <- factor(targets.tGFR.tGF$Experiment)
Strain <- factor(targets.tGFR.tGF$Strain, levels=c("GF","GFR"))
Treatment <- factor(targets.tGFR.tGF$Treatment, levels=c("control","drug"))

#this compares all drug vs no drug samples
design <- model.matrix(~Exp+Strain)

#I'm just redefining the eset name so taht Ch.17 is easier to copy
eset <- project.bgcorrect.norm.avg.filt

fit <- lmFit(eset, design)
fit <- eBayes(fit, trend=T, robust=T)
results <- decideTests(fit)
summary(results)
write.table(summary(results), "treated.summary", quote=F)

#this generates a simple table of all differentially expressed genes
diff.tGFR.tGF <- topTable(fit, coef="StrainGFR", n=2337)
write.table(diff.tGFR.tGF, "diff.tGFR.tGF", quote = F)

######################  tGF vs. cGF #########################################################

dat <- read.celfiles(targets.tGF.cGF$FileName)
project.bgcorrect.norm.avg <- rma(dat, background=TRUE, normalize=TRUE)
project.bgcorrect.norm.avg.filt <- exprs(project.bgcorrect.norm.avg)[which(rownames(project.bgcorrect.norm.avg) %in% Annotation[,1]),]

######################  LOADING IN METADATA 
Exp <- factor(targets.tGF.cGF$Experiment)
Strain <- factor(targets.tGF.cGF$Strain, levels=c("GF","GFR"))
Treatment <- factor(targets.tGF.cGF$Treatment, levels=c("control","drug"))

#this compares all drug vs no drug samples
design <- model.matrix(~Exp+Treatment)

#I'm just redefining the eset name so that Ch.17 is easier to copy
eset <- project.bgcorrect.norm.avg.filt

fit <- lmFit(eset, design)
fit <- eBayes(fit, trend=T, robust=T)
results <- decideTests(fit)
summary(results)
write.table(summary(results), "GF.summary", quote=F)

#this generates a simple table of all differentially expressed genes
diff.tGF.cGF <- topTable(fit, coef="Treatmentdrug", n=996)
write.table(diff.tGF.cGF, "diff.tGF.cGF", quote = F)

######################  tGFR vs. cGFR #########################################################

dat <- read.celfiles(targets.tGFR.cGFR$FileName)
project.bgcorrect.norm.avg <- rma(dat, background=TRUE, normalize=TRUE)
project.bgcorrect.norm.avg.filt <- exprs(project.bgcorrect.norm.avg)[which(rownames(project.bgcorrect.norm.avg) %in% Annotation[,1]),]

######################  LOADING IN METADATA 
Exp <- factor(targets.tGFR.cGFR$Experiment)
Strain <- factor(targets.tGFR.cGFR$Strain, levels=c("GF","GFR"))
Treatment <- factor(targets.tGFR.cGFR$Treatment, levels=c("control","drug"))

#this compares all drug vs no drug samples
design <- model.matrix(~Exp+Treatment)

#I'm just redefining the eset name so that Ch.17 is easier to copy
eset <- project.bgcorrect.norm.avg.filt

fit <- lmFit(eset, design)
fit <- eBayes(fit, trend=T, robust=T)
results <- decideTests(fit)
summary(results)
write.table(summary(results), "GFR.summary", quote=F)

#this generates a simple table of all differentially expressed genes
diff.tGFR.cGFR <- topTable(fit, coef="Treatmentdrug", n=106)
write.table(diff.tGFR.cGFR, "diff.tGFR.cGFR", quote = F)


######################  tGFR vs. tGF controlling for cGFR vs. cGF #################################

dat <- read.celfiles(targets$FileName)
project.bgcorrect.norm.avg <- rma(dat, background=TRUE, normalize=TRUE)
project.bgcorrect.norm.avg.filt <- exprs(project.bgcorrect.norm.avg)[which(rownames(project.bgcorrect.norm.avg) %in% Annotation[,1]),]

######################  LOADING IN METADATA 
Exp <- factor(targets$Experiment)
Strain <- factor(targets$Strain, levels=c("GF","GFR"))
Treatment <- factor(targets$Treatment, levels=c("control","drug"))

#this attempts to compare effect of drug on GFR while controlling for strain differences and experimental differences
design <- model.matrix(~Exp+~Strain+Treatment)

#I'm just redefining the eset name so that Ch.17 is easier to copy
eset <- project.bgcorrect.norm.avg.filt

fit <- lmFit(eset, design)
fit <- eBayes(fit, trend=T, robust=T)
results <- decideTests(fit)
summary(results)
write.table(summary(results), "tGFR.cGFR.strain.control.summary", quote=F)

#this generates a simple table of all differentially expressed genes
diff.tGFR.cGFR <- topTable(fit, coef="Treatmentdrug", n=1536)
write.table(diff.tGFR.cGFR, "diff.tGFR.cGFR", quote = F)
