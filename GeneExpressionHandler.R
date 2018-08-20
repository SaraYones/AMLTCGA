#Read Gene Expression data and correct for Batch, effects
#-----------------Preprocessing-----------------------------------------------
#Read one file from the GE files
#temp=read.table(file = '/Users/sarayones/Downloads/TCGA_LAML/GE/bcgsc.ca_LAML.IlluminaGA_RNASeq.Level_3/TCGA-AB-3008-03A-01T-0736-13.gene.quantification.txt',header = T,sep="\t")
#dim(temp)
#GeneExpressionMatrix=matrix(0,nrow=dim(temp)[1],ncol=length(significant_genes))
#Get the list of samples number in the GE folder and see return the number of matching samples with CNV and SNP files 

library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
library(sva)
library(limma)
library(ggplot2)
require(reshape2)

files=list.files("/Users/sarayones/Downloads/TCGA_LAML/GE/bcgsc.ca_LAML.IlluminaGA_RNASeq.Level_3/")
SampleGEBarcode=NULL
SamplesGE=NULL
for(i in 1:length(files))
{
  
  if(grepl("*.gene.quantification.txt", files[i]) == TRUE)
{
  #  SampleGEBarcode==gsub("TCGA-AB-([0-9]{4})-03{A,B}-01T-[0-9]{4}-13.gene.quantification.txt","\\1",files)
    #SamplesCNV=gsub("TCGA_AB_([0-9]{4})_03A_01T_[0-9]{4}-13","\\1",SamplesCNVbarcodes)
   SampleGEBarcode=append(SampleGEBarcode,files[i])
   
}
  
}
#Sample Numbers
SampleGE=gsub("(TCGA-AB-([0-9]{4})-(03A)-01T-[0-9]{4}-13.gene.quantification.txt)","\\2",SampleGEBarcode)
SampleGE=gsub("(TCGA-AB-([0-9]{4})-(03B)-01T-[0-9]{4}-13.gene.quantification.txt)","\\2",SampleGE)
#Samples Barcode
SampleGEBarcode=gsub("(TCGA-AB-[0-9]{4}-03A-01T-[0-9]{4}-13).gene.quantification.txt","\\1",SampleGEBarcode)
SampleGEBarcode=gsub("(TCGA-AB-[0-9]{4}-03B-01T-[0-9]{4}-13).gene.quantification.txt","\\1",SampleGEBarcode)
#Remove samples from GE data that are not part of CNV or mutation
IndexGEnotCNV=which(!(SampleGE %in% SamplesCNV))
SampleGE=SampleGE[-IndexGEnotCNV]
SampleGEBarcode=SampleGEBarcode[-IndexGEnotCNV]
#update Binary Matrix to only  have the samples which go with Gene expression samples and also make new FilteredGeneList and Ones list
GEBinaryMatrix=BinaryMatrix[SampleGE,]
GEdetailedBinaryMatrix=detailedBinaryMatrix[SampleGE,]
Oneslist<-createOneslist(GEBinaryMatrix,SampleGEBarcode,significant_genes)
hist(Oneslist,main = "Frequency of Aberrations", xlab="Frequency of Aberrations", ylab="Number of genes")
#Filter Ones list (if you want)
Oneslist[sapply(Oneslist, function(x) x>=quantile(Oneslist)[4])]
#Check the number of genes in each mutation profile after filteration
table(Oneslist[sapply(Oneslist, function(x) x>=quantile(Oneslist)[4])])
FilteredGeneList<-FilterGenes(significant_genes,quantile(Oneslist)[4])
#--------------#Create Gene/Samples matrix----------------------------------------------
temp=NULL
temp=read.table(file = '/Users/sarayones/Downloads/TCGA_LAML/GE/bcgsc.ca_LAML.IlluminaGA_RNASeq.Level_3/TCGA-AB-3008-03A-01T-0736-13.gene.quantification.txt',header = T,sep="\t")
names(FilteredGeneList)
namesFiteredGeneList=names(FilteredGeneList)
#************************Extract the indices of the genes that match with the regex i create for the Filtered Gene list but  first put them in a form that can be matched by this regex*******************************************
tempGenes=temp[,"gene"]
#tempGenes=tempGenes[which(grepl(RegexGE,tempGenes))]
#Transform the gene names into something i can match a regex with 
#Extract the genes names from this format AMY1A|276_calculated
tempGenes=gsub("((*)(\\|)([0-9]+_calculated))","\\2",tempGenes)
#Extract the genes names from this format AMY1A|276|3of3_calculated
tempGenes=gsub("((*)(\\|)([0-9]+(\\|)([0-9]of[0-9]_calculated)))","\\2",tempGenes)
#Find the indices that match with the names of the Filtered Gene List

for(i in 1:length(namesFiteredGeneList))
{
  
  namesFiteredGeneList[i]=paste("^",as.character(namesFiteredGeneList[i]),sep="")
  namesFiteredGeneList[i]=paste(as.character(namesFiteredGeneList[i]),"$",sep="")
}
#there are 2 genes that do not exist in the file of the quantification so they are not included in the analysis these are:
#ZNF783 and KATA6A

RegexGE=paste("",as.character(namesFiteredGeneList),collapse="|",sep="")
#IndexGenes=tempGenes[which(grepl(RegexGE,tempGenes))]
IndexGenes=which(grepl(RegexGE,tempGenes))
#---------------------------------------------------- Create the Matrix----------------------------------------------------------------------
#GeneExpressionMatrix=matrix(0,ncol=length(namesFiteredGeneList),nrow=length(SampleGEBarcode))
#colnames(GeneExpressionMatrix)=tempGenes[IndexGenes]
#row.names(GeneExpressionMatrix)=SampleGEBarcode
##########Old version
#ncol=length(namesFiteredGeneList)
#nrow=length(SampleGEBarcode)
#######################################
#I did that because some genes in the FilteredGeneList do no exist in the file for Gene Expression like this gene ZNF783
ncol=length(IndexGenes)
nrow=length(SampleGEBarcode)

GeneExpressionMatrix<-createGEMatrix(tempGenes,ncol,nrow,SampleGEBarcode,IndexGenes)
#-------------------------------Read the Quantification after creating the matrix------------------------------------------
#for (i in i:length(SampleGEBarcode))
GeneExpressionMatrix<-readQuantification('/Users/sarayones/Downloads/TCGA_LAML/GE/bcgsc.ca_LAML.IlluminaGA_RNASeq.Level_3/',GeneExpressionMatrix,SampleGEBarcode,IndexGenes)
readQuantification= function(filepath,GEMatrixlocal,SampleBarcodeLocal,IndexGenes){

for (i in  1:length(SampleGEBarcode))
{
  tempfile=paste(as.character(SampleGEBarcode[i]),as.character(".gene.quantification.txt"),sep="")
  tempfile=read.table(file=paste('/Users/sarayones/Downloads/TCGA_LAML/GE/bcgsc.ca_LAML.IlluminaGA_RNASeq.Level_3/',as.character(tempfile),sep=""),header = T,sep="\t")
 # names(FilteredGeneList)
  if(!is.null(IndexGenes))
  { 
    RPKMValues=tempfile[IndexGenes,"RPKM"]
  }
  else
  {
    RPKMValues=tempfile[,"RPKM"]
  }
  #View(tempfile)
  GEMatrixlocal[as.character(SampleGEBarcode[i]),]=RPKMValues

}
  return(GEMatrixlocal)
}
#--------------------------------------------------------------------------------------------------------------------
#Draw denisty plots for the genes to check the distribution if it's not normal apply a transformation to make it look kind of normalish
plotGE= function(filepath,GEMatrix,genestoplot){
  pdf(filepath)
  for(i in 1:length(genestoplot))
  {
    d <- density(GEMatrix[,genestoplot[i]])
    plot(d,main=genestoplot[i])
  }
  dev.off()
}
#-----------------------------------------------------------------------------------------------------------
#Transform the Matrix using log2
#Check this link https://www.researchgate.net/post/What_is_Log_transformation_and_why_do_we_do_it_in_gene_expression_analysis very 
#useful for understanding why we do log transformation
plotGE('/Users/sarayones/Desktop/Projects/AMLProject/Code/GEv1.1.pdf',GeneExpressionMatrix,tempGenes[IndexGenes])
logGeneExpressionMatrix=GeneExpressionMatrix+0.00001
logGeneExpressionMatrix=log2(logGeneExpressionMatrix)
plotGE('/Users/sarayones/Desktop/Projects/AMLProject/Code/GETransformationv1.1.pdf',logGeneExpressionMatrix,tempGenes[IndexGenes])

#-------------------------------------------------------------------------------------------------------------
#Check For batch effects by plotting a PCA for all samples with all the genes in the Gene expression files
row.names(GEBinaryMatrix)=SampleGEBarcode
MutantWildtype=getMutantWildtypeSamples("IDH1",GEBinaryMatrix)
MutantBatches=extractBatchesClasses(MutantWildtype$Mutant)
WildtypeBatches=extractBatchesClasses(MutantWildtype$Wildtype)
Batches=extractBatchesClasses(SampleGEBarcode)
extractBatchesClasses= function(MutantWildtypelocal){
Batches=gsub("TCGA-AB-[0-9]{4}-03A-01T-([0-9]{4})-13","\\1",MutantWildtypelocal)
Batches=gsub("TCGA-AB-[0-9]{4}-03B-01T-([0-9]{4})-13","\\1",Batches)
  return(Batches)
}
#Create a gene expression matrix for all the genes available to be able to catch if there is a batch effect or not

ncol=length(temp[,"gene"])
nrow=length(SampleGEBarcode)
GeneExpressionMatrixBatch<-createGEMatrix(temp[,"gene"],ncol,nrow,SampleGEBarcode,NULL)
GeneExpressionMatrixBatch<-readQuantification('/Users/sarayones/Downloads/TCGA_LAML/GE/bcgsc.ca_LAML.IlluminaGA_RNASeq.Level_3/',GeneExpressionMatrixBatch,SampleGEBarcode,NULL)
logGeneExpressionMatrixBatch=GeneExpressionMatrixBatch+0.00001
logGeneExpressionMatrixBatch=log2(logGeneExpressionMatrixBatch)
#----------------------------------------------------------------------------------------------------------------
plotPCA(logGeneExpressionMatrix,Batches)
plotPCA(logGeneExpressionMatrixBatch,Batches)
plotPCA= function(GeneExpressionMatrixlocal,Groups){

# log transform 
# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 

#ir.pca <- prcomp(GeneExpressionMatrixlocal,
 #                center = TRUE,
  #               scale. = TRUE) 
  ir.pca <- prcomp(GeneExpressionMatrixlocal,
                                 center = TRUE) 


g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
              groups =Groups, ellipse = FALSE, 
              circle = FALSE,var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
}
#--------------------------------------------------------------------------
createGEMatrix= function(tempGeneslocal,ncollocal,nrowlocal,SampleBarcodeslocal,IndexGenes){
  GeneExpressionMatrix=matrix(0,ncol=ncollocal,nrow=nrowlocal)

  if(!is.null(IndexGenes))
  {
  print("hello")
 
   colnames(GeneExpressionMatrix)=tempGeneslocal[IndexGenes]
   print("hello")
  }
  else
  {
  colnames(GeneExpressionMatrix)=tempGeneslocal
  }
  row.names(GeneExpressionMatrix)=SampleBarcodeslocal
  return(GeneExpressionMatrix)
}

#---------------------------Adjust for batch effect using SVA package----------------------------------------------------
clinicalfile='/Users/sarayones/Downloads/TCGA_LAML/metadataTCGA/laml_metadata/Clinical/Biotab/nationwidechildrens.org_clinical_patient_laml.txt'

requiredPheno="cyto_risk_group"
RegexSamplesGE=Regex=paste("",as.character(SampleGE),collapse="|",sep="")
matchedcol="bcr_patient_barcode"
pheno<-createPheno(SampleGE,clinicalfile,RegexSamplesGE,requiredPheno,matchedcol,Batches,"batch")
logGeneExpressionMatrixWithoutBatch<-removeBatchEffect(logGeneExpressionMatrixBatch,pheno,requiredPheno)
#For survival classifier
AllGenesGE=logGeneExpressionMatrixWithoutBatch
#-------------------------------------------------------------------------------------------------------
#Keep only genes that  have been filtered before after removing the batch effects from the matrix
#Plot Gene expression for the updated genes and see if anything has changed
logGeneExpressionMatrixWithoutBatch=t(logGeneExpressionMatrixWithoutBatch[IndexGenes,])
colnames(logGeneExpressionMatrixWithoutBatch)=colnames(logGeneExpressionMatrix)
plotGEboxplots('/Users/sarayones/Desktop/Projects/AMLProject/Code/GEboxplotwithBatcheffectggplotv1.1.pdf',logGeneExpressionMatrix)
plotGEboxplots('/Users/sarayones/Desktop/Projects/AMLProject/Code/GEboxplotwithoutBatcheffectggplotv1.1.pdf',logGeneExpressionMatrixWithoutBatch)
plotGE('/Users/sarayones/Desktop/Projects/AMLProject/Code/GERemovedBatcheffect.pdf',logGeneExpressionMatrixWithoutBatch,tempGenes[IndexGenes])

#-------------------------------------------------------------------------------------------------------
createPheno= function(SampleGE,clinicalfile,RegexSamplesGE,requiredPheno,matchedcol,Batches,Batchescol){
  pheno=matrix(0,ncol=2,nrow=length(SampleGE))
  rownames(pheno)=SampleGE
   clinical=read.table(file=clinicalfile,header = T,sep="\t")
  #Remove first two rows
  clinical=clinical[-c(1,2),]
  colnames(pheno)=list(Batchescol,requiredPheno)
  SamplesClinicalbarcode=clinical[,matchedcol]
  IndicesClinical=which(grepl(RegexSamplesGE, SamplesClinicalbarcode))  
  pheno[,requiredPheno]=as.character(clinical[IndicesClinical,requiredPheno])
  pheno[,Batchescol]=as.character(Batches)
  return(pheno)
}


removeBatchEffect= function(logGeneExpressionMatrixBatch,pheno,requiredPheno){
  
 #Ask Mateusz how to deal with this because I dont want an exact name for the factor i want to substitute the value of requiredPheno
  #pheno2=as.data.frame(pheno) 
  
  #Creating a model for sva 
 # mod = model.matrix(~cyto_risk_group, data=pheno)
  # mod0 = model.matrix(~1,as.data.frame(pheno))
#   n.sv = num.sv(t(logGeneExpressionMatrixBatch),mod,method="leek")
 #  svobj = sva(t(logGeneExpressionMatrixBatch),mod,mod0,n.sv=n.sv)
#   modSv = cbind(mod,svobj$sv)
 #  mod0Sv = cbind(mod0,svobj$sv)
#fit = lmFit(t(logGeneExpressionMatrixBatch),modSv)
  #Creating a model for COMBAT with known batches
   batch = pheno[,"batch"]
   modcombat = model.matrix(~1, data=as.data.frame(pheno))
   combat_edata = ComBat(dat=as.data.frame(t(logGeneExpressionMatrixBatch)), batch=batch, mod=modcombat)
   return(combat_edata)
}
#---------------------------------------------------------------------------------------------------
plotGEboxplots= function(filepath,GEMatrix){
  pdf(filepath)
  #boxplot(GEMatrix)
  sample <- rownames(GEMatrix)
  d.f <- data.frame(sample,GEMatrix)
  d.f2 <- melt(d.f, id.vars = "sample")
  geneNamesplot=colnames(GEMatrix)
  
  print(ggplot(d.f2, aes(x=as.factor(variable), y=value)) + geom_boxplot(fill="slateblue", alpha=0.2)+xlab("genes"))
  
  
  dev.off()
}
#---------------------------Compare Expression Levels between Mutant and Wildtype samples---------------------------
#Prepare Matrix before comparing Gene Expression (Make the row names to carry only the number in the barcode)
rownames(logGeneExpressionMatrixWithoutBatch)=SampleGE
rownames(GEBinaryMatrix)=SampleGE
GeneExpressionSignificance=computeSignificanceWilcoxon(logGeneExpressionMatrixWithoutBatch,significance_level,GEBinaryMatrix)
#Some times wilcoxon gives NA neither true nor falso so you have to find a way to deal with these genes which have no significane or they do not change like CCL11 and POM121
#----------------------------Find the genes that are significantly altered then access the GEdetailedBinaryMatrix to check their mutation status---------------------
SignificantGenes=GeneExpressionSignificance<=0.05
MutationStatus=GEdetailedBinaryMatrix[,names(SignificantGenes)]
table(MutationStatus[,"TET2"])
#Plot details of each of the genes (CNV,SNV and CNV-SNV mutational status)
plotdetails(MutationStatus,SignificantGenes,'/Users/sarayones/Desktop/Projects/AMLProject/Code/detailsMutationStatusv1.1.pdf')
#------------------------------plot boxplot for gene expression for wildtype and mutatnt for each of the genes---------------------------------------------
FilteredSignificantGenes=SignificantGenes[sapply(SignificantGenes, function(x) x == TRUE)]
plotGEMutatantWildtype(GEBinaryMatrix,logGeneExpressionMatrixWithoutBatch,FilteredSignificantGenes,'/Users/sarayones/Desktop/Projects/AMLProject/Code/boxplotGEMutantWildtypev1.1.pdf')
#-----------------------------------------------------------------------------------------------------------------------------------------
computeSignificanceWilcoxon= function(logGeneExpressionMatrixWithoutBatch,significance_level,GEBinaryMatrix){
geneNamesLocal=colnames(logGeneExpressionMatrixWithoutBatch)
z=NULL
for(i in 1:length(geneNamesLocal)){
  MutantWildtype=getMutantWildtypeSamples(geneNamesLocal[i],GEBinaryMatrix)
 # MutantBatches=extractBatchesClasses(MutantWildtype$Mutant)
 # WildtypeBatches=extractBatchesClasses(MutantWildtype$Wildtype)
  x=logGeneExpressionMatrixWithoutBatch[MutantWildtype$Mutant,geneNamesLocal[i]]
  y=logGeneExpressionMatrixWithoutBatch[MutantWildtype$Wildtype,geneNamesLocal[i]]
  
  names(x)<-NULL
  names(y)<-NULL
 # newList<- list("x"=as.data.frame(x),"y"=as.data.frame(y))
# wilcox.test(x ~ y, data=as.data.frame(newList)) 
 # t.test(x,y)
  z=append(z,wilcox.test(x,y)$p.value)
  if(geneNamesLocal[i]=="CCL11")
  {
    print("this is x")
    print(x)
    print("this is y")
    print(y)
    print("this is the result of test")
    print(wilcox.test(x,y)$p.value)
    
  }
}
p.adjust(z, method = "BH")
names(z)<-geneNamesLocal
return(z)
}
#---------------------------------------------------------------------------------------
plotdetails= function(MutationStatus,SignificantGenes,filepath)
{
  namesSignificantGenes=names(SignificantGenes)
  pdf(filepath)
  for(i in 1:length(SignificantGenes))
  {
    x=table(MutationStatus[,namesSignificantGenes[i]])
    x=x[-c(1)]
    y=length(names(x))
    if(y==1)
    {barplot(x, col=c(1),main=namesSignificantGenes[i])}
    else if(y==2)
    {barplot(x, col=c(1,2),main=namesSignificantGenes[i])}
    else
    {barplot(x, col=c(1,2,3),main=namesSignificantGenes[i])}  
  }
  dev.off()
  
}
#--------------------------------------------------------------------------------------
plotGEMutatantWildtype=function(GEBinaryMatrixlocal,logGeneExpressionMatrixWithoutBatchlocal,SignificantGeneslocal,filepath)
{
namesSignificantGenes=names(SignificantGeneslocal)
pdf(filepath)
for(i in 1:length(SignificantGeneslocal))
{
  MutantWildtypelocal=getMutantWildtypeSamples(namesSignificantGenes[i],GEBinaryMatrixlocal)
  # MutantBatches=extractBatchesClasses(MutantWildtype$Mutant)
  # WildtypeBatches=extractBatchesClasses(MutantWildtype$Wildtype)
  xlocal=logGeneExpressionMatrixWithoutBatchlocal[MutantWildtypelocal$Mutant,namesSignificantGenes[i]]
  ylocal=logGeneExpressionMatrixWithoutBatchlocal[MutantWildtypelocal$Wildtype,namesSignificantGenes[i]]
  
  mutationstatuslocal1=rep(c("mutant"),length(xlocal))
  mutationstatuslocal2=rep(c("wildtype"),length(ylocal))
  MutationStatus=c(mutationstatuslocal1,mutationstatuslocal2)
  
  GeneExpression=c(xlocal,ylocal)
  names(GeneExpression)<-NULL
  GeneName=rep(c(namesSignificantGenes[i]),length(MutationStatus))
  data=data.frame(GeneName,MutationStatus,GeneExpression)
  
  print(ggplot(data, aes(x=GeneName, y=GeneExpression, fill=MutationStatus)) + 
    geom_boxplot())
}
dev.off()
}