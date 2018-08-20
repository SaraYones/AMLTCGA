library(survminer)
library(survival)
require(Ckmeans.1d.dp)
library(CoRegNet)
library(Boruta)
library("data.table")
library(biomaRt)
library(HGNChelper)
install.packages("refGenome")
library(refGenome)
library(data.table)
library('TCGAbiolinks')
library(edgeR)
library('maftools')
library(devtools)
install_github("mategarb/R.ROSETTA")

library(R.ROSETTA)

install.packages("xlsx")
install.packages("glmnet", repos = "http://cran.us.r-project.org")
library(glmnet)
library(xlsx)
library(arules)
#biocLite("RTCGA.clinical")
#Second Pipline

LAMLClinical=read.table(file="/Users/sarayones/Downloads/laml_tcga_pub/data_clinical_patient.txt",header = T,sep="\t")

#LAMLClinical=read.table(file="TCGA_AML_NEJM_2013_SupplementalTable02.xls",quote = "\"",skipNul = TRUE)


LAMLClinical=read.xlsx(file="TCGA_AML_NEJM_2013_SupplementalTable02.xlsx", 1)

#LAMLSurvival=LAMLClinical[,"OS_MONTHS"]
#event Free survival
LAMLSurvival=LAMLClinical[which(as.character(LAMLClinical$TCGA.Patient.ID) %in% SamplesCNV),"EFS.Months"]
LAMLGE=LAMLClinical[which(as.character(LAMLClinical$TCGA.Patient.ID) %in% SamplesGE),"TCGA.Patient.ID"]

#If I want to  plot Survival Analysis using Kaplan meir plots and not assuming that an event is censored
SurvObj<-Surv(LAMLSurvival)

km.as.one <- survfit(SurvObj ~ 1, data = LAMLClinical)
plot(km.as.one)

Clusters=Ckmeans.1d.dp(LAMLSurvival)
medianLAML=median(LAMLSurvival)
decision=NULL
decision=sapply(LAMLSurvival, function(x) if(x<=6) append(decision,c("Low")) else append(decision,c("High")))
decisionMedian=sapply(LAMLSurvival, function(x) if(x<=median(LAMLSurvival)) append(decision,c("Low")) else append(decision,c("High")))
#Create matrix for GE for all genes
#---------------------------------------------------------------------------------------
#features=paste(rep(significant_genes,4),c("-mutation","-cnv","-GE","-fusion"),sep="")
#featuresGE=paste(rep(colnames(t(AllGenesGE)),"-GE"),sep="")
#row.names(AllGenesGE)<-tempGenes
AllGenesGE=t(AllGenesGE[(which(!grepl("\\?", tempGenes))),])
#Remove all the genes that begin with ? because these are unkown probes
colnames(AllGenesGE)=tempGenes[which(!grepl("\\?", tempGenes))]
colnames(AllGenesGE)=paste(colnames(AllGenesGE),"-GE",sep="")
row.names(AllGenesGE)<-SampleGE
#rownames(AllGenesGE)=tempGenes
#Some of the gene have multiple probes so the here i calculate the mean for all the genes sharing a common name
AllGenesGE=t(apply(AllGenesGE, 1, function(x) tapply(x, colnames(AllGenesGE), mean)))
AllGenesGEDiscretized=discretizeExpressionData(AllGenesGE,standardDeviationThreshold=1)

#---------------------------------------------------------------------------------------
#Create A CNV matrix for the classifier

CNVMatrixClassifier=returnCNVMatrixClassifier(TRUE)
returnCNVMatrixClassifier=function(focal)
{
  if(focal==TRUE)
  {
CNVMatrixClassifier=as.data.frame(genesCNVlocal$FocalGenes)
  }
  else
  {
    CNVMatrixClassifier=as.data.frame(genesCNVlocal$State)
  }
    
CNVMatrixClassifier=CNVMatrixClassifier[,c(-1,-2,-3)]
rownamesCNVMatrixClassifier=gsub("TCGA_AB_([0-9]{4})_03[A-Z]_01[A-Z]_[0-9]{4}_[0-9]{2}","\\1",colnames(CNVMatrixClassifier))
#CNVMatrixClassifier=CNVMatrixClassifier[,SampleGE]
colnames(CNVMatrixClassifier)=genesCNVlocal$State$Gene.Symbol
CNVMatrixClassifier=t(CNVMatrixClassifier)
CNVMatrixClassifier=as.data.frame(CNVMatrixClassifier)
#CNVMatrixClassifier=as.data.frame.numeric(CNVMatrixClassifier)
CNVMatrixClassifier=as.data.frame(lapply(CNVMatrixClassifier, as.numeric))
rownames(CNVMatrixClassifier)=rownamesCNVMatrixClassifier
colnames(CNVMatrixClassifier)=colnamesCNVMatrixClassifier
if(focal==TRUE)
{
  CNVMatrixClassifier[CNVMatrixClassifier<0]<-1
  CNVMatrixClassifier[CNVMatrixClassifier>0]<-1
  
}
return(CNVMatrixClassifier)
}
#CNVMatrixClassifier=CNVMatrixClassifier[,SampleGE]

#colnames(CNVMatrixClassifier)=paste(colnames(CNVMatrixClassifier),"-CNV",sep="")
#---------------------------------------------------------------------------------------
#Create A mutation matrix for the classifier
MutationMatrixClassifier=returnMutationMatrixClassifer(mutationStatePerSample,AllGenesMutation,SampleGE,"TCGA-AB-([0-9]{4})-03[A-Z]-01[A-Z]-[0-9]{4}-[0-9]{2}")
returnMutationMatrixClassifer=function(mutationStatePerSample,AllGenesMutation,SampleGE,regexsamples)
{
MutationMatrixClassifier<-matrix(0,nrow=length(SampleGE),ncol=length(AllGenesMutation))
row.names(MutationMatrixClassifier)<-SampleGE
colnames(MutationMatrixClassifier)<-AllGenesMutation

#Create A Mutation Matrix for the classifier
for (i in 1: length(AllGenesMutation))
{
  
  #transform the result of Mutated Samples Per genes into matrix so that I can access it easily
  
  #MutatedSamplesPerGeneResultsloop=as.matrix(MutatedSamplesPerGeneResultslocal[i])
  #print(MutatedSamplesPerGeneResultslocal[[i]])
  #print(i)
  MutatedSamplesPerGeneResultslocal=genesToBarcodes(mutationStatePerSample,genes=AllGenesMutation)
  MutatedSamplesPerGeneResultslocal=as.character(MutatedSamplesPerGeneResultslocal[[i]]$Tumor_Sample_Barcode)
 # print(MutatedSamplesPerGeneResultslocal[[i]]$Tumor_Sample_Barcode)
  MutatedSamplesPerGeneResultslocal=gsub(regexsamples,"\\1",  MutatedSamplesPerGeneResultslocal)
   
    # print(i)
     #Apply an Or operation with the already existing value inside the matrix
 for(j in 1:length( MutatedSamplesPerGeneResultslocal))
   {
#      print("Hello ")
      #Because we are interested of samples which have mutation information and Structural variants information
    if(MutatedSamplesPerGeneResultslocal[j] %in% SampleGE )
     {
        
       MutationMatrixClassifier[j,AllGenesMutation[i]]=1
    }
  }
    
}
return(MutationMatrixClassifier)

}
#-----------------------------Building Data structures for expirements------------------------------------------
#Find significant genes 
colnames(MutationMatrixClassifier)=paste(colnames(MutationMatrixClassifier),"-Mutation",sep="")
BinaryMatrixlocal=matrix(0,nrow=length(SamplesCNV),ncol=length(significant_genes))
colnames(BinaryMatrixlocal)<-significant_genes
row.names(BinaryMatrixlocal)=SamplesCNVbarcodes
#Create a Detailed Binary Matrix for the reason of 1 values
detailedBinaryMatrixlocal=matrix(0,nrow=length(SamplesCNV),ncol=length(significant_genes))
row.names(detailedBinaryMatrixlocal)=SamplesCNVbarcodes
colnames(detailedBinaryMatrixlocal)<-significant_genes
newList=createBinaryMatrix(significant_genes,BinaryMatrixlocal,detailedBinaryMatrixlocal)
rownames(BinaryMatrixlocal)<-SamplesCNV
rownames(detailedBinaryMatrixlocal)<-SamplesCNV
#We need only the genes that have CNV those that don't have we just remove
significant_genes_cnv_classifier=names(which((colSums(CNVMatrixClassifier != 0) > 0 ) == TRUE))
significant_genes_mutation_classifier=getGeneSummary(mutationStatePerSample)$Hugo_Symbol
significant_genes_classifier=unique(append(significant_genes_cnv_classifier,significant_genes_mutation_classifier))
significant_genes_classifier=unique(append(significant_genes_classifier,significant_genes_Fusion))
#CNVMatrixClassifier[,apply(CNVMatrixClassifier,2,function(x) table(as.numeric(x)[names(table(as.numeric(x)))=="0"])), drop=FALSE]
BinaryMatrixlocal=matrix(0,nrow=length(SamplesCNV),ncol=length(significant_genes_classifier))
colnames(BinaryMatrixlocal)<-significant_genes_classifier
row.names(BinaryMatrixlocal)=SamplesCNVbarcodes
#Create a Detailed Binary Matrix for the reason of 1 values
detailedBinaryMatrixlocal=matrix(0,nrow=length(SamplesCNV),ncol=length(significant_genes_classifier))
row.names(detailedBinaryMatrixlocal)=SamplesCNVbarcodes
colnames(detailedBinaryMatrixlocal)<-significant_genes_classifier
newList=createClassifierBinaryMatrix(significant_genes_classifier,BinaryMatrixlocal,detailedBinaryMatrixlocal)
#Filter any coloumns or any genes in the matrix that are all 0s 
BinaryMatrixlocal=newList$BinaryMatrix[,which(colSums(newList$BinaryMatrix)!=0)]
detailedBinaryMatrixlocal=newList$detailedBinaryMatrix[,which(colSums(newList$BinaryMatrix)!=0)]
BinaryMatrixlocalFiltered=BinaryMatrix[,which(colSums(BinaryMatrix)!=0)]
detailedBinaryMatrixlocalFiltered=detailedBinaryMatrix[,which(colSums(BinaryMatrix)!=0)]
#DataStructs for expirement 7
BinaryMatrixlocalFilteredRows=BinaryMatrixlocalFiltered[which(rowSums(BinaryMatrixlocalFiltered)>5),]
detailedBinaryMatrixlocalFilteredRows=detailedBinaryMatrixlocal[which(rowSums(BinaryMatrixlocalFiltered)>5),]
decisionRowFiltered=decision[which(rowSums(BinaryMatrixlocalFiltered)>5)]

#All Genes filtered Rows 50 % quartile
BinaryMatrixlocalAllGenesFilteredRows=BinaryMatrixlocal[which(rowSums(BinaryMatrixlocal)> as.numeric(quantile(rowSums(BinaryMatrixlocal))[2])),]
BinaryMatrixlocalAllGenesFilteredRows=apply(as.matrix(BinaryMatrixlocalAllGenesFilteredRows), 2, as.character)
detailedBinaryMatrixAllGenesFilteredRows=detailedBinaryMatrixlocal[which(rowSums(BinaryMatrixlocal)> as.numeric(quantile(rowSums(BinaryMatrixlocal))[2])),]
decisionAllGenesFilteredRows=decision[which(rowSums(BinaryMatrixlocal)>as.numeric(quantile(rowSums(BinaryMatrixlocal))[2]))]
#DataStructs for Expirement 9
BinaryMatrixlocal9=BinaryMatrixlocalFiltered[which(rowSums(BinaryMatrixlocalFiltered)>5),]
detailedBinaryMatrixlocal9=detailedBinaryMatrixlocalFiltered[which(rowSums(BinaryMatrixlocalFiltered)>5),]
decisionMedian9=decisionMedian[which(rowSums(BinaryMatrixlocalFiltered)>5)]
BinaryMatrixlocal9=apply(as.matrix(BinaryMatrixlocal9), 2, as.character)

#Compute weights using Generalized Linear Regression model
lassoModel=computeWeights(BinaryMatrixlocalFiltered,detailedBinaryMatrixlocalFiltered,decision)
lassoModeltrial=computeWeights(BinaryMatrixlocalFiltered,detailedBinaryMatrixlocalFiltered,decision)
#To extract only the GE samples from Binary Matrix and detailed Binary Matrix
newList$BinaryMatrix=newList$BinaryMatrix[SampleGE,]
tempMatrix=newList$detailedBinaryMatrix
row.names(tempMatrix)=SamplesCNV
newList$detailedBinaryMatrix=tempMatrix
newList$detailedBinaryMatrix=newList$detailedBinaryMatrix[SampleGE,]

#---------------------------------------Expirements after the modification of the pipeline to include Focal Events-----------------------------------------------
#----Expirement 1--------------------------------- 
BinaryMatrixlocalFiltered=BinaryMatrix[,which(colSums(BinaryMatrix)!=0)]
detailedBinaryMatrixlocalFiltered=detailedBinaryMatrix[,which(colSums(BinaryMatrix)!=0)]
BinaryMatrixlocalFilteredMod1=BinaryMatrixlocalFiltered[which(rowSums(BinaryMatrixlocalFiltered)> as.numeric(quantile(rowSums(BinaryMatrixlocalFiltered))[2])),]
BinaryMatrixlocalFilteredMod1=apply(as.matrix(BinaryMatrixlocalFilteredMod1), 2, as.character)
detailedBinaryMatrixlocalFilteredMod1=apply(as.matrix(BinaryMatrixlocalFilteredMod1), 2, as.character)
detailedBinaryFilteredMod1=detailedBinaryMatrixlocal[which(rowSums(BinaryMatrixlocalFiltered)> as.numeric(quantile(rowSums(BinaryMatrixlocalFiltered))[2])),]
decisionRowFilteredMod1=decision[which(rowSums(BinaryMatrixlocalFiltered)>as.numeric(quantile(rowSums(BinaryMatrixlocalFiltered))[2]))]
detailedBinaryMatrixlocalFiltered=apply(as.matrix(detailedBinaryMatrixlocalFiltered), 2, as.character)
#--------------------------------------------Important --------------------------------------------------
#Before running Rosetta or RMCFS consider filtering the Decision table to include only the genes which have colSums >0
##newList$BinaryMatrix[,(which(colSums(newList$BinaryMatrix)!=0))]
#Also the detailed Binary Matrix should be filtered according to the Binary Matrix
#For any Binary matrix that I am gonna run classifier on remove the coloumns with all 0s 
#If you are going to run each classifer individually (Binary and GE alone) do not remove the samples that are not included in GE (use SamplesCNV)
#--------------------------------------------------------------------------------------------------------
#-------Expirements----------------Check the Expirements log---------------------------------------------
#For Rosetta we need to change the data type of binary values to charachter and Gene expression to integers or numeric
#newList$BinaryMatrix=as.data.frame.character(newList$BinaryMatrix)
#newList$BinaryMatrix=lapply(newList$BinaryMatrix, as.character)
#AllGenesGE=as.numeric(AllGenesGE)
Decisiontable=cbind(BinaryMatrixlocal,AllGenesGE,decision)
DecisiontableBinary=cbind(BinaryMatrixlocal,decision)
detailedDecisiontableBinary=cbind(detailedBinaryMatrixlocal,decision)
#Remove all "LOC" from the DecisiontableBinary---Those are the ones used in the expirements after modication
DecisiontableBinary=DecisiontableBinary[,-which(grepl("^LOC.*",colnames(DecisiontableBinary)))]
detailedDecisiontableBinary=detailedDecisiontableBinary[,-which(grepl("^LOC.*",colnames(detailedDecisiontableBinary)))]
#-------------------------------------------------------------------------------------------------
DecisiontableBinary=apply(as.matrix(DecisiontableBinary), 2, as.character)
DecisiontableGE=cbind(AllGenesGE,decision)
DecisiontableBinaryFiltered=cbind(BinaryMatrixlocalFiltered,decision)
DecisiontableBinaryFilteredRows=cbind.data.frame(BinaryMatrixlocalFilteredRows,decisionRowFiltered)
DecisiontableBinaryFilteredWeights=cbind(lassoModel$BinaryMatrix,decision)
DecisiontableBinaryAllGenesFilteredRows=cbind.data.frame(as.data.frame(BinaryMatrixlocalAllGenesFilteredRows),decisionRowFilteredAllGenesFilteredRows)
DecisiontableBinary9=cbind.data.frame(as.data.frame(BinaryMatrixlocal9),decisionMedian9)
#Trying expirement 4 (Lasso couldnt found that all genes are correlated to the output)
cvfit=cvfit=cv.glmnet(BinaryMatrixlocalFilteredRows,y=as.factor(decisionRowFiltered),family="binomial",alpha=0,lambda=exp(seq(log(0.001), log(5), length.out=100))) 
x=coef(cvfit, s = "lambda.min",exact=FALSE) 
which(x[,1]==0)
#------------------Expirements After Modification---------------------------------------
DecisiontableBinaryMod1=cbind.data.frame(as.data.frame(BinaryMatrixlocalFilteredMod1),decisionRowFilteredMod1)
DecisiontabledetailedBinaryMod1=cbind.data.frame(detailedBinaryFilteredMod1,decisionRowFilteredMod1)
DecisiontabledetailedBinaryMod3=cbind.data.frame(as.data.frame(detailedBinaryMatrixlocalFiltered),as.character(decision))
DecisiontabledetailedBinaryMod3=apply(as.matrix(DecisiontabledetailedBinaryMod3), 2, as.character)
#----------------------------------------------------------------------------------------
#To write the decision tables in files to be able to run it on server
write.csv(Decisiontable, file ="SurvivalClassifierMCFS.csv" )
write.table(DecisiontableBinary, file ="SurvivalClassifierBinaryMCFS" )
write.table(DecisiontableGE, file ="SurvivalClassifierGEMCFS" )
write.table(DecisiontableBinaryFiltered, file ="SurvivalClassifierBinaryFilteredGEMCFS" )
write.table(DecisiontableBinaryFilteredRows, file ="SurvivalClassifierBinaryFilteredRows" )
write.table(DecisiontableBinaryFilteredWeights, file ="SurvivalClassifierBinaryFilteredWeights" )
write.table(DecisiontableBinaryAllGenesFilteredRows, file ="SurvivalClassifierBinaryAllGenesFilteredRows" )
write.table(DecisiontableBinary9, file ="SurvivalClassifierBinaryExp9" )


#------------------Expirements After Modification---------------------------------------
write.table(DecisiontableBinaryMod1, file ="SurvivalClassifierBinaryFilteredMod1" )
write.table(DecisiontabledetailedBinaryMod1, file ="SurvivalClassifierdetailedBinaryFilteredMod1" )
write.table(DecisiontableBinary9, file ="SurvivalClassifierBinaryExp9" )




boruta.train=Boruta(decision~., data = as.data.frame(DecisiontableBinaryFilteredRowa), doTrace = 2)
#This is because MCFS doesnt like the | so we replaced it with -


orignialColNames=updatedcolnames$orinigalColNames
updateColNames=function(DecisiontableBinary)
{
originalColNames=colnames(DecisiontableBinary)
colnames=gsub("\\|", ".",colnames(DecisiontableBinary))
#colnames=gsub("\\-", ".",colnames(DecisiontableBinary))
colnames(DecisiontableBinary)<-colnames

newList<-list("Decisiontable"=DecisiontableBinary,"originialColNames"=originalColNames)
return(newList)
}
#Extract features for Decision table with weights
Features=FilterFeatures("/Users/sarayones/Desktop/Projects/AMLProject/AMLTCGA/outputclassifier/outputBinaryFilteredWeights/outputBinaryFilteredWeights__RI.csv",100)
transform=apply(DecisiontableBinaryFilteredWeights[,Features], 2, as.numeric)
resultRosettaWeights=rosetta(cbind.data.frame(transform,decision),classifier="StandardVoter",discrete=FALSE , discreteMethod="EqualFrequencyScaler",clroc="decision")
recalculatedRules=recalculateRules(cbind.data.frame(transform,decision),resultRosettaWeights$main)

#------Running Expirement 5-----------------------------------------------------
FeaturesFilteredRows=FilterFeatures("/Users/sarayones/Desktop/Projects/AMLProject/AMLTCGA/outputclassifier/outputBinaryFilteredRows/outputBinaryFilteredRows__RI.csv",40)
resultRosettaFilteredRows=rosetta(cbind.data.frame(DecisiontableBinaryFilteredRows[,FeaturesFilteredRows],decisionRowFiltered),classifier="StandardVoter",discrete=TRUE)

updatedcolnames=updateColNames(DecisiontableBinaryFilteredRows) 
DecisionBinaryFilteredRows=updatedcolnames$Decisiontable

recalculatedRules5=recalculateRules(cbind.data.frame(DecisiontableBinaryFilteredRows[,FeaturesFilteredRows],decisionRowFiltered),resultRosettaFilteredRows$main,discrete=TRUE)

#---------------------Running Expirement 6-------------------------------------------------------------
resultRosettaFilteredRowsPip1=rosetta(cbind.data.frame(DecisiontableBinaryFilteredRows[,names(FilteredGeneList)],decisionRowFiltered),classifier="StandardVoter",discrete=TRUE)
resultRosettaFilteredRowsPip1=rosetta(cbind.data.frame(DecisiontableBinaryFilteredRows[,names(sort(FilteredGeneList,decreasing=TRUE)[1:100])],decisionRowFiltered),classifier="StandardVoter",discrete=TRUE)
#-----------------------Running Expirement 8----------------------
FeaturesAllGenesFilteredRows=FilterFeatures("/Users/sarayones/Desktop/Projects/AMLProject/AMLTCGA/outputclassifier/outputBinaryAllGenesFilteredRows/outputBinaryAllGenesFilteredRows__RI.csv",100)

updatedcolnames=updateColNames(DecisiontableBinaryAllGenesFilteredRows) 
DecisiontableBinaryAllGenesFilteredRows=updatedcolnames$Decisiontable

resultRosettaFilteredAllGenesFilteredRows=rosetta(as.data.frame(apply(cbind.data.frame(as.data.frame(DecisiontableBinaryAllGenesFilteredRows[,FeaturesAllGenesFilteredRows]),decisionAllGenesFilteredRows ),2,as.character)) ,classifier="StandardVoter",discrete=TRUE)
recalculatedRules5=recalculateRules(as.data.frame(apply(cbind.data.frame(as.data.frame(DecisiontableBinaryAllGenesFilteredRows[,FeaturesAllGenesFilteredRows]),decisionAllGenesFilteredRows ),2,as.character)),resultRosettaFilteredAllGenesFilteredRows$main,discrete=TRUE)

#----------------------Running Expirement 9-----------------------------------------------------------
FeaturesFilteredExp9=FilterFeatures("/Users/sarayones/Desktop/Projects/AMLProject/AMLTCGA/outputclassifier/outputBinaryExp9/outputBinaryExp9__RI.csv",10)

updatedcolnames=updateColNames(DecisiontableBinary9) 
DecisiontableBinary9=updatedcolnames$Decisiontable
orignialColNames=updatedcolnames$orinigalColNames

resultRosettaFilteredExp9=rosetta(as.data.frame(apply(cbind.data.frame(as.data.frame(DecisiontableBinary9[,FeaturesFilteredExp9]),decisionMedian9 ),2,as.character)) ,classifier="StandardVoter",discrete=TRUE)
recalculatedRules9=recalculateRules(cbind.data.frame(apply(cbind.data.frame(as.data.frame(DecisiontableBinary9[,FeaturesFilteredExp9]),decisionMedian9 ),2,as.character)),resultRosettaFilteredExp9$main,discrete=TRUE)


#Extract features with decision table BinaryFiltered
Features=FilterFeatures("/Users/sarayones/Desktop/Projects/AMLProject/AMLTCGA/outputclassifier/outputBinaryFiltered/outputBinaryFiltered__RI.csv",25)
hello=apply(DecisiontableBinaryFiltered[,Features], 2, as.character)
resultRosetta=rosetta(as.data.frame(hello,clroc="decision"))
#-----------------------Expirements After Modification------------------------------
#--------------------------1--------------------------------------------------------
FeaturesFilteredMod1=FilterFeatures("/Users/sarayones/Desktop/Projects/AMLProject/AMLTCGA/outputclassifier/Modification/outputBinaryModExp1/outputBinaryModExp1__RI.csv",5)

updatedcolnames=updateColNames(DecisiontableBinaryMod1) 
DecisiontableBinaryMod1=updatedcolnames$Decisiontable

resultRosettaFilteredMod1=rosetta(as.data.frame(apply(cbind.data.frame(as.data.frame(DecisiontableBinaryMod1[,FeaturesFilteredMod1]),decisionRowFilteredMod1 ),2,as.character)) ,classifier="StandardVoter",underSample = TRUE,discrete=TRUE)
recalculatedRulesMod1=recalculateRules(as.data.frame(apply(cbind.data.frame(as.data.frame(DecisiontableBinaryMod1[,FeaturesFilteredMod1]),decisionRowFilteredMod1 ),2,as.character)),resultRosettaFilteredMod1$main,discrete=TRUE)

#------------------------2---------------------------------------------------------
FeaturesFilteredMod2=FilterFeatures("/Users/sarayones/Desktop/Projects/AMLProject/AMLTCGA/outputclassifier/Modification/outputdetailedBinaryModExp2/outputdetailedBinaryModExp2__RI.csv",50)

updatedcolnames=updateColNames(DecisiontabledetailedBinaryMod1) 
DecisiontabledetailedBinaryMod1=updatedcolnames$Decisiontable

resultRosettaFilteredMod2=rosetta(as.data.frame(apply(cbind.data.frame(as.data.frame(DecisiontabledetailedBinaryMod1[,FeaturesFilteredMod2]),decisionRowFilteredMod1 ),2,as.character)) ,classifier="StandardVoter",discrete=TRUE)
recalculatedRulesMod2=recalculateRules(as.data.frame(apply(cbind.data.frame(as.data.frame(DecisiontableBinaryMod1[,FeaturesFilteredMod1]),decisionRowFilteredMod1 ),2,as.character)),resultRosettaFilteredMod1$main,discrete=TRUE)

#---------------------------------Running Borouta on Expirement 2------------------------------------------------
resultRosettaFilteredMod2=rosetta(as.data.frame(apply(as.data.frame(DecisiontableBinaryMod1 ),2,as.character) ),classifier="StandardVoter",discrete=TRUE)

boruta.train=Boruta(decision~., data = as.data.frame(DecisiontableBinaryFilteredRowa), doTrace = 2)

#---------------------3-----------------------------------------------------------
FeaturesFilteredMod3=FilterFeatures("/Users/sarayones/Desktop/Projects/AMLProject/AMLTCGA/outputclassifier/Modification/outputBinaryModExp3/outputBinaryModExp3__RI.csv",40)

updatedcolnames=updateColNames(DecisiontableBinaryFiltered) 
DecisiontableBinaryFiltered=updatedcolnames$Decisiontable

resultRosettaFilteredMod3=rosetta(as.data.frame(apply(cbind.data.frame(as.data.frame(DecisiontableBinaryFiltered[,FeaturesFilteredMod3]),decision ),2,as.character)) ,classifier="StandardVoter",discrete=TRUE)
resultRosettaFilteredModGenetic3=rosetta(as.data.frame(apply(cbind.data.frame(as.data.frame(DecisiontableBinaryFiltered[,FeaturesFilteredMod3]),decision ),2,as.character)) ,classifier="StandardVoter",discrete=TRUE,reducer='Genetic')

#---------------------4-----------------------------------------------------------
FeaturesFilteredMod4=FilterFeatures("/Users/sarayones/Desktop/Projects/AMLProject/AMLTCGA/outputclassifier/Modification/outputBinaryModExp4/outputBinaryModExp4__RI.csv",400)

updatedcolnames=updateColNames(DecisiontableBinary) 
DecisiontableBinary=updatedcolnames$Decisiontable

DecisiontableBinarytemp=as.data.frame((apply(cbind.data.frame(as.data.frame(DecisiontableBinary[,FeaturesFilteredMod4]),decision ),2,as.character)))

#resultRosettaFilteredMod5=rosetta(as.data.frame(DecisiontableBinarytemp[,-which(colnames(DecisiontableBinarytemp)=="TP53")]),classifier="StandardVoter",discrete=TRUE,reducer = 'Genetic')
resultRosettaWithUSMod5=rosetta(DecisiontableBinarytemp,classifier="StandardVoter",discrete=TRUE,reducer = 'Genetic',underSample = TRUE ,underSampleNum = 100)
resultRosettaWithoutUSMod5=rosetta(DecisiontableBinarytemp,classifier="StandardVoter",discrete=TRUE,reducer = 'Genetic')

resultRosettaWithoutUSMod4=rosetta(as.data.frame(DecisiontableBinarytemp),classifier="StandardVoter",discrete=TRUE)
resultRosettaWithUSMod4=rosetta(as.data.frame(DecisiontableBinarytemp),classifier="StandardVoter",discrete=TRUE,underSample = TRUE, underSampleNum = 500)

#-----------------------Reporting-------------------------------Filter Rules-------------------------
filterResultRosettaWithoutUSMod5=filterRules(resultRosettaWithoutUSMod5)
filterResultRosettaWithoutUSMod4=filterRules(resultRosettaWithoutUSMod4)


filterRules=function(resultRosetta)
{
  #Choose only low EFS
  filterResultRosetta=resultRosetta$main[resultRosetta$main$DECISION=="Low",]
  
  #choose only rules where length >1 (to filter the rules with interactions)
  
  filterResultRosetta=filterResultRosetta[which(lapply(as.list(strsplit(filterResultRosetta$CUTS_COND, ",")),length)>1),]
  
  #Keep only those rules that have atleast 2 1s
  filterResultRosetta=filterResultRosetta[which(lapply(lapply(as.list(strsplit(filterResultRosetta$CUTS_COND, ",")),table),function(x){x[["1"]]})>=2),]
  
  return(filterResultRosetta)
  
  
}

genes=NULL
genes=readGeneDB("gencode.v20.annotation.gtf")


filterResultRosettaWithoutUSMod5=AnnotateRules(filterResultRosettaWithoutUSMod5)
filterResultRosettaWithoutUSMod4=AnnotateRules(filterResultRosettaWithoutUSMod4)

#Get only the chromosome part of the rules in a list of lists 
temp=lapply(filterResultRosettaWithoutUSMod5$CHROMOSOMES,function(x) {return (unlist(strsplit(as.character(x),split=",")))})

#Distinguish between the rules that have features with all chromosomes == chr21 and not
filter=which(unlist(lapply(temp,function(x) {if(length(x) ==length(which(x %in% "chr21"))) return (FALSE) else return (TRUE)})))
View(filterResultRosettaWithoutUSMod5[filter,])
ReportResults(filterResultRosetta[filter,],DecisiontableBinary,detailedDecisiontableBinary,"outputclassifier/ReportedResultsLinda/GeneticAlgorithmResults/")
ReportResults(filterResultRosettaWithoutUSMod4,DecisiontableBinary,detailedDecisiontableBinary,"outputclassifier/ReportedResultsLinda/JohnsonAlgorithmResults/AllRules/")
ReportResults(filterResultRosettaWithoutUSMod4,DecisiontableBinary,detailedDecisiontableBinary,"outputclassifier/ReportedResultsLinda/JohnsonAlgorithmResults/CombinationRules/")
#trial
ReportResults=function(filterResultRosetta,DecisiontableBinary,detailedDecisiontableBinary,path)
{
write.csv(filterResultRosetta, file =paste(path,"outputRules.csv",sep=""))
write.csv(DecisiontableBinary,file =paste(path,"DecisiontableBinary.csv",sep=""))
write.csv(detailedDecisiontableBinary,file =paste(path,"detailedDecisiontableBinary.csv",sep=""))
}


AnnotateRules=function(filterResultRosetta)
{
  
  geneslocal=genes[genes$source=="HAVANA"]
  chromosomes=matrix(0,nrow=dim(filterResultRosetta)[1],ncol=1)
  chr=NULL
for(i in 1: length(filterResultRosetta$FEATURES))
{ 
  features=as.list(strsplit(as.character(filterResultRosetta$FEATURES[[i]]), ","))
  features=features[[1]]
  # features=gsub("\\.", "-",features)
  features=unlist(lapply(features, function(x) gsub("\\.", "-",x)))
 temp=as.matrix(as.data.frame(strsplit(features,"\\|"))[1,])
 colnames(temp)=NULL
 features=temp[1,]
 temp=cbind(genes$chr[which(genes$gene_name %in% features)],genes$gene_name[which(genes$gene_name %in% features)])
 print(i)
 if(dim(temp)[1]==1)
   {
   print("hi")
   temp=temp[match(features,temp[2])] 
   chr=paste(temp[1],collapse = ",")
 }else if(dim(temp)[1]>1){
     print("hello")
   if(length(unique(temp[,2]))==1){
     temp=temp[match(features,unique(temp[,2]))] 
     chr=paste(temp[1],collapse = ",")
   }else{
   temp=temp[match(features,unique(temp[,2])),]
   chr=paste(temp[,1],collapse = ",")
   }
 }else {chr=NA}

 #chr=paste(genes$chr[which(genes$gene_name %in% features[[1]])],collapse = ",")

 chromosomes[i]=chr
}
  filterResultRosetta=cbind(FEATURES=filterResultRosetta$FEATURES,CUTS_COND=filterResultRosetta$CUTS_COND,DECISION=filterResultRosetta$DECISION,SUPP_LHS=filterResultRosetta$SUPP_LHS,SUPP_RHS=filterResultRosetta$SUPP_RHS,CHROMOSOMES=chromosomes)
  filterResultRosetta=as.data.frame(filterResultRosetta)
  #filterResultRosetta=filterResultRosetta[length]
  colnames(filterResultRosetta)=c("FEATURES", "CUTS_COND", "DECISION","SUPP_LHS","SUPP_RHS","CHROMOSOMES")
  return(filterResultRosetta)
  
  
}
#lapply(lapply(as.list(strsplit(filterResultRosetta$CUTS_COND, ",")),table),function(x){x[["1"]]})

readGeneDB=function(file)
{
#genes <- fread("gencode.v20.annotation.gtf")
genes <- fread(file)
setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )
genes$gene_name <- unlist(lapply(genes$attributes, extract_attributes, "gene_name"))

# [optional] focus, for example, only on entries of type "gene", 
# which will drastically reduce the file size
genes <- genes[type == "gene"]
return(genes)

}



FilterFeatures=function(file,numberOfFeatures)
{
  
  features=fread(file)
  features=features[seq(1,numberOfFeatures),"attribute"]
  return(features$attribute)
  
}


#GEBinaryMatrixlocal=newList$BinaryMatrix[SampleGE,]
#GEdetailedBinaryMatrixlocal=newList$detailedBinaryMatrix[SampleGE,]
#GeneExpressionClassifier=t(logGeneExpressionMatrixWithoutBatch)
#GeneExpressionClassifier=GeneExpressionClassifier[,-which(tempGenes == "?")]

#------------------------------------------------------------------------------------------------
#Decisiontable=cbind(AllGenesGE,CNVMatrixClassifier,MutationMatrixClassifier,decision)

#-------------------------------Calculating the weights for each gene and updating weights---------------------------
computeWeights=function(BinaryMatrixlocal,detailedBinaryMatrixlocal,decision)
{
  removeindices=NULL;
  #For debugging 
 # BinaryMatrixlocal1=BinaryMatrixlocalFiltered
#detailedBinaryMatrixlocal1=detailedBinaryMatrixlocalFiltered
 for (i in 1: dim(BinaryMatrixlocal)[2])
  # for (i in 1: 1)
  { 
    # 7 because CNV,SNV and Fusion ,CNV-Fusion SNV-CNV SNV-Fusion ,CNV-SNV and Fusion
    weightsMatrix=matrix(0,nrow=dim(BinaryMatrixlocal)[1],ncol=7)
    print("hii")
    colnames(weightsMatrix)<-c("CNV","SNV","Fusion","CNVFusion","CNVSNV","SNVFusion" ,"CNVSNVFusion")
    rownames(weightsMatrix)<-rownames(BinaryMatrixlocal1)
    count=0
    #Number of samples in the matrix
   for (j in 1:dim(BinaryMatrixlocal)[1])
     
    {
      #print("hello")
      values=as.character(detailedBinaryMatrixlocal1[j,i])
      if (values!="0")
      {
        
        #StrSplit
        if(grepl(paste("^",values,"$",sep = ""),"CNV"))
        {
         # print("CNV")
        #  print(j)
          weightsMatrix[j,"CNV"]=1
          count=count+1
        }
        else if(grepl(paste("^",values,"$",sep = ""),"SNV"))
        {# print("SNV")
        #  print(j)
          weightsMatrix[j,"SNV"]=1
          count=count+1
        }
        else if(grepl(paste("^",values,"$",sep = ""),"Fusion"))
        {print("Fusion")
         # print(j)
          weightsMatrix[j,"Fusion"]=1
          count=count+1
        }
        else if(grepl(paste("^",values,"$",sep = ""),"CNV -Fusion"))
        {
          #print("CNV-Fusion")
          #print(j)
          weightsMatrix[j,"CNVFusion"]=1
          print("CNV-FUSION finish")
          count=count+1
        }
        else if(grepl(paste("^",values,"$",sep = ""),"SNV -Fusion"))
        {#print("SNV-Fusion")
          weightsMatrix[j,"SNVFusion"]=1
          count=count+1
        }
        else if(grepl(paste("^",values,"$",sep = ""),"CNV-SNV"))
        {#print("CNV-SNV")
          #print(j)
          #print("CNV-SNV")
          weightsMatrix[j,"CNVSNV"]=1
          count=count+1
        }
        
        else if(grepl(paste("^",values,"$",sep = ""),"CNV-SNV -Fusion"))
        {#print("CNV-SNV -Fusion")
          weightsMatrix[j,"CNVSNVFusion"]=1
          count=count+1
        }
        
        
        
      }
   }
  
  print("Which gene is it")
  print(colnames(BinaryMatrixlocal)[i])
  
  print("HEYYYYYYYYYYYYYYYYYYYY")
  if(count>3)
  {
  print(colnames(BinaryMatrixlocal)[i])
  
  print(weightsMatrix)
  
  weightsMatrix=Matrix(weightsMatrix, sparse=TRUE)
  #View(weightsMatrix)
    
 #print(weightsMatrix)  
#cvfit = cv.glmnet(as.matrix(as.data.frame(lapply(as.data.frame(weightsMatrix), as.double))),y=as.factor(decision),family="binomial" ,alpha=0,
 #                 lambda=exp(seq(log(0.001), log(5), length.out=100)))
  
  cvfit=cv.glmnet(weightsMatrix,y=as.factor(decision),family="binomial",alpha=0,lambda=exp(seq(log(0.001), log(5), length.out=100))) 
  #cvfit=cv.glmnet(weightsMatrix,y=as.factor(decision),family="binomial")
  
  x=coef(cvfit, s = "lambda.min",exact=FALSE) 
print("Coefficents are ")
print(x)
#coef(cvfit,)
#  x=coef(cvfit)
  
for (j in 1:dim(BinaryMatrixlocal)[1])
  
{
  temp=0
  values=as.character(detailedBinaryMatrixlocal[j,i])
  print("detailedBinaryMatrix is")
  print(values)
  if (values!="0")
  {
    
    #StrSplit
    #if(grepl(paste("^",values,"$",sep = ""),"CNV"))
    if(grepl("CNV",values))
    {
      print("CNV")
      temp=temp+x["CNV",]
    }
    #else if(grepl(paste("^",values,"$",sep = ""),"SNV"))
    if(grepl("SNV",values))
    { print("SNV")
      temp=temp+x["SNV",]
    }
   if(grepl(paste("^",values,"$",sep = ""),"Fusion"))
    {print("Fusion")
      temp=temp+x["Fusion",]
    }
    if(grepl(paste("^",values,"$",sep = ""),"CNV -Fusion"))
    {
      print("CNV-Fusion")
      print("hello")
      temp=temp+x["CNVFusion",]
   
    }
     if(grepl(paste("^",values,"$",sep = ""),"SNV -Fusion"))
    {print("SNV-Fusion")
      temp=temp+x["SNVFusion",]
    }
    if(grepl(paste("^",values,"$",sep = ""),"CNV-SNV"))
    {print("CNV-SNV")
      temp=temp+x["CNVSNV",]
    }
    
    if(grepl(paste("^",values,"$",sep = ""),"CNV-SNV -Fusion"))
    {print("CNV-SNV -Fusion")
      temp=temp+x["CNVSNVFusion",]
    }
    print("temp is ")
    print(temp)
    BinaryMatrixlocal[j,i]=temp+x["(Intercept)",]
}
 
}
  }
  
  else
  {
    
    removeindices=append(removeindices,i)
  }

    }  
  if(!is.null(removeindices))
  newList<-list("BinaryMatrix"=BinaryMatrixlocal[,-removeindices],"detailedBinaryMatrix"=detailedBinaryMatrix[,-removeindices])
  else
    newList<-list("BinaryMatrix"=BinaryMatrixlocal,"detailedBinaryMatrix"=detailedBinaryMatrix)
  
return(newList)
  
}


#-------------------------If we normalize using GTEX we use this script-----------------------------------------------
#---------------------------------------------------------------------------------
#Filter Samples that are whole blood and have batch number not "not _reported" those are the ones that we are going to use for GE normalization

SampleAttributeGtex=fread(file="SampleAttributesGtex.txt",header = T)
SampleAttributeGtex=as.data.frame(SampleAttributeGtex)
SampleAttributeGtex=SampleAttributeGtex[which(SampleAttributeGtex$SMTSD=="Whole Blood"),]
SampleAttributeGtex=SampleAttributeGtex[which(SampleAttributeGtex$SMNABTCH!="not reported"),]

GTEXData=fread(file="All_Tissue_Site_Details.combined.reads.gct",header = T)
GTEXData=as.data.frame(GTEXData)
#GTEXData=as.matrix(GTEXData)
GeneExpressionGTEX=GTEXData[,c(1,2,which(colnames(GTEXData) %in% SampleAttributeGtex$SAMPID))]
GeneExpressionGTEX=t(GeneExpressionGTEX)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#--------------------------------------------------------------------------------------------------------
#We will first create a quantification matrix for all the genes in gtex then remove batch effects
View(head(GTEXData))
#I also need to return the genes matrix (the one that has all the annotation) for this function
significant_genes_Filtered=updateSignificantGenes(file ,significant_genes_classifier,tempGenes)
#Filter out genes more according to the ones in GTEX
significant_genes_Filtered_GTEX=significant_genes_Filtered[which(!(significant_genes_Filtered %in% GeneExpressionGTEX["Description",]))]
View(head(SampleAttributeGtex))
BatchesGTEX=unique(SampleAttributeGtex$SMNABTCH)
logGeneExpressionGTEX=GeneExpressionGTEX
colnames(logGeneExpressionGTEX)=GeneExpressionGTEX["Description",]
logGeneExpressionGTEX=as.data.frame(logGeneExpressionGTEX)
#View(head(logGeneExpressionGTEX))
logGeneExpressionGTEX=logGeneExpressionGTEX[-c(1,2),]
rpkm(t(logGeneExpressionGTEX))

#logGeneExpressionGTEX=logGeneExpressionGTEX[which(rownames(logGeneExpressionGTEX) %in% SampleAttributeGtex$SAMPID),]


#We will then 
#---------------Some Helper script used in updateSignificantGenes-----------------------------------------------------

GeneExpressionGTEX["Description",which(grepl(unique(tempGenes[which(!(grepl("\\?",tempGenes)))]),GeneExpressionGTEX["Description",]))]
hello2=(unlist(lapply(unique(expressionlaml$Hugo_Symbol), function(x) x[!is.na(x)])))

hgnc_swissprot <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','hgnc_id'),filters = 'hgnc_symbol', 
 values = significant_genes, mart = ensembl)          
                        #values = hello2[(which(!(hello2 %in% temp[,3])))]
hgnc_swissprot <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','hgnc_id'),filters = 'hgnc_symbol', 
                        values = checkGeneSymbols(hello[which (!(hello %in% temp[,3]))], unmapped.as.na=TRUE, hgnc.table=NULL)$Suggested.Symbol, mart = ensembl)          

ens <- ensemblGenome()
RefAnnotation= fread(file="gencode.v20.annotation.gtf",header = F,sep="")

RefAnnotation$V1=as.character(RefAnnotation$V1)
RefAnnotation$V2=as.integer(RefAnnotation$V2)
RefAnnotation$V3=as.integer(RefAnnotation$V3)
temp=as.matrix(strsplit( as.character(RefAnnotation$V9) , "\\;" ))
temp=t(do.call("cbind",as.list(temp)))
temp=temp[,5]
temp=do.call(rbind,as.matrix(strsplit(temp , " " )))
Data  <- gsub("[^0-9A-Za-z///' ]","'" ,temp[1,3] ,ignore.case = TRUE)

#--------------------------------------------------------------------------------------------------

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}
genes <- fread("gencode.v20.annotation.gtf")
setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )

# [optional] focus, for example, only on entries of type "gene", 
# which will drastically reduce the file size
genes <- genes[type == "gene"]


#check each gene if it has duplicate using grepl 
#save the name of each gene and when finished remove the gene and

genes$gene_name <- unlist(lapply(genes$attributes, extract_attributes, "gene_name"))

##significant_genes[which(!(significant_genes %in% genes$gene_name))]
##hello=checkGeneSymbols(significant_genes[which(!(significant_genes %in% genes$gene_name))])


# Update Significant Genes function according to TCGA and RefAnnnotation
#Update SignifcantGenes According to Annotation in the REF genome file (To get GE from GTEX later on using ENSMBL IDs) and 
#Those that are actually in the GE files of TCGA

updateSignificantGenes=function(file ,significant_genes_local,tempGenes)
{

  genes <- fread("gencode.v18.genes.bed")
  
  setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )
  
  # [optional] focus, for example, only on entries of type "gene", 
  # which will drastically reduce the file size
  genes <- genes[type == "gene"]
  genes$gene_name <- unlist(lapply(genes$attributes, extract_attributes, "gene_name"))
  genesNotAnnotated=significant_genes[which(!(significant_genes %in% genes$gene_name))]
  genesAnnotated=significant_genes[which(significant_genes %in% genes$gene_name)]
  #AnnotationSignificantGenes$gene_symbols=genesAnnotated
  #AnnotationSignificantGenes$Ensmbl_ID=genes$gene_id[which(genes$gene_name %in% genesAnnotated)]
  #geneAlternativeSymbol=checkGeneSymbols(genesNotAnnotated)
  #We need ENSMBL geneIDs for normalizations
  significantNotInGEAlternative=checkGeneSymbols(genesNotAnnotated)
  #significantNotInGEAlternative=checkGeneSymbols(significant_genes[which(!(significant_genes %in% tempGenes[which(!(grepl("\\?",tempGenes)))]))])
  #Originial Data structure with that has false values of approved
  genesChangedSymbols=significantNotInGEAlternative[which(significantNotInGEAlternative$Approved== 'FALSE'),]
 #get suggested symbols of changed genes
   genesChangedSymbolsChanged=genesChangedSymbols$Suggested.Symbol[which(!is.na(genesChangedSymbols$Suggested.Symbol))]
  #genesChangedSymbolsOriginal=genesChangedSymbols$x[which(!is.na(genesChangedSymbols$Suggested.Symbol))]
 
   #Updated Symbols here to found and not found
   
  genesChangedSymbolsFound=genesChangedSymbolsChanged[which(genesChangedSymbolsChanged %in% genes$gene_name)]
  genesChangedSymbolsNotFound=genesChangedSymbols[which(!(genesChangedSymbols$Suggested.Symbol %in% genesChangedSymbolsFound )),]
  #Originals of the found genes
  genesChangedSymbolsFoundOriginal=genesChangedSymbols[which(genesChangedSymbols$Suggested.Symbol %in% genesChangedSymbolsFound ),]
   #These are the ones to look for in the Annotation for Normal Gene Expression Values
   GenesAnnotatedFinal=append(genesAnnotated,genesChangedSymbolsFound)
#Remove the rest from significant genes   
  significant_genes_updated_annotation=(unique((genesAnnotated),genesChangedSymbolsFoundOriginal$x))
  
  #End of updating significant genes according to annotation. GenesAnnotated Final is the One we are going to use to get the 
  #GE because these are the ones annotated
  genesinGE=significant_genes[which(significant_genes %in% tempGenes[which(!(grepl("\\?",tempGenes)))])]
 
 genesNotinGE=significant_genes[which(!(significant_genes %in% tempGenes[which(!(grepl("\\?",tempGenes)))]))]
 genesNotinGEState=checkGeneSymbols(genesNotinGE)
 genesNotinGEChanged=genesNotinGEState[which(genesNotinGEState$Approved== 'FALSE'),]
 genesNotinGEChangedSymbols= genesNotinGEChanged$Suggested.Symbol[which(!is.na(genesNotinGEChanged$Suggested.Symbol))]
 
 genesNotinGEChangedSymbolsFound=genesNotinGEChangedSymbols[which(genesNotinGEChangedSymbols %in% tempGenes[which(!(grepl("\\?",tempGenes)))])]
 genesNotinGEChangedSymbolsNotFound=genesNotinGEChanged[which(!( genesNotinGEChanged$Suggested.Symbol %in%  genesNotinGEChangedSymbolsFound )),]
 genesNotinGEChangedFoundOriginal=genesNotinGEChanged[which(genesNotinGEChanged$Suggested.Symbol %in% genesNotinGEChangedSymbolsFound ),]
 
 significant_genes_updated_GE=(unique((genesinGE),genesNotinGEChangedFoundOriginal$x))
 significant_genes_Filtered=(unique(significant_genes_updated_annotation, significant_genes_updated_GE))
 
 
  #if(length(genesChangedSymbols[which(genesChangedSymbols %in% tempGenes[which(!(grepl("\\?",tempGenes)))])])==0)
  #{
    
    
  #}
  #hello=append(genesAnnotated,geneAlternativeSymbol$Suggested.Symbol[which(!(is.na(geneAlternativeSymbol$Suggested.Symbol)))])
  
  #what i am going to do is make a map of genes in the significant gene list that are annotated and that are not annotated.
  #Try first with gene names the ones that cannot be found in gene expression files then change to alternative names and search again 
  #whatever is found  will get it from gene expression files and whatever is not found will search in the map if it has annotation or something then we can get 
  #its ENSMBL IDs and search in the other TCGA files
 
 return (significant_genes_Filtered)
}


#x1=significant_genes[which(!(significant_genes %in% tempGenes[which(!(grepl("\\?",tempGenes)))]))]
#x2=significant_genes[which(!(significant_genes %in% cbioportal$Hugo_Symbol[which(!(is.na(cbioportal$Hugo_Symbol)))]))]
#x3=checkGeneSymbols(significant_genes[which(!(significant_genes %in% tempGenes[which(!(grepl("\\?",tempGenes)))]))])
#x4=x3$Suggested.Symbol
#x5=x4[which(!(is.na(x4[which(x3$Approved== 'FALSE')])))]
#x5=

query.exp.hg19 <- GDCquery(project = "TCGA-LAML",
                           data.category = "Gene expression",
                           data.type = "Gene expression quantification",
                           workflow.type="HTSeq - FPKM",
                           legacy= TRUE)

query.exp.hg38 <- GDCquery(project = "TCGA-LAML",
                           data.category = "Transcriptome Profiling",
                           data.type = "Gene Expression Quantification",
                           workflow.type="HTSeq - FPKM")
GDCdownload(query.exp.hg38, method = "api")
expdat <- GDCprepare(query = query.exp.hg38,
                     save = TRUE, save.filename = "exp.rda")
View(expdat)
x1=load(file = "exp.rda")

#trial hello hello hello hello