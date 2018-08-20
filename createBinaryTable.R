  #This script creates the binary table for driver genes from mutation file from MutsigCV and driver genes from CNV files from GISTIC 2.0
#Run the functions first to include them in the memory or enviroment
library(maftools)
readRequiredCNV= function(significantgenesAmp,significantgenesDel,geneState,genefocal){
  allsegmentslocal=NULL
  #Trying github  trying versioning this time
  #For debugging 
  #significantgenesAmp='/Users/sarayones/Desktop/CNVAMLAnalysis/amp_genes.conf_90.txt'
  #significantgenesDel='/Users/sarayones/Desktop/CNVAMLAnalysis/del_genes.conf_90.txt'
  #geneState='/Users/sarayones/Desktop/CNVAMLAnalysis/all_thresholded.by_genes.txt'
  
 genes.Amp=read.table(file = significantgenesAmp,header = T,sep="\t")
 genes.Del=read.table(file = significantgenesDel,header = T,sep="\t")
 
 genes.State=read.table(file = geneState,header = T,sep="\t")
 genes.focal=read.table(file = genefocal,header = T,sep="\t")
 
 newList<- list("Amp"=as.data.frame(genes.Amp),"Del"=as.data.frame(genes.Del),"State"=as.data.frame(genes.State),"FocalGenes"=as.data.frame(genes.focal))
 return(newList)
}


#Keep only the gene symbols or gene names remove everything else    
returnReleventgenes<-function(CNVfile)
{
alteredGenes=NULL;
#-1 because there is an empty coloumn at the end
for(i in 1: dim(CNVfile)[2]-1)
{
  indices=which(grepl("[0-9]{1}\\.[0-9]", levels(CNVfile[,i])))
  indices=append(indices,which(grepl("chr*", levels(CNVfile[,i]))))
  alteredGenes=append(alteredGenes,levels(CNVfile[,i])[-indices])
}
alteredGenes=alteredGenes[alteredGenes != ""]
return(alteredGenes)
}

returnBinaryMatrixCNV<-function(genesCNV)
{
  #For debugging
#genesCNV=genesCNVlocal
amplifiedGenes<-returnReleventgenes(genesCNV$Amp)
deletedGenes<-returnReleventgenes(genesCNV$Del)
#Append the significant genes from both files together
significant_genes=append(amplifiedGenes,deletedGenes)
##StatePerSample=genesCNV$State[(which( genesCNV$State[,1] %in% significant_genes)),]

#extract the non matching pattern from the significant_Genes list that do not match with all the genes in the State matrix
##nonmatching_significantgenes=significant_genes[which(!(significant_genes %in% genesCNV$State[,1]))]
#extract the indices in the State matrix that haave the same "pattern" of the first non matching element in the significant genes
##if(length(nonmatching_significantgenes)!=0)
##{
##for(i in 1: length(nonmatching_significantgenes))
##{
##indices_nonmatching_state_matrix=which((grepl(nonmatching_significantgenes[i], genesCNV$State[,1])))
#append these rows in the State matrix to the StatePerSample that we have already extracted
##StatePerSample=rbind(StatePerSample,genesCNV$State[indices_nonmatching_state_matrix,])
#}
#}
##return(StatePerSample)
return(significant_genes)
}

returnBinaryMatrixMutation<-function(genesCNV)
{
  significantMutatedGenes=read.table(file = '/Users/sarayones/Desktop/MutsigCVFiles/MutSigOutput.sig_genes.txt',header = T,sep="\t")
 mutationStatePerSample=read.table(file = '/Users/sarayones/Desktop/MutsigCVFiles/genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2/genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2.2.13.0.somatic.maf',header = T,as.is=TRUE,sep="\t",fill = TRUE)
  significantMutatedGenes=significantMutatedGenes[significantMutatedGenes$p<=0.05,]
   #which((mutationStatePerSample[,1] %in% significantMutatedGenes[,1]))
 ## mutationStatePerSample=read.maf('/Users/sarayones/Desktop/MutsigCVFiles/genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2.2.13.0.somatic.maf', useAll = TRUE)
  #SamplesMutationbarcode=getSampleSummary(mutationStatePerSample)
 ## mutationStatePerSample1=genesToBarcodes(mutationStatePerSample,genes=x)
   ##return(StatePerSample)
  return(as.character(significantMutatedGenes$gene))
}

getFusiongenes=function(FusionFilepath,RegexSamplesFusion)
{
  FusionGenes=read.table(file = '/Users/sarayones/Desktop/Projects/AMLProject/AMLTCGA/pancanfus-2.txt',header = T,sep="\t")
  FusionGenes=FusionGenes[which(FusionGenes$Cancer=='LAML'),]
  SamplesFusion=sub("AB-([0-9]{4})-03[A-Z]","\\1",FusionGenes$TCGA_barcode)
  newList<- list("FusionGenes"=as.data.frame(FusionGenes),"SamplesFusion"=as.data.frame(SamplesFusion))
  return(newList)
  
  
}
#------------------------Preprocesss the data--------------------------------------------------
genesCNVlocal<-readRequiredCNV('/Users/sarayones/Desktop/CNVAMLAnalysis/amp_genes.conf_90.txt','/Users/sarayones/Desktop/CNVAMLAnalysis/del_genes.conf_90.txt','/Users/sarayones/Desktop/CNVAMLAnalysis/all_thresholded.by_genes.txt','/Users/sarayones/Desktop/CNVAMLAnalysis/focal_data_by_genes.txt')


##StatePerSample<-returnBinaryMatrixCNV(genesCNVlocal)
mutationStatePerSample=read.maf('/Users/sarayones/Desktop/MutsigCVFiles/genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2/genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2.2.13.0.somatic.maf', useAll = TRUE)
                                   
#mutationStatePerSample=genesToBarcodes(mutationStatePerSample,genes=x)

significant_genes_cnv=returnBinaryMatrixCNV(genesCNVlocal)
significant_genes_mutation=returnBinaryMatrixMutation()
#significant_genes_common=getcommongenes("/Users/sarayones/Desktop/MutsigCVFiles/genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2/genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2.2.13.0.somatic.maf","/Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3")
#update that later
#Get fusion genes from Tumor Fusion Gene data portal and add them to the list of significant genes
FusionGenes=as.data.frame(getFusiongenes("/Users/sarayones/Desktop/Projects/AMLProject/AMLTCGA/pancanfus-2.txt"))
significant_genes_Fusion=unique(append(as.character(FusionGenes$FusionGenes.Gene_A),as.character(FusionGenes$FusionGenes.Gene_B)))
significant_genes=unique(append(significant_genes_cnv,significant_genes_mutation))
significant_genes=unique(append(significant_genes,significant_genes_common))
significant_genes=unique(append(significant_genes,significant_genes_Fusion))


#significant_genes=significant_genes_cnv
#Find common samples between CNV data and mutation data
SamplesCNVbarcodes=colnames(genesCNVlocal$State)
indicesCNV=c(1,2,3)
SamplesCNVbarcodes=SamplesCNVbarcodes[-indicesCNV]
SamplesCNV=gsub("TCGA_AB_([0-9]{4})_03A_01D_0756_21","\\1",SamplesCNVbarcodes)

SamplesMutationbarcode=getSampleSummary(mutationStatePerSample)
SamplesMutationbarcode=list(SamplesMutationbarcode[[1]])
SamplesMutationbarcode=unlist(SamplesMutationbarcode)
SamplesMutationbarcode=levels(SamplesMutationbarcode)
SamplesMutation=gsub("TCGA.AB.([0-9]{4}).03[A-Z].01[A-Z].[0-9]{4}.[0-9]{2}","\\1",SamplesMutationbarcode)
IndexCNVnotMutation=which(!(SamplesCNV %in% SamplesMutation))
IndexMutationnotCNV=which(!(SamplesMutation %in% SamplesCNV))
#Remove the indices of the elements in CNV but not in Mutation and vice versa
SamplesMutation=SamplesMutation[-IndexMutationnotCNV]
SamplesCNV=SamplesCNV[-IndexCNVnotMutation]
SamplesFusion=as.character(FusionGenes$SamplesFusion)
FusionGenes=cbind(FusionGenes,SamplesFusion)


#----------------Building the Binary Matrix--------------------------------

#Build the matrix such that col are significant gene names and rows are samples
BinaryMatrix=matrix(0,nrow=length(SamplesCNV),ncol=length(significant_genes))
colnames(BinaryMatrix)<-significant_genes
Regex=paste("",as.character(SamplesCNV),collapse="|",sep="")
#This is if I want to use the row names as the barcodes of CNV or mutation for examples
SamplesCNVbarcodes=SamplesCNVbarcodes[which(grepl(Regex,SamplesCNVbarcodes))]
SamplesMutationbarcode=SamplesMutationbarcode[which(grepl(Regex,SamplesMutationbarcode))]
row.names(BinaryMatrix)=SamplesCNVbarcodes
#Create a Detailed Binary Matrix for the reason of 1 values
detailedBinaryMatrix=matrix(0,nrow=length(SamplesCNV),ncol=length(significant_genes))
row.names(detailedBinaryMatrix)=SamplesCNVbarcodes
colnames(detailedBinaryMatrix)<-significant_genes
#####################################################
#This is if i want to use rownames as cases
#row.names(BinaryMatrix)=SamplesCNV
#Check CNV state and Mutation state then fill the matrix
#Relevent
#Loop for CNV
#Note:there is some kind of a bug here because I found a gene that is in the significant gene list =MIR1244 which has multiple entries in different chromosomes
#in the enesCNVlocal$State$Gene.Symbol matrix so when I run this routine it cannot find this gene in the enesCNVlocal$State$Gene.Symbol
# I suggest using regular expression which(grepl("MIR1244.*",genesCNVlocal$State$Gene.Symbol)) to look for this gene in the state matrix but how about the values there ? it is different on each chromosome
#I suggest first processing  the genesCNVlocal$State then if it finds a gene with mutliple indices we make an or operation for all the values in all the chromosomes and put it in one single entry 
#then use it here instead of as.matrix(genesCNVlocal$State)
#-------------------------------------------------------------------

#I want to check if the gene has one entry in the CNV state matrix or multiple ones like MIR1244  
#This part is for CNV only so we add ones to the locations where there are CNVs
matrices=createBinaryMatrix(significant_genes,BinaryMatrix,detailedBinaryMatrix)
createBinaryMatrix=function(significant_genes,BinaryMatrix,detailedBinaryMatrix)
{
  significant_genes_regex=NULL
  for(i in 1:length(significant_genes))
  {
    
    significant_genes_regex[i]=paste("^",as.character(significant_genes[i]),sep="")
    significant_genes_regex[i]=paste(as.character(significant_genes_regex[i]),"$",sep="")
    significant_genes_regex[i]=paste(as.character(significant_genes_regex[i]),".*",sep="")
  }
  
  #---------------------------------------------------------------------
 #If i want to use the broad CNVs i use the state matrix else if i want to use the focal ones i use the Focal Genes  
  CNVstateMatrix=as.matrix(genesCNVlocal$State)
  #CNVstateMatrix=as.matrix(genesCNVlocal$FocalGenes)
  indexgene=NULL
  for (i in 1: length(significant_genes))
  {
    #Old version   
    # if (significant_genes[i] %in% genesCNVlocal$State$Gene.Symbol)
    #We can keep it this way if("TRUE" %in%  grepl(significant_genes_regex[i],genesCNVlocal$State$Gene.Symbol))
    #Or better check the length of the return value
    #This says
    
    temp=length(which(grepl(significant_genes_regex[i],genesCNVlocal$State$Gene.Symbol)))
    #this is in case we find the gene in the table of threshold file 
    if(temp>=0)
    {
      #This is in case we find it only once
      if(temp==1)
      {
        
        indexgene=which(CNVstateMatrix[,1]==significant_genes[i])
        "ALL single significant genes"
        #  print(significant_genes[i])
      }
      else #This is in case we find it multiple time the exception case (As i Explained above the exception case has a "|chrx" at the end but sometimes the genes
        #That are not found from the first time is not an exception case but they are not even in the table for example a gene 
        #like MED12 if i apply the below regex it will find it MED12L in the file but it is actually not the same gene so I make this else for only the exception case )
      {
        print("Exceptioncase")
        print(significant_genes[i])
        temp1=paste("^",as.character(significant_genes[i]),sep="")
        temp1=paste(as.character(temp1),".*",sep="")
        temp1=paste(as.character(temp1),"$",sep="")
        temp=length(which(grepl(temp1,genesCNVlocal$State$Gene.Symbol)))
        #For debugging
        # if(significant_genes[i]=='MIR1244-2')
        # {
        #  print("Iam in the expceptional case")
        #   print(significant_genes[i])
        # }
        
        if(temp>=1) #This  is in case it matches but it is not an exception case
        {
          # temp2=genesCNVlocal$State$Gene.Symbol[which(grepl(temp,genesCNVlocal$State$Gene.Symbol))]
          #For debugging
          #  if(significant_genes[i]=='MIR1244-2')
          # {
          #   print("Iam in the expceptional case")
          #   print(significant_genes[i])
          # }
          #I check if this gene gotten from the regular expression has a | in it or not . if it has then its an exception case if not then we should leave it with 0's in the Binary table
          temp2= grepl("\\|",genesCNVlocal$State$Gene.Symbol[which(grepl(temp1,genesCNVlocal$State$Gene.Symbol))][1])
          #    print("This is temp2")
          #   print(temp2)
          
        }
        #if(significant_genes[i] %in% significant_genes_cnv)
        # if(significant_genes[i] %in% genesCNVlocal$State$Gene.Symbol)
        #if it is an exception case (temp2 marks that it is definitely an exceptional case)
        if(temp>=1&&temp2==TRUE)
        { print("Exceptioncase")
          print(significant_genes[i])
          indexChromosome=which(as.matrix((genesCNVlocal$Del)==significant_genes[i]),arr.ind=TRUE)
          relevantChromosome=as.matrix(genesCNVlocal$Del)[3,indexChromosome[1,2]]
          relevantChromosome=gsub("(chr[0-9]{1}[0-9]{0,1}):.*","\\1",relevantChromosome)
          indexgene=which(CNVstateMatrix[,1]==paste(as.character(significant_genes[i]),"|",as.character(relevantChromosome),sep = ""))
          
          #print(indexgene)
          # print("Exception case")
        }
      }
      #  if(temp>=1&&significant_genes[i] %in% significant_genes_cnv)
      # if(temp>=1&&significant_genes[i] %in% genesCNVlocal$State$Gene.Symbol)
      if(temp>=1&&!is.null(indexgene))
      {
        #print(significant_genes[i])
        for (j in 1:length(SamplesCNVbarcodes))
          # for (j in 1:length(SamplesCNVbarcodes))
        {
          #The index of the gene inside the CNV matrix to access it with respect to the sample
          #gsub("TCGA_AB_([0-9]{4})_03A_01D_0756_21","\\1",SamplesCNVbarcodes)
          #   print("hi")
          ValueCNV=CNVstateMatrix[indexgene,SamplesCNVbarcodes[j]]
          #If value is -1 or -2 put =1 we will consider any type of CNV at the moment without deffrentiation
          ValueCNV=as.numeric(ValueCNV)
          #     print("hello")
          #For debugging
          #     if(significant_genes[i]=='DNAH10')
          #    { print('this is DNAH10')
          #     print("this is valueCNV")
          #    print(ValueCNV)
          #   }
          # print(ValueCNV)
          #print(SamplesCNVbarcodes[j])
          if(ValueCNV==2| ValueCNV==-2)
          {
            
            #  print(significant_genes[i])
            ValueCNV=1
            detailedBinaryMatrix[SamplesCNVbarcodes[j],significant_genes[i]]='CNV'
            BinaryMatrix[SamplesCNVbarcodes[j],significant_genes[i]]=ValueCNV
          }
          #   print("yes")
          #print(SamplesCNVbarcodes[j])
          #print(significant_genes[i])
          #print(ValueCNV)
          #print("-------------")
          #print("SampleCNVbatcodes[j]")
          #print(SamplesCNVbarcodes[j])
          #print("significant_genes")
          #print(significant_genes[i])
         ## BinaryMatrix[SamplesCNVbarcodes[j],significant_genes[i]]=ValueCNV
          
          
        }
        
        
      }
      indexgene=NULL
    }
    
    #  else
    # {
    #   temp=paste("^",as.character(significant_genes[i]),sep="")
    #    significant_genes_regex[i]=paste(as.character(temp),".*",sep="")
    #   significant_genes_regex[i]=paste(as.character(temp),"$",sep="")
    #    if(length(which(grepl(temp,genesCNVlocal$State$Gene.Symbol))))
    #   {
    #      indexgene=which(grepl(temp,genesCNVlocal$State$Gene.Symbol))
    #   }
    
    
    
    
    
  }
  
  newList=AddMutationToBinaryMatrix(significant_genes,BinaryMatrix,detailedBinaryMatrix)
  newList=AddFusionToBinaryMatrix(newList$BinaryMatrix,newList$detailedBinaryMatrix)
  
  return(newList)  
}
#------------------------------------------------------------------------------------
#Create simple classifier binary matrix and this I will use only with the survival classifier

createClassifierBinaryMatrix=function(significant_genes,BinaryMatrix,detailedBinaryMatrix)
{
  
  
  CNVstateMatrix=as.matrix(genesCNVlocal$State)
  indexgene=NULL
  for (i in 1: length(significant_genes))
  {
    
    
    indexgene=which(CNVstateMatrix[,1]==significant_genes[i])
    
    if(length(indexgene)!=0)
    {
    for (j in 1:length(SamplesCNVbarcodes))
      # for (j in 1:length(SamplesCNVbarcodes))
    {
      print(significant_genes[i])
      ValueCNV=CNVstateMatrix[indexgene,SamplesCNVbarcodes[j]]
      #If value is -1 or -2 put =1 we will consider any type of CNV at the moment without deffrentiation
      ValueCNV=as.numeric(ValueCNV)
      #     print("hello")
      #For debugging
      #     if(significant_genes[i]=='DNAH10')
      #    { print('this is DNAH10')
      #     print("this is valueCNV")
      #    print(ValueCNV)
      #   }
      # print(ValueCNV)
      #print(SamplesCNVbarcodes[j])
      if(ValueCNV==2|ValueCNV==-2)
      {
        #  print(significant_genes[i])
        ValueCNV=1
        BinaryMatrix[SamplesCNVbarcodes[j],significant_genes[i]]=ValueCNV
        detailedBinaryMatrix[SamplesCNVbarcodes[j],significant_genes[i]]='CNV'
      }
      #   print("yes")
      #print(SamplesCNVbarcodes[j])
      #print(significant_genes[i])
      #print(ValueCNV)
      #print("-------------")
     
      
      
    }
    
  indexgene=NULL
  
    }
    
  }









newList=AddMutationToBinaryMatrix(significant_genes,BinaryMatrix,detailedBinaryMatrix)
newList=AddFusionToBinaryMatrix(newList$BinaryMatrix,newList$detailedBinaryMatrix)

return(newList)  
}
#-------------------------------------------------------------------------------------
#Adding Mutations to the Binary Matrix and Detailed Binary Matrix

#BinaryMatrix=as.numeric(BinaryMatrix)
#Loop for mutation
#Make the row names the same as the samples for mutations
AddMutationToBinaryMatrix=function(significant_genes,BinaryMatrix,detailedBinaryMatrix)
{
row.names(BinaryMatrix)=SamplesMutationbarcode
row.names(detailedBinaryMatrix)=SamplesMutationbarcode
mutationState=getSampleSummary(mutationStatePerSample)
#For SurvivalClassifier Script
AllGenesMutation= getGeneSummary(mutationStatePerSample)$Hugo_Symbol

#significant_genes=unique(append(significant_genes_cnv,significant_genes_mutation))

for (i in 1: length(significant_genes))
{
  if (significant_genes[i] %in% AllGenesMutation)
  { 
 print(i)
MutatedSamplesPerGeneResults=genesToBarcodes(mutationStatePerSample,genes=significant_genes[i])
#transform the result of Mutated Samples Per genes into matrix so that I can access it easily
MutatedSamplesPerGeneResults=as.matrix(MutatedSamplesPerGeneResults[1])
MutatedSamplesPerGene=as.character(MutatedSamplesPerGeneResults[[1]]$Tumor_Sample_Barcode)
#Apply an Or operation with the already existing value inside the matrix
if(length(MutatedSamplesPerGene)!=0)
   {
for(j in 1:length(MutatedSamplesPerGene))
{
  print("Hello ")
  #Because we are interested of samples which have mutation information and Structural variants information
  if(MutatedSamplesPerGene[j] %in% rownames(BinaryMatrix) )
  {
    print(significant_genes[i])
    CNVvalue=as.integer(BinaryMatrix[MutatedSamplesPerGene[j],significant_genes[i]]) 
    Mutationvalue=1
   if (CNVvalue==1 && Mutationvalue==1)
   {detailedBinaryMatrix[MutatedSamplesPerGene[j],significant_genes[i]]='CNV-SNV'}
   else if(CNVvalue==0 && Mutationvalue==1)
   {detailedBinaryMatrix[MutatedSamplesPerGene[j],significant_genes[i]]='SNV'}
    
    temp=as.integer(BinaryMatrix[MutatedSamplesPerGene[j],significant_genes[i]])  | 1 
    if(temp==TRUE)
   {
      BinaryMatrix[MutatedSamplesPerGene[j],significant_genes[i]]=1
    }
      
  }
}
}
  
  }
}
newList<- list("BinaryMatrix"=BinaryMatrix,"detailedBinaryMatrix"=detailedBinaryMatrix)

return(newList)  
}
#-------------------------------------------------------------------------------------------------
#I want to return the matrix as it was (Row names) and debug by checking the genes that have only fusions are the same as the ones in the main file
AddFusionToBinaryMatrix=function(BinaryMatrix,detailedBinaryMatrix)
{
row.names(BinaryMatrix)=SamplesMutation
SamplesFusion=unique(SamplesFusion)
for (i in 1: length(SamplesFusion))
{
 #The location of the sample in the BinaryMatrix 
  x=which(row.names(BinaryMatrix)%in%SamplesFusion[i])
  #get the genes that have fusion for this sample from GeneA and append to this list the genes related to this samples from Gene B
  y=as.character(FusionGenes[FusionGenes$SamplesFusion==SamplesFusion[i],]$FusionGenes.Gene_A)
  y=unique(append(y,as.character(FusionGenes[FusionGenes$SamplesFusion==SamplesFusion[i],]$FusionGenes.Gene_B)))

  #Update the Values inside the BinaryMatrix for each gene that has a fusion for this sample
  if(length(x)!=0)
  {    
    for(j in 1:length(y)) 
  
     #Gene fusion that exists for this sample
 {
  
  print(as.integer(BinaryMatrix[x,y[j]])|1)
  temp=as.integer(BinaryMatrix[x,y[j]])|1
  if(temp==TRUE)
  {BinaryMatrix[x,y[j]]=1}
  if(detailedBinaryMatrix[x,y[j]]=="0")
  {
    detailedBinaryMatrix[x,y[j]]="Fusion"
  }
  else
  {
    detailedBinaryMatrix[x,y[j]]=paste(detailedBinaryMatrix[x,y[j]],"-Fusion")
  }
  
    }
  }
}
newList<- list("BinaryMatrix"=BinaryMatrix,"detailedBinaryMatrix"=detailedBinaryMatrix)

return(newList) 
}

#------------------------Choose a Threshold for picking up the interesting genes---------------------------------------
#Create a list of all the percentages for all the 1's in all the genes (the percentage of samples that are mutated for this gene out of all the samples)
BinaryMatrix=matrices$BinaryMatrix
detailedBinaryMatrix=matrices$detailedBinaryMatrix
Oneslist<-createOneslist(BinaryMatrix,SamplesCNVbarcodes,significant_genes)
createOneslist<-function(BinaryMatrixlocal,barcodes,significant_geneslocal)
{
max=0
count=0
PercentageList=matrix(0,nrow=1,ncol=length(significant_geneslocal))
Oneslist=matrix(0,nrow=1,ncol=length(significant_geneslocal))
colnames(PercentageList)=significant_geneslocal
colnames(Oneslist)=significant_geneslocal
for(i in 1: length(significant_geneslocal) )
{
  print(i)
  temp=table(BinaryMatrixlocal[,i])
  if(dim(temp)>1)
 {
    print(temp)
    percentage=(temp[["1"]][1]/length(barcodes))*100
   Oneslist[1,i]= temp[["1"]][1] 
    PercentageList[1,i]=percentage
 if (percentage>=max)
  {
    max=percentage
 }
  }
  else
  {
    PercentageList[1,i]=0
  }
  
}
return(Oneslist)
}
#from this plot we can see the frequency of the number of mutated samples we can decide upon a cutoff for the genes of interest
#from here
hist(Oneslist,main = "Frequency of Abberations", xlab="Frequency of Abberations", ylab="Number of genes")
#Filter Ones list (if you want)
quantile(Oneslist)
Oneslist[sapply(Oneslist, function(x) x>=as.numeric(quantile(Oneslist)[4]))]
##Oneslist[sapply(Oneslist, function(x) x>=20)]
#Check the number of genes in each mutation profile after filteration
##table(Oneslist[sapply(Oneslist, function(x) x>=20)])
#I choose my cutoff for the TCGA data to be 20 and up
#Select genes which have mutated number of samples 20 and up
#FilteredGeneList<-FilterGenes(significant_genes,20)
#Updating the cutoff point to be the data above the 4th quartile of the data
FilteredGeneList<-FilterGenes(significant_genes,as.numeric(quantile(Oneslist)[4]))
##FilteredGeneList<-FilterGenes(significant_genes,20)
FilterGenes<-function(significant_geneslocal,cutoff)
{
FilteredGeneList=NULL
IncludedIndices=NULL;
for(i in 1:length(significant_geneslocal))
{
  #print(i)
  if(Oneslist[1,i]>=cutoff)
   {
    print(i)
   IncludedIndices=append(IncludedIndices,i)
  }
}
for(i in 1:length(IncludedIndices))
{
  #print(i)
  
    FilteredGeneList=append(FilteredGeneList,Oneslist[1,IncludedIndices[i]])
 
}
return(FilteredGeneList)
}
#----------------------------Determine Mutant and Wildtype samples for each gene-------------------------------------
row.names(BinaryMatrix)=SamplesCNV
row.names(detailedBinaryMatrix)=SamplesCNV
MutantWildtype=getMutantWildtypeSamples("NPM1",BinaryMatrix)
getMutantWildtypeSamples<-function(genelocal,BinaryMatrixLocal)
{
  Mutant=BinaryMatrixLocal[BinaryMatrixLocal[, genelocal] == 1, ]
  Wildtype=BinaryMatrixLocal[BinaryMatrixLocal[, genelocal] == 0, ]
  Mutant=row.names(Mutant)
  Wildtype=row.names(Wildtype)
  newList<- list("Mutant"=Mutant,"Wildtype"=Wildtype)
  return(newList)
}

#----------------------------------------------------------------------------------------------------------------------

