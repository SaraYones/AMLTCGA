library(maftools)
#This script is to find the genes colliding within a somatic CNV that also has mutations in them.

#Find all the genes in the CNV file by merging each CNV file with the annotation file
requiredgene=NULL
RefAnnotation=read.table(file="/Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/gencode.v18.genes.bed",header = F,sep="\t")
RefAnnotation$V1=as.character(RefAnnotation$V1)
RefAnnotation$V2=as.integer(RefAnnotation$V2)
RefAnnotation$V3=as.integer(RefAnnotation$V3)
temp=as.matrix(strsplit( as.character(RefAnnotation$V10) , "\\;" ))
#from this statement you can find all the information about a certain gene
temp=t(do.call("cbind",as.list(temp)))
#Here you only want to get the gene name
temp=temp[,5]
temp=do.call(rbind,as.matrix(strsplit(temp , " " )))
#--------------------------------------------------------------------------------------------------------
returnSegmentFile= function(segments_file){
  allsegmentslocal=NULL
  
  seg.spl=read.table(file = segments_file,header = T)
  return(seg.spl)
}
#This function to make bedinrsect between the CNVfile and the annotation file
returnMergedMatrix=function(bedfile,annotationfile,coltosplit,pathtobed,outputfile)
{
  #x=system("/Users/sarayones/Downloads/bedtools2/bin/bedtools intersect -a /Users/sarayones/Desktop/CNVfile.bed -b /Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/gencode.v18.genes.bed -wa -wb > /Users/sarayones/Desktop/CNVfileintersect.bed",intern=TRUE)
  x=system(paste(pathtobed ,"intersect -a ",bedfile, " -b",annotationfile," -wa -wb >",outputfile),intern=TRUE)
                    
 #onekb=system("/Users/sarayones/Downloads/bedtools2/bin/bedtools intersect -a /Users/sarayones/Desktop/onekb.bed -b /Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/gencode.v18.genes.bed -wa -wb > /Users/sarayones/Desktop/onekbfileintersect.bed",intern=TRUE)
  x=as.matrix(read.table(outputfile, sep="\t", header=TRUE))
 # onekb=as.matrix(read.table("/Users/sarayones/Desktop/onekbfileintersect.bed", sep="\t", header=FALSE))   
  y=as.matrix(strsplit( as.character(x[,16]) , "\\;" ))
  z=do.call(rbind, y)
  j=do.call(rbind,as.matrix(strsplit( as.character(z[,5]) , " " )))
  merged=cbind(x[,c(1,2,3,4,6)],j[,3])
  merged=as.data.frame(merged)

  colnames(merged)<-c("chr","start","end","sample","log2","gene")
 
 # colnames(onekb)<-c("chr","start","end","sample","log2","gene")
  #colnames(merged)<-c("chr","start","end")
  #merged=merged[,c(1,2,3,6,5,4)]
  #bedr::bed2vcf(merged,file = "/Users/sarayones/Desktop/output.bed",zero.based = TRUE, header = NULL, fasta ="/Users/sarayones/Desktop/Projects/AMLProject/hg38.fa")
  #merged[,'chr']=as.character(merged[,'chr'])
  merged$chr=as.character(merged$chr)
  merged$start=strtoi(merged$start)
  merged$end=strtoi(merged$end)
  merged=bedr::bedr.sort.region(
    merged,
    method = "lexicographical",
    engine = "R",
    chr.to.num = c("X" = 23, "Y" = 24, "M" = 25),
    check.zero.based = TRUE,
    check.chr = FALSE,
    check.valid = TRUE,
    check.merge = FALSE,
    verbose = TRUE
  )
  
  print('I finished this function')
  return(merged)
  
}
# returns those genes that are common between CNVgenes and between the genes that are mutated for a certain  sample
getSigGenesPersample = function (CNVgenes,trial,mutationStatePerSample,regex,CNVsample)
{
#CNVgenes=unique(merged[,6])
#trial=genesToBarcodes(mutationStatePerSample,genes=as.character(CNVgenes))  
#trial=as.matrix(trial)
IntersectedGenes=rownames(trial)
significantgenesCNVmutation=NULL
for(j in 1:length(trial))
{
  Tumor_Sample_Barcode= trial[[j]]$Tumor_Sample_Barcode
 # barcode=gsub("TCGA.AB.([0-9]{4}).03[A-Z].01[A-Z].[0-9]{4}.[0-9]{2}","\\1",Tumor_Sample_Barcode)
  #print("Tumor_sample_barcode")
 #print(Tumor_Sample_Barcode)
  barcode=gsub(regex,"\\1",Tumor_Sample_Barcode)
 # print("This is the CNVsample")
#  print(CNVsample)
#  print("barcode after regex")
 # print(barcode)
  if( CNVsample %in% barcode)
  {
    print('I am here')
    significantgenesCNVmutation=append( significantgenesCNVmutation,IntersectedGenes[j])
  }
}
print('end of getSigGenesPersample function')
return(significantgenesCNVmutation)
}
#--------------------------------------------------------------------------------------------------------------------
#This function is used to get common genes that have both mutation and CNV at the same time coinciding together so I first find 
#The  genes that fall in the coordinates of the CNV segment file  (I add +1000 to the end and -1000 to the CNV coordinates)for each case then I check the genes that are mutated
#and if this gene is mutated and has a CNV  for a sample put it in a list to be added to the binary matrix later 
getcommongenes= function(mutationStateFile,filespath)
{
#Read all mutation data from maf file then start processing each gene to find if it exists in the sample CNV file 
#mutationStatePerSample1=read.maf('/Users/sarayones/Desktop/MutsigCVFiles/genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2/genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2.2.13.0.somatic.maf', useAll = TRUE)
mutationStatePerSample1=read.maf(mutationStateFile, useAll = TRUE)
#TCGA
#files=list.files("/Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3")
files=list.files(filespath)
#AMLCohort
#files=list.files("/Users/sarayones/Desktop/Projects/AMLProject/AMLCohortData/CNV/Diagnosis/Processed")
#temp=NULL
#Put them in a function later
#-----------------------------------------------------------------------------
#get the samples IDs for all CNVs
CNVbarcodes=files[which(grepl("TCGA_AB_[0-9]{4}_03A_01D_0756_21.nocnv_hg19.seg.txt", files))]
CNVbarcodes=gsub("(TCGA_AB_[0-9]{4}_03A_01D_0756_21).nocnv_hg19.seg.txt","\\1", CNVbarcodes)
CNVsamples=gsub("TCGA_AB_([0-9]{4})_03A_01D_0756_21","\\1", CNVbarcodes)
#---------------------------------------------------------------------------
#get the samples IDs for all Mutations 
Mutationbarcode=getSampleSummary(mutationStatePerSample1)
Mutationbarcode=list(Mutationbarcode[[1]])
Mutationbarcode=unlist(Mutationbarcode)
Mutationbarcode=levels(Mutationbarcode)
MutationIDSamples=gsub("TCGA.AB.([0-9]{4}).03[A-Z].01[A-Z].[0-9]{4}.[0-9]{2}","\\1",SamplesMutationbarcode)
#------------------------------------------------------------------------
#get a list of all files for CNV files
requiredFiles=files[which(grepl("TCGA_AB_[0-9]{4}_03A_01D_0756_21.nocnv_hg19.seg.txt", files))]
significantgenesCNVmutation=NULL
#onekb=NULL
#for(i in 1:length(files))
for(i in 1:length(requiredFiles))
{
  #Only use this If condition with TCGA data
 # if(grepl("TCGA_AB_[0-9]{4}_03A_01D_0756_21.nocnv_hg19.seg.txt", files[i]) == TRUE)
  #{
   
  
   #barcode=gsub("(TCGA_AB_[0-9]{4}_03A_01D_0756_21).nocnv_hg19.seg.txt","\\1",requiredFiles)
    
    #print('hello')
    #tcga_sample_number=gsub("TCGA_AB_([0-9]{4})_03A_01D_0756_21.nocnv_hg19.seg.txt","\\1",files[i])
    #normal=paste("TCGA_AB_",tcga_sample_number,"_11A_01D_0756_21.hg19.seg.txt",sep = "")
    #indexlocal=which(files %in% normal) 
    #matchedNormal=files[indexlocal]
    allsegments=NULL
    #TCGA
    path=paste("/Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/",requiredFiles[i],sep="")
    #AMLCohort
    #path=paste("/Users/sarayones/Desktop/Projects/AMLProject/AMLCohortData/CNV/Diagnosis/Processed/",files[i],sep="")
    #Read the CNV file at this path and transform it into bed (Add +1000 and -1000) and then write it to a file
    CNVfile<-returnSegmentFile(path)
    colnames(CNVfile)<-c("sample","chr","start","end","Num_Probes","Segment_Mean")
    CNVfile$chr=as.character(CNVfile$chr)
    CNVfile$start=as.integer(CNVfile$start)
    CNVfile$end=as.integer(CNVfile$end)
    CNVfile=CNVfile[,c(2,3,4,1,5,6)]
    CNVfile[,1]=paste("chr",CNVfile[,1], sep="")
    CNVfile[,"start"]=CNVfile[,"start"]-1000
    CNVfile[,"end"]=CNVfile[,"end"]+1000
    write.table(CNVfile,file = "/Users/sarayones/Desktop/CNVfile.bed",sep = "\t",quote = F,col.names = F ,row.names = F)
    #CNVfile=cbind(CNVfile[,1],CNVfile[,"start"]-1000,CNVfile[,"end"]+1000,CNVfile[,"Num_Probes"],CNVfile[,"Segment_Mean"])
   # colnames(CNVfile)<-c("chr","start","end","Num_Probes","Segment_Mean")
     #write.table(CNVfile,file = "/Users/sarayones/Desktop/onekb.bed",sep = "\t",quote = F,col.names = F ,row.names = F)
     #   n <- "foo.txt"
    #  if (file.exists(fn)) file.remove(fn)  
  #  View(CNVfile)
  #  print('what is sthe minimum')
   # print(min(CNVfile[,"start"]))
   print('before sorted region')
     bedr::is.sorted.region(
      CNVfile,
      method = "lex",
      engine = "R",
      check.zero.based = TRUE,
      check.chr = FALSE,
      check.valid = TRUE,
      check.merge = FALSE,
      verbose =TRUE
    )
    print('MAAAAAN')
  #   x=system("/Users/sarayones/Downloads/bedtools2/bin/bedtools intersect -a /Users/sarayones/Desktop/CNVfile.bed -b /Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/gencode.v18.genes.bed -wa -wb > /Users/sarayones/Desktop/CNVfileintersect.bed",intern=TRUE)
  #  onekb=system("/Users/sarayones/Downloads/bedtools2/bin/bedtools intersect -a /Users/sarayones/Desktop/onekb.bed -b /Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/gencode.v18.genes.bed -wa -wb > /Users/sarayones/Desktop/onekbfileintersect.bed",intern=TRUE)
   # x=as.matrix(read.table("/Users/sarayones/Desktop/CNVfileintersect.bed", sep="\t", header=TRUE))
    #  onekb=as.matrix(read.table("/Users/sarayones/Desktop/onekbfileintersect.bed", sep="\t", header=FALSE))   
    #  y=as.matrix(strsplit( as.character(x$V16) , "\\;" ))
    #z=do.call(rbind, y)
    #j=do.call(rbind,as.matrix(strsplit( as.character(z[,5]) , " " )))
    #merged=cbind(x[,c(1,2,3,4,6)],j[,3])
    #y=as.matrix(strsplit( as.character(onekb[,16]) , "\\;" ))
    #z=do.call(rbind, y)
    #j=do.call(rbind,as.matrix(strsplit( as.character(z[,5]) , " " )))
    #mergedonekb=cbind(x[,c(1,2,3,4,6)],j[,3])
    
    #merged=j
    #colnames(merged)<-c("chr","start","end","sample","log2","gene")
    #colnames(onekb)<-c("chr","start","end","sample","log2","gene")
    #colnames(merged)<-c("chr","start","end")
    #merged=merged[,c(1,2,3,6,5,4)]
    #bedr::bed2vcf(merged,file = "/Users/sarayones/Desktop/output.bed",zero.based = TRUE, header = NULL, fasta ="/Users/sarayones/Desktop/Projects/AMLProject/hg38.fa")
    #merged$chr=as.character(merged$chr)
    #merged$start=as.integer(merged$start)
    #merged$end=as.integer(merged$end)
    #merged=bedr::bedr.sort.region(
    # merged,
    # method = "lexicographical",
    # engine = "R",
    # chr.to.num = c("X" = 23, "Y" = 24, "M" = 25),
    # check.zero.based = TRUE,
    # check.chr = FALSE,
    #  check.valid = TRUE,
    # check.merge = FALSE,
    #  verbose = TRUE
    #)
    #return the Merged matrix (make intersect bed between the annotation file from gencode and the CNVfile) and return it into a matrix
    merged=returnMergedMatrix(pathtobed="/Users/sarayones/Downloads/bedtools2/bin/bedtools",bedfile="/Users/sarayones/Desktop/CNVfile.bed",annotationfile="/Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/gencode.v18.genes.bed",outputfile="/Users/sarayones/Desktop/CNVfileintersect.bed")
   # mergedonekb=returnMergedMatrix(pathtobed="Users/sarayones/Downloads/bedtools2/bin/bedtools",bedfile="/Users/sarayones/Desktop/onekbfileintersect.bed",annotationfile="/Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/gencode.v18.genes.bed",outputfile="/Users/sarayones/Desktop/onekbfileintersect.bed")
 #  print(CNVsamples[i])
    #Get the genes after this intersection 
    CNVgenes=unique(merged[,6])
    #get the genes that are mutated and their barcodes
    trial=genesToBarcodes(mutationStatePerSample1,genes=as.character(CNVgenes))  
    trial=as.matrix(trial)
    
    #Call this function to return those genes that are common between CNVgenes and between the genes that are mutated for each sample
   temp=getSigGenesPersample(CNVgenes,trial,mutationStatePerSample1,"TCGA.AB.([0-9]{4}).03[A-Z].01[A-Z].[0-9]{4}.[0-9]{2}",CNVsamples[i])
   #Append the returned genes to significantgenesCNVmutation list
    significantgenesCNVmutation=unique(append(significantgenesCNVmutation,temp))
    }
return(unique(significantgenesCNVmutation))
    
}


