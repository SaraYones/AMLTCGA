#update MAF with genes that have copy number alterations
#Read Annnotation file and prepare it
library(bedr)

requiredgene=NULL
RefAnnotation=read.table(file="/Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/gencode.v18.genes.bed",header = F,sep="\t")
RefAnnotation$V1=as.character(RefAnnotation$V1)
RefAnnotation$V2=as.integer(RefAnnotation$V2)
RefAnnotation$V3=as.integer(RefAnnotation$V3)
temp=as.matrix(strsplit( as.character(RefAnnotation$V10) , "\\;" ))
temp=t(do.call("cbind",as.list(temp)))
temp=temp[,5]
temp=do.call(rbind,as.matrix(strsplit(temp , " " )))
#which(temp[,3]=='ETS2') # the gene in question
chr=RefAnnotation[which(temp[,3]=='ETS2'),][[1]]
start=RefAnnotation[which(temp[,3]=='ETS2'),][[2]]
end=RefAnnotation[which(temp[,3]=='ETS2'),][[3]]
requiredgene=rbind(requiredgene,c(chr,start,end))
colnames(requiredgene)<-c("chr","start","end")
write.table(requiredgene,file = "/Users/sarayones/Desktop/CNVAnnotated.bed",sep = "\t",quote = F,col.names = F, row.names = F)
x=system("/Users/sarayones/Downloads/bedtools2/bin/bedtools getfasta -fi /Users/sarayones/Desktop/Projects/AMLProject/hg38.fa -bed /Users/sarayones/Desktop/CNVAnnotated.bed -fo  /Users/sarayones/Desktop/test.fa.out  ",intern=TRUE)
reference_allele=read.table(file="/Users/sarayones/Desktop/test.fa.out",header = F,sep="\t")
trialMAF=matrix(NA,nrow=0,ncol=70)

returnSegmentFile= function(segments_file){
  allsegmentslocal=NULL
  
  seg.spl=read.table(file = segments_file,header = T)
  return(seg.spl)
}
#TCGA
files=list.files("/Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3")
#AMLCohort
#files=list.files("/Users/sarayones/Desktop/Projects/AMLProject/AMLCohortData/CNV/Diagnosis/Processed")
#temp=NULL
trialMAF=NULL
for(i in 1:length(files))
  {
  #Only use this If condition with TCGA data
  if(grepl("TCGA_AB_[0-9]{4}_03A_01D_0756_21.nocnv_hg19.seg.txt", files[i]) == TRUE)
  {
    barcode=gsub("(TCGA_AB_[0-9]{4}_03A_01D_0756_21).nocnv_hg19.seg.txt","\\1",files[i])
  

    #tcga_sample_number=gsub("TCGA_AB_([0-9]{4})_03A_01D_0756_21.nocnv_hg19.seg.txt","\\1",files[i])
    #normal=paste("TCGA_AB_",tcga_sample_number,"_11A_01D_0756_21.hg19.seg.txt",sep = "")
    #indexlocal=which(files %in% normal) 
    #matchedNormal=files[indexlocal]
    allsegments=NULL
    #TCGA
    path=paste("/Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/",files[i],sep="")
    #AMLCohort
    #path=paste("/Users/sarayones/Desktop/Projects/AMLProject/AMLCohortData/CNV/Diagnosis/Processed/",files[i],sep="")
    CNVfile<-returnSegmentFile(path)
    colnames(CNVfile)<-c("sample","chr","start","end","Num_Probes","Segment_Mean")
    CNVfile$chr=as.character(CNVfile$chr)
    CNVfile$start=as.integer(CNVfile$start)
    CNVfile$end=as.integer(CNVfile$end)
    CNVfile=CNVfile[,c(2,3,4,1,5,6)]
    CNVfile[,1]=paste("chr",CNVfile[,1], sep="")
    write.table(CNVfile,file = "/Users/sarayones/Desktop/CNVfile.bed",sep = "\t",quote = F,col.names = F ,row.names = F)
 #   n <- "foo.txt"
  #  if (file.exists(fn)) file.remove(fn)  
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
 x=system("/Users/sarayones/Downloads/bedtools2/bin/bedtools intersect -a /Users/sarayones/Desktop/CNVfile.bed -b /Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/gencode.v18.genes.bed -wa -wb > /Users/sarayones/Desktop/CNVfileintersect.bed",intern=TRUE)
 x=as.matrix(read.table("/Users/sarayones/Desktop/CNVfileintersect.bed", sep="\t", header=TRUE))
 x=read.table(file="/Users/sarayones/Desktop/CNVfileintersect.bed",header = F,sep="\t")
 y=as.matrix(strsplit( as.character(x$V16) , "\\;" ))
 z=do.call(rbind, y)
 j=do.call(rbind,as.matrix(strsplit( as.character(z[,5]) , " " )))
 
 merged=cbind(x[,c(1,2,3,4,6)],j[,3])
 #merged=j
 colnames(merged)<-c("chr","start","end","sample","log2","gene")
 #colnames(merged)<-c("chr","start","end")
 #merged=merged[,c(1,2,3,6,5,4)]
 #bedr::bed2vcf(merged,file = "/Users/sarayones/Desktop/output.bed",zero.based = TRUE, header = NULL, fasta ="/Users/sarayones/Desktop/Projects/AMLProject/hg38.fa")
 merged$chr=as.character(merged$chr)
 merged$start=as.integer(merged$start)
 merged$end=as.integer(merged$end)
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

#colnames(trialMAF)<-c("Hugo_Symbol"	Entrez_Gene_Id	Center	NCBI_Build	Chromosome	Start_Position	End_Position	Strand	Variant_Classification	Variant_Type	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2	dbSNP_RS	dbSNP_Val_Status	Tumor_Sample_Barcode	Matched_Norm_Sample_Barcode	Match_Norm_Seq_Allele1	Match_Norm_Seq_Allele2	Tumor_Validation_Allele1	Tumor_Validation_Allele2	Match_Norm_Validation_Allele1	Match_Norm_Validation_Allele2	Verification_Status	Validation_Status	Mutation_Status	Sequencing_Phase	Sequence_Source	Validation_Method	Score	BAM_File	Sequencer	Tumor_Sample_UUID	Matched_Norm_Sample_UUID	chromosome_name	start	stop	reference	variant	type	gene_name	transcript_name	transcript_species	transcript_source	transcript_version	strand	transcript_status	trv_type	c_position	amino_acid_change	ucsc_cons	domain	all_domains	deletion_substructures	transcript_error	NormalRefReads_WU	NormalVarReads_WU	NormalVAF_WU	TumorRefReads_WU	TumorVarReads_WU	TumorVAF_WU	RNARefReads_WU	RNAVarReads_WU	RNAVAF_WU)
requiredgene=merged[which(merged$gene%in% 'ETS2'),] #required genes
if(dim(requiredgene)[1]!=0)
{
  print('hello')
 # write.table(requiredgene,file = "/Users/sarayones/Desktop/CNVAnnotated.bed",sep = "\t",quote = F,col.names = F, row.names = F)
 #x=system("/Users/sarayones/Downloads/bedtools2/bin/bedtools getfasta -fi /Users/sarayones/Desktop/Projects/AMLProject/hg38.fa -bed /Users/sarayones/Desktop/CNVAnnotated.bed -fo  /Users/sarayones/Desktop/test.fa.out  ",intern=TRUE)
 #x=read.table(file="/Users/sarayones/Desktop/test.fa.out",header = F,sep="\t")

 state=2^(requiredgene$log2)
 refrenceAllele=as.character(reference_allele[2,1])
 if(state[1]<2)
 {
   #Deletion
   Tumor_Seq_Allele2=as.character(reference_allele[2,1])
   trialMAF=rbind(trialMAF,c("ETS2","0","genome.wustl.edu","36",requiredgene$chr,requiredgene$start[1],requiredgene$start[1],"Frame_Shift_Del","DEL",Tumor_Seq_Allele2,Tumor_Seq_Allele2,"-","","",barcode,barcode,Tumor_Seq_Allele2,Tumor_Seq_Allele2,Tumor_Seq_Allele2,"-",Tumor_Seq_Allele2,Tumor_Seq_Allele2,"Verified"
                    ,"Valid","Somatic","Phase_IV","WXS","Hybrid_Capture_Illumina_Seq","1","dbGAP","Illumina","HiSeq","166d3943-4006-4434-9af3-4ba8e86b6ac6","a109df8d-1667-48a5-a2bc-91d25851124a",requiredgene$chr,requiredgene$start[1],requiredgene$start[1],Tumor_Seq_Allele2,"-","DEL",
                    "ETS2","ENSG00000157557","human","genbank","54_36p","-1","reviewed","frame_shift_del","c.1812_1811","p.605in_frame_insYDLKWE","1.000:1.000","superfamily_Kinase_like","HMMPfam_Pkinase_Tyr","HMMSmart_TyrKc","PatternScan_RECEPTOR_TYR_KIN_III","PatternScan_PROTEIN_KINASE_TYR","superfamily_Kinase_like","HMMPfam_ig","PatternScan_PROTEIN_KINASE_ATP,superfamily_SSF48726","-","no_errors",
                    "NA","NA","NA","NA","NA","NA","NA","NA","NA"))
   
 }
   else
   {
     #Insertion
     Tumor_Seq_Allele2=as.character(reference_allele[2,1])
     for(j in as.integer(state))
     {Tumor_Seq_Allele2=paste(Tumor_Seq_Allele2,as.character(x[2,1]))}
     trialMAF=rbind(trialMAF,c("ETS2","0","genome.wustl.edu","36",requiredgene$chr,requiredgene$start[1],as.character(as.integer(requiredgene$start[1]+1)),"In_Frame_Ins","INS","-","-",Tumor_Seq_Allele2,"","",barcode,barcode,"-","-","-",Tumor_Seq_Allele2,"-","-","Verified"
                      ,"Valid","Somatic","Phase_IV","WXS","Hybrid_Capture_Illumina_Seq","1","dbGAP","Illumina","HiSeq","166d3943-4006-4434-9af3-4ba8e86b6ac6","a109df8d-1667-48a5-a2bc-91d25851124a",requiredgene$chr,requiredgene$start[1],as.character(as.integer(requiredgene$start[1]+1)),"-",Tumor_Seq_Allele2,"INS",
                      "ETS2","ENSG00000157557","human","genbank","54_36p","-1","reviewed","in_frame_ins","c.1812_1811","p.605in_frame_insYDLKWE","1.000:1.000","superfamily_Kinase_like","HMMPfam_Pkinase_Tyr","HMMSmart_TyrKc","PatternScan_RECEPTOR_TYR_KIN_III","PatternScan_PROTEIN_KINASE_TYR","superfamily_Kinase_like","HMMPfam_ig","PatternScan_PROTEIN_KINASE_ATP,superfamily_SSF48726","-","no_errors",
                      "NA","NA","NA","NA","NA","NA","NA","NA","NA"))
   }
}
#see if it is deleted or amplified
if (file.exists("/Users/sarayones/Desktop/CNVAnnotated.bed")) file.remove("/Users/sarayones/Desktop/CNVAnnotated.bed")  
if (file.exists("/Users/sarayones/Desktop/CNVfile.bed")) file.remove("/Users/sarayones/Desktop/CNVfile.bed")  
if (file.exists("/Users/sarayones/Desktop/test.fa.out")) file.remove("/Users/sarayones/Desktop/test.fa.out")  
 
   
  }
  
}

write.table(trialMAF,file = "/Users/sarayones/Desktop/trialMAF.maf",sep = "\t",quote = F,col.names = F ,row.names = F)
trialMAF=trialMAF[trialMAF[,9]=="DEL",]
