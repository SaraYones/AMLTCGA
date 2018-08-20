#This script is for merging all CNV files for TCGA into one where first coloumn represents the sample name 


returnSegmentFile= function(segments_file){
allsegmentslocal=NULL

seg.spl=read.table(file = segments_file,header = T)
return(seg.spl)
}
#TCGA
#files=list.files("/Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3")
#AMLCohort
files=list.files("/Users/sarayones/Desktop/Projects/AMLProject/AMLCohortData/CNV/Diagnosis/Processed")
temp=NULL
for(i in 1:length(files)){
  #Only use this If condition with TCGA data
  #if(grepl("TCGA_AB_[0-9]{4}_03A_01D_0756_21.nocnv_hg19.seg.txt", files[i]) == TRUE)
  {
    
    #tcga_sample_number=gsub("TCGA_AB_([0-9]{4})_03A_01D_0756_21.nocnv_hg19.seg.txt","\\1",files[i])
    #normal=paste("TCGA_AB_",tcga_sample_number,"_11A_01D_0756_21.hg19.seg.txt",sep = "")
    #indexlocal=which(files %in% normal) 
    #matchedNormal=files[indexlocal]
    allsegments=NULL
    #TCGA
    #path=paste("/Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/",files[i],sep="")
    #AMLCohort
    path=paste("/Users/sarayones/Desktop/Projects/AMLProject/AMLCohortData/CNV/Diagnosis/Processed/",files[i],sep="")
    segmentFile<-returnSegmentFile(path)
     temp=rbind(temp,segmentFile)
  }
  
  }

length(files)
    
write.table(temp, file = "/Users/sarayones/Desktop/Projects/AMLProject/AMLCohortData/CNV/Diagnosis/Processed/MergedCNVAMLCohort.txt", sep = "\t", col.names = FALSE, row.names = FALSE,quote = FALSE)