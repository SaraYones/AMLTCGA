#Read each CNV file and compute segment mean then save it again in the /Processed Folder with the same name
#Move files with names following a certain pattern to another folder
#find . -name '*-D.clean.dedup.bam_CNVs' -exec mv {} Diagnosis/ \;
Inputpath="/Users/sarayones/Desktop/Projects/AMLProject/AMLCohortData/CNV/Diagnosis"
Outputpath="/Users/sarayones/Desktop/Projects/AMLProject/AMLCohortData/CNV/Diagnosis/Processed"
Regex="QC-[0-9]{4}-(AML[0-9]{3})-D.clean.dedup.bam_CNVs"

processAMLCohortCNV(Inputpath,Outputpath,Regex)
processAMLCohortCNV=function(Inputpath,Outputpath,Regex){

  
  files=list.files(Inputpath,recursive = TRUE)
  print(files)
  samples=gsub(Regex,"\\1",files)

for(i in 1:length(files))
{
 #These samples we will not work on them because they might have variants of healthy donors
 if(!(samples[i] %in% c("AML033","AML050","AML041")))
 {
  print(samples[i])

  AMLfile=read.table(paste(Inputpath,"/",files[i],sep=""),header = T,sep="\t")

  sampleid=rep(c(samples[i]),dim(AMLfile)[1])
  AMLfile=cbind(sampleid,AMLfile)
  colnames(AMLfile)<-c("Sample","Chromosome","Start","End","predicted_copy_number","type_of_alteration")
  AMLfile[AMLfile[,"predicted_copy_number"]== 0,"predicted_copy_number"]=AMLfile[AMLfile[,"predicted_copy_number"]== 0,"predicted_copy_number"]+0.00000000001
  Segment_Mean=log2(AMLfile[,"predicted_copy_number"]) -1
  AMLfile=cbind(AMLfile,Segment_Mean)
  AMLfile=AMLfile[,-c(5,6)]
  write.table(AMLfile,file = paste(Outputpath,"/",files[i],".processed",sep = ""),sep = "\t",quote = F)
 }
}


}