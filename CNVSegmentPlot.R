#hello
#This script is mainly for plotting segments for each fragment and finding out which segments are significantly amplified or deleted.
files=list.files("/Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3")
#files=files[8:399]
plot_segments = function(segments_file){
  
  
  allsegmentslocal=NULL
  print(segments_file)
  seg = read.delim(segments_file, sep="\t")
  
  seg.spl = split(seg,as.factor(as.character(seg$Chromosome)))
  #print(seg.spl)
  print(length(seg.spl))
  
 ##pdf(file=paste(segments_file,"pdf",sep="."),paper="special",width=12,onefile=T,pointsize=8)
  
  #dev.new()
  
  par(mfrow=c(4,1))
  
  print('before')
  #print(seg.spl)
  length(seg.spl)
  for(i in 1:length(seg.spl)){
    
    
    
    x = seg.spl[[i]]    
    
    
   # print(x$Start)
    #print(x$Segment_Mean)
   # #plot(x$Start,x$Segment_Mean,xlim = c(x[1,3],x[nrow(x),4]),pch = "",ylim = c(-3,3),xlab = paste("chr",names(seg.spl[i]),sep="_"),ylab = "log2 ratio")
#  hello=x$Segment_Mean
  # # print("SegmentMean")
   # #print(x$Segment_Mean)
allsegmentslocal=append(allsegmentslocal,x$Segment_Mean)
  #  dev.new()
   ## hist(x$Segment_Mean)
    
    #print(x$Segment_Mean)
  ##  points(x$End,x$Segment_Mean,pch = "")
    
    
    
  ##  segments(x0 = x$Start , y0 = x$Segment_Mean , x1 = x$End, y1 = x$Segment_Mean, lwd = 2, col = "maroon")
    
    #Calculate 
    
  ##  abline(h = 0, lty = 1,lwd = 0.5)    
    
  }
#hist(x$Segment_Mean)
#dev.new()
 ## dev.off()
  print("printing all segments")
  print(allsegmentslocal)
  newList<- c("segments"=as.data.frame(allsegmentslocal),"segment_file"=seg.spl)
  return(newList)
}

#calculate the mean and variance for cnv in normal tissue (Germline+tumor(maybe)) and compare each seg_mean in the nocnv for each sample check if it is 
#gain, loss or normal according to the statistics you calclulated
compute_SD = function(){
meanSegmentslocal=NULL


  for(i in 1:length(files)){
  if(grepl("TCGA_AB_[0-9]{4}_03A_01D_0756_21.nocnv_hg19.seg.txt", files[i]) == TRUE)
  {
    tcga_sample_number=gsub("TCGA_AB_([0-9]{4})_03A_01D_0756_21.nocnv_hg19.seg.txt","\\1",files[i])
    normal=paste("TCGA_AB_",tcga_sample_number,"_11A_01D_0756_21.hg19.seg.txt",sep = "")
    indexlocal=which(files %in% normal) 
    matchedNormal=files[indexlocal]
    allsegments=NULL
    path=paste("/Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/",files[indexlocal],sep="")
    allsegmentslist<-plot_segments(path)
    NormalmeanSegmentslocal=mean(allsegmentslist$segments)
    
    NormalstandardDeviationlocal=sd(allsegments$segments)
    allsegments=NULL
    path=paste("/Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/",files[7],sep="")
   # print(path)
  allsegments<-plot_segments(path)
  #print("I am here")
  #print(allsegments)
    CancermeanSegmentslocal=append(meanSegmentslocal,c(mean(allsegments)))
    standardDeviationlocal=sd(allsegments)
    print(standardDeviationlocal)
  ## print(meanSegmentslocal)
  }
#  meanSegmentslocal=append(meanSegmentslocal,c(mean(allsegments)))
# print(meanSegmentslocal)
  }
print(meanSegmentslocal)
  return(meanSegmentslocal)
  }
meanSegments=NULL
meanSegments<-compute_SD()
print(meanSegments)