    #library(devtools)
  #install_github('davetang/bedr')
library(bedr)
  
  CNVfile=read.table(file="/home/saryou/AMLProject/AMLTCGA/TCGA_AB_2877_03A_01D_0756_21.hg19.seg.txt",header = T,sep="\t")
colnames(CNVfile)<-c("sample","chr","start","end","Num_Probes","Segment_Mean")
CNVfile$chr=as.character(CNVfile$chr)
CNVfile$start=as.integer(CNVfile$start)
CNVfile$end=as.integer(CNVfile$end)
CNVfile=CNVfile[,c(2,3,4,1,5,6)]
#CNVfile=CNVfile[,c(1,2,3)]
#drops <- c("Num_Probes")
#CNVfile=CNVfile[ , !(names(CNVfile) %in% drops)]
#drops <- c("sample")
#CNVfile=CNVfile[ , !(names(CNVfile) %in% drops)]
#drops <- c("Segment_Mean")
#CNVfile=CNVfile[ , !(names(CNVfile) %in% drops)]
RefAnnotation=read.table(file="/home/saryou/AMLProject/AMLTCGA/gencode.v18.genes.bed",header = F,sep="\t")
RefAnnotation$V1=as.character(RefAnnotation$V1)
RefAnnotation$V2=as.integer(RefAnnotation$V2)
RefAnnotation$V3=as.integer(RefAnnotation$V3)
hello2=bedr::convert2bed( RefAnnotation,
                         set.type = TRUE,
                         check.zero.based = TRUE,
                         check.chr = FALSE,
                         check.valid = TRUE,
                         check.sort = FALSE,
                         check.merge = FALSE,
                         verbose = TRUE
)
#if you want to merge two bed files both have to follow the same naming convention
CNVfile[,1]=paste("chr",CNVfile[,1], sep="")
RefAnnotation=bedr::bedr.sort.region(
  RefAnnotation,
  method = "lexicographical",
  engine = "R",
  chr.to.num = c("X" = 23, "Y" = 24, "M" = 25),
  check.zero.based = TRUE,
  check.chr = FALSE,
  check.valid = TRUE,
  check.merge = FALSE,
  verbose = TRUE
)
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
hello=bedr::convert2bed( CNVfile,
  set.type = TRUE,
  check.zero.based = TRUE,
  check.chr = FALSE,
  check.valid = TRUE,
  check.sort = FALSE,
  check.merge = FALSE,
  verbose = TRUE
)
CNVfile=bedr::bedr.sort.region(
 CNVfile,
  method = "lexicographical",
  engine = "R",
  chr.to.num = c("X" = 23, "Y" = 24, "M" = 25),
  check.zero.based = TRUE,
  check.chr = FALSE,
  check.valid = TRUE,
  check.merge = FALSE,
  verbose = TRUE
)

write.table(CNVfile,file = "/home/saryou/AMLProject/AMLTCGA/CNVfile.bed",sep = "\t",quote = F,col.names = F ,row.names = F)

a.int1 <- bedr::bedr(engine="bedtools2",input = list(a = RefAnnotation, b = CNVfile), method = "intersect", params = "-loj")

#system("/Users/sarayones/Downloads/bedtools2/bin/bedtools intersect -a /Users/sarayones/Desktop/CNVfile.bed -b /Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/gencode.v18.genes.bed /Users/sarayones/Desktop/exons.bed")
#system("/Users/sarayones/Downloads/bedtools2/bin/bedtools intersect -a /Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/gencode.v18.genes.bed -b /Users/sarayones/Downloads/TCGA_LAML/CNV/broad.mit.edu_LAML.Genome_Wide_SNP_6.Level_3/gencode.v18.genes.bed -wa -wb")
#system("/Users/sarayones/Downloads/bedtools2/bin/bedtools intersect -a /Users/sarayones/Desktop/CNVfile.bed -b /Users/sarayones/Desktop/CNVfile.bed -wa -wb")
x=system("bedtools intersect -a /home/saryou/AMLProject/AMLTCGA/CNVfile.bed -b /home/saryou/AMLProject/AMLTCGA/gencode.v18.genes.bed -wa -wb > /home/saryou/AMLProject/AMLTCGA/CNVfileintersect.bed",intern=TRUE)

x=as.matrix(read.table("/home/saryou/AMLProject/AMLTCGA/CNVfileintersect.bed", sep="\t", header=TRUE))
x=read.table(file="/home/saryou/AMLProject/AMLTCGA/CNVfileintersect.bed",header = F,sep="\t")
#extract gene name from the coloumn by transforming the required coloumn into matrix then  extracting the gene name coloumn then transforming it into matrix
#again to extract only the required gene name
y=as.matrix(strsplit( as.character(x$V16) , "\\;" ))
z=do.call(rbind, y)
j=do.call(rbind,as.matrix(strsplit( as.character(z[,5]) , " " )))
merged=cbind(x[,c(1,2,3,4,6)],j[,3])
#merged=j
colnames(merged)<-c("chr","start","end","sample","log2","gene")
#colnames(merged)<-c("chr","start","end")
#merged=merged[,c(1,2,3,6,5,4)]
write.table(merged,file = "/home/saryou/AMLProject/AMLTCGA/CNVAnnotated.bed",sep = "\t",quote = F,col.names = T, row.names = F)
merged$chr=as.character(merged$chr)
merged$start=as.integer(merged$start)
merged$end=as.integer(merged$end)
#bedr::bed2vcf(merged,file = "/home/saryou/AMLProject/AMLTCGA/output.bed",zero.based = TRUE, header = NULL, fasta ="/home/saryou/indexedChromosomes/hg38.fa")

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

#bedr::bed2vcf(CNVfile,file = "/Users/sarayones/Desktop/output.bed",zero.based = FALSE, header = NULL, fasta = "/Users/sarayones/Desktop/Projects/AMLProject/hg38.fa")
#is.valid.region(
 # CNVfile,
#  check.zero.based = TRUE,
#  check.chr = TRUE,
#  throw.error = FALSE,
#  verbose = TRUE
#)