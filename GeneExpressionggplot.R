# library
library(ggplot2)
require(reshape2)

#Drawing a regular box plot for Gene expression

sample <- rownames(logGeneExpressionMatrix)
d.f <- data.frame(sample, logGeneExpressionMatrix)
d.f2 <- melt(d.f, id.vars = "sample")
geneNamesplot=colnames(logGeneExpressionMatrix)
#ggplot(d.f2, aes(sample, value, group = variable, colour = variable)) + geom_line()


ggplot(d.f2, aes(x=as.factor(variable), y=value)) + geom_boxplot(fill="slateblue", alpha=0.2)+xlab("genes")




# library
library(ggplot2)

# create a data frame
variety=rep(LETTERS[1:7], each=40)
treatment=rep(c("high","low"),each=20)
note=seq(1:280)+sample(1:150, 280, replace=T)
data=data.frame(variety, treatment ,  note)




# grouped boxplot
ggplot(data, aes(x=variety, y=note, fill=treatment)) + 
  geom_boxplot()

#Drawing a grouped box plot for Gene expression for each gene beside eachother     

MutantWildtypelocal=getMutantWildtypeSamples(geneNamesLocal[1],GEBinaryMatrix)
# MutantBatches=extractBatchesClasses(MutantWildtype$Mutant)
# WildtypeBatches=extractBatchesClasses(MutantWildtype$Wildtype)
xlocal=logGeneExpressionMatrixWithoutBatch[MutantWildtypelocal$Mutant,geneNamesLocal[1]]
ylocal=logGeneExpressionMatrixWithoutBatch[MutantWildtypelocal$Wildtype,geneNamesLocal[1]]
      
mutationstatuslocal1=rep(c("mutant"),length(xlocal))
mutationstatuslocal2=rep(c("wildtype"),length(ylocal))
mutationstatuslocal=c(mutationstatuslocal1,mutationstatuslocal2)

note=c(xlocal,ylocal)
names(note)<-NULL
variety=rep(c("FLT3"),length(SampleGE))
data=data.frame(variety,mutationstatuslocal,note)

ggplot(data, aes(x=variety, y=note, fill=mutationstatuslocal)) + 
  geom_boxplot()
