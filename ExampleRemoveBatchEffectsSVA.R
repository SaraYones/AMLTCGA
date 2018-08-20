#Trial to remove batch effects

library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)
pheno2 = pData(bladderEset)
edata = exprs(bladderEset)
batch2= pheno2$batch
modcombat2 = model.matrix(~1, pheno2)
combat_edata = ComBat(dat=t(edata), batch=batch2, mod=modcombat2)