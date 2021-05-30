##
library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm_SNPs_and_XY_removed.rda.rda')

# keep only Control samples (drop cases)
GRset=GRset[,which(GRset$Dx=='Control')] # 64 samples; 456928 probes

# Enable parallelization
library(doParallel)
registerDoParallel(cores = 2)

# Find blocks
pd=pData(GRset)
beta=getBeta(GRset)

sex=factor(as.vector(GRset$sex), levels=c('Female','Male')) # Female=0; Male=1

mod=model.matrix(~sex)
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/MaleVsFemale_Syed_svaobj.rda')
mod=cbind(mod,svaobj$sv)

cobj=cpgCollapse(GRset, what="Beta")

blocks=blockFinder(cobj$object, design=mod, coef = 2, what = 'Beta', cluster=NULL, nullMethod='bootstrap',
cutoff = 0.1, pickCutoff = FALSE, smooth = TRUE, smoothFunction = locfitByCluster,
B = 1000, verbose = TRUE, bpSpan = 2.5 * 10^5)

dat=data.frame(chr=blocks$tab$chr,start=blocks$tab$start,end=blocks$tab$end,p.value=blocks$tab$p.value, FWER=blocks$tab$fwer)

save(blocks, dat, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/MaleVsFemale_Syed_blocks.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/MaleVsFemale_Syed_blocks.rda')

sig_blocks=blocks$tab[which(blocks$tab$fwer<=0.1),]
nrow(sig_blocks) # 0