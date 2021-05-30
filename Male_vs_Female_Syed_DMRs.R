##
library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm_SNPs_and_XY_removed.rda.rda')

# keep only Control samples (drop cases)
GRset=GRset[,which(GRset$Dx=='Control')] # 64 samples; 456928 probes

# Enable parallelization
library(doParallel)
registerDoParallel(cores = 4)

# Find bumps
library(bumphunter)

pd=pData(GRset)
sex=factor(as.vector(GRset$sex), levels=c('Female','Male')) # Female=0; Male=1

#Arguments for bumphunter
mod=model.matrix(~sex)
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/MaleVsFemale_Syed_svaobj.rda')
mod=cbind(mod,svaobj$sv)
p=getBeta(GRset)

anno=getAnnotation(GRset)
chr=as.vector(anno$chr)
pos=as.vector(anno$pos)

bumps = bumphunterEngine(p, mod, chr = chr, 
pos = pos, cutoff= 0.1, nullMethod = "bootstrap",
smooth=TRUE, B=1000)

dat=data.frame(chr=bumps$tab$chr,start=bumps$tab$start,end=bumps$tab$end,p.value=bumps$tab$p.value, FWER=bumps$tab$fwer)

save(bumps, dat, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/MaleVsFemale_Syed_DMRs.rda')

num_dmrs=length(which(dat$FWER<=0.1))
num_dmrs
