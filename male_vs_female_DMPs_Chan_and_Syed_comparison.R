load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/MaleVsFemale_Syed_all_DMPs.rda')
control_dmps=all_dmps

load("/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/MaleVsFemale_Chan_all_DMPs.rda")
PNET_dmps=all_dmps

mm=match(PNET_dmps$probe, control_dmps$probe)
mm
# [1] 13  1 NA  2

#### demonstatrate loss of gene imprinting (in IGF2)
library(minfi)

## IGF2 methylation in controls
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm_SNPs_and_XY_removed.rda.rda')

# keep only Control samples (drop PNET cases)
GRset=GRset[,which(GRset$Dx=='Control')] # 64 samples; 456928 probes
pd=pData(GRset)
beta=getBeta(GRset)

load("/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/MaleVsFemale_Syed_DMPs_proximal_to_genes.rda")
Syed_DMPs_proximal_to_genes=out2
Syed_DMPs_overlapping_IGF2=Syed_DMPs_proximal_to_genes[which(Syed_DMPs_proximal_to_genes$Ensembl=='ENSG00000167244'),]

mm=match(Syed_DMPs_overlapping_IGF2$probe,rownames(beta))
beta=beta[mm,]

pdf(file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/pdfs/MaleVsFemale_Syed_DMPs_overlapping_IGF2.pdf')

par(mar=c(5.1,5.3,4.1,2.1))

n=nrow(Syed_DMPs_overlapping_IGF2)

mycol=as.vector(pd$sex)
mycol=gsub('Male','blue',mycol)
mycol=gsub('Female','red',mycol)

x=pd$sex
x=factor(x, levels=c('Female','Male'))

for (i in 1:n){
	zag1=paste0('probe ', Syed_DMPs_overlapping_IGF2$probe[i])
	#zag2=paste0('p=',signif(Syed_DMPs_overlapping_IGF2$p.value[i],3), '; t=',signif(Syed_DMPs_overlapping_IGF2$t[i],3))
	zag2=paste0('FDR=',signif(Syed_DMPs_overlapping_IGF2$fdr[i],3))
	boxplot(as.vector(beta[i,])~x, ylab='DNAm (Beta)', xlab='Sex', main=paste0(zag1,' (',zag2,')'), cex.main=2, cex.lab=2, cex.axis=2, outline=FALSE, col='lightgrey', ylim=c(min(as.vector(beta[i,])),max(as.vector(beta[i,]))) )
	#stripchart(as.vector(beta[i,])~x, vertical = TRUE, method = "jitter", add = TRUE, pch = 21, col = 'black', bg=factor(mycol,levels=c('red','blue'))) #bg='gray'
	points(	as.vector(beta[i,]) ~ jitter(as.numeric(x), amount=0.2), pch =21, col='black', bg=mycol, cex=1.2)
}

dev.off()

## IGF2 methylation in PNET samples

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm_SNPs_and_XY_removed.rda.rda')

# keep only PNET samples (drop controls)
GRset=GRset[,which(GRset$Dx=='PNET')] # 23 samples; 456928 probes

pd=pData(GRset)
beta=getBeta(GRset)

mm=match(Syed_DMPs_overlapping_IGF2$probe,rownames(beta))
beta=beta[mm,]

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/MaleVsFemale_Chan_probeFit.rda')

mm=match(Syed_DMPs_overlapping_IGF2$probe, probeFit$probe)
probeFit=probeFit[mm,]

all(probeFit$probe==rownames(beta)) #TRUE

pdf(file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/pdfs/MaleVsFemale_Chan_DMPs_overlapping_IGF2.pdf')

par(mar=c(5.1,5.3,4.1,2.1))

n=nrow(Syed_DMPs_overlapping_IGF2)

mycol=as.vector(pd$sex)
mycol=gsub('Male','blue',mycol)
mycol=gsub('Female','red',mycol)

x=pd$sex
x=factor(x, levels=c('Female','Male'))

for (i in 1:n){
	zag1=paste0('probe ', probeFit$probe[i])
	#zag2=paste0('p=',signif(Syed_DMPs_overlapping_IGF2$p.value[i],3), '; t=',signif(Syed_DMPs_overlapping_IGF2$t[i],3))
	zag2=paste0('FDR=',signif(probeFit$fdr[i],3))
	boxplot(as.vector(beta[i,])~x, ylab='DNAm (Beta)', xlab='Sex', main=paste0(zag1,' (',zag2,')'), cex.main=2, cex.lab=2, cex.axis=2, outline=FALSE, col='lightgrey', ylim=c(min(as.vector(beta[i,])),max(as.vector(beta[i,]))) )
	#stripchart(as.vector(beta[i,])~x, vertical = TRUE, method = "jitter", add = TRUE, pch = 21, col = 'black', bg=factor(mycol,levels=c('red','blue'))) #bg='gray'
	points(	as.vector(beta[i,]) ~ jitter(as.numeric(x), amount=0.2), pch =21, col='black', bg=mycol, cex=1.2)
}

dev.off()



