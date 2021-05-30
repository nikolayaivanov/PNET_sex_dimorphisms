##
library(minfi)
library(limma)
library(sva)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm_SNPs_and_XY_removed.rda.rda')

# keep only PNET samples (drop controls)
GRset=GRset[,which(GRset$Dx=='PNET')] # 23 samples; 456928 probes

library(AnnotationHub)
hub = AnnotationHub()
dm = query(hub, c("gencode", "sapiens", "v31"))
genes=dm[['AH75183']] # AH75183 | Annotated genes for Gencode v31 on hg19 coordinates

# find DMPs
all_DMPs_rda_filename="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/MaleVsFemale_Chan_all_DMPs.rda"
all_DMPs_csv_filename="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/data_tables/MaleVsFemale_Chan_all_DMPs.csv"
DMPs_proximal_to_genes_rda_filename="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/MaleVsFemale_Chan_DMPs_proximal_to_genes.rda"
DMPs_proximal_to_genes_csv_filename="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/data_tables/MaleVsFemale_Chan_DMPs_proximal_to_genes.csv"
pdf_filename='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/pdfs/MaleVsFemale_Chan_DMPs.pdf'

pd=pData(GRset)
beta=getBeta(GRset)

sex=factor(as.vector(GRset$sex), levels=c('Female','Male')) # Female=0; Male=1
table(sex)
# Female   Male
#      9     14
     
mod = model.matrix(~sex)

svaobj = sva(beta, mod)
# Number of significant surrogate variables is:  8

save(svaobj, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/MaleVsFemale_Chan_svaobj.rda')
#load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/MaleVsFemale_Chan_svaobj.rda')

mod=cbind(mod,svaobj$sv)

probe_fit=lmFit(object=beta,design=mod)
eb=eBayes(probe_fit)

probeFit=data.frame(probe=rownames(eb$p.value),
	intercept=probe_fit$coefficients[,1],
	slope=probe_fit$coefficients[,2], 
	p.value=eb$p.value[,2],
	fdr=p.adjust(as.vector(eb$p.value[,2]),method='fdr'),
	t=eb$t[,2])

rownames(probeFit)=NULL

o=order(probeFit$fdr)
probeFit=probeFit[o,]
beta=beta[o,]

save(probeFit, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/MaleVsFemale_Chan_probeFit.rda')

num_dmps=length(which(probeFit$fdr<=0.1)) # 4

percent_dmps=length(which(probeFit$fdr<=0.1))/nrow(probeFit) #8.754114e-06

dmps=probeFit[which(probeFit$fdr<=0.1),]

# Determine nearest genes to the DMPs

anno=getAnnotation(GRset)
mm=match(rownames(beta),rownames(anno))
anno=anno[mm,]

rti=as.data.frame(table(factor(anno$Relation_to_Island)))
num_shore_probes=sum(rti$Freq[grep('Shore', rti$Var1, ignore.case=TRUE)])
num_shelf_probes=sum(rti$Freq[grep('Shelf', rti$Var1, ignore.case=TRUE)])
num_island_probes=rti$Freq[grep('Island', rti$Var1, ignore.case=TRUE)]
num_opensea_probes=rti$Freq[grep('OpenSea', rti$Var1, ignore.case=TRUE)]
num_promoter_probes=length(grep('Promoter',anno$Regulatory_Feature_Group, , ignore.case=TRUE))
num_enhancer_probes=length(which(anno$Enhancer=="TRUE"))

anno_DMPs=anno[1:nrow(dmps),]

ad_df=as.data.frame(anno_DMPs)

gr=GRanges(seqnames=ad_df$chr, 
	ranges=IRanges(ad_df$pos, ad_df$pos),
	strand=ad_df$strand)

nn=nearest(gr, genes, select='arbitrary', ignore.strand=TRUE)
hits=genes[nn,]
dist=distance(gr,hits,ignore.strand=TRUE)

Enhancer=ad_df$Enhancer
Enhancer[which(Enhancer=="")]=FALSE

Promoter=rep('FALSE', times=nrow(dmps))
Promoter[grep('Promoter', ad_df$Regulatory_Feature_Group, ignore.case=TRUE)]='TRUE'

out=as.data.frame(cbind(dmps, anno_DMPs[,1:2], distance_to_nearest_gene=dist, Enhancer=Enhancer, Promoter=Promoter))
out$distance_to_nearest_gene=as.numeric(out$distance_to_nearest_gene)

all_dmps=as.data.frame(out)

save(all_dmps, file=all_DMPs_rda_filename) # all DMPs
write.csv(all_dmps,file=all_DMPs_csv_filename, row.names=FALSE) # all_DMPs

out1=out[which(out$distance_to_nearest_gene <= 5000 | out$Enhancer == 'TRUE' | out$Promoter == 'TRUE'),]

# Find genes that overlap or are <= 5kb from genes
oo=as.data.frame(findOverlaps(gr, genes, maxgap=5000, select='all', ignore.strand=TRUE))
hits=genes[oo$subjectHits,]
dist=distance(gr[oo$queryHits],hits,ignore.strand=TRUE)

pt1=dmps[oo$queryHits,]
pt2=anno_DMPs[,1:2]
pt2=pt2[oo$queryHits,]

out2=cbind(pt1, pt2, gene=as.vector(hits$Gene), Ensembl=as.vector(hits$Geneid), distance=dist)
mm=match(out2$probe,all_dmps$probe)
out2$Enhancer=all_dmps$Enhancer[mm]
out2$Promoter=all_dmps$Promoter[mm]

# label which genes are imprinted
load('/athena/masonlab/scratch/users/nai2008/Imprinted_genes/imprinted_genes.rda') #ig
ig=ig[-which(is.na(ig$Ensembl.ID))]
out2$Imprinted_Gene='FALSE'
out2$ExpressedAllele=NA
mm=match(out2$Ensembl,ig$Ensembl.ID)
out2$Imprinted_Gene[which(!is.na(mm))]=TRUE
out2$ExpressedAllele=as.vector(ig$ExpressedAllele)[mm]

save(out2, file=DMPs_proximal_to_genes_rda_filename)

write.csv(out2,file=DMPs_proximal_to_genes_csv_filename, row.names=FALSE) # only the DMPs that overlap or are within 5kb of genes

return=c(num_dmps=nrow(dmps),
	dmps_that_overlap_promoters_or_enhancers=length(which(out$Enhancer=='TRUE' | out$Promoter=='TRUE')),
	dmps_proximal_to_genes=length(unique(out2$probe)),
	dmps_that_overlap_promoters_or_enhancers_or_genes=nrow(out1),
	dmps_that_overlap_imprinted_genes=length(unique(out2$probe[which(out2$Imprinted_Gene==TRUE)])),
	number_of_imprinted_genes=length(unique(out2$Ensembl[which(out2$Imprinted_Gene==TRUE)])),
	range_of_abs_values_of_DNAm_differences_in_DMPs_MIN=signif(min(abs(out$slope)),2),
	range_of_abs_values_of_DNAm_differences_in_DMPs_MAX=signif(max(abs(out$slope)),2) )
#                                            num_dmps
#                                                4.00
#            dmps_that_overlap_promoters_or_enhancers
#                                                1.00
#                              dmps_proximal_to_genes
#                                                3.00
#   dmps_that_overlap_promoters_or_enhancers_or_genes
#                                                3.00
#                   dmps_that_overlap_imprinted_genes
#                                                   0
#                           number_of_imprinted_genes
#                                                   0
# range_of_abs_values_of_DNAm_differences_in_DMPs_MIN
#                                                0.15
# range_of_abs_values_of_DNAm_differences_in_DMPs_MAX
#                                                0.39

# output DMP pdf

pdf(file=pdf_filename)

center_hist='TRUE'

if (center_hist=='TRUE'){
	hist(dmps$slope, xlim = c(-max(abs(dmps$slope)),max(abs(dmps$slope))), xlab='Male vs. Female\n change in methylation', ylab='Frequency', main='DMPs', col='grey')
		} else {
	hist(dmps$slope, xlab='Male vs. Female\n change in methylation', ylab='Frequency', main='DMPs', col='grey')
}

rti_dmps=as.data.frame(table(factor(anno_DMPs$Relation_to_Island)))
Island = rti_dmps$Freq[which(rti_dmps$Var1=='Island')]/num_island_probes
Shore= sum(rti_dmps$Freq[grep('Shore', rti_dmps$Var1, ignore.case=TRUE)])/num_shore_probes
Shelf= sum(rti_dmps$Freq[grep('Shelf', rti_dmps$Var1, ignore.case=TRUE)])/num_shelf_probes
OpenSea= rti_dmps$Freq[which(rti_dmps$Var1=='OpenSea')]/num_opensea_probes
Promoter=length(grep('Promoter',anno_DMPs$Regulatory_Feature_Group, , ignore.case=TRUE))/num_promoter_probes
Enhancer=length(which(anno_DMPs$Enhancer=="TRUE"))/num_enhancer_probes

normalized_Relation_to_Island=data.frame(Island=Island, Shore=Shore, Shelf=Shelf, OpenSea=OpenSea, Promoter=Promoter, Enhancer=Enhancer)

bp=barplot(height=as.matrix(normalized_Relation_to_Island), xlab='Genomic Position', ylab='Normalized Frequency', main='DMPs', col='lightgrey')
labs=signif(as.vector(normalized_Relation_to_Island),3)
text(bp, 0, labs, cex=1, pos=3)

par(mar=c(5.1,5.3,4.1,2.1))

if (num_dmps<100) { n = num_dmps} else {n =100 }

mycol=as.vector(pd$sex)
mycol=gsub('Male','blue',mycol)
mycol=gsub('Female','red',mycol)

x=pd$sex
x=factor(x, levels=c('Female','Male'))

for (i in 1:n){
	zag1=paste0('probe ', probeFit$probe[i])
	#zag2=paste0('p=',signif(probeFit$p.value[i],3), '; t=',signif(probeFit$t[i],3))
	zag2=paste0('FDR=',signif(probeFit$fdr[i],3))
	boxplot(as.vector(beta[i,])~x,ylab='DNAm (Beta)', xlab='Sex', main=paste0(zag1,' (',zag2,')'), cex.main=2, cex.lab=2, cex.axis=2, outline=FALSE, col='lightgrey', ylim=c(min(as.vector(beta[i,])),max(as.vector(beta[i,]))) )
	#stripchart(as.vector(beta[i,])~x, vertical = TRUE, method = "jitter", add = TRUE, pch = 21, col = 'black', bg=factor(mycol,levels=c('red','blue'))) #bg='gray'
	points(	as.vector(beta[i,]) ~ jitter(as.numeric(x), amount=0.2), pch =21, col='black', bg=mycol, cex=1.2)
}

dev.off()






















