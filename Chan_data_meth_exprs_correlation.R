## read in DNAm data
library(minfi)

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm_SNPs_and_XY_removed.rda.rda')
colnames(GRset)=GRset$my_unique_ID

# keep only PNET samples (drop controls)
GRset=GRset[,which(GRset$Dx=='PNET')] # 23 samples; 456928 probes

pd=pData(GRset)
beta=getBeta(GRset)

## load sex-DMP data
dmp_table=read.csv("/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/data_tables/MaleVsFemale_Chan_DMPs_proximal_to_genes.csv")

dmp_table$meth_expres_corr_coeff=NA
dmp_table$meth_expres_p_value=NA

## read in RNAseq data
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_Chan_dataset_DESeq2_DEbySex_BUNDLE.rda') # se, gse, tpm, dds, ddsClean, rr, vsd, rld

log2TPMplus1=log2(tpm+1)

geneExps_metadata=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Chan_et_al_GSE118014_PNET_RNAseq_data/full_pheno_data_with_all_IDs.csv')

id_conversion_table=geneExps_metadata[,c(1,28)]

mm=match(colnames(log2TPMplus1),id_conversion_table$SRA_ID)
colnames(log2TPMplus1)=id_conversion_table$Paper_ID[mm]

mm=match(colnames(beta),colnames(log2TPMplus1))
log2TPMplus1=log2TPMplus1[,mm]

# correlate DNAm w/ gene expression

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/pdfs/Chan_data_meth_exprs_correlation.pdf')

for(i in 1: nrow(dmp_table)){

	m_exprs=match(dmp_table$Ensembl[i], rownames(log2TPMplus1))

	if(!is.na(m_exprs)){
		m_beta=match(dmp_table$probe[i],rownames(beta))
		cc=cor.test(as.numeric(beta[m_beta,]), as.numeric(log2TPMplus1[m_exprs,]), alternative='two.sided', method='spearman', exact=FALSE)
		dmp_table$meth_expres_corr_coeff[i]=cc$estimate
		dmp_table$meth_expres_p_value[i]=cc$p.value
		
		if(!is.na(as.vector(cc$estimate))){

			legend_label=paste0('rho = ',signif(cc$estimate,2),'; p = ',signif(cc$p.value,2))
			
			fit=lm(as.numeric(log2TPMplus1[m_exprs,])~as.numeric(beta[m_beta,]))
			eb=summary(fit)
			t=eb$coef[2,3]
			p=eb$coef[2,4]
			intercept=eb$coef[1,1]
			slope=eb$coef[2,1]

			title=paste0('CpG: ', dmp_table$probe[i],'; Proximal gene: ',dmp_table$gene[i],'; Distance: ', dmp_table$distance[i],' bps')

			par(mar=c(5.1,5.2,4.1,2.1))

			plot(as.numeric(beta[m_beta,]), as.numeric(log2TPMplus1[m_exprs,]),xlab='DNAm(Beta)',ylab=as.expression(bquote(log[2](TPM+1))),main=title, pch=21, col='black', bg='grey', cex=2, cex.axis=2, cex.lab=2, cex.main=1)
			legend('topright', legend_label, cex=.8, bty = "n")

			#ylab='log2(TPM+1)'

			#as.expression(bquote(log[2](TPM+1)))

			abline(intercept,slope,lty='solid',col='black',lwd=3)

		}
		
	}

}

dev.off()



