#
library(org.Hs.eg.db)

ENSEMBL2EG=unlist(as.list(org.Hs.egENSEMBL2EG)) # Map Ensembl gene accession numbers with Entrez Gene identifiers
ENSEMBL2EG_df=data.frame(ensembl=names(ENSEMBL2EG), entrez=as.vector(ENSEMBL2EG))

load("/athena/masonlab/scratch/users/nai2008/genes.rda")

ig=read.csv('/athena/masonlab/scratch/users/nai2008/Imprinted_genes/imprinted_human_genes.csv')
ig$Gene=as.vector(ig$Gene)
ig$Gene[which(ig$Gene=='TNDM')]='DMTN'
ig$Gene[which(ig$Gene=='INPP5F V2')]='INPP5F'
ig$Gene[which(ig$Gene=='PSIMCT-1')]='MCTS2P'

mm=match(ig$Gene, genes$Gene.Name)
ig_pt1=ig[which(!is.na(mm)),]
ig_pt2=ig[which(is.na(mm)),]

mm=match(ig_pt2$Gene, genes$Gene.Synonym)
ig_pt2A=ig_pt2[which(!is.na(mm)),]
ig_pt2B=ig_pt2[which(is.na(mm)),]

for(i in 1:nrow(ig_pt2A)){
	m=match(ig_pt2A$Gene[i], genes$Gene.Synonym)
	ig_pt2A$Gene[i]=as.vector(genes$Gene.Name[m])
}

ig=rbind(ig_pt1,ig_pt2A)

mm=match(ig$Gene, genes$Gene.Name)
ig$Ensembl.ID=genes$Ensembl.ID[mm]

mm=match(ig$Ensembl.ID,ENSEMBL2EG_df$ensembl)
ig$Entrez=ENSEMBL2EG_df$entrez[mm]

# remove duplicates
ig=ig[-which(duplicated(ig$Ensembl.ID)=='TRUE'),]

# lets manually look up and add the Entrez IDs for the genes that didn't map

ig_pt2B$Ensembl.ID=NA
ig_pt2B$Entrez=NA

# ig_pt2B:

# COPG2IT1 = 53844 (no Ensembl available on GeneCards)
which(ENSEMBL2EG_df$entrez == 53844) # NA
ig_pt2B$Entrez[which(ig_pt2B$Gene=='COPG2IT1')]='53844'

# PWCR1 = 692236 (no Ensembl available on GeneCards)
which(ENSEMBL2EG_df$entrez == 692236) # NA
ig_pt2B$Entrez[which(ig_pt2B$Gene=='PWCR1')]='692236'

ig=rbind(ig, ig_pt2B)

save(ig, file='/athena/masonlab/scratch/users/nai2008/Imprinted_genes/imprinted_genes.rda')

























