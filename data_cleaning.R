#FILTERING, CLEANING, ORGANIZING DATA

#making updated GWAS
single_complete_disease_gwas<-c()
for(i in 1:nrow(complete_disease_gwas)){
  split_data<-c()
  if(grepl(",",complete_disease_gwas$SNP_GENE_IDS[i])==FALSE){
    single_complete_disease_gwas<-rbind(single_complete_disease_gwas,complete_disease_gwas[i,])
  }else{
    genes<-unlist(strsplit(as.character(complete_disease_gwas$SNP_GENE_IDS[i]),","))
    genes<-trimws(genes)
    for(gene in genes){
      single_complete_disease_gwas<-rbind(single_complete_disease_gwas,complete_disease_gwas[i,])
      single_complete_disease_gwas$SNP_GENE_IDS[nrow(single_complete_disease_gwas)]<-gene
    }
  }
}

#Updating Qualified Genes:
snp_ids<-c()
for(i in 1:nrow(qualified_genes)){
  snp<-single_complete_disease_gwas$SNPS[which(single_complete_disease_gwas$CHR_ID==qualified_genes$Chr.Number[i]&single_complete_disease_gwas$CHR_POS==qualified_genes$SNP.Position[i]&single_complete_disease_gwas$SNP_GENE_IDS==qualified_genes$Genes[i])]
  if(length(snp)>1){
    snp_ids<-c(snp_ids,snp[1])
  }else{
    snp_ids<-c(snp_ids,snp)
  }
}
#adding traits:
traits<-c()
for(i in 1:nrow(qualified_genes)){
  trait<-single_complete_disease_gwas$DISEASE.TRAIT[which(single_complete_disease_gwas$SNPS==qualified_genes$snp_ids[i])]
  if(length(trait)>1){
    traits<-c(traits,trait[1])
  }else{
    traits<-c(traits,trait)
  }
}