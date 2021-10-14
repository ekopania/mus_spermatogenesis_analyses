#PURPOSE: Test if interesting modules are over or under-enriched for sex chromosomes

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line arguments.\nArgument 1: dataset")
}

dataset<-args[1]

print(paste0("Running sex chromosome enrichment analysis for ", dataset))

library(WGCNA)
library(biomaRt)

options(stringsAsFactors = FALSE)

#Load data from previous script
lnames<-load(file=paste("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/03_WGCNA_dataChecks",dataset,"RData",sep="."))
lnames<-c(lnames, load(file=paste0("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/",dataset,"Data-block.1.RData")))
lnames<-c(lnames, load(file=paste("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/05-WGCNA",dataset,"RData",sep=".")))
allTraits<-read.table("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/traits.txt",header=TRUE)
if(dataset=="all"){
        interesting_mods<-c("blue","turquoise","green","yellow","red","brown")
#} else if(dataset=="LZ"){
#        interesting_mods<-c()
#} else if(dataset=="RS"){
#        interesting_mods<-c("")
#} else if(dataset=="LZind"){
#        interesting_mods<-c("")
#} else if(dataset=="RSind"){
#        interesting_mods<-c("")
} else{
        stop("Invalid dataset: must be 'all', 'LZ', 'RS', 'LZind', or 'RSind'")
}

#Get Entrez IDs
mouse_mart<-useMart("ENSEMBL_MART_ENSEMBL", "mmusculus_gene_ensembl")
myIDs<-getBM(attributes=c("ensembl_gene_id","entrezgene_id","chromosome_name"), mart=mouse_mart)
entrez<-c()
for(gene in names(moduleLabels)){
	e<-myIDs$entrezgene_id[which(myIDs$ensembl_gene_id==gene)]
	entrez<-c(entrez, e)
}
table(is.finite(entrez))

#Get sex chromosome genes
XchrGenes_df<-myIDs[which(myIDs$chromosome_name=="X"), ]
XchrGenes_names<-unique(XchrGenes_df$ensembl_gene_id)
XchrGenes_inDataset<-XchrGenes_names[which(XchrGenes_names %in% names(moduleLabels))]
YchrGenes_df<-myIDs[which(myIDs$chromosome_name=="Y"), ]
YchrGenes_names<-unique(YchrGenes_df$ensembl_gene_id)
YchrGenes_inDataset<-YchrGenes_names[which(YchrGenes_names %in% names(moduleLabels))]

#Binomial Test
Xbinom_ps<-c()
Ybinom_ps<-c()
for(m in interesting_mods){
        print(m)
        mod_genes<-names(moduleLabels)[which(moduleColors==m)]
        print("X chromosome:")
        X_only<-XchrGenes_inDataset[which(!(XchrGenes_inDataset %in% mod_genes))]
        mod_only<-mod_genes[which(!(mod_genes %in% XchrGenes_inDataset))]
        both_genes<-intersect(mod_genes, XchrGenes_inDataset)
        X_prop<-length(which(names(moduleLabels) %in% XchrGenes_inDataset)) / length(names(moduleLabels))
        result<-binom.test(length(both_genes), length(mod_genes), p=X_prop, alternative="two.sided", conf.level=0.95)
        print(result)
        Xbinom_ps<-c(Xbinom_ps, result$p.value)
	print("Y chromosome:")
        Y_only<-YchrGenes_inDataset[which(!(YchrGenes_inDataset %in% mod_genes))]
        mod_only<-mod_genes[which(!(mod_genes %in% YchrGenes_inDataset))]
        both_genes<-intersect(mod_genes, YchrGenes_inDataset)
	Y_prop<-length(which(names(moduleLabels) %in% YchrGenes_inDataset)) / length(names(moduleLabels))
        result<-binom.test(length(both_genes), length(mod_genes), p=Y_prop, alternative="two.sided", conf.level=0.95)
        print(result)
        Ybinom_ps<-c(Ybinom_ps, result$p.value)
}

Xcor_ps_b<-p.adjust(Xbinom_ps)
names(Xcor_ps_b)<-interesting_mods
print("FDR-corrected p-values for Binomial Tests, X chromosome:")
print(Xcor_ps_b)
Ycor_ps_b<-p.adjust(Ybinom_ps)
names(Ycor_ps_b)<-interesting_mods
print("FDR-corrected p-values for Binomial Tests, Y chromosome:")
print(Ycor_ps_b)

print("Done with 08_WGCNA_sex_chr_enrichment.r")
