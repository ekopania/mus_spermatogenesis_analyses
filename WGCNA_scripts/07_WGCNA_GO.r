#PURPOSE: Perform GO enrichment analysis for interesting modules

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line arguments.\nArgument 1: dataset")
}

dataset<-args[1] #all #LZ #RS

print(paste0("Running GO enrichment analysis for ", dataset))

library(WGCNA)
library(anRichment)
library(biomaRt)

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads=8)

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

#GO analysis
GOcollection<-buildGOcollection(organism = "mouse")
GOenrichment<-enrichmentAnalysis(classLabels = moduleColors, identifiers = entrez, refCollection = GOcollection, useBackground = "given", threshold = 1e-4, thresholdType = "Bonferroni", getOverlapEntrez = TRUE, getOverlapSymbols = TRUE, ignoreLabels = "grey")
collectGarbage()
write.csv(GOenrichment$enrichmentTable, file=paste("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/GOenrichment-enrichmentTable",dataset,"csv", sep="."), row.names=FALSE)

print("Done with 07_WGCNA_GO.r") 
