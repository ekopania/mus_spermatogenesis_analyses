#PURPOSE: Run automatic clustering and network analysis with WGCNA

args<-commandArgs(TRUE)
if(length(args) != 2){
        stop("Missing command line arguments.\nArgument 1: dataset (options: all, LZ, RS, LZind, RSind)\nArgument 2: signed or unsigned? (options: TRUE, FALSE)")
}

dataset<-args[1] #all #LZ #RS #LZind #RSind
signed<-args[2] #FALSE #Perform signed network analysis? SIGN RECOMMENDED

print(paste("Performing WGCNA analysis for",dataset,"; signed =",signed, sep="."))

library(WGCNA)

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads=8)

lnames<-load(file=paste("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/03_WGCNA_dataChecks",dataset,"RData",sep=".")) #Load data from previous script
allTraits<-read.table("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/traits.txt",header=TRUE)
if(dataset=="all"){
	myTraits<-allTraits
	stPower<-10
#} else if(dataset=="LZ"){
#	myTraits<-allTraits[grepl("LZ", rownames(allTraits)),]
#	stPower<-
#} else if(dataset=="RS"){
#	myTraits<-allTraits[grepl("RS", rownames(allTraits)),]
#        stPower<-
#} else if(dataset=="LZind"){
#	myTraits<-allTraits[grepl("LZ", rownames(allTraits)),]
#	stPower<-
#} else if(dataset=="RSind"){
#        myTraits<-allTraits[grepl("RS", rownames(allTraits)),]
#	stPower<-
} else{
	stop("Invalid dataset: must be 'all', 'LZ', 'RS', 'LZind', 'RSind'")
}

print("Performing network analysis...")
#power=stPower (where R^2 levels off, different for each dataset; see softThresholding_power_curves.<dataset>.pdf from 03_WGCNA_datachecks)
#maxBlockSize = number of genes; can make smaller if taking too much memory
if(signed){
	myNetwork<-blockwiseModules(final_expData, corType="pearson", power = stPower, networkType = "signed", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = paste0(dataset,"Data"), verbose = 3, maxBlockSize = 12951)
} else{
	myNetwork<-blockwiseModules(final_expData, corType="pearson", power = stPower, networkType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = paste0(dataset,"Data"), verbose = 3, maxBlockSize = 12951)
}

if(signed){
	pdf(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/",dataset,"Data.network.signed.pdf"))
} else{
	pdf(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/",dataset,"Data.network.pdf"))
}
mergedColors<-labels2colors(myNetwork$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(myNetwork$dendrograms[[1]], mergedColors[myNetwork$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

#Get some info from network analysis
moduleLabels<-myNetwork$colors
moduleColors<-labels2colors(myNetwork$colors)
MEs<-myNetwork$MEs
geneTree<-myNetwork$dendrograms[[1]]
#Get relationship between label (number) and color
n_mods<-length(unique(moduleLabels))
labelOrder<-c(0:(n_mods-1))
colorOrder<-c()
for(i in labelOrder){
	colorOrder<-c(colorOrder, unique(moduleColors[which(moduleLabels==i)]))
}
label_to_color<-as.data.frame(cbind(lab=labelOrder, col=colorOrder))

#Get per-gene module membership (also known as kME or signed eigengene-based connectivity)
myKME<-signedKME(final_expData, MEs, outputColumnName = "kME", corFnc = "cor", corOptions = "use = 'p'")
#Get column names as module colors for kME
colname_nums<-sapply(colnames(myKME), function(x) gsub("kME","",x))
colname_cols<-sapply(colname_nums, function(x) label_to_color$col[which(label_to_color$lab==x)])
colnames(myKME)<-colname_cols

#Save for future scripts/analyses
save(MEs, moduleLabels, moduleColors, geneTree, label_to_color, myKME, file=paste("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/05-WGCNA",dataset,"RData",sep="."))

print("Testing for associations with trait data...")
# Define numbers of genes and samples
nGenes<-ncol(final_expData)
nSamples<-nrow(final_expData)
# Recalculate MEs with color labels
MEs0<-moduleEigengenes(final_expData, moduleColors)$eigengenes
MEs<-orderMEs(MEs0)
#Convert categorical traits to binary
myTraitsBin<-binarizeCategoricalColumns.forRegression(myTraits)
print(head(myTraitsBin))
#Calculate significant associations between trait and expression module
moduleTraitCor<-cor(MEs, myTraitsBin, use = "p")
moduleTraitPvalue_raw<-corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitPvalue<-p.adjust(moduleTraitPvalue_raw, method="fdr")

if(signed){
	pdf(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/",dataset,"Data.traitHeatmap.signed.pdf"), height=11, width=8.5)
} else{
	pdf(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/",dataset,"Data.traitHeatmap.pdf"), height=11, width=8.5)
}
if( (dataset=="LZ") || (dataset=="RS") || (dataset=="LZind") || (dataset=="RSind") ){
        myTraitsBin<-myTraitsBin[,which(colnames(myTraitsBin) != "cell_type.RS.vs.all")]
}
# Will display correlations and their p-values
textMatrix<-paste(signif(moduleTraitCor, 2), " (", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix)<-dim(moduleTraitCor)
par(mar = c(9, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(myTraitsBin), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()

print("Done with 05_WGCNA.r")
