#PURPOSE: Reviewer comments seemed most interested in questions related to connectivity, so addressing some of those questions here.

#Reviewer questions:
	#Do lineage-specific genes have fewer connections than other genes?
	#Is there a difference in how central genes on the X vs. autosomes are?

#Additional questions that might be interesting:
	#Something about connectivity and/or module membership of genes evolving rapidly across the whole phylogeny (EVE) vs those evolving rapidly in just one or a few lineages (pairwise comparisons)
	#Connectivity of cis vs trans genes - important not to confound co-expression network with actual regulatory networks, but I would still predict trans to have higher connectivity in general, as reg networks should be one predictor of co-expression
	#Changes in overall connectivity between cell types? Might predict less connectivity late? Not sure, need to think about specific predictions a bit more

#Combine this with ERC? Not in this dataset probably (underpowered, but something to keep in mind for murines and working with Nathan)

#Approach
	#I think the best way to get at most of these questions is to compare median KME values across groups of genes (ex: within a module, is the median KME for X genes significantly higher than that of auto genes?)
	#Probably violin plots or box plots and Wilcoxon rank sum test?
	#Don't think EM's linear model approach is necessary here because we are comparing across genes, not individuals

library(biomaRt)
library(ggplot2)

#Loading WGCNA data
dataset<-"all"
lnames<-load(file=paste("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/03_WGCNA_dataChecks",dataset,"RData",sep="."))
lnames<-c(lnames, load(file=paste0("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/",dataset,"Data-block.1.RData")))
lnames<-c(lnames, load(file=paste("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/05-WGCNA",dataset,"RData",sep=".")))
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

#Getting gene info from biomaRt
mouse_mart<-useMart("ENSEMBL_MART_ENSEMBL", "mmusculus_gene_ensembl")
myIDs<-getBM(attributes=c("ensembl_gene_id","chromosome_name"), mart=mouse_mart)

#Assigning chromosome ID to each gene
myChrs<-c()
chrType<-c()
for(i in rownames(myKME)){
	if(i %in% myIDs$ensembl_gene_id){
		new_chr<-myIDs$chromosome_name[which(myIDs$ensembl_gene_id==i)]
		if(length(new_chr) > 1){
			print(myIDs[which(myIDs$ensembl_gene_id==i), ])
		}
	} else{
		new_chr<-NA
	}
	if(new_chr %in% c(1:19)){
		chrType<-c(chrType, "auto")
	} else{
		chrType<-c(chrType, new_chr)
	}
	myChrs<-c(myChrs, new_chr)
}

KME_wChr<-as.data.frame(cbind(myKME, chr=myChrs, chr_type=chrType))

pdf("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/gene_eigenvalues.XvsAuto.all.pdf", onefile=TRUE)
wilcox.ps<-c()
for(m in interesting_mods){
	print(paste("Working on module", m))
	print(paste("X median:", median(KME_wChr[which(KME_wChr$chr_type=="X"),m])))
	print(paste("auto median:", median(KME_wChr[which(KME_wChr$chr_type=="auto"),m])))
	print("Wilcoxon test:")
	result<-wilcox.test(KME_wChr[which(KME_wChr$chr_type=="X"),m], KME_wChr[which(KME_wChr$chr_type=="auto"),m])
	print(result)
	wilcox.ps<-c(wilcox.ps, result$p.value)
	temp_df<-as.data.frame(KME_wChr[which(KME_wChr$chr_type %in% c("X", "auto")), c("chr_type",m)])
	colnames(temp_df)<-c("chr_type","mod_eig")
	p<-ggplot(temp_df, aes(x=chr_type, y=mod_eig)) + geom_boxplot()
	p<-p + labs(title=paste("Per-gene eigengene module membership: X vs auto,", m), y="Eigengene value")
	p<-p + theme(axis.text = element_text(size=18), axis.title = element_text(size=21), plot.title = element_text(size=36))
	p<-p + theme_minimal()
	print(p)
	#print(paste("X median, absolute value:", median(abs(KME_wChr[which(KME_wChr$chr_type=="X"),m]))))
        #print(paste("auto median, absolute value:", median(abs(KME_wChr[which(KME_wChr$chr_type=="auto"),m]))))
        #print("Wilcoxon test, absolute value:")
        #print(wilcox.test(abs(KME_wChr[which(KME_wChr$chr_type=="X"),m]), abs(KME_wChr[which(KME_wChr$chr_type=="auto"),m])))
}
dev.off()

cor.ps<-p.adjust(wilcox.ps, method="fdr")
names(cor.ps)<-interesting_mods
print(cor.ps)

print("Done with 09_gene_connectivity.r")
