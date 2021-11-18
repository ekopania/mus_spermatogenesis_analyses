#PURPOSE: Plot network connectivity (KME values) against EVE divergence and dN/dS

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line arguments.\nArgument 1: dataset")
}

dataset<-args[1]

print(paste0("Running gene connectivity analysis for ", dataset))

library(biomaRt)
library(ggplot2)

#Loading WGCNA data
lnames<-load(file=paste("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/03_WGCNA_dataChecks",dataset,"RData",sep="."))
lnames<-c(lnames, load(file=paste0("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/",dataset,"Data-block.1.RData")))
lnames<-c(lnames, load(file=paste("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/05-WGCNA",dataset,"RData",sep=".")))
if(dataset=="all"){
        interesting_mods<-c("blue","turquoise","green","yellow","red","brown")
} else if(dataset=="LZ"){
        interesting_mods<-c("darkturquoise","grey60","salmon","brown","turquoise","tan","greenyellow","magenta","black","lightyellow","darkred","lightgreen","midnightblue","blue","pink","yellow","darkgreen","orange","royalblue","darkgrey","green","cyan","darkorange","lightcyan","purple","red")
} else if(dataset=="RS"){
        interesting_mods<-c("turquoise","black","magenta","midnightblue","red","cyan","greenyellow","lightcyan","yellow","green","tan","lightyellow","salmon","grey60","lightgreen","purple","blue","brown","pink")
} else{
        stop("Invalid dataset: must be 'all', 'LZ', 'RS', 'LZind', or 'RSind'")
}

#Plot EVE divergence by eigengene value (KME) (only works for cell types separate)
if(dataset != "all"){
	print("EVE divergence comparison")
	eveData<-read.table(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/results/RPKM1/indivBetaMLparams",dataset,"expressed_RPKM1.res"))
	eveGenes<-scan(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/gene_list_eve_",dataset,"_edgeR_wholeGenome.ensemblOrthos.rpkm1.txt"),what=character())
	print(head(eveData))
	print(head(eveGenes))
	eveDiv0<-as.data.frame(cbind(betai=0-log(eveData[,4]),gene=eveGenes))
	eveDiv<-eveDiv0[which(as.numeric(as.character(eveDiv0$betai))>-5),]
	
	eveDiv_filtered<-eveDiv[which(eveDiv$gene %in% rownames(myKME)),]
	eveDiv_sorted<-eveDiv_filtered[order(match(eveDiv_filtered$gene, rownames(myKME))),]
	
	myKME_filtered<-myKME[which(rownames(myKME) %in% eveDiv_sorted$gene), ]
	
	stopifnot(all.equal(eveDiv_sorted$gene, rownames(myKME_filtered)))
	KME_wDiv<-as.data.frame(cbind(eveDiv_sorted, myKME_filtered))
	print(head(KME_wDiv))
		
	pdf(paste("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/gene_eigenvaluesVSeveDiv",dataset,"pdf", sep="."), onefile=TRUE)
	spearman.ps<-c()
	for(m in interesting_mods){
		print(paste("Working on module", m))
	        print("Spearman's rank correlation test:")
	        result<-cor.test(as.numeric(as.character(KME_wDiv[,"betai"])), as.numeric(as.character(KME_wDiv[,m])), method="spearman")
	        print(result)
	        spearman.ps<-c(spearman.ps, result$p.value)
	        temp_df<-as.data.frame(KME_wDiv[, c("betai",m)])
	        colnames(temp_df)<-c("betai","mod_eig")
	        p<-ggplot(temp_df, aes(x=as.numeric(as.character(betai)), y=as.numeric(as.character(mod_eig)))) + geom_point()
	        p<-p + labs(title=paste("Per-gene eigengene module membership vs EVE divergence,", m), y="Eigengene value", x="EVE divergence")
	        p<-p + theme(axis.text = element_text(size=18), axis.title = element_text(size=21), plot.title = element_text(size=36))
	        p<-p + theme_minimal()
	        print(p)
	}
	dev.off()
	
	print("FDR-corrected p-values from Spearman's rank correlation tests, EVE divergence:")
	cor.ps<-p.adjust(spearman.ps, method="fdr")
	names(cor.ps)<-interesting_mods
	print(cor.ps)
}

print("Repeating for dN/dS")
dnds0<-read.table("/mnt/beegfs/ek112884/PAML_runs/WHOLE_GENOMES/omega_list.txt", col.names=c("gene","dN","dS","dN.dS"))
dnds<-dnds0[which(as.numeric(as.character(dnds0$dN.dS)) < 1.5),]

dnds_filtered<-dnds[which(dnds$gene %in% rownames(myKME)),]
dnds_sorted<-dnds_filtered[order(match(dnds_filtered$gene, rownames(myKME))),]

myKME_filtered<-myKME[which(rownames(myKME) %in% dnds_sorted$gene), ]

stopifnot(all.equal(dnds_sorted$gene, rownames(myKME_filtered)))
KME_wDnds<-as.data.frame(cbind(dnds_sorted, myKME_filtered))
print(head(KME_wDnds))

pdf(paste("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/gene_eigenvaluesVSdnds",dataset,"pdf", sep="."), onefile=TRUE)
spearman.ps<-c()
for(m in interesting_mods){
        print(paste("Working on module", m))
        print("Spearman's rank correlation test:")
        result<-cor.test(as.numeric(as.character(KME_wDnds[,"dN.dS"])), as.numeric(as.character(KME_wDnds[,m])), method="spearman")
        print(result)
        spearman.ps<-c(spearman.ps, result$p.value)
        temp_df<-as.data.frame(KME_wDnds[, c("dN.dS",m)])
        colnames(temp_df)<-c("dN.dS","mod_eig")
        p<-ggplot(temp_df, aes(x=as.numeric(as.character(dN.dS)), y=as.numeric(as.character(mod_eig)))) + geom_point()
        p<-p + labs(title=paste("Per-gene eigengene module membership vs EVE divergence,", m), y="Eigengene value", x="dN/dS")
        p<-p + theme(axis.text = element_text(size=18), axis.title = element_text(size=21), plot.title = element_text(size=36))
        p<-p + theme_minimal()
        print(p)
}
dev.off()

print("FDR-corrected p-values from Spearman's rank correlation tests, dN/dS:")
cor.ps<-p.adjust(spearman.ps, method="fdr")
names(cor.ps)<-interesting_mods
print(cor.ps)
