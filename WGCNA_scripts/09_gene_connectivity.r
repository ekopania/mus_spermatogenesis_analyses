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

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line arguments.\nArgument 1: dataset")
}

dataset<-args[1]

print(paste0("Running gene connectivity analysis for ", dataset))

library(biomaRt)
library(ggbeeswarm)
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

pdf(paste("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/gene_eigenvalues.XvsAuto",dataset,"pdf", sep="."), onefile=TRUE)
wilcox.ps<-c()
for(m in interesting_mods){
	print(paste("Working on module", m))
	print(paste("X median:", median(KME_wChr[which(KME_wChr$chr_type=="X"),m], na.rm=TRUE)))
	print(paste("auto median:", median(KME_wChr[which(KME_wChr$chr_type=="auto"),m], na.rm=TRUE)))
	print("Wilcoxon test:")
	result<-wilcox.test(KME_wChr[which(KME_wChr$chr_type=="X"),m], KME_wChr[which(KME_wChr$chr_type=="auto"),m])
	print(result)
	wilcox.ps<-c(wilcox.ps, result$p.value)
	temp_df<-as.data.frame(KME_wChr[which(KME_wChr$chr_type %in% c("X", "auto")), c("chr_type",m)])
	colnames(temp_df)<-c("chr_type","mod_eig")
	p<-ggplot(temp_df, aes(x=chr_type, y=mod_eig)) + geom_violin() + geom_boxplot(width=0.1) + geom_quasirandom()
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

print("FDR-corrected p-values from Wilcoxon rank sum tests, X vs auto:")
cor.ps<-p.adjust(wilcox.ps, method="fdr")
names(cor.ps)<-interesting_mods
print(cor.ps)

print("Lineage-specific genes:")
#Get lineage-specific genes
ls_path<-"/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/LINEAGE_SPECIFIC_ZIGZAG/"
if(dataset=="all"){
	domAutoLZ<-scan(paste0(ls_path,"lineage_specific_genes.spec.dom_only_auto.LZ.txt"), what=character())
	musAutoLZ<-scan(paste0(ls_path,"lineage_specific_genes.spec.mus_only_auto.LZ.txt"), what=character())
	mus_domAutoLZ<-scan(paste0(ls_path,"lineage_specific_genes.spec.mus_dom_only_auto.LZ.txt"), what=character())
	sprAutoLZ<-scan(paste0(ls_path,"lineage_specific_genes.spec.spr_only_auto.LZ.txt"), what=character())
	domXLZ<-scan(paste0(ls_path,"lineage_specific_genes.spec.dom_only_x.LZ.txt"), what=character())
	musXLZ<-scan(paste0(ls_path,"lineage_specific_genes.spec.mus_only_x.LZ.txt"), what=character())
	mus_domXLZ<-scan(paste0(ls_path,"lineage_specific_genes.spec.mus_dom_only_x.LZ.txt"), what=character())
	sprXLZ<-scan(paste0(ls_path,"lineage_specific_genes.spec.spr_only_x.LZ.txt"), what=character())
	domAutoRS<-scan(paste0(ls_path,"lineage_specific_genes.spec.dom_only_auto.RS.txt"), what=character())
        musAutoRS<-scan(paste0(ls_path,"lineage_specific_genes.spec.mus_only_auto.RS.txt"), what=character())
        mus_domAutoRS<-scan(paste0(ls_path,"lineage_specific_genes.spec.mus_dom_only_auto.RS.txt"), what=character())
        sprAutoRS<-scan(paste0(ls_path,"lineage_specific_genes.spec.spr_only_auto.RS.txt"), what=character())
        domXRS<-scan(paste0(ls_path,"lineage_specific_genes.spec.dom_only_x.RS.txt"), what=character())
        musXRS<-scan(paste0(ls_path,"lineage_specific_genes.spec.mus_only_x.RS.txt"), what=character())
        mus_domXRS<-scan(paste0(ls_path,"lineage_specific_genes.spec.mus_dom_only_x.RS.txt"), what=character())
        sprXRS<-scan(paste0(ls_path,"lineage_specific_genes.spec.spr_only_x.RS.txt"), what=character())
	domAllLZ<-c(domAutoLZ,domXLZ)
        musAllLZ<-c(musAutoLZ,musXLZ)
        mus_domAllLZ<-c(mus_domAutoLZ,mus_domXLZ)
        sprAllLZ<-c(sprAutoLZ,sprXLZ)
	domAllRS<-c(domAutoRS,domXRS)
        musAllRS<-c(musAutoRS,musXRS)
        mus_domAllRS<-c(mus_domAutoRS,mus_domXRS)
        sprAllRS<-c(sprAutoRS,sprXRS)

	#Assigning lineage-specificity to each gene
        myLS<-c()
        lsType<-c()
        for(i in rownames(myKME)){
                if(i %in% domAllLZ){
                	if(i %in% domAllRS){
				myLS<-c(myLS, "domBoth")
                        	lsType<-c(lsType, "lin-specBoth")
                        } else{
				myLS<-c(myLS, "domLZ")
				lsType<-c(lsType, "lin-specLZ")
                	}
		} else if(i %in% musAllLZ){
                        if(i %in% musAllRS){
				myLS<-c(myLS, "musBoth")
				lsType<-c(lsType, "lin-specBoth")
			} else{
				myLS<-c(myLS, "musLZ")
				lsType<-c(lsType, "lin-specLZ")
			}
                } else if(i %in% mus_domAllLZ){
                        if(i %in% mus_domAllRS){
				 myLS<-c(myLS, "mus_domBoth")
				lsType<-c(lsType, "lin-specBoth")
			} else{
				myLS<-c(myLS, "mus_domLZ")
				lsType<-c(lsType, "lin-specLZ")
			}
                } else if(i %in% sprAllLZ){
                        if(i %in% sprAllRS){
				myLS<-c(myLS, "sprBoth")
				lsType<-c(lsType, "lin-specBoth")
			} else{
				myLS<-c(myLS, "sprLZ")
				lsType<-c(lsType, "lin-specLZ")
			}
                } else if(i %in% domAllRS){
                        myLS<-c(myLS, "domRS")
                        lsType<-c(lsType, "lin-specRS")
                } else if(i %in% musAllRS){
                        myLS<-c(myLS, "musRS")
                        lsType<-c(lsType, "lin-specRS")
                } else if(i %in% mus_domAllRS){
                        myLS<-c(myLS, "mus_domRS")
                        lsType<-c(lsType, "lin-specRS")
                } else if(i %in% sprAllRS){
                        myLS<-c(myLS, "sprRS")
                        lsType<-c(lsType, "lin-specRS")
                } else{
                        myLS<-c(myLS, NA)
                        lsType<-c(lsType, "NOT_lin-spec")
                }
        }

        KME_wLS<-as.data.frame(cbind(myKME, LS=myLS, lsType))
        print(head(KME_wLS))

} else{
	domAuto<-scan(paste0(ls_path,"lineage_specific_genes.spec.dom_only_auto.",dataset,".txt"), what=character())
	musAuto<-scan(paste0(ls_path,"lineage_specific_genes.spec.mus_only_auto.",dataset,".txt"), what=character())
	mus_domAuto<-scan(paste0(ls_path,"lineage_specific_genes.spec.mus_dom_only_auto.",dataset,".txt"), what=character())
	sprAuto<-scan(paste0(ls_path,"lineage_specific_genes.spec.spr_only_auto.",dataset,".txt"), what=character())
	domX<-scan(paste0(ls_path,"lineage_specific_genes.spec.dom_only_x.",dataset,".txt"), what=character())
	musX<-scan(paste0(ls_path,"lineage_specific_genes.spec.mus_only_x.",dataset,".txt"), what=character())
	mus_domX<-scan(paste0(ls_path,"lineage_specific_genes.spec.mus_dom_only_x.",dataset,".txt"), what=character())
	sprX<-scan(paste0(ls_path,"lineage_specific_genes.spec.spr_only_x.",dataset,".txt"), what=character())
	domAll<-c(domAuto,domX)
	musAll<-c(musAuto,musX)
	mus_domAll<-c(mus_domAuto,mus_domX)
	sprAll<-c(sprAuto,sprX)
	
	#Assigning lineage-specificity to each gene
	myLS<-c()
	lsType<-c()
	for(i in rownames(myKME)){
		if(i %in% domAll){
			myLS<-c(myLS, "dom")
			lsType<-c(lsType, "lin-spec")
		} else if(i %in% musAll){
			myLS<-c(myLS, "mus")
	                lsType<-c(lsType, "lin-spec")
		} else if(i %in% mus_domAll){
			myLS<-c(myLS, "mus_dom")
			lsType<-c(lsType, "lin-spec")
		} else if(i %in% sprAll){
			myLS<-c(myLS, "spr")
	                lsType<-c(lsType, "lin-spec")
		} else{
			myLS<-c(myLS, NA)
	                lsType<-c(lsType, "NOT_lin-spec")
		}
	}
	
	KME_wLS<-as.data.frame(cbind(myKME, LS=myLS, lsType))
	print(head(KME_wLS))
}

pdf(paste("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/gene_eigenvalues.lineage-specific",dataset,"pdf", sep="."), onefile=TRUE)
wilcox.ps<-c()
for(m in interesting_mods){
        print(paste("Working on module", m))
        temp_df<-as.data.frame(KME_wLS[, c("lsType",m)])
	colnames(temp_df)<-c("lsType","mod_eig")
	temp_df$lsType<-as.factor(temp_df$lsType)
	temp_df$mod_eig<-as.numeric(as.character(temp_df$mod_eig))
	print(head(temp_df))
	if(dataset=="all"){
		print(paste("lineage-specific in both median:", median(temp_df$mod_eig[which(temp_df$lsType=="lin-specBoth")], na.rm=TRUE)))
		print(paste("lineage-specific early median:", median(temp_df$mod_eig[which(temp_df$lsType=="lin-specLZ")], na.rm=TRUE)))
		print(paste("lineage-specific late median:", median(temp_df$mod_eig[which(temp_df$lsType=="lin-specRS")], na.rm=TRUE)))
		print(paste("not lineage-specific median:", median(temp_df$mod_eig[which(temp_df$lsType=="NOT_lin-spec")], na.rm=TRUE)))
		print("Wilcox tests")
		result_LZvRS<-wilcox.test(temp_df$mod_eig[which(temp_df$lsType=="lin-specLZ")], temp_df$mod_eig[which(temp_df$lsType=="lin-specRS")])
		result_LZ<-wilcox.test(temp_df$mod_eig[which(temp_df$lsType=="lin-specLZ")], temp_df$mod_eig[which(temp_df$lsType=="NOT_lin-spec")])
		result_RS<-wilcox.test(temp_df$mod_eig[which(temp_df$lsType=="lin-specRS")], temp_df$mod_eig[which(temp_df$lsType=="NOT_lin-spec")])
		print(result_LZvRS)
		print(result_LZ)
		print(result_RS)
		result<-pairwise.wilcox.test(x=temp_df$mod_eig, g=temp_df$lsType, p.adjust.method="none")
		print(result)
		wilcox.ps<-c(wilcox.ps, result$p.value[1], result$p.value[2], result$p.value[3], result$p.value[5], result$p.value[6], result$p.value[9])
	} else{
		print(paste("lineage-specific median:", median(temp_df$mod_eig[which(temp_df$lsType=="lin-spec")], na.rm=TRUE)))
		print(paste("not lineage-specific median:", median(temp_df$mod_eig[which(temp_df$lsType=="NOT_lin-spec")], na.rm=TRUE)))
		result<-wilcox.test(temp_df$mod_eig[which(temp_df$lsType=="lin-spec")], temp_df$mod_eig[which(temp_df$lsType=="NOT_lin-spec")])
        	print(result)
        	wilcox.ps<-c(wilcox.ps, result$p.value)
        }
        p<-ggplot(temp_df, aes(x=lsType, y=mod_eig)) + geom_violin() + geom_boxplot(width=0.1) + geom_quasirandom()
        p<-p + labs(title=paste("Per-gene eigengene module membership: lineage-specific vs not,", m), y="Eigengene value")
        p<-p + theme(axis.text = element_text(size=18), axis.title = element_text(size=21), plot.title = element_text(size=36))
        p<-p + theme_minimal()
        print(p)
}
dev.off()

print("FDR-corrected p-values from Wilcoxon rank sum tests, lineage-specific:")
cor.ps<-p.adjust(wilcox.ps, method="fdr")
if(dataset=="all"){
	interesting_mods_rep<-sapply(interesting_mods, function(x) rep(x, 6))
	myComparisons<-rep(c("bothVlz","bothVrs","bothVnone","lzVrs","lsVnone","rsVnone"), length(interesting_mods))
	myLabels<-paste(interesting_mods_rep, myComparisons, sep="_")
	names(cor.ps)<-myLabels
} else{
	names(cor.ps)<-interesting_mods
}
print(cor.ps)

#Lineage specificity by strain
print("Plotting lineage-specificity by strains")
pdf(paste("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/gene_eigenvalues.lineage-specific.byLineage",dataset,"pdf", sep="."), onefile=TRUE)
wilcox.ps<-c()
allCategories<-unique(KME_wLS$LS)
for(m in interesting_mods){
        print(paste("Working on module", m))
        temp_df<-as.data.frame(KME_wLS[, c("LS",m)])
        colnames(temp_df)<-c("LS","mod_eig")
        temp_df$LS<-as.factor(temp_df$LS)
        temp_df$mod_eig<-as.numeric(as.character(temp_df$mod_eig))
        print(head(temp_df))
	for(i in allCategories){
		print(paste("lineage-specific median for", i, median(temp_df$mod_eig[which(temp_df$LS==i)])))
	}
	result<-pairwise.wilcox.test(x=temp_df$mod_eig, g=temp_df$LS, p.adjust.method="fdr")
	print(result)
	p<-ggplot(temp_df, aes(x=LS, y=mod_eig)) + geom_violin() + geom_boxplot(width=0.1) + geom_quasirandom()
	p<-p + labs(title=paste("Per-gene eigengene module membership: lineage-specific vs not by lineage,", m), y="Eigengene value")
	p<-p + theme(axis.text = element_text(size=18), axis.title = element_text(size=21), plot.title = element_text(size=36))
	p<-p + theme_minimal()
	print(p)
}
dev.off()

print("Done with 09_gene_connectivity.r")
