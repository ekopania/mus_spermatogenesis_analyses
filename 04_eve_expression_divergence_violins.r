#PURPOSE: make expression divergence violin plots as in 02 BUT use EVEmodel beta-i parameter instead of pairwise divergence calculation; EVE beta-i should be a measure of population variation/evolutioanry divergence controlling for phylogeny (instead of pairwise comparison)

library(ggplot2)
library(EnsDb.Mmusculus.v79)

rpkm<-1
myDir<-"WHOLE_GENOME_MULTIMAP_FRACTIONAL"
induced<-FALSE
testis_spec<-FALSE

#Read in data
print("Reading in files...")
if(induced){
	#INDUCED ONLY
	mydata<-read.table(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/",myDir,"/results/RPKM",rpkm,"/indivBetaMLparams_ensemblOrthos_BOTHinduced_RPKM",rpkm,".res"))
	#sharedData<-read.table("results/sharedBetaMLparams_ensemblOrthos_BOTHinduced.res")
	#LZinduced<-read.table("results/indivBetaMLparams_ensemblOrthos_LZinduced.res")
	#RSinduced<-read.table("results/indivBetaMLparams_ensemblOrthos_RSinduced.res")
	#print(head(LZinduced))
	#print(head(RSinduced))
	LZgenes<-scan(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/",myDir,"/gene_list_LZinduced_edgeR_wholeGenome.ensemblOrthos.rpkm",rpkm,".txt"),what=character())
	RSgenes<-scan(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/",myDir,"/gene_list_RSinduced_edgeR_wholeGenome.ensemblOrthos.rpkm",rpkm,".txt"),what=character())
	betai<-cbind(betai=0-log(mydata[,4]),group=c(rep("LZ",3375),rep("RS",2769)),gene=c(LZgenes,RSgenes))
} else{
	#ALL GENES
	mydataLZ<-read.table(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/",myDir,"/results/RPKM",rpkm,"/indivBetaMLparamsLZexpressed_RPKM",rpkm,".res"))
	mydataRS<-read.table(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/",myDir,"/results/RPKM",rpkm,"/indivBetaMLparamsRSexpressed_RPKM",rpkm,".res"))
	mydata<-as.data.frame(rbind(mydataLZ,mydataRS))
	#mydata<-read.table(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/",myDir,"/results/RPKM",rpkm,"/indivBetaMLparamsBOTHexpressed_RPKM",rpkm,".res"))
	LZgenes<-scan(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/",myDir,"/gene_list_eve_LZ_edgeR_wholeGenome.ensemblOrthos.rpkm",rpkm,".txt"),what=character())
	RSgenes<-scan(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/",myDir,"/gene_list_eve_RS_edgeR_wholeGenome.ensemblOrthos.rpkm",rpkm,".txt"),what=character())
	betai<-cbind(betai=0-log(mydata[,4]),group=c(rep("LZ",nrow(mydataLZ)),rep("RS",nrow(mydataRS))),gene=c(LZgenes,RSgenes))
}
if(testis_spec){
	testis_genes<-scan(gzfile("/mnt/beegfs/ek112884/mus_expression_analysis/chalmel_testis_specific.txt.gz"),what=character())
        testis_genes<-unique(testis_genes)
	betai<-betai[which(as.character(betai[,"gene"]) %in% testis_genes),]
}
#print(head(mydata))
print(head(betai))
print(dim(betai))
#print(betaShared)

#combo<-as.data.frame(rbind(LZinduced_betai,RSinduced_betai))
#combo<-combo[which(as.numeric(as.character(combo$betai))<10),]
combo<-as.data.frame(betai)
combo<-combo[which(as.numeric(as.character(combo$betai))>-5),]
print(dim(combo))
print(head(combo))
print("Number of unique genes in dataset:")
print(length(unique(combo$gene)))

print("Plotting...")
p<-ggplot(combo, aes(x=group, y=as.numeric(as.character(betai)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Gene expression divergence (EVEmodel) in early and late spermatogenesis", x="Cell type", y="Expression Divergence")
p<-p + theme(axis.text.y = element_text(size=20))
p<-p + theme_minimal()
p<-p + scale_fill_manual(values=c("chocolate1","lightsteelblue1"))
if(induced){
	ggsave("EVEmodel_divergence_violins_wholeGenome_multiMap.ensemblOrthos.noOutliers.pdf",p,width=11,height=8.5,units="in")
} else if(testis_spec){
	ggsave("EVEmodel_divergence_violins_wholeGenome_multiMap.ensemblOrthos.noOutliers.testisSpecificGenes.pdf",p,width=11,height=8.5,units="in")
} else{
	ggsave("EVEmodel_divergence_violins_wholeGenome_multiMap.ensemblOrthos.noOutliers.allGenes.pdf",p,width=11,height=8.5,units="in")
}

print("Median EVE divergence, autos early:")
print(median(as.numeric(as.character(combo$betai[which(combo$group=="LZ.auto")]))))
print("# genes, autos early ")
print(length(combo$betai[which(combo$group=="LZ.auto")]))

print("Here are the Wilcox test results:")
wilcox.test(as.numeric(as.character(combo[which(combo$group=="LZ"),]$betai)),as.numeric(as.character(combo[which(combo$group=="RS"),]$betai)))#, alternative="less")

#Repeat for alpha
#LZinduced_alpha<-cbind(alpha=log(LZinduced[,3]),group=rep("LZ",nrow(LZinduced)))
#RSinduced_alpha<-cbind(alpha=log(RSinduced[,3]),group=rep("RS",nrow(RSinduced)))
#print(head(LZinduced_alpha))
#print(head(RSinduced_alpha))

#alpha<-cbind(alpha=0-log(mydata[,3]),group=c(rep("LZ",3375),rep("RS",2769)))
#print(head(alpha))

#combo<-as.data.frame(rbind(LZinduced_alpha,RSinduced_alpha))
#combo<-combo[which(as.numeric(as.character(combo$alpha))<10),]
#print(dim(combo))
#print(head(combo))
#combo<-as.data.frame(alpha)
#combo<-combo[which(as.numeric(as.character(combo$alpha))>-5),]
#combo<-combo[which(as.numeric(as.character(combo$alpha)) < exp(5)),]
#print(dim(combo))
#print(head(combo))

#print("Plotting...")
#p<-ggplot(combo, aes(x=group, y=as.numeric(as.character(alpha)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Rate of adaptation (EVEmodel) in early and late spermatogenesis", x="Cell type", y="Rate of adaptation")
#p<-p + theme_minimal()
#p<-p + theme(axis.text.y = element_text(size=20))
#p<-p + scale_fill_manual(values=c("chocolate1","lightsteelblue1"))
#if(induced){
#	ggsave("EVEmodel_adaptation_rate_violins_wholeGenome_multiMap.ensemblOrthos.noOutliers.pdf",p,width=11,height=8.5,units="in")
#} else{
#	ggsave("EVEmodel_adaptation_rate_violins_wholeGenome_multiMap.ensemblOrthos.noOutliers.allGenes.pdf",p,width=11,height=8.5,units="in")
#}
#
#print("Here are the Wilcox test results:")
#wilcox.test(as.numeric(as.character(combo[which(combo$group=="LZ"),]$alpha)),as.numeric(as.character(combo[which(combo$group=="RS"),]$alpha)))#, alternative="less")

#Repeat for X chr
edb<-EnsDb.Mmusculus.v79
edb_y<-addFilter(edb, SeqNameFilter("Y"))
y_genes<-genes(edb_y)
edb_x<-addFilter(edb, SeqNameFilter("X"))
x_genes<-genes(edb_x)
edb_auto<-addFilter(edb, SeqNameFilter(c(1:19)))
auto_genes<-genes(edb_auto)

#append gene names
#if(induced){
#	LZgenes<-scan("gene_list_LZinduced_edgeR_wholeGenome.ensemblOrthos.txt",what=character())
#	RSgenes<-scan("gene_list_RSinduced_edgeR_wholeGenome.ensemblOrthos.txt",what=character())
#} else if(testis_spec){
#	LZgenes<-as.character(combo$gene[which(combo$group=="LZ")])
#	RSgenes<-as.character(combo$gene[which(combo$group=="RS")])
#} else{
#	LZgenes<-scan("gene_list_eve_LZ_edgeR_wholeGenome.ensemblOrthos.txt",what=character())
#	RSgenes<-scan("gene_list_eve_RS_edgeR_wholeGenome.ensemblOrthos.txt",what=character())
#}
if(testis_spec){
	LZgenes<-as.character(combo$gene[which(combo$group=="LZ")])
	RSgenes<-as.character(combo$gene[which(combo$group=="RS")])
}
gene_names<-c(LZgenes,RSgenes)
print(gene_names)

chr_type<-c()
for(g in gene_names){
	if(g %in% y_genes$gene_id){
		chr_type<-c(chr_type,"Y")
	} else if(g %in% x_genes$gene_id){
		chr_type<-c(chr_type,"X")
	} else if(g %in% auto_genes$gene_id){
		chr_type<-c(chr_type,"auto")
	} else{
		chr_type<-c(chr_type,"NA")
	}
}

betai_final<-as.data.frame(cbind(betai,gene_names,chr_type,ct_chr=paste(betai[,"group"],chr_type,sep=".")))
#alpha_final<-as.data.frame(cbind(alpha,gene_names,chr_type,ct_chr=paste(alpha[,"group"],chr_type,sep=".")))
betai_final<-betai_final[c(which(betai_final$chr_type=="X"),which(betai_final$chr_type=="auto")),]
#alpha_final<-alpha_final[c(which(alpha_final$chr_type=="X"),which(alpha_final$chr_type=="auto")),]

betai_X<-betai_final[which(betai_final$chr_type=="X"),]
betai_auto<-betai_final[which(betai_final$chr_type=="auto"),]
#alpha_X<-alpha_final[which(alpha_final$chr_type=="X"),]
#alpha_auto<-alpha_final[which(alpha_final$chr_type=="auto"),]

#NO OUTLIERS - after neg log
betai_final<-betai_final[which(as.numeric(as.character(betai_final$betai)) > -5),]
betai_X<-betai_X[which(as.numeric(as.character(betai_X$betai)) > -5),]
betai_auto<-betai_auto[which(as.numeric(as.character(betai_auto$betai)) > -5),]
#alpha_final<-alpha_final[which(as.numeric(as.character(alpha_final$alpha)) > -5),]
#alpha_X<-alpha_X[which(as.numeric(as.character(alpha_X$alpha)) > -5),]
#alpha_auto<-alpha_auto[which(as.numeric(as.character(alpha_auto$alpha)) > -5),]
print(head(betai_final))

#NO OUTLIERS - no neg log
#betai_final<-betai_final[which(as.numeric(as.character(betai_final$betai)) < exp(5)),]
#betai_X<-betai_X[which(as.numeric(as.character(betai_X$betai)) < exp(5)),]
#betai_auto<-betai_auto[which(as.numeric(as.character(betai_auto$betai)) < exp(5)),]
#alpha_final<-alpha_final[which(as.numeric(as.character(alpha_final$alpha)) < exp(5)),]
#alpha_X<-alpha_X[which(as.numeric(as.character(alpha_X$alpha)) < exp(5)),]
#alpha_auto<-alpha_auto[which(as.numeric(as.character(alpha_auto$alpha)) < exp(5)),]

print("Plotting...")
#divergence X
p<-ggplot(betai_X, aes(x=group, y=as.numeric(as.character(betai)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Expression divergence (EVEmodel) in early and late spermatogenesis - X chromosome", x="Cell type", y="Expression divergence")
p<-p + theme(axis.text.y = element_text(size=20))
p<-p + theme_minimal()
p<-p + scale_fill_manual(values=c("chocolate1","lightsteelblue1"))
if(induced){
	ggsave("EVEmodel_divergence_violins_wholeGenome_multiMap.ensemblOrthos.Xchr.noOutliers.pdf",p,width=11,height=8.5,units="in")
} else if(testis_spec){
	ggsave("EVEmodel_divergence_violins_wholeGenome_multiMap.ensemblOrthos.Xchr.noOutliers.testisSpecific.Genes.pdf",p,width=11,height=8.5,units="in")
} else{
	ggsave("EVEmodel_divergence_violins_wholeGenome_multiMap.ensemblOrthos.Xchr.noOutliers.allGenes.pdf",p,width=11,height=8.5,units="in")
}
#divergence autos
p<-ggplot(betai_auto, aes(x=group, y=as.numeric(as.character(betai)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Expression divergence (EVEmodel) in early and late spermatogenesis - autosomes", x="Cell type", y="Expression divergence")
p<-p + theme(axis.text.y = element_text(size=20))
p<-p + theme_minimal()
p<-p + scale_fill_manual(values=c("chocolate1","lightsteelblue1"))
if(induced){
	ggsave("EVEmodel_divergence_violins_wholeGenome_multiMap.ensemblOrthos.autos.noOutliers.pdf",p,width=11,height=8.5,units="in")
} else if(testis_spec){
	ggsave("EVEmodel_divergence_violins_wholeGenome_multiMap.ensemblOrthos.autos.noOutliers.testisSpecificGenes.pdf",p,width=11,height=8.5,units="in")
} else{
	ggsave("EVEmodel_divergence_violins_wholeGenome_multiMap.ensemblOrthos.autos.noOutliers.allGenes.pdf",p,width=11,height=8.5,units="in")
}
#divergence combined
p<-ggplot(betai_final, aes(x=ct_chr, y=as.numeric(as.character(betai)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Expression divergence (EVEmodel) in early and late spermatogenesis", x="Cell type", y="Expression divergence")
p<-p + theme(axis.text.y = element_text(size=20))
p<-p + theme_minimal()
#p<-p + scale_x_discrete(limits=c("auto.LZ","X.LZ","auto.RS","X.RS"))
p<-p + scale_fill_manual(values=c("chocolate1","lightsteelblue1"))
#p<-p + geom_hline(yintercept=betaShared, linetype="dashed")
if(induced){
	ggsave("EVEmodel_divergence_violins_wholeGenome_multiMap.ensemblOrthos.XandAutos.noOutliers.pdf",p,width=11,height=8.5,units="in")
} else if(testis_spec){
	ggsave("EVEmodel_divergence_violins_wholeGenome_multiMap.ensemblOrthos.XandAutos.noOutliers.testisSpecificGenes.pdf",p,width=11,height=8.5,units="in")
} else{
	ggsave("EVEmodel_divergence_violins_wholeGenome_multiMap.ensemblOrthos.XandAutos.noOutliers.allGenes.pdf",p,width=11,height=8.5,units="in")
}

#adaptation rate X
#p<-ggplot(alpha_X, aes(x=group, y=as.numeric(as.character(alpha)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Rate of adaptation (EVEmodel) in early and late spermatogenesis - X chromosome", x="Cell type", y="Adaptation rate")
#p<-p + theme(axis.text.y = element_text(size=20))
#p<-p + theme_minimal()
#p<-p + scale_fill_manual(values=c("chocolate1","lightsteelblue1"))
#ggsave("EVEmodel_adaptation_rate_violins_wholeGenome_multiMap.ensemblOrthos.Xchr.noOutliers.pdf",p,width=11,height=8.5,units="in")
#adaptation rate autos
#p<-ggplot(alpha_auto, aes(x=group, y=as.numeric(as.character(alpha)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Rate of adaptation (EVEmodel) in early and late spermatogenesis - autosomes", x="Cell type", y="Adaptation rate")
#p<-p + theme(axis.text.y = element_text(size=20))
#p<-p + theme_minimal()
#p<-p + scale_fill_manual(values=c("chocolate1","lightsteelblue1"))
#ggsave("EVEmodel_adaptation_rate_violins_wholeGenome_multiMap.ensemblOrthos.autos.noOutliers.pdf",p,width=11,height=8.5,units="in")
#adaptation rate combined
#p<-ggplot(alpha_final, aes(x=ct_chr, y=as.numeric(as.character(alpha)))) + geom_violin(aes(fill=group)) + geom_boxplot(width=0.1) + labs(title="Rate of adaptation (EVEmodel) in early and late spermatogenesis", x="Cell type", y="Adaptation rate")
#p<-p + theme(axis.text.y = element_text(size=20))
#p<-p + theme_minimal()
#p<-p + scale_x_discrete(limits=c("auto.LZ","X.LZ","auto.RS","X.RS"))
#p<-p + scale_fill_manual(values=c("chocolate1","lightsteelblue1"))
#ggsave("EVEmodel_adaptation_rate_violins_wholeGenome_multiMap.ensemblOrthos.XandAutos.noOutliers.pdf",p,width=11,height=8.5,units="in")

print("Median EVE divergence, autos early:")
print(median(as.numeric(as.character(betai_final$betai[which(betai_final$ct_chr=="LZ.auto")]))))
print("# genes, autos early:")
print(length(betai_final$betai[which(betai_final$ct_chr=="LZ.auto")]))
print("Median EVE divergence, autos late:")
print(median(as.numeric(as.character(betai_final$betai[which(betai_final$ct_chr=="RS.auto")]))))
print("# genes, autos late:")
print(length(betai_final$betai[which(betai_final$ct_chr=="RS.auto")]))
print("Median EVE divergence, X early:")
print(median(as.numeric(as.character(betai_final$betai[which(betai_final$ct_chr=="LZ.X")]))))
print("# genes, X early:")
print(length(betai_final$betai[which(betai_final$ct_chr=="LZ.X")]))
print("Median EVE divergence, X late:")
print(median(as.numeric(as.character(betai_final$betai[which(betai_final$ct_chr=="RS.X")]))))
print("# genes, X late:")
print(length(betai_final$betai[which(betai_final$ct_chr=="RS.X")]))

print("Here are the Wilcox test results:")
print("Expression divergence early vs late autos:")
wilcox.test(as.numeric(as.character(betai_auto[which(betai_auto$group=="LZ"),]$betai)),as.numeric(as.character(betai_auto[which(betai_auto$group=="RS"),]$betai)))#, alternative="less")
print("Expression divergence early vs late X:")
wilcox.test(as.numeric(as.character(betai_X[which(betai_X$group=="LZ"),]$betai)),as.numeric(as.character(betai_X[which(betai_X$group=="RS"),]$betai)))#, alternative="less")
print("Expression divergence early X vs autos:")
wilcox.test(as.numeric(as.character(betai_X[which(betai_X$group=="LZ"),]$betai)),as.numeric(as.character(betai_auto[which(betai_auto$group=="LZ"),]$betai)))
print("Expression divergence late X vs autos:")
wilcox.test(as.numeric(as.character(betai_X[which(betai_X$group=="RS"),]$betai)),as.numeric(as.character(betai_auto[which(betai_auto$group=="RS"),]$betai)))
pairwise.wilcox.test(x=as.numeric(as.character(betai_final$betai)), g=betai_final$ct_chr, p.adjust.method="fdr")

#print("Adaptation rate early vs late autos:")
#wilcox.test(as.numeric(as.character(alpha_auto[which(alpha_auto$group=="LZ"),]$alpha)),as.numeric(as.character(alpha_auto[which(alpha_auto$group=="RS"),]$alpha)))#, alternative="less")
#print("Adaptation rate early vs late X:")
#wilcox.test(as.numeric(as.character(alpha_X[which(alpha_X$group=="LZ"),]$alpha)),as.numeric(as.character(alpha_X[which(alpha_X$group=="RS"),]$alpha)))#, alternative="less")
#print("Adaptation rate early X vs autos:")
#wilcox.test(as.numeric(as.character(alpha_X[which(alpha_X$group=="LZ"),]$alpha)),as.numeric(as.character(alpha_auto[which(alpha_auto$group=="LZ"),]$alpha)))
#print("Adaptation rate late X vs autos:")
#wilcox.test(as.numeric(as.character(alpha_X[which(alpha_X$group=="RS"),]$alpha)),as.numeric(as.character(alpha_auto[which(alpha_auto$group=="RS"),]$alpha)))

print("Done!")
