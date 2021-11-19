#PURPOSEL generate violin plots to compare dN/dS values for genes induced early. those induced late, and all genes

#library(EnsDb.Mmusculus.v79)
suppressPackageStartupMessages(library(EnsDb.Mmusculus.v75))
library(ggplot2)

mapping<-"MULTI_MAP"
rpkm<-10
induced_cutoff<-2
testis_spec<-FALSE

print("Reading in data...")
#LZinduced<-read.table(paste0("omega_list_LZinduced",induced_cutoff,".ensemblOrthos.",mapping,".txt"))
#RSinduced<-read.table(paste0("omega_list_RSinduced",induced_cutoff,".ensemblOrthos.",mapping,".txt"))
LZinduced<-read.table(paste0("omega_list_LZ.ensemblOrthos.",mapping,".rpkm",rpkm,".txt")) #NOT induced but keeping variable names for simplicity
RSinduced<-read.table(paste0("omega_list_RS.ensemblOrthos.",mapping,".rpkm",rpkm,".txt")) #NOT induced but keeping variable names for simplicity
all<-read.table("omega_list.ensemblOrthos.txt")

LZinduced<-LZinduced[which(as.numeric(as.character(LZinduced$V4)) < 1.5),]
RSinduced<-RSinduced[which(as.numeric(as.character(RSinduced$V4)) < 1.5),]
all<-all[which(as.numeric(as.character(all$V4)) < 1.5),]

if(testis_spec){
	testis_genes<-scan(gzfile("../../mus_expression_analysis/chalmel_testis_specific.txt.gz"),what=character())
	testis_genes<-unique(testis_genes)
	LZinduced<-LZinduced[which(LZinduced$V1 %in% testis_genes),]
	RSinduced<-RSinduced[which(RSinduced$V1 %in% testis_genes),]
	all<-all[which(all$V1 %in% testis_genes),]
}

print("LZ dimensions:")
print(dim(LZinduced))
print("RS dimensions:")
print(dim(RSinduced))
print("all dimensions:")
print(dim(all))

#append chr type
#edb<-EnsDb.Mmusculus.v79
edb<-EnsDb.Mmusculus.v75
edb_y<-addFilter(edb, SeqNameFilter("Y"))
y_genes<-genes(edb_y)
edb_x<-addFilter(edb, SeqNameFilter("X"))
x_genes<-genes(edb_x)
edb_auto<-addFilter(edb, SeqNameFilter(c(1:19)))
auto_genes<-genes(edb_auto)

LZ_chr_type<-c()
for(g in LZinduced$V1){
        if(g %in% y_genes$gene_id){
                LZ_chr_type<-c(LZ_chr_type,"Y")
        } else if(g %in% x_genes$gene_id){
                LZ_chr_type<-c(LZ_chr_type,"X")
        } else if(g %in% auto_genes$gene_id){
                LZ_chr_type<-c(LZ_chr_type,"auto")
        } else{
                LZ_chr_type<-c(LZ_chr_type,"NA")
        }
}
RS_chr_type<-c()
for(g in RSinduced$V1){
        if(g %in% y_genes$gene_id){
                RS_chr_type<-c(RS_chr_type,"Y")
        } else if(g %in% x_genes$gene_id){
                RS_chr_type<-c(RS_chr_type,"X")
        } else if(g %in% auto_genes$gene_id){
                RS_chr_type<-c(RS_chr_type,"auto")
        } else{
                RS_chr_type<-c(RS_chr_type,"NA")
        }
}
all_chr_type<-c()
for(g in all$V1){
        if(g %in% y_genes$gene_id){
                all_chr_type<-c(all_chr_type,"Y")
        } else if(g %in% x_genes$gene_id){
                all_chr_type<-c(all_chr_type,"X")
        } else if(g %in% auto_genes$gene_id){
                all_chr_type<-c(all_chr_type,"auto")
        } else{
                all_chr_type<-c(all_chr_type,"NA")
        }
}

LZinduced_final<-as.data.frame(cbind(LZinduced,chr_type=LZ_chr_type,ct=rep("LZ",nrow(LZinduced)),group=paste("LZ",LZ_chr_type,sep=".")))
RSinduced_final<-as.data.frame(cbind(RSinduced,chr_type=RS_chr_type,ct=rep("RS",nrow(RSinduced)),group=paste("RS",RS_chr_type,sep=".")))
all_final<-as.data.frame(cbind(all,chr_type=all_chr_type,ct=rep("all",nrow(all)),group=paste("all",all_chr_type,sep=".")))

#combo<-as.data.frame(rbind(LZinduced_final,RSinduced_final,all_final))
#combo<-combo[c(which(combo$group=="LZ.auto"),which(combo$group=="LZ.X"),which(combo$group=="RS.auto"),which(combo$group=="RS.X"),which(combo$group=="all.auto"),which(combo$group=="all.X")),]
combo<-as.data.frame(rbind(LZinduced_final,RSinduced_final))
combo<-combo[c(which(combo$group=="LZ.auto"),which(combo$group=="LZ.X"),which(combo$group=="RS.auto"),which(combo$group=="RS.X")),]

#print("Plotting...")
#p<-ggplot(combo, aes(x=group, y=as.numeric(as.character(V4)))) + geom_violin(aes(fill=ct)) + geom_boxplot(width=0.1) + labs(title="dN/dS in early and late spermatogenesis", x="Cell type", y="dN/dS")
#p<-p + theme(axis.text.y = element_text(size=20))
#p<-p + theme_minimal()
#p<-p + scale_fill_manual(values=c("chocolate1","lightsteelblue1")) #,"darkgray"))
#ggsave("dnds_violins_wholeGenome_multiMap.ensemblOrthos.cellTypesOnly.pdf",p,width=11,height=8.5,units="in")

print("# genes and Median dN/dS autos early:")
print(dim(combo[which(combo$group=="LZ.auto"),]))
print(median(as.numeric(as.character(combo$V4[which(combo$group=="LZ.auto")]))))
print("# genes and Median dN/dS X early:")
print(dim(combo[which(combo$group=="LZ.X"),]))
print(median(as.numeric(as.character(combo$V4[which(combo$group=="LZ.X")]))))
print("# genes and Median dN/dS autos late:")
print(dim(combo[which(combo$group=="RS.auto"),]))
print(median(as.numeric(as.character(combo$V4[which(combo$group=="RS.auto")]))))
print("# genes and Median dN/dS X late:")
print(dim(combo[which(combo$group=="RS.X"),]))
print(median(as.numeric(as.character(combo$V4[which(combo$group=="RS.X")]))))

print("Here are the Wilcoxon test results:")
print("Early vs late autosomes:")
wilcox.test(as.numeric(as.character(combo[which(combo$group=="LZ.auto"),]$V4)),as.numeric(as.character(combo[which(combo$group=="RS.auto"),]$V4))) #, alternative="less")
#print("Early vs all autos:")
#wilcox.test(as.numeric(as.character(combo[which(combo$group=="LZ.auto"),]$V4)),as.numeric(as.character(combo[which(combo$group=="all.auto"),]$V4))) #, alternative="less")
#print("Late vs all autos:")
#wilcox.test(as.numeric(as.character(combo[which(combo$group=="RS.auto"),]$V4)),as.numeric(as.character(combo[which(combo$group=="all.auto"),]$V4))) #, alternative="greater")
print("Early vs late X:")
wilcox.test(as.numeric(as.character(combo[which(combo$group=="LZ.X"),]$V4)),as.numeric(as.character(combo[which(combo$group=="RS.X"),]$V4)))
#print("Early vs all X:")
#wilcox.test(as.numeric(as.character(combo[which(combo$group=="LZ.X"),]$V4)),as.numeric(as.character(combo[which(combo$group=="all.X"),]$V4)))
#print("Late vs all X:")
#wilcox.test(as.numeric(as.character(combo[which(combo$group=="RS.X"),]$V4)),as.numeric(as.character(combo[which(combo$group=="all.X"),]$V4)))
print("X vs autos early:")
wilcox.test(as.numeric(as.character(combo[which(combo$group=="LZ.X"),]$V4)),as.numeric(as.character(combo[which(combo$group=="LZ.auto"),]$V4)))
print("X vs autos late:")
wilcox.test(as.numeric(as.character(combo[which(combo$group=="RS.X"),]$V4)),as.numeric(as.character(combo[which(combo$group=="RS.auto"),]$V4)))
#print("X vs autos all:")
#wilcox.test(as.numeric(as.character(combo[which(combo$group=="all.X"),]$V4)),as.numeric(as.character(combo[which(combo$group=="all.auto"),]$V4)))
print("X early vs auto late:")
wilcox.test(as.numeric(as.character(combo[which(combo$group=="LZ.X"),]$V4)),as.numeric(as.character(combo[which(combo$group=="RS.auto"),]$V4)))

print("Pairwise Wilcoxon Test:")
pairwise.wilcox.test(x=as.numeric(as.character(combo$V4)), g=combo$group, method="fdr")

print("Done!")
