#PURPOSE: Test if X chromosome or autosomes tend to have more protein-protein interactions

suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))
suppressPackageStartupMessages(library(EnsDb.Mmusculus.v75))

print("Reading in data...")

#Get interaction counts
LZinteractions<-read.table("interaction_counts.LZinduced.txt", header=TRUE)
RSinteractions<-read.table("interaction_counts.RSinduced.txt", header=TRUE)

#edb<-EnsDb.Mmusculus.v79
edb<-EnsDb.Mmusculus.v75
myGenes<-genes(edb)
print(head(myGenes))
myMapLZ<-mapIds(edb, keys=as.character(LZinteractions$Entrez), keytype="ENTREZID", column="GENEID", multiVals="first")
myMapRS<-mapIds(edb, keys=as.character(RSinteractions$Entrez), keytype="ENTREZID", column="GENEID", multiVals="first")
edb_y<-addFilter(edb, SeqNameFilter("Y"))
y_genes<-genes(edb_y)
edb_x<-addFilter(edb, SeqNameFilter("X"))
x_genes<-genes(edb_x)
edb_auto<-addFilter(edb, SeqNameFilter(c(1:19)))
auto_genes<-genes(edb_auto)

#Append Ensembl IDs and chromosome types
LZ_ensIDs<-c()
LZ_chrTypes<-c()
for(i in 1:nrow(LZinteractions)){
        myEntrez<-as.character(LZinteractions$Entrez[i])
        myID<-myMapLZ[myEntrez]
	LZ_ensIDs<-c(LZ_ensIDs, myID)
	if(myID %in% y_genes$gene_id){
                LZ_chrTypes<-c(LZ_chrTypes,"Y")
        } else if(myID %in% x_genes$gene_id){
                LZ_chrTypes<-c(LZ_chrTypes,"X")
        } else if(myID %in% auto_genes$gene_id){
                LZ_chrTypes<-c(LZ_chrTypes,"auto")
        } else{
                LZ_chrTypes<-c(LZ_chrTypes,"NA")
        }
}

RS_ensIDs<-c()
RS_chrTypes<-c()
for(i in 1:nrow(RSinteractions)){
        myEntrez<-as.character(RSinteractions$Entrez[i])
        myID<-myMapRS[myEntrez]
        RS_ensIDs<-c(RS_ensIDs, myID)
        if(myID %in% y_genes$gene_id){
                RS_chrTypes<-c(RS_chrTypes,"Y")
        } else if(myID %in% x_genes$gene_id){
                RS_chrTypes<-c(RS_chrTypes,"X")
        } else if(myID %in% auto_genes$gene_id){
                RS_chrTypes<-c(RS_chrTypes,"auto")
        } else{
                RS_chrTypes<-c(RS_chrTypes,"NA")
        }
}

LZdf0<-as.data.frame(cbind(ct=rep("LZ", nrow(LZinteractions)), ens_id=LZ_ensIDs, chrType=LZ_chrTypes, LZinteractions))
LZdf<-LZdf0[which(LZdf0$chrType %in% c("X", "auto")),]
RSdf0<-as.data.frame(cbind(ct=rep("RS", nrow(RSinteractions)), ens_id=RS_ensIDs, chrType=RS_chrTypes, RSinteractions))
RSdf<-RSdf0[which(RSdf0$chrType %in% c("X", "auto")),]
ALLdf<-as.data.frame(rbind(LZdf, RSdf))

#Print medians
print("Median number of interactions:")
print(paste("Early, auto, all:", median(LZdf$num_interactions[which(LZdf$chrType=="auto")], na.rm=TRUE)))
print(paste("Early, auto, mid:", median(LZdf$num_midScore_interactions[which(LZdf$chrType=="auto")], na.rm=TRUE)))
print(paste("Early, auto, high:", median(LZdf$num_highScore_interactions[which(LZdf$chrType=="auto")], na.rm=TRUE)))
print(paste("Early, X, all:", median(LZdf$num_interactions[which(LZdf$chrType=="X")], na.rm=TRUE)))
print(paste("Early, X, mid:", median(LZdf$num_midScore_interactions[which(LZdf$chrType=="X")], na.rm=TRUE)))
print(paste("Early, X, high:", median(LZdf$num_highScore_interactions[which(LZdf$chrType=="X")], na.rm=TRUE)))
print(paste("Late, auto, all:", median(RSdf$num_interactions[which(RSdf$chrType=="auto")], na.rm=TRUE)))
print(paste("Late, auto, mid:", median(RSdf$num_midScore_interactions[which(RSdf$chrType=="auto")], na.rm=TRUE)))
print(paste("Late, auto, high:", median(RSdf$num_highScore_interactions[which(RSdf$chrType=="auto")], na.rm=TRUE)))
print(paste("Late, X, all:", median(RSdf$num_interactions[which(RSdf$chrType=="X")], na.rm=TRUE)))
print(paste("Late, X, mid:", median(RSdf$num_midScore_interactions[which(RSdf$chrType=="X")], na.rm=TRUE)))
print(paste("Late, X, high:", median(RSdf$num_highScore_interactions[which(RSdf$chrType=="X")], na.rm=TRUE)))

#Wilcoxon rank sum tests to compare early vs late
print("Wilcoxon rank sum tests for X vs auto, early (all, mid score, high score):")
wilcox_all<-wilcox.test(LZdf$num_interactions[which(LZdf$chrType=="auto")], LZdf$num_interactions[which(LZdf$chrType=="X")])
wilcox_mid<-wilcox.test(LZdf$num_midScore_interactions[which(LZdf$chrType=="auto")], LZdf$num_midScore_interactions[which(LZdf$chrType=="X")])
wilcox_high<-wilcox.test(LZdf$num_highScore_interactions[which(LZdf$chrType=="auto")], LZdf$num_highScore_interactions[which(LZdf$chrType=="X")])
print(wilcox_all)
print(wilcox_mid)
print(wilcox_high)
p.values<-c(wilcox_all$p.value, wilcox_mid$p.value, wilcox_high$p.value)
print("FDR-corrected p-values")
p.adjust(p.values, method="fdr")
print("Bonferroni-corrected p-values")
p.adjust(p.values, method="bonferroni")
print("Wilcoxon rank sum tests for X vs auto, late (all, mid score, high score):")
wilcox_all<-wilcox.test(RSdf$num_interactions[which(RSdf$chrType=="auto")], RSdf$num_interactions[which(RSdf$chrType=="X")])
wilcox_mid<-wilcox.test(RSdf$num_midScore_interactions[which(RSdf$chrType=="auto")], RSdf$num_midScore_interactions[which(RSdf$chrType=="X")])
wilcox_high<-wilcox.test(RSdf$num_highScore_interactions[which(RSdf$chrType=="auto")], RSdf$num_highScore_interactions[which(RSdf$chrType=="X")])
print(wilcox_all)
print(wilcox_mid)
print(wilcox_high)
p.values<-c(wilcox_all$p.value, wilcox_mid$p.value, wilcox_high$p.value)
print("FDR-corrected p-values")
p.adjust(p.values, method="fdr")
print("Bonferroni-corrected p-values")
p.adjust(p.values, method="bonferroni")

#Plot
print("Plotting...")
pdf("number_interactions.XVSautos.boxplots.pdf", onefile=TRUE)
#High-score interactions - early
p<-ggplot(LZdf, aes(x=chrType, y=log(as.numeric(as.character(num_highScore_interactions))), fill=ct)) + geom_boxplot()
p<-p + labs(title="Number of high-score interactions; X vs auto; early", x="Chr type", y="log(Number of high-score interactions)")
p<-p + theme(axis.text.y = element_text(size=20))
p<-p + theme_minimal()
p<-p + scale_fill_manual(values=c("chocolate1"))
print(p)
#High-score interactions - late
p<-ggplot(RSdf, aes(x=chrType, y=log(as.numeric(as.character(num_highScore_interactions))), fill=ct)) + geom_boxplot()
p<-p + labs(title="Number of high-score interactions; X vs auto; late", x="Chr type", y="log(Number of high-score interactions)")
p<-p + theme(axis.text.y = element_text(size=20))
p<-p + theme_minimal()
p<-p + scale_fill_manual(values=c("lightsteelblue1"))
print(p)

dev.off()

print("Done with X vs auto")
