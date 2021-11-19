#PURPOSE: Test if number of interactions is correlated with EVE divergence or dN/dS

#suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))
suppressPackageStartupMessages(library(EnsDb.Mmusculus.v75))
suppressPackageStartupMessages(library(ggplot2))

print("Reading in data...")

#Get interaction counts
LZinteractions<-read.table("interaction_counts.LZinduced.txt", header=TRUE)
RSinteractions<-read.table("interaction_counts.RSinduced.txt", header=TRUE)

#Get EVE data and filter
mydataLZ<-read.table(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/results/RPKM1_INDUCED2/indivBetaMLparamsLZinduced2_RPKM1.res"))
mydataRS<-read.table(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/results/RPKM1_INDUCED2/indivBetaMLparamsRSinduced2_RPKM1.res"))
#mydata<-as.data.frame(rbind(mydataLZ,mydataRS))
LZgenes<-scan(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/INDUCED_CUTOFF2/gene_list_LZinduced_edgeR_wholeGenome.ensemblOrthos.rpkm1.txt"),what=character())
RSgenes<-scan(paste0("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/INDUCED_CUTOFF2/gene_list_RSinduced_edgeR_wholeGenome.ensemblOrthos.rpkm1.txt"),what=character())
betai_LZ<-cbind(betai=0-log(mydataLZ[,4]), gene=LZgenes)
betai_RS<-cbind(betai=0-log(mydataRS[,4]), gene=RSgenes)

EVEdata_LZ<-as.data.frame(betai_LZ)
EVEdata_LZ<-EVEdata_LZ[which(as.numeric(as.character(EVEdata_LZ$betai))>-5),]
print(head(EVEdata_LZ))
EVEdata_RS<-as.data.frame(betai_RS)
EVEdata_RS<-EVEdata_RS[which(as.numeric(as.character(EVEdata_RS$betai))>-5),]
print(head(EVEdata_RS))

#Get dN/dS values and filter
dnds0<-read.table("/mnt/beegfs/ek112884/PAML_runs/WHOLE_GENOMES/omega_list.txt", col.names=c("gene","dN","dS","dN.dS"))
dnds<-dnds0[which(as.numeric(as.character(dnds0$dN.dS)) < 1.5),]

#For converting Entrez to Ensembl ID
#edb<-EnsDb.Mmusculus.v79
edb<-EnsDb.Mmusculus.v75
myGenes<-genes(edb)
print(head(myGenes))
myMap<-mapIds(edb, keys=as.character(LZinteractions$Entrez), keytype="ENTREZID", column="GENEID", multiVals="first")
print(head(mapIds))

#Get final dataframes for ggplot
LZall<-c()
for(i in 1:nrow(LZinteractions)){
	myEntrez<-as.character(LZinteractions$Entrez[i])
	#myID<-myGenes$gene_id[which(unlist(myGenes$entrezid)==i)]
	myID<-myMap[myEntrez]
	#print(myID)
	if(myID %in% EVEdata_LZ$gene){
		myEVE<-EVEdata_LZ$betai[which(EVEdata_LZ$gene==myID)]
	} else{
		myEVE<-NA
	}
	if(myID %in% dnds$gene){
		myDnds<-dnds$dN.dS[which(dnds$gene==myID)]
	} else{
		myDnds<-NA
	}
	LZall<-rbind(LZall, c(myID, myEVE, myDnds, LZinteractions[i,]))
}
LZdf<-as.data.frame(LZall)
print(head(LZdf))
colnames(LZdf)<-c("ensemblID", "eve_betai", "dnds", "entrez", "interaction_counts", "midInteraction_counts", "highInteraction_counts")
print(head(LZdf))

RSall<-c()
for(i in 1:nrow(RSinteractions)){
        myEntrez<-as.character(RSinteractions$Entrez[i])
        #myID<-myGenes$gene_id[which(unlist(myGenes$entrezid)==i)]
        myID<-myMap[myEntrez]
        #print(myID)
        if(myID %in% EVEdata_RS$gene){
                myEVE<-EVEdata_RS$betai[which(EVEdata_RS$gene==myID)]
        } else{
                myEVE<-NA
        }
        if(myID %in% dnds$gene){
                myDnds<-dnds$dN.dS[which(dnds$gene==myID)]
        } else{
                myDnds<-NA
        }
        RSall<-rbind(RSall, c(myID, myEVE, myDnds, RSinteractions[i,]))
}
RSdf<-as.data.frame(RSall)
print(head(RSdf))
colnames(RSdf)<-c("ensemblID", "eve_betai", "dnds", "entrez", "interaction_counts", "midInteraction_counts", "highInteraction_counts")
print(head(RSdf))

rhos<-c()
p.values<-c()

#Get correlations and plot
pdf("interactions_vs_EVE.LZinduced.pdf", onefile=TRUE)
#LZ interactions vs EVE
result<-cor.test(as.numeric(as.character(LZdf$interaction_counts)), as.numeric(as.character(LZdf$eve_betai)), method="spearman")
print(result)
rhos<-c(rhos, result$estimate)
p.values<-c(p.values, result$p.value)
p<-ggplot(LZdf, aes(x=log(as.numeric(as.character(interaction_counts))), y=as.numeric(as.character(eve_betai)))) + geom_point()
p<-p+labs(title="EVE divergence vs number of interactions", x="log(Number of Interactions)", y="EVE divergence")
print(p)

#LZ mid-score interactions vs EVE
tempDF<-LZdf[which(LZdf$midInteraction_counts != 0), ]
result<-cor.test(as.numeric(as.character(tempDF$midInteraction_counts)), as.numeric(as.character(tempDF$eve_betai)), method="spearman")
print(result)
rhos<-c(rhos, result$estimate)
p.values<-c(p.values, result$p.value)
p<-ggplot(tempDF, aes(x=log(as.numeric(as.character(midInteraction_counts))), y=as.numeric(as.character(eve_betai)))) + geom_point()
p<-p+labs(title="EVE divergence vs number of midInteractions", x="log(Number of mid-score or higher Interactions)", y="EVE divergence")
print(p)

#LZ high-score interactions vs EVE
tempDF<-LZdf[which(LZdf$highInteraction_counts != 0), ]
result<-cor.test(as.numeric(as.character(tempDF$highInteraction_counts)), as.numeric(as.character(tempDF$eve_betai)), method="spearman")
print(result)
rhos<-c(rhos, result$estimate)
p.values<-c(p.values, result$p.value)
p<-ggplot(tempDF, aes(x=log(as.numeric(as.character(highInteraction_counts))), y=as.numeric(as.character(eve_betai)))) + geom_point()
p<-p+labs(title="EVE divergence vs number of highInteractions", x="log(Number of high-score Interactions)", y="EVE divergence")
print(p)
dev.off()

pdf("interactions_vs_dnds.LZinduced.pdf", onefile=TRUE)
#LZ interactions vs dnds
result<-cor.test(as.numeric(as.character(LZdf$interaction_counts)), as.numeric(as.character(LZdf$dnds)), method="spearman")
print(result)
rhos<-c(rhos, result$estimate)
p.values<-c(p.values, result$p.value)
p<-ggplot(LZdf, aes(x=log(as.numeric(as.character(interaction_counts))), y=as.numeric(as.character(dnds)))) + geom_point()
p<-p+labs(title="dN/dS vs number of interactions", x="log(Number of Interactions)", y="dN/dS")
print(p)

#LZ mid-score interactions vs dnds
tempDF<-LZdf[which(LZdf$midInteraction_counts != 0), ]
result<-cor.test(as.numeric(as.character(tempDF$midInteraction_counts)), as.numeric(as.character(tempDF$dnds)), method="spearman")
print(result)
rhos<-c(rhos, result$estimate)
p.values<-c(p.values, result$p.value)
p<-ggplot(tempDF, aes(x=log(as.numeric(as.character(midInteraction_counts))), y=as.numeric(as.character(dnds)))) + geom_point()
p<-p+labs(title="dN/ds divergence vs number of midInteractions", x="log(Number of mid-score or higher Interactions)", y="dN/dS divergence")
print(p)

#LZ high-score interactions vs dnds
tempDF<-LZdf[which(LZdf$highInteraction_counts != 0), ]
result<-cor.test(as.numeric(as.character(tempDF$highInteraction_counts)), as.numeric(as.character(tempDF$dnds)), method="spearman")
print(result)
rhos<-c(rhos, result$estimate)
p.values<-c(p.values, result$p.value)
p<-ggplot(tempDF, aes(x=log(as.numeric(as.character(highInteraction_counts))), y=as.numeric(as.character(dnds)))) + geom_point()
p<-p+labs(title="dN/dS divergence vs number of highInteractions", x="log(Number of high-score Interactions)", y="dN/dS divergence")
print(p)
dev.off()

pdf("interactions_vs_EVE.RSinduced.pdf", onefile=TRUE)
#RS interactions vs EVE
result<-cor.test(as.numeric(as.character(RSdf$interaction_counts)), as.numeric(as.character(RSdf$eve_betai)), method="spearman")
print(result)
rhos<-c(rhos, result$estimate)
p.values<-c(p.values, result$p.value)
p<-ggplot(RSdf, aes(x=log(as.numeric(as.character(interaction_counts))), y=as.numeric(as.character(eve_betai)))) + geom_point()
p<-p+labs(title="EVE divergence vs number of interactions", x="log(Number of Interactions)", y="EVE divergence")
print(p)

#RS mid-score interactions vs EVE
tempDF<-RSdf[which(RSdf$midInteraction_counts != 0), ]
result<-cor.test(as.numeric(as.character(tempDF$midInteraction_counts)), as.numeric(as.character(tempDF$eve_betai)), method="spearman")
print(result)
rhos<-c(rhos, result$estimate)
p.values<-c(p.values, result$p.value)
p<-ggplot(tempDF, aes(x=log(as.numeric(as.character(midInteraction_counts))), y=as.numeric(as.character(eve_betai)))) + geom_point()
p<-p+labs(title="EVE divergence vs number of midInteractions", x="log(Number of mid-score or higher Interactions)", y="EVE divergence")
print(p)

#RS high-score interactions vs EVE
tempDF<-RSdf[which(RSdf$highInteraction_counts != 0), ]
result<-cor.test(as.numeric(as.character(tempDF$highInteraction_counts)), as.numeric(as.character(tempDF$eve_betai)), method="spearman")
print(result)
rhos<-c(rhos, result$estimate)
p.values<-c(p.values, result$p.value)
p<-ggplot(tempDF, aes(x=log(as.numeric(as.character(highInteraction_counts))), y=as.numeric(as.character(eve_betai)))) + geom_point()
p<-p+labs(title="EVE divergence vs number of highInteractions", x="log(Number of high-score Interactions)", y="EVE divergence")
print(p)
dev.off()

pdf("interactions_vs_dnds.RSinduced.pdf", onefile=TRUE)
#RS interactions vs dnds
result<-cor.test(as.numeric(as.character(RSdf$interaction_counts)), as.numeric(as.character(RSdf$dnds)), method="spearman")
print(result)
rhos<-c(rhos, result$estimate)
p.values<-c(p.values, result$p.value)
p<-ggplot(RSdf, aes(x=log(as.numeric(as.character(interaction_counts))), y=as.numeric(as.character(dnds)))) + geom_point()
p<-p+labs(title="dN/dS vs number of interactions", x="log(Number of Interactions)", y="dN/dS")
print(p)

#RS mid-score interactions vs dnds
tempDF<-RSdf[which(RSdf$midInteraction_counts != 0), ]
result<-cor.test(as.numeric(as.character(tempDF$midInteraction_counts)), as.numeric(as.character(tempDF$dnds)), method="spearman")
print(result)
rhos<-c(rhos, result$estimate)
p.values<-c(p.values, result$p.value)
p<-ggplot(tempDF, aes(x=log(as.numeric(as.character(midInteraction_counts))), y=as.numeric(as.character(dnds)))) + geom_point()
p<-p+labs(title="dN/ds divergence vs number of midInteractions", x="log(Number of mid-score or higher Interactions)", y="dN/dS divergence")
print(p)

#RS high-score interactions vs dnds
tempDF<-RSdf[which(RSdf$highInteraction_counts != 0), ]
result<-cor.test(as.numeric(as.character(tempDF$highInteraction_counts)), as.numeric(as.character(tempDF$dnds)), method="spearman")
print(result)
rhos<-c(rhos, result$estimate)
p.values<-c(p.values, result$p.value)
p<-ggplot(tempDF, aes(x=log(as.numeric(as.character(highInteraction_counts))), y=as.numeric(as.character(dnds)))) + geom_point()
p<-p+labs(title="dN/dS divergence vs number of highInteractions", x="log(Number of high-score Interactions)", y="dN/dS divergence")
print(p)
dev.off()

myNames<-c("early_allVSeve","early_midVSeve","early_highVSeve","early_allVSdnds","early_midVSdnds","early_highVSdnds","late_allVSeve","late_midVSeve","late_highVSeve","late_allVSdnds","late_midVSdnds","late_highVSdnds")
names(rhos)<-myNames
print("Rho values:")
print(rhos)

fdr.ps<-p.adjust(p.values, method="fdr")
names(fdr.ps)<-myNames
print("FDR-corrected Ps:")
print(fdr.ps)

print("Done with correlations")
