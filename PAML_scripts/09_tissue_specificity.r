#PURPOSE: Test for a relationship between dN/dS and tsi (tissue specificity index)
#Analagous to test for relationship between exp div and tsi
#Also test for significant association between pos selection and high tsi? - I think we WOULDN'T expect this necessarily b/c rapid evo (high dN/dS) in high tsi genes probably due to relaxed functional constraint

library(mouse4302.db)
#TISSUE SPECIFICITY using the  method from Yanai et al. 2005; formula more clear in Liao and Zhang 2006
#t = sum(1-(log2(Sh(i,j))/log2(Sh(i,max))))/(#tissues-1)
#Where sum is over all tissues, i is particular gene, Sh(i,j) is expression of gene i in tissue j, Sh(i,max) is maximum expression of gene i in any tissue
#Get expression levels across tissues from mouse expression atlas (MOE430 2.0); available from BioGPS

#Mouse expression atlas has Affymetrix microarray IDs; change these to ensembl IDs
print("Reading in expression atlas and converting to Ensembl IDs...")
ensembl_data<-toTable(mouse4302ENSEMBL)
print(dim(ensembl_data))
ensembl_data_nodup<-ensembl_data[which(!(duplicated(ensembl_data$probe_id))),]
print(dim(ensembl_data_nodup))
head(ensembl_data_nodup)
ensembl_data<-ensembl_data_nodup

#BioGPS has two different atlasses - and "averaged" one and a not averaged one
#The averaged one has one entry per tissue while the other has two entries per tissue
#BUT the averaged file is not the mean of the two tissue values in the other file
#It is close though; maybe they added samples and updated one file but not other? Or maybe their average is nota simple mean?
#Regardless, I think it's best to use the "average" one
#It looks like that's what Erica did too, based on her Rscript
#atlas csv available from BioGPS: http://biogps.org/downloads/
atlas<-read.csv("../../mus_expression_analysis/geneatlas_MOE430_20090327.raw.avg.csv")
print(dim(atlas))
head(atlas$X)

ensembl_id<-c()
for(i in 1:nrow(atlas)){
        if(atlas$X[i] %in% ensembl_data$probe_id){
                id<-ensembl_data$ensembl_id[which(ensembl_data$probe_id==atlas$X[i])]
        } else{
                id<-NA
        }
        ensembl_id<-c(ensembl_id,id)
}
atlas<-cbind(ensembl_id,atlas)

#I think the first several columns are mouse cell lines, not actual tissues; getting rid of those
atlas_noCellLines<-atlas[,c(1,2,18:ncol(atlas))]
atlas<-atlas_noCellLines

max_exp<-rowMax(as.matrix(atlas[,3:ncol(atlas)]))

tsi<-rowSums(1-(log(atlas[,3:ncol(atlas)])/log(max_exp)))/(ncol(atlas)-3) #denominator should be #tissues-1; here the #tissues is the  #columns-2 because the first two columns are gene IDs

tsi_table<-as.data.frame(cbind(as.character(atlas$ensembl_id),tsi))
tsi_table_unique<-tsi_table[which(!(duplicated(tsi_table$V1))),]

#Test for correlation between dN/dS and tsi
#Spearman's Rho?
#Also plot to visualize?
print("Pairing dN/dS with tsi for each gene...")
pdf("dnds_vs_tissues_specificity.pdf",height=11,width=8.5,onefile=TRUE)

dnds<-read.table("omega_list.txt")
print(paste("dnds dim:",dim(dnds)))
dnds<-dnds[which(dnds$V4 != 999),]
print(paste("dnds dim after filtering >>1 dnds:",dim(dnds)))
ordered_tsi<-c()
for(i in dnds$V1){
        if(i %in% as.character(tsi_table_unique$V1)){
#               print(tsi_table_unique[which(as.character(tsi_table_unique$V1)==i),])
                ordered_tsi<-c(ordered_tsi,as.numeric(as.character(tsi_table_unique$tsi[which(as.character(tsi_table_unique$V1)==i)])))
        } else{
                ordered_tsi<-c(ordered_tsi,NA)
        }
}
dnds_tsi<-as.data.frame(cbind(dnds,tsi=ordered_tsi))
print("dnds-tsi table:")
print(head(dnds_tsi))

print("Calculating spearman's rho and plotting...")
sr<-cor.test(dnds_tsi$tsi,dnds_tsi$V4,method="spearman")
plot(x=dnds_tsi$tsi,y=dnds_tsi$V4,xlab="Tissue Specificity Index",ylab="dN/dS",main="Relationship between dN/dS and tissue specificity - genome-wide")
legend("topleft",legend=c(paste("Spearman's Rho:",formatC(sr$estimate,digits=3)),paste("p-value:",formatC(sr$p.value,format="e",digits=3))))

print("Repeating for LZ...")
LZ_dnds<-read.table("omega_list_LZ.txt")
print(paste("LZ_dnds dim:",dim(LZ_dnds)))
LZ_dnds<-LZ_dnds[which(LZ_dnds$V4 != 999),]
print(paste("LZ_dnds dim after filtering >>1 LZ_dnds:",dim(LZ_dnds)))
ordered_tsi<-c()
for(i in LZ_dnds$V1){
        if(i %in% as.character(tsi_table_unique$V1)){
#               print(tsi_table_unique[which(as.character(tsi_table_unique$V1)==i),])
                ordered_tsi<-c(ordered_tsi,as.numeric(as.character(tsi_table_unique$tsi[which(as.character(tsi_table_unique$V1)==i)])))
        } else{
                ordered_tsi<-c(ordered_tsi,NA)
        }
}
LZ_dnds_tsi<-as.data.frame(cbind(LZ_dnds,tsi=ordered_tsi))
print("LZ_dnds-tsi table:")
print(head(LZ_dnds_tsi))

print("Calculating spearman's rho and plotting...")
sr<-cor.test(LZ_dnds_tsi$tsi,LZ_dnds_tsi$V4,method="spearman")
plot(x=LZ_dnds_tsi$tsi,y=LZ_dnds_tsi$V4,xlab="Tissue Specificity Index",ylab="dN/dS",main="Relationship between dN/dS and tissue specificity - expressed early",col="chocolate1")
legend("topleft",legend=c(paste("Spearman's Rho:",formatC(sr$estimate,digits=3)),paste("p-value:",formatC(sr$p.value,format="e",digits=3))))

print("Repeating for RS...")
RS_dnds<-read.table("omega_list_RS.txt")
print(paste("RS_dnds dim:",dim(RS_dnds)))
RS_dnds<-RS_dnds[which(RS_dnds$V4 != 999),]
print(paste("RS_dnds dim after filtering >>1 RS_dnds:",dim(RS_dnds)))
ordered_tsi<-c()
for(i in RS_dnds$V1){
        if(i %in% as.character(tsi_table_unique$V1)){
#               print(tsi_table_unique[which(as.character(tsi_table_unique$V1)==i),])
                ordered_tsi<-c(ordered_tsi,as.numeric(as.character(tsi_table_unique$tsi[which(as.character(tsi_table_unique$V1)==i)])))
        } else{
                ordered_tsi<-c(ordered_tsi,NA)
        }
}
RS_dnds_tsi<-as.data.frame(cbind(RS_dnds,tsi=ordered_tsi))
print("RS_dnds-tsi table:")
print(head(RS_dnds_tsi))

print("Calculating spearman's rho and plotting...")
sr<-cor.test(RS_dnds_tsi$tsi,RS_dnds_tsi$V4,method="spearman")
plot(x=RS_dnds_tsi$tsi,y=RS_dnds_tsi$V4,xlab="Tissue Specificity Index",ylab="dN/dS",main="Relationship between dN/dS and tissue specificity - expressed late",col="lightsteelblue1")
legend("topleft",legend=c(paste("Spearman's Rho:",formatC(sr$estimate,digits=3)),paste("p-value:",formatC(sr$p.value,format="e",digits=3))))

print("Repeating for LZinduced...")
LZinduced_dnds<-read.table("omega_list_LZinduced.txt")
print(paste("LZinduced_dnds dim:",dim(LZinduced_dnds)))
LZinduced_dnds<-LZinduced_dnds[which(LZinduced_dnds$V4 != 999),]
print(paste("LZinduced_dnds dim after filtering >>1 LZinduced_dnds:",dim(LZinduced_dnds)))
ordered_tsi<-c()
for(i in LZinduced_dnds$V1){
        if(i %in% as.character(tsi_table_unique$V1)){
#               print(tsi_table_unique[which(as.character(tsi_table_unique$V1)==i),])
                ordered_tsi<-c(ordered_tsi,as.numeric(as.character(tsi_table_unique$tsi[which(as.character(tsi_table_unique$V1)==i)])))
        } else{
                ordered_tsi<-c(ordered_tsi,NA)
        }
}
LZinduced_dnds_tsi<-as.data.frame(cbind(LZinduced_dnds,tsi=ordered_tsi))
print("LZinduced_dnds-tsi table:")
print(head(LZinduced_dnds_tsi))

print("Calculating spearman's rho and plotting...")
sr<-cor.test(LZinduced_dnds_tsi$tsi,LZinduced_dnds_tsi$V4,method="spearman")
plot(x=LZinduced_dnds_tsi$tsi,y=LZinduced_dnds_tsi$V4,xlab="Tissue Specificity Index",ylab="dN/dS",main="Relationship between dN/dS and tissue specificity - induced early",col="chocolate1")
legend("topleft",legend=c(paste("Spearman's Rho:",formatC(sr$estimate,digits=3)),paste("p-value:",formatC(sr$p.value,format="e",digits=3))))

print("Repeating for RSinduced...")
RSinduced_dnds<-read.table("omega_list_RSinduced.txt")
print(paste("RSinduced_dnds dim:",dim(RSinduced_dnds)))
RSinduced_dnds<-RSinduced_dnds[which(RSinduced_dnds$V4 != 999),]
print(paste("RSinduced_dnds dim after filtering >>1 RSinduced_dnds:",dim(RSinduced_dnds)))
ordered_tsi<-c()
for(i in RSinduced_dnds$V1){
        if(i %in% as.character(tsi_table_unique$V1)){
#               print(tsi_table_unique[which(as.character(tsi_table_unique$V1)==i),])
                ordered_tsi<-c(ordered_tsi,as.numeric(as.character(tsi_table_unique$tsi[which(as.character(tsi_table_unique$V1)==i)])))
        } else{
                ordered_tsi<-c(ordered_tsi,NA)
        }
}
RSinduced_dnds_tsi<-as.data.frame(cbind(RSinduced_dnds,tsi=ordered_tsi))
print("RSinduced_dnds-tsi table:")
print(head(RSinduced_dnds_tsi))

print("Calculating spearman's rho and plotting...")
sr<-cor.test(RSinduced_dnds_tsi$tsi,RSinduced_dnds_tsi$V4,method="spearman")
plot(x=RSinduced_dnds_tsi$tsi,y=RSinduced_dnds_tsi$V4,xlab="Tissue Specificity Index",ylab="dN/dS",main="Relationship between dN/dS and tissue specificity - induced late",col="lightsteelblue1")
legend("topleft",legend=c(paste("Spearman's Rho:",formatC(sr$estimate,digits=3)),paste("p-value:",formatC(sr$p.value,format="e",digits=3))))

dev.off()

print("Done!")
