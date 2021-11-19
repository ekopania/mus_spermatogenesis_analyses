#PURPOSE: Use EdgeR to calculate RPKM and filter by RPKM, then plot density using Zigzag (and also maybe EdgeR, just to make sure they look the same?)

#change these to test different parameters
min_rpkm<-1
min_samples<-8 #8

library(coda)
library(zigzag)

print("Loading EdgeR and setting up DGE object...")
library(edgeR)
#ALL FOUR SPECIES
#Mus ref geneid from ensembl pairwise 1:1 orthos appended:
files<-list.files(path="/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/",pattern=".orthoIDappended.txt$")
myfiles<-c()
for(f in files) (myfiles<-c(myfiles,paste0("../",f)))
head(myfiles)
mygroups<-c(1,2,1,2,3,4,3,4,1,2,1,2,1,2,1,2,1,2,1,2,3,4,3,4,3,4,7,8,7,8,7,8,3,4,3,4,3,4,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,1,2,1,2,1,2,1,2)
myDGE<-readDGE(myfiles,columns=c(2,8),group=mygroups)
dim(myDGE)

print("Filtering by expression level...")
#Get gene lengths for each species and for mouse reference
#Getting from geature counts output b/c these are the correct lengths to use for fpkm calc (no introns, etc)
#fc_dom<-read.table("BIK4665-1M-LZ_counts.refIDappended.compara.txt",header=TRUE)
fc_dom<-read.table("../BIK4665-1M-LZ_counts.orthoIDappended.txt",header=TRUE)
fc_domLengths<-cbind(as.numeric(fc_dom$Length),as.character(fc_dom$Geneid))
head(fc_domLengths)
fc_mus<-read.table("../PPPP98-3M-LZ_counts.orthoIDappended.txt",header=TRUE)
fc_musLengths<-cbind(as.numeric(fc_mus$Length),as.character(fc_mus$Geneid))
head(fc_musLengths)
fc_spr<-read.table("../SEG4156-LZ_counts.orthoIDappended.txt",header=TRUE)
fc_sprLengths<-cbind(as.numeric(fc_spr$Length),as.character(fc_spr$Geneid))
head(fc_sprLengths)
fc_pah<-read.table("../PAHNew-1M-LZ_counts.orthoIDappended.txt",header=TRUE)
fc_pahLengths<-cbind(as.numeric(fc_pah$Length),as.character(fc_pah$Geneid))
head(fc_pahLengths)
fc_ref<-read.table("../../PSEUDOREF_MAPPINGS/labeled_BIK4665-1M-LZ.counts.txt",header=TRUE,skip=1)
fc_refLengths<-cbind(as.numeric(fc_ref$Length),as.character(fc_ref$Geneid))
head(fc_refLengths)
#overlap between genes included in species reference mapping and mouse reference mapping
overlap_fc_domLengths<-fc_domLengths[which(fc_domLengths[,2] %in% fc_refLengths[,2]),]
overlap_fc_musLengths<-fc_musLengths[which(fc_musLengths[,2] %in% fc_refLengths[,2]),]
overlap_fc_sprLengths<-fc_sprLengths[which(fc_sprLengths[,2] %in% fc_refLengths[,2]),]
overlap_fc_pahLengths<-fc_pahLengths[which(fc_pahLengths[,2] %in% fc_refLengths[,2]),]
overlap_fc_refLengths<-fc_refLengths[which(fc_refLengths[,2] %in% fc_domLengths[,2]),]
#sort by gene name
ordered_fc_domLengths<-overlap_fc_domLengths[match(rownames(myDGE$counts),overlap_fc_domLengths[,2]),]
ordered_fc_musLengths<-overlap_fc_musLengths[match(rownames(myDGE$counts),overlap_fc_musLengths[,2]),]
ordered_fc_sprLengths<-overlap_fc_sprLengths[match(rownames(myDGE$counts),overlap_fc_sprLengths[,2]),]
ordered_fc_pahLengths<-overlap_fc_pahLengths[match(rownames(myDGE$counts),overlap_fc_pahLengths[,2]),]
ordered_fc_refLengths<-overlap_fc_refLengths[match(rownames(myDGE$counts),overlap_fc_refLengths[,2]),]
#Figure out FPKM for each species separately, using proper gene lengths for each
myDGE<-calcNormFactors(myDGE)
rpkm_dom<-rpkm(myDGE[,which(myDGE$samples$group==c(1,2))],gene.length=as.numeric(ordered_fc_domLengths[,1]))
rpkm_mus<-rpkm(myDGE[,which(myDGE$samples$group==c(3,4))],gene.length=as.numeric(ordered_fc_musLengths[,1]))
rpkm_spr<-rpkm(myDGE[,which(myDGE$samples$group==c(5,6))],gene.length=as.numeric(ordered_fc_sprLengths[,1]))
rpkm_pah<-rpkm(myDGE[,which(myDGE$samples$group==c(7,8))],gene.length=as.numeric(ordered_fc_pahLengths[,1]))
rpkm_all<-cbind(rpkm_dom,rpkm_mus,rpkm_spr,rpkm_pah)
#rpkm_all<-cbind(rpkm_dom,rpkm_mus)
new_rpkm_all<-rpkm_all[which(rownames(rpkm_all) %in% rownames(myDGE$counts)),]
head(rowSums(rpkm_all) > min_rpkm)
keep<-rowSums(rpkm_all > min_rpkm) >= min_samples
keep[is.na(keep)]<-FALSE #NAs show up for gene iDs that aren't in myDGE; replace then with false otherwise myDGE gets confused

print(head(keep))
myDGE<-myDGE[keep, , keep.lib.sizes=FALSE]

rpkm_LZ<-new_rpkm_all[keep,which(grepl("LZ",colnames(new_rpkm_all)))]
rpkm_RS<-new_rpkm_all[keep,which(grepl("RS",colnames(new_rpkm_all)))]

#head(round(rpkm_LZ, digits = 3))
#head(round(rpkm_RS, digits = 3))
dim(rpkm_LZ)
dim(rpkm_RS)

pdf("LZ_rpkm.log_exp_distribution.pdf")
plot(density(log(rpkm_LZ)), main ="", xlab = "log Expression", ylim = c(0, 0.5))
for(i in seq(ncol(rpkm_LZ))) lines(density(log(rpkm_LZ[,i])))
dev.off()

pdf("RS_rpkm.log_exp_distribution.pdf")
plot(density(log(rpkm_RS)), main ="", xlab = "log Expression", ylim = c(0, 0.5))
for(i in seq(ncol(rpkm_RS))) lines(density(log(rpkm_RS[,i])))
dev.off()

rpkm_combo<-cbind(rpkm_LZ,rpkm_RS)
pdf("all_rpkm.log_exp_distribution.pdf")
plot(density(log(rpkm_combo)), main ="", xlab = "log Expression", xlim = c(-10, 10), ylim = c(0, 0.5))
legend(-8, y=0.4, legend=c("LZ","RS"), fill=c("chocolate1","lightsteelblue1"))
for(i in seq(ncol(rpkm_LZ))) lines(density(log(rpkm_LZ[,i])),col="chocolate1")
par(new=TRUE)
for(i in seq(ncol(rpkm_RS))) lines(density(log(rpkm_RS[,i])),col="lightsteelblue1")
dev.off()


print("Done!")
