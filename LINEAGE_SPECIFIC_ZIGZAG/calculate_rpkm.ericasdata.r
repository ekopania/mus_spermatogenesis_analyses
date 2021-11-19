#PURPOSE: Use EdgeR to calculate RPKM and filter by RPKM, then plot density using Zigzag (and also maybe EdgeR, just to make sure they look the same?)

#change these to test different parameters
min_rpkm<-1
min_samples<-4 #48 - all samples  #4 - what Erica used

library(coda)
library(zigzag)

print("Loading EdgeR and setting up DGE object...")
library(edgeR)
#READ IN ERICA'S F1 DATA
mydata<-read.table("/mnt/beegfs/ek112884/RNAseqData_4cellType_fromErica/paired/CZIIhybrid_paired_suspenders_default_11Mar16_M_C_featureCounts.csv",header=TRUE)
dim(mydata)
#Get groups
#CCPP21.1M.DIP.bam
mygroups<-c()
for(sample in colnames(mydata)[7:ncol(mydata)]){
	ct<-unlist(strsplit(sample,"\\."))[3]
	#print(ct)
	if(ct=="SP"){
		mygroups<-c(mygroups,1)
	} else if(ct=="LZ"){
		mygroups<-c(mygroups,2)
	} else if(ct=="DIP"){
		mygroups<-c(mygroups,3)
	} else if(ct=="RS"){
		mygroups<-c(mygroups,4)
	} else{
		stop(paste("ERROR: invalid cell type",ct))
	}
}
myDGE<-DGEList(counts=mydata[,7:ncol(mydata)],group=mygroups)
dim(myDGE)
print(head(myDGE$samples))

print("Filtering by expression level...")
#Get gene lengths 
lengths<-cbind(as.numeric(mydata$Length),as.character(mydata$Geneid))
head(lengths)
#sort by gene name
ordered_lengths<-lengths[match(rownames(myDGE$counts),lengths[,2]),]
#Figure out FPKM for each species separately, using proper gene lengths for each
myDGE<-calcNormFactors(myDGE)
myrpkm<-rpkm(myDGE,gene.length=as.numeric(lengths[,1]))
new_rpkm_all<-myrpkm[which(rownames(myrpkm) %in% rownames(myDGE$counts)),]
keep<-rowSums(new_rpkm_all > min_rpkm) >= min_samples
keep[is.na(keep)]<-FALSE #NAs show up for gene iDs that aren't in myDGE; replace then with false otherwise myDGE gets confused

print(head(keep))

rpkm_SP<-new_rpkm_all[keep,which(grepl("SP",colnames(new_rpkm_all)))]
rpkm_LZ<-new_rpkm_all[keep,which(grepl("LZ",colnames(new_rpkm_all)))]
rpkm_DIP<-new_rpkm_all[keep,which(grepl("DIP",colnames(new_rpkm_all)))]
rpkm_RS<-new_rpkm_all[keep,which(grepl("RS",colnames(new_rpkm_all)))]

#head(round(rpkm_LZ, digits = 3))
#head(round(rpkm_RS, digits = 3))
dim(rpkm_SP)
dim(rpkm_LZ)
dim(rpkm_DIP)
dim(rpkm_RS)

pdf("SP_rpkm.ericasdata.log_exp_distribution.pdf")
plot(density(log(rpkm_SP)), main ="", xlab = "log Expression", ylim = c(0, 0.5))
for(i in seq(ncol(rpkm_SP))) lines(density(log(rpkm_SP[,i])))
dev.off()

pdf("LZ_rpkm.ericasdata.log_exp_distribution.pdf")
plot(density(log(rpkm_LZ)), main ="", xlab = "log Expression", ylim = c(0, 0.5))
for(i in seq(ncol(rpkm_LZ))) lines(density(log(rpkm_LZ[,i])))
dev.off()

pdf("DIP_rpkm.ericasdata.log_exp_distribution.pdf")
plot(density(log(rpkm_DIP)), main ="", xlab = "log Expression", ylim = c(0, 0.5))
for(i in seq(ncol(rpkm_DIP))) lines(density(log(rpkm_DIP[,i])))
dev.off()

pdf("RS_rpkm.ericasdata.log_exp_distribution.pdf")
plot(density(log(rpkm_RS)), main ="", xlab = "log Expression", ylim = c(0, 0.5))
for(i in seq(ncol(rpkm_RS))) lines(density(log(rpkm_RS[,i])))
dev.off()

#Plot all cell types on one graph
#rpkm_SP_labeled<-as.data.frame(cbind(rpkm_SP,ct=rep("SP",nrow(rpkm_SP))))
#rpkm_LZ_labeled<-as.data.frame(cbind(rpkm_LZ,ct=rep("LZ",nrow(rpkm_LZ))))
#rpkm_DIP_labeled<-as.data.frame(cbind(rpkm_DIP,ct=rep("DIP",nrow(rpkm_DIP))))
#rpkm_RS_labeled<-as.data.frame(cbind(rpkm_RS,ct=rep("RS",nrow(rpkm_RS))))
#rpkm_combo<-as.data.frame(cbind(rpkm_SP_labeled,rpkm_LZ_labeled,rpkm_DIP_labeled,rpkm_RS_labeled))
rpkm_combo<-cbind(rpkm_SP,rpkm_LZ,rpkm_DIP,rpkm_RS)

print(rpkm_combo[1:10,1:10])
pdf("all_rpkm.ericasdata.log_exp_distribution.pdf")
plot(density(log(rpkm_combo)), main ="", xlab = "log Expression", xlim = c(-10, 10), ylim = c(0, 0.5))
legend(-8, y=0.4, legend=c("SP","LZ","DIP","RS"), fill=c("red","chocolate1","grey","lightsteelblue1"))
for(i in seq(ncol(rpkm_SP))) lines(density(log(rpkm_SP[,i])),col="red")
par(new=TRUE)
for(i in seq(ncol(rpkm_LZ))) lines(density(log(rpkm_LZ[,i])),col="chocolate1")
par(new=TRUE)
for(i in seq(ncol(rpkm_DIP))) lines(density(log(rpkm_DIP[,i])),col="grey")
par(new=TRUE)
for(i in seq(ncol(rpkm_RS))) lines(density(log(rpkm_RS[,i])),col="lightsteelblue1")
dev.off()


print("Done!")
