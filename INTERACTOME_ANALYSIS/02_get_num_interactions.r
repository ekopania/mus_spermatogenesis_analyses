#PURPOSE: parse mippie output and count number of interactions for each gene within the query set (induced early or induced late); do for all interactions, mid-scoring, and high-scoring

print("Parsing LZ data...")
LZdata<-read.table("mippie.LZinduced.txt", header=TRUE)
print(dim(LZdata))
print(head(LZdata))

print("Get mid- and high-scoring interactions...")
LZdata_midScore<-LZdata[which(LZdata$MIPPIE_score > 0.53),]
LZdata_highScore<-LZdata[which(LZdata$MIPPIE_score > 0.6),]
print(dim(LZdata_midScore))
print(dim(LZdata_highScore))
print(head(LZdata_midScore))
print(head(LZdata_highScore))

print("Loop through all unique genes in LZ dataset...")
allLZgenes<-union(unique(LZdata$entrezA), unique(LZdata$entrezB))
LZcounts<-c()
for(g in allLZgenes){
	allCount<-length(union(which(LZdata$entrezA==g), which(LZdata$entrezB==g)))
	midCount<-length(union(which(LZdata_midScore$entrezA==g), which(LZdata_midScore$entrezB==g)))
	highCount<-length(union(which(LZdata_highScore$entrezA==g), which(LZdata_highScore$entrezB==g)))
	LZcounts<-rbind(LZcounts, c(g, allCount, midCount, highCount))
}
LZdf<-as.data.frame(LZcounts)
colnames(LZdf)<-c("Entrez","num_interactions","num_midScore_interactions","num_highScore_interactions")
write.table(LZdf, file="interaction_counts.LZinduced.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

print("Parsing RS data...")
RSdata<-read.table("mippie.RSinduced.txt", header=TRUE)
print(dim(RSdata))
print(head(RSdata))

print("Get mid- and high-scoring interactions...")
RSdata_midScore<-RSdata[which(RSdata$MIPPIE_score > 0.53),]
RSdata_highScore<-RSdata[which(RSdata$MIPPIE_score > 0.6),]
print(dim(RSdata_midScore))
print(dim(RSdata_highScore))
print(head(RSdata_midScore))
print(head(RSdata_highScore))

print("Loop through all unique genes in RS dataset...")
allRSgenes<-union(unique(RSdata$entrezA), unique(RSdata$entrezB))
RScounts<-c()
for(g in allRSgenes){
        allCount<-length(union(which(RSdata$entrezA==g), which(RSdata$entrezB==g)))
        midCount<-length(union(which(RSdata_midScore$entrezA==g), which(RSdata_midScore$entrezB==g)))
        highCount<-length(union(which(RSdata_highScore$entrezA==g), which(RSdata_highScore$entrezB==g)))
        RScounts<-rbind(RScounts, c(g, allCount, midCount, highCount))
}
RSdf<-as.data.frame(RScounts)
colnames(RSdf)<-c("Entrez","num_interactions","num_midScore_interactions","num_highScore_interactions")
write.table(RSdf, file="interaction_counts.RSinduced.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

print("Done counting interactions")
