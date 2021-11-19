#PURPOSE: count X-X interactions
#After re-reading the reviewer comment, it was specifically about X-X trans regulatory interactions, so need to seprate X-X interactions from X-auto interactions; counting each of these separately
#Only counting high-scoring interactions

#suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))
suppressPackageStartupMessages(library(EnsDb.Mmusculus.v75))

print("Reading in data...")
#Read in data
LZmippie<-read.table("mippie.LZinduced.txt", header=TRUE)
RSmippie<-read.table("mippie.RSinduced.txt", header=TRUE)

#Only counting high-scoring interactions
LZmippie_high<-LZmippie[which(LZmippie$MIPPIE_score > 0.6),]
RSmippie_high<-RSmippie[which(RSmippie$MIPPIE_score > 0.6),]

#edb<-EnsDb.Mmusculus.v79
edb<-EnsDb.Mmusculus.v75
myGenes<-genes(edb)
print(head(myGenes))
myMapLZ<-mapIds(edb, keys=as.character(union(LZmippie_high$entrezA, LZmippie_high$entrezB)), keytype="ENTREZID", column="GENEID", multiVals="first")
myMapRS<-mapIds(edb, keys=as.character(union(RSmippie_high$entrezA, RSmippie_high$entrezB)), keytype="ENTREZID", column="GENEID", multiVals="first")
edb_y<-addFilter(edb, SeqNameFilter("Y"))
y_genes<-genes(edb_y)
edb_x<-addFilter(edb, SeqNameFilter("X"))
x_genes<-genes(edb_x)
edb_auto<-addFilter(edb, SeqNameFilter(c(1:19)))
auto_genes<-genes(edb_auto)

#Get Ensembl IDs and chromosome types for each gene in each interaction
print("Appending Ensembl IDs and chromosomes")
LZmippie_wChrs<-c()
for(i in 1:nrow(LZmippie_high)){
        ea<-as.character(LZmippie_high$entrezA[i])
        eb<-as.character(LZmippie_high$entrezB[i])
        IDa<-myMapLZ[ea]
        IDb<-myMapLZ[eb]
        if(IDa %in% x_genes$gene_id){
                chrA<-"X"
        } else if(IDa %in% auto_genes$gene_id){
                chrA<-"auto"
        } else{
                chrA<-NA
        }
        if(IDb %in% x_genes$gene_id){
                chrB<-"X"
        } else if(IDb %in% auto_genes$gene_id){
                chrB<-"auto"
        } else{
                chrB<-NA
        }
        LZmippie_wChrs<-rbind(LZmippie_wChrs, c(chrA, chrB, IDa, IDb, LZmippie_high[i,]))
}
LZmippie_df<-as.data.frame(LZmippie_wChrs)
colnames(LZmippie_df)<-c("chrA","chrB","ensA","ensB",colnames(LZmippie_high))

RSmippie_wChrs<-c()
for(i in 1:nrow(RSmippie_high)){
        ea<-as.character(RSmippie_high$entrezA[i])
        eb<-as.character(RSmippie_high$entrezB[i])
        IDa<-myMapRS[ea]
        IDb<-myMapRS[eb]
        if(IDa %in% x_genes$gene_id){
                chrA<-"X"
        } else if(IDa %in% auto_genes$gene_id){
                chrA<-"auto"
        } else{
                chrA<-NA
        }
        if(IDb %in% x_genes$gene_id){
                chrB<-"X"
        } else if(IDb %in% auto_genes$gene_id){
                chrB<-"auto"
        } else{
                chrB<-NA
        }
        RSmippie_wChrs<-rbind(RSmippie_wChrs, c(chrA, chrB, IDa, IDb, RSmippie_high[i,]))
}
RSmippie_df<-as.data.frame(RSmippie_wChrs)
colnames(RSmippie_df)<-c("chrA","chrB","ensA","ensB",colnames(RSmippie_high))

print(head(LZmippie_df))
print(head(RSmippie_df))

#Separate by chromosome types involved in interactions
print("Separating by chromosome types...")
LZmippie_XX<-LZmippie_df[intersect(which(LZmippie_df$chrA=="X"), which(LZmippie_df$chrB=="X")),]
LZmippie_auto<-LZmippie_df[intersect(which(LZmippie_df$chrA=="auto"), which(LZmippie_df$chrB=="auto")),]
#Make sure you get both directions for X-auto (can be X-auto or auto-X)
LZmippie_Xauto<-LZmippie_df[intersect(which(LZmippie_df$chrA=="X"), which(LZmippie_df$chrB=="auto")),]
LZmippie_autoX<-LZmippie_df[intersect(which(LZmippie_df$chrA=="auto"), which(LZmippie_df$chrB=="X")),]
LZmippie_XA<-as.data.frame(rbind(LZmippie_Xauto, LZmippie_autoX))
print("LZ, X-X:")
print(dim(LZmippie_XX))
print(head(LZmippie_XX))
print("LZ, auto-auto:")
print(dim(LZmippie_auto))
print(head(LZmippie_auto))
print("LZ, X-auto or auto-X:")
print(dim(LZmippie_XA))
print(head(LZmippie_XA))

RSmippie_XX<-RSmippie_df[intersect(which(RSmippie_df$chrA=="X"), which(RSmippie_df$chrB=="X")),]
RSmippie_auto<-RSmippie_df[intersect(which(RSmippie_df$chrA=="auto"), which(RSmippie_df$chrB=="auto")),]
#Make sure you get both directions for X-auto (can be X-auto or auto-X)
RSmippie_Xauto<-RSmippie_df[intersect(which(RSmippie_df$chrA=="X"), which(RSmippie_df$chrB=="auto")),]
RSmippie_autoX<-RSmippie_df[intersect(which(RSmippie_df$chrA=="auto"), which(RSmippie_df$chrB=="X")),]
RSmippie_XA<-as.data.frame(rbind(RSmippie_Xauto, RSmippie_autoX))
print("RS, X-X:")
print(dim(RSmippie_XX))
print(head(RSmippie_XX))
print("RS, auto-auto:")
print(dim(RSmippie_auto))
print(head(RSmippie_auto))
print("RS, X-auto or auto-X:")
print(dim(RSmippie_XA))
print(head(RSmippie_XA))

#Count interactions
print("Counting interactions...")
#LZ
#X-X
LZ_XX_genes<-union(unique(LZmippie_XX$entrezA), unique(LZmippie_XX$entrezB))
LZ_XX_counts<-c()
for(g in LZ_XX_genes){
        myCount<-length(union(which(LZmippie_XX$entrezA==g), which(LZmippie_XX$entrezB==g)))
        LZ_XX_counts<-rbind(LZ_XX_counts, c(g, myCount))
}
LZ_XX_df<-as.data.frame(LZ_XX_counts)
colnames(LZ_XX_df)<-c("Entrez","num_highScore_XX_interactions")
write.table(LZ_XX_df, file="interaction_counts.LZinduced.XX.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#auto-auto
LZ_auto_genes<-union(unique(LZmippie_auto$entrezA), unique(LZmippie_auto$entrezB))
LZ_auto_counts<-c()
for(g in LZ_auto_genes){
        myCount<-length(union(which(LZmippie_auto$entrezA==g), which(LZmippie_auto$entrezB==g)))
        LZ_auto_counts<-rbind(LZ_auto_counts, c(g, myCount))
}
LZ_auto_df<-as.data.frame(LZ_auto_counts)
colnames(LZ_auto_df)<-c("Entrez","num_highScore_auto_interactions")
write.table(LZ_auto_df, file="interaction_counts.LZinduced.auto.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#X-auto
LZ_XA_genes<-union(unique(LZmippie_XA$entrezA), unique(LZmippie_XA$entrezB))
LZ_XA_counts<-c()
for(g in LZ_XA_genes){
        myCount<-length(union(which(LZmippie_XA$entrezA==g), which(LZmippie_XA$entrezB==g)))
        LZ_XA_counts<-rbind(LZ_XA_counts, c(g, myCount))
}
LZ_XA_df<-as.data.frame(LZ_XA_counts)
colnames(LZ_XA_df)<-c("Entrez","num_highScore_XA_interactions")
write.table(LZ_XA_df, file="interaction_counts.LZinduced.XA.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#RS
#X-X
RS_XX_genes<-union(unique(RSmippie_XX$entrezA), unique(RSmippie_XX$entrezB))
RS_XX_counts<-c()
for(g in RS_XX_genes){
        myCount<-length(union(which(RSmippie_XX$entrezA==g), which(RSmippie_XX$entrezB==g)))
        RS_XX_counts<-rbind(RS_XX_counts, c(g, myCount))
}
RS_XX_df<-as.data.frame(RS_XX_counts)
colnames(RS_XX_df)<-c("Entrez","num_highScore_XX_interactions")
write.table(RS_XX_df, file="interaction_counts.RSinduced.XX.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#auto-auto
RS_auto_genes<-union(unique(RSmippie_auto$entrezA), unique(RSmippie_auto$entrezB))
RS_auto_counts<-c()
for(g in RS_auto_genes){
        myCount<-length(union(which(RSmippie_auto$entrezA==g), which(RSmippie_auto$entrezB==g)))
        RS_auto_counts<-rbind(RS_auto_counts, c(g, myCount))
}
RS_auto_df<-as.data.frame(RS_auto_counts)
colnames(RS_auto_df)<-c("Entrez","num_highScore_auto_interactions")
write.table(RS_auto_df, file="interaction_counts.RSinduced.auto.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#X-auto
RS_XA_genes<-union(unique(RSmippie_XA$entrezA), unique(RSmippie_XA$entrezB))
RS_XA_counts<-c()
for(g in RS_XA_genes){
        myCount<-length(union(which(RSmippie_XA$entrezA==g), which(RSmippie_XA$entrezB==g)))
        RS_XA_counts<-rbind(RS_XA_counts, c(g, myCount))
}
RS_XA_df<-as.data.frame(RS_XA_counts)
colnames(RS_XA_df)<-c("Entrez","num_highScore_XA_interactions")
write.table(RS_XA_df, file="interaction_counts.RSinduced.XA.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Get median counts
print("Median counts:")
print(paste("X-X interactions, early:", median(LZ_XX_df$num_highScore_XX_interactions, na.rm=TRUE)))
print(paste("auto-auto interactions, early:", median(LZ_auto_df$num_highScore_auto_interactions, na.rm=TRUE)))
print(paste("X-auto interactions, early:", median(LZ_XA_df$num_highScore_XA_interactions, na.rm=TRUE)))
print(paste("X-X interactions, late:", median(RS_XX_df$num_highScore_XX_interactions, na.rm=TRUE)))
print(paste("auto-auto interactions, late:", median(RS_auto_df$num_highScore_auto_interactions, na.rm=TRUE)))
print(paste("X-auto interactions, late:", median(RS_XA_df$num_highScore_XA_interactions, na.rm=TRUE)))

print("Done with X-X interactions")
