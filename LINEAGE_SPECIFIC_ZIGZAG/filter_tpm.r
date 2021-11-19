#PURPOSE: Remove all genes that have median tpm < 1
#NOTE: This is for visualization purposes and should NOT be the data that goes into zigzag

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line argument. Enter the file name of a *combined.tpm.txt file")
}
sample<-gsub(".combined.tpm.txt$","",args[1])
print(sample)
mydata<-read.table(args[1],header=TRUE)
print(dim(mydata))
filtered_data<-c()
for(i in 1:nrow(mydata)){
	if(median(as.numeric(as.character(mydata[i,])))>1){
		filtered_data<-rbind(filtered_data,mydata[i,])
	}
}
print(dim(filtered_data))
print(head(filtered_data))
write.table(filtered_data,file=paste(sample,"filtered.combined.tpm.txt",sep="."),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
lengthdata<-read.table(paste(sample,"length.txt",sep="."),header=FALSE)
newlengths<-lengthdata[which(lengthdata$V1 %in% rownames(filtered_data)),]
print(dim(newlengths))
print(head(newlengths))
write.table(newlengths,file=paste(sample,"filtered.length.txt",sep="."),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
print("Done!")
