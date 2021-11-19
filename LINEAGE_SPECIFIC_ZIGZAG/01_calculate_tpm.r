#PURPOSE: Read in a featurecounts output and convert raw counts to tpm; format for zigzag

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line argument. Enter the file name of a feature counts output file")
}

mycounts<-read.table(args[1],header=TRUE)
print(dim(mycounts))
#print(head(mycounts))

kb_len<-mycounts$Length/1000
print(head(kb_len))
rpk<-mycounts$Counts/kb_len
print(head(rpk))
scaling_factor<-sum(rpk)/1000000
print(scaling_factor)
tpm<-rpk/scaling_factor
print(length(tpm))
print(head(tpm))
print(tail(tpm))
print(tpm[100:110])
print(length(which(tpm>0)))

sample_name<-gsub("_.*","",gsub("\\.\\.\\/","",args[1]))
names(tpm)<-mycounts$Geneid
#out_table<-as.data.frame(cbind(mycounts$Geneid,tpm))
write.table(tpm,file=paste(sample_name,"tpm.txt",sep="."),row.names=TRUE,col.names=sample_name,quote=FALSE,sep="\t")
write.table(mycounts$Length,file=paste(sample_name,"length.txt",sep="."),row.names=mycounts$Geneid,col.names=FALSE,quote=FALSE,sep="\t")

print("Done with 01_calculate_tpm")
