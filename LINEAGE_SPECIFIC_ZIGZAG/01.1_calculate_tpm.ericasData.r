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
tpm_all<-c()
for(i in 7:ncol(mycounts)){
	print(colnames(mycounts))[i]
	rpk<-mycounts[,i]/kb_len
	print(head(rpk))
	scaling_factor<-sum(rpk)/1000000
	print(scaling_factor)
	tpm<-rpk/scaling_factor
	print(length(tpm))
	print(head(tpm))
	print(tail(tpm))
#	print(tpm[100:110])
	print(length(which(tpm>0)))
	tpm_all<-cbind(tpm_all,tpm)
}
tpm_all<-as.data.frame(tpm_all)
#sample_name<-gsub("_.*","",gsub("\\.\\.\\/","",args[1]))
colnames(tpm_all)<-colnames(mycounts)[7:ncol(mycounts)]
rownames(tpm_all)<-mycounts$Geneid
print(head(tpm_all))
tpm_SP<-c()
tpm_LZ<-c()
tpm_DIP<-c()
tpm_RS<-c()
samps_SP<-c()
samps_LZ<-c()
samps_DIP<-c()
samps_RS<-c()
for(i in 1:ncol(tpm_all)){
	sample<-colnames(tpm_all)[i]
	if("SP" %in% unlist(strsplit(sample,"\\."))){
		tpm_SP<-cbind(tpm_SP,tpm_all[,i])
		samps_SP<-c(samps_SP,sample)
	} else if("LZ" %in% unlist(strsplit(sample,"\\."))){
		tpm_LZ<-cbind(tpm_LZ,tpm_all[,i])
		samps_LZ<-c(samps_LZ,sample)
	} else if("DIP" %in% unlist(strsplit(sample,"\\."))){
		tpm_DIP<-cbind(tpm_DIP,tpm_all[,i])
		samps_DIP<-c(samps_DIP,sample)
	} else if("RS" %in% unlist(strsplit(sample,"\\."))){
		tpm_RS<-cbind(tpm_RS,tpm_all[,i])
		samps_RS<-c(samps_RS,sample)
	} else{
		print(paste("Problem sample:",colnames(tpm_all)[i]))
	}
}
tpm_SP<-as.data.frame(tpm_SP)
colnames(tpm_SP)<-samps_SP
rownames(tpm_SP)<-rownames(tpm_all)
tpm_LZ<-as.data.frame(tpm_LZ)
colnames(tpm_LZ)<-samps_LZ
rownames(tpm_LZ)<-rownames(tpm_all)
tpm_DIP<-as.data.frame(tpm_DIP)
colnames(tpm_DIP)<-samps_DIP
rownames(tpm_DIP)<-rownames(tpm_all)
tpm_RS<-as.data.frame(tpm_RS)
colnames(tpm_RS)<-samps_RS
rownames(tpm_RS)<-rownames(tpm_all)
print(dim(tpm_SP))
print(head(tpm_SP))
print(dim(tpm_LZ))
print(head(tpm_LZ))
print(dim(tpm_DIP))
print(head(tpm_DIP))
print(dim(tpm_RS))
print(head(tpm_RS))
#out_table<-as.data.frame(cbind(mycounts$Geneid,tpm))
write.table(tpm_all,file="erica_samples.combined.tpm.txt",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")
write.table(mycounts$Length,file="erica_samples.length.txt",row.names=mycounts$Geneid,col.names=FALSE,quote=FALSE,sep="\t")
write.table(tpm_SP,file="erica_samples_SP.combined.tpm.txt",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")
write.table(tpm_LZ,file="erica_samples_LZ.combined.tpm.txt",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")
write.table(tpm_DIP,file="erica_samples_DIP.combined.tpm.txt",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")
write.table(tpm_RS,file="erica_samples_RS.combined.tpm.txt",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")

print("Done with 01.1_calculate_tpm_ericasData")
