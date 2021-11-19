#PURPOSE: extract CDS in bedfile format for longest txpt of all protein coding genes

args<-commandArgs(TRUE)
if(length(args)<1){
	stop("Enter the name of a CDS bedfile as a command line argument")
}
cdsfile<-as.character(args[1])
cdsdata<-read.table(cdsfile)
cdskeep<-NULL
current_gene<-""
current_txpt<-""
for(i in 1:nrow(cdsdata)){
	if(cdsdata$V4[i]!=current_gene){
		if(exists("longest_txpt")){
			cdskeep<-rbind(cdskeep,strsplit(longest_txpt,split=" ")[[1]])
		}
		current_gene<-cdsdata$V4[i]
		current_txpt<-do.call(paste,cdsdata[i,])
		longest_txpt<-do.call(paste,cdsdata[i,])
		current_len<-as.numeric(abs(cdsdata$V3[i]-cdsdata$V2[i]))
		longest_len<-current_len
		print(current_txpt)
		#Make sure last gene in file gets included if it only has one cds
		if(i==nrow(cdsdata)){
			cdskeep<-rbind(cdskeep,strsplit(longest_txpt,split=" ")[[1]])
		}
	}
	else if((cdsdata$V4[i]==current_gene) & (do.call(paste,cdsdata[i,])!=current_txpt)){
		if(current_len>longest_len){
			longest_len<-current_len
			longest_txpt<-current_txpt
		}
		current_txpt<-do.call(paste,cdsdata[i,])
		current_len<-as.numeric(abs(cdsdata$V3[i]-cdsdata$V2[i]))
	}
	else if((cdsdata$V4[i]==current_gene) & (do.call(paste,cdsdata[i,])==current_txpt)){
		current_len<-current_len+as.numeric(abs(cdsdata$V3[i]-cdsdata$V2[i]))
		if(current_len>longest_len){
			longest_len<-current_len
			longest_txpt<-current_txpt
		}
	}
	else{
		stop("ERROR: incremented txpt w/o incrementing gene")
	}
	print(paste(i,current_gene,current_txpt,current_len,longest_txpt,longest_len,sep=", "))
}
print(head(cdskeep))
cdskeep<-as.data.frame(cdskeep)
print(head(cdskeep))
###IF NEED TO CONVERT START TO ZERO_BASED
#zbase_start<-as.numeric(cdskeep$V2)-1
#zbase_cdskeep<-cbind(as.character(paste("chr",cdskeep$V1,sep="")),zbase_start,cdskeep$V3,as.character(cdskeep$V4),as.character(cdskeep$V5), as.character(cdskeep$V6))
#write.table(zbase_cdskeep, file="cds_in_longest_txpts2.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE, sep="\t")

#ALREADY CONVERTED START TO ZERO_BASED
write.table(cdskeep, file="cds_in_longest_txpts2.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE, sep="\t")
