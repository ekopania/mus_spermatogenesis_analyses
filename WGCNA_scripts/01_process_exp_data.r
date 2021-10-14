#PURPOSE: Log transform normalize expression data; generate table of "traits" to associate with expression modules

print("Reading in normalized expressiond data...")
#These are filtered by expression level
myExp<-read.table("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/fpkm_filtered_table.rpkm1.txt")
print(myExp[1:5,1:5])
print(dim(myExp))

print("Performing log transformation: log2(x+1)...")
myLogData<-apply(myExp, c(1,2), function(x) log2(as.numeric(x)+1))
print(dim(myLogData))
write.table(myLogData, file="/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/WGCNA/logTransformed_expData.txt", quote=FALSE, sep="\t")

print("Done with 01_process_exp_data.r")
