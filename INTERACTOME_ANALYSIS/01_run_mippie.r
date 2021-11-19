args<-commandArgs(TRUE)
if(length(args) != 2){
        stop("Missing command line arguments.\nArgument 1: query file\nArgument 2: output file")
}

#suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))
suppressPackageStartupMessages(library(EnsDb.Mmusculus.v75))
source("mippie_nc.R")

ppi<-"/mnt/beegfs/ek112884/REFERENCE_DIR/mippie_ppi_v1_0.tsv"
prots<-"/mnt/beegfs/ek112884/REFERENCE_DIR/mippie_proteins_v1_0.tsv"
queryList<-scan(args[1], what=character())

#edb<-EnsDb.Mmusculus.v79
edb<-EnsDb.Mmusculus.v75
myGenes<-genes(edb)
print(head(myGenes))

ppi_tab<-read.table(ppi, header=TRUE, sep="\t")
inMippie<-ppi_tab$entrezA

myEntrez<-unlist(myGenes$entrezid[which(myGenes$gene_id %in% queryList)])
print(head(myEntrez))
filtEntrez<-myEntrez[which(myEntrez %in% inMippie)]
write(unlist(filtEntrez), file="temp_entrez_query.txt", ncolumns=1)

#myNames<-myGenes$gene_name[which(myGenes$gene_id %in% queryList)]
#write(unlist(myNames), file="temp_name_query.txt", ncolumns=1)

mippie_nc(query.file="temp_entrez_query.txt", path.to.mippie=ppi, path.to.proteins=prots, output.file=args[2])

print(paste("Done constructing mippie network. Output is in file", args[2]))
