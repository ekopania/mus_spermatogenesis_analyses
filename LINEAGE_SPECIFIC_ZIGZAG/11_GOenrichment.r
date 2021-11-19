#PURPOSE: Perform GO enrichment test

library(topGO)

cell_type<-"RS"
print(paste("Performing GO enrichment analysis for lineage specific genes in", cell_type))

#Custom funtion for converting active genes list into boolean (TRUE if lineage specific)
isLinSpec<-function(allActive) {return(allActive==1)}

#Get all lineage-specific genes
lineages<-c("dom","mus","spr","mus_dom")
for(i in lineages){
	print(paste("Working on lineage:", i))
	assign(paste0(i,"Auto"), scan(paste0("lineage_specific_genes.spec.",i,"_only_auto.",cell_type,".txt"), what=character()))
	assign(paste0(i,"X"), scan(paste0("lineage_specific_genes.spec.",i,"_only_x.",cell_type,".txt"), what=character()))
	assign(paste0(i,"LinSpecAll"), c(get(paste0(i,"Auto")), get(paste0(i,"X"))))
	if(i=="mus_dom"){
		temp1<-scan(paste0("active_genes.spec.mus.",cell_type,".txt"), what=character())
		temp2<-scan(paste0("active_genes.spec.dom.",cell_type,".txt"), what=character())
		assign(paste0(i,"Active"), union(temp1, temp2))
		print(head(mus_domActive))
	} else{
		assign(paste0(i,"Active"), scan(paste0("active_genes.spec.",i,".",cell_type,".txt"), what=character()))
	}
	#GO enrichment for this lineage only
	isLS<-sapply(get(paste0(i,"Active")), function(x) if(x %in% get(paste0(i,"LinSpecAll"))) 1 else 0)
	sampleGOdata<-new("topGOdata", description=paste("Lineage specific genes -", i), ontology="BP", allGenes=isLS, geneSel=isLinSpec, nodeSize=10, annot=annFUN.org, mapping="org.Mm.eg.db", ID="ensembl")
	resultFisher<-runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
	allRes<-GenTable(sampleGOdata, classicFisher=resultFisher, topNodes = 10)
	write.table(allRes, paste("GOresults", i, cell_type, "txt", sep="."), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

#GO enrichment for all lineage-specific genes, regardless of which lineage they are l-s in
all_ls<-c(domLinSpecAll, musLinSpecAll, sprLinSpecAll, mus_domLinSpecAll)
all_active<-Reduce(union, list(domActive, musActive, sprActive, mus_domActive))

isLS<-sapply(all_active, function(x) if(x %in% all_ls) 1 else 0)
sampleGOdata<-new("topGOdata", description="Lineage specific genes - any lineage", ontology="BP", allGenes=isLS, geneSel=isLinSpec, nodeSize=10, annot=annFUN.org, mapping="org.Mm.eg.db", ID="ensembl")
resultFisher<-runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
allRes<-GenTable(sampleGOdata, classicFisher=resultFisher, topNodes = 10)
write.table(allRes, paste("GOresults.anyLineage", cell_type, "txt", sep="."), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

print("Done with GO enrichment test")
