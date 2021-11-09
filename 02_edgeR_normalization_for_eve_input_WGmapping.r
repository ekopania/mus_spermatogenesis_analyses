#PURPOSE: Normalize gene expression data, filter by fpkm, using EdgeR to generate data fr EVE input

#NOTE: Run this from the directory containing the count files you want to input into EdgeR

#change these to test different parameters
min_rpkm<-1
min_samples<-8
min_logFC<-0
induced_cutoff<-2

print("Loading EdgeR and setting up DGE object...")
library(edgeR)
#Mus ref geneid from compara appended:
#files<-list.files(pattern=".compara.txt$")
#Mus ref geneid from ensembl pairwise 1:1 orthos appended:
files<-list.files(pattern=".orthoIDappended.txt$")
mygroups<-c(1,2,1,2,3,4,3,4,1,2,1,2,1,2,1,2,1,2,1,2,3,4,3,4,3,4,7,8,7,8,7,8,3,4,3,4,3,4,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,1,2,1,2,1,2,1,2)
myDGE<-readDGE(files,columns=c(2,8),group=mygroups)

print("Filtering by expression level...")
#Get gene lengths for each species and for mouse reference
#Getting from feature counts output b/c these are the correct lengths to use for fpkm calc (no introns, etc)
fc_dom<-read.table("BIK4665-1M-LZ_counts.orthoIDappended.txt",header=TRUE)
fc_domLengths<-cbind(as.numeric(fc_dom$Length),as.character(fc_dom$Geneid))
head(fc_domLengths)
fc_mus<-read.table("PPPP98-3M-LZ_counts.orthoIDappended.txt",header=TRUE)
fc_musLengths<-cbind(as.numeric(fc_mus$Length),as.character(fc_mus$Geneid))
head(fc_musLengths)
fc_spr<-read.table("SEG4156-LZ_counts.orthoIDappended.txt",header=TRUE)
fc_sprLengths<-cbind(as.numeric(fc_spr$Length),as.character(fc_spr$Geneid))
head(fc_sprLengths)
fc_pah<-read.table("PAHNew-1M-LZ_counts.orthoIDappended.txt",header=TRUE)
fc_pahLengths<-cbind(as.numeric(fc_pah$Length),as.character(fc_pah$Geneid))
head(fc_pahLengths)
fc_ref<-read.table("../TEST_TOPHAT_MAPPING/MUS_REF_MAPPINGS/PPPP98-3M-LZ.counts.txt",header=TRUE,skip=1)
fc_refLengths<-cbind(as.numeric(fc_ref$Length),as.character(fc_ref$Geneid))
head(fc_refLengths)
#overlap between genes included in species reference mapping and mouse reference mapping
overlap_fc_domLengths<-fc_domLengths[which(fc_domLengths[,2] %in% fc_refLengths[,2]),]
overlap_fc_musLengths<-fc_musLengths[which(fc_musLengths[,2] %in% fc_refLengths[,2]),]
overlap_fc_sprLengths<-fc_sprLengths[which(fc_sprLengths[,2] %in% fc_refLengths[,2]),]
overlap_fc_pahLengths<-fc_pahLengths[which(fc_pahLengths[,2] %in% fc_refLengths[,2]),]
overlap_fc_refLengths<-fc_refLengths[which(fc_refLengths[,2] %in% fc_domLengths[,2]),]
#Sort by gene name
ordered_fc_domLengths<-overlap_fc_domLengths[match(rownames(myDGE$counts),overlap_fc_domLengths[,2]),]
ordered_fc_musLengths<-overlap_fc_musLengths[match(rownames(myDGE$counts),overlap_fc_musLengths[,2]),]
ordered_fc_sprLengths<-overlap_fc_sprLengths[match(rownames(myDGE$counts),overlap_fc_sprLengths[,2]),]
ordered_fc_pahLengths<-overlap_fc_pahLengths[match(rownames(myDGE$counts),overlap_fc_pahLengths[,2]),]
ordered_fc_refLengths<-overlap_fc_refLengths[match(rownames(myDGE$counts),overlap_fc_refLengths[,2]),]
#Figure out FPKM for each species separately, using proper gene lengths for each
rpkm_dom<-rpkm(myDGE[,which(myDGE$samples$group==c(1,2))],gene.length=as.numeric(ordered_fc_domLengths[,1]))
rpkm_mus<-rpkm(myDGE[,which(myDGE$samples$group==c(3,4))],gene.length=as.numeric(ordered_fc_musLengths[,1]))
rpkm_spr<-rpkm(myDGE[,which(myDGE$samples$group==c(5,6))],gene.length=as.numeric(ordered_fc_sprLengths[,1]))
rpkm_pah<-rpkm(myDGE[,which(myDGE$samples$group==c(7,8))],gene.length=as.numeric(ordered_fc_pahLengths[,1]))
rpkm_all<-cbind(rpkm_dom,rpkm_mus,rpkm_spr,rpkm_pah)
new_rpkm_all<-rpkm_all[which(rownames(rpkm_all) %in% rownames(myDGE$counts)),]
write.table(new_rpkm_all,paste0("fpkm_full_table.rpkm",min_rpkm,".txt"),quote=FALSE,sep="\t")

keep<-rowSums(rpkm_all > min_rpkm) >= min_samples
keep[is.na(keep)]<-FALSE #NAs show up for gene iDs that aren't in myDGE; replace then with false otherwise myDGE gets confused
#ALTERNATIVELY figure out FPKM for all species based on mouse reference gene lengths
#keep<-rowSums(rpkm(myDGE, gene.length=as.numeric(geneLengths[,1]))>min_rpkm) >= min_samples

myDGE<-myDGE[keep, , keep.lib.sizes=FALSE]
myDGE<-calcNormFactors(myDGE)

spec<-c()
for(i in myDGE$samples$group){
	if((i==1) | (i==2)){
		spec<-c(spec,"dom")
	}
	else if((i==3) | (i==4)){
		spec<-c(spec,"mus")
	}
	else if((i==5) | (i==6)){
		spec<-c(spec,"spr")
	}
	else{
		spec<-c(spec,"pah")
	}
}

ct<-c()
for(i in myDGE$samples$files){
	temp<-gsub(".*-","",i)
	#ct<-c(ct,sub(".counts.refIDappended.compara.txt","",temp))
	ct<-c(ct,sub(".counts.orthoIDappended.txt","",temp))
}

design_table<-as.data.frame(cbind(rownames(myDGE$samples), spec, ct, paste(spec,ct,sep=".")))
colnames(design_table)<-c("file","species","cell.type","group")
design_matrix<-model.matrix(~0+design_table$group, data=myDGE$samples)
colnames(design_matrix)<-levels(design_table$group)
myDGE<-estimateDisp(myDGE,design_matrix)

print("Setting up contrasts...")
myContrasts<-makeContrasts(MusLZvsDomLZ="mus.LZ-dom.LZ", SprLZvsDomLZ="spr.LZ-dom.LZ", PahLZvsDomLZ="pah.LZ-dom.LZ", SprLZvsMusLZ="spr.LZ-mus.LZ", PahLZvsMusLZ="pah.LZ-mus.LZ", PahLZvsSprLZ="pah.LZ-spr.LZ", MusRSvsDomRS="mus.RS-dom.RS", SprRSvsDomRS="spr.RS-dom.RS", PahRSvsDomRS="pah.RS-dom.RS", SprRSvsMusRS="spr.RS-mus.RS", PahRSvsMusRS="pah.RS-mus.RS", PahRSvsSprRS="pah.RS-spr.RS", DomRSvsDomLZ="dom.RS-dom.LZ", MusRSvsMusLZ="mus.RS-mus.LZ", SprRSvsSprLZ="spr.RS-spr.LZ", PahRSvsPahLZ="pah.RS-pah.LZ", levels=design_table$group)
#save lengths for all genes in myDGE; will need later to calculate FPKM
domLengths<-ordered_fc_domLengths[ordered_fc_domLengths[,2] %in% rownames(myDGE$counts),]
musLengths<-ordered_fc_musLengths[ordered_fc_musLengths[,2] %in% rownames(myDGE$counts),]
sprLengths<-ordered_fc_sprLengths[ordered_fc_sprLengths[,2] %in% rownames(myDGE$counts),]
pahLengths<-ordered_fc_pahLengths[ordered_fc_pahLengths[,2] %in% rownames(myDGE$counts),]
refLengths<-ordered_fc_refLengths[ordered_fc_refLengths[,2] %in% rownames(myDGE$counts),]

for(i in colnames(myDGE$design)){
	result<-rownames(myDGE$design)[which(myDGE$design[,i]==1)]
	assign(paste("samples",i,sep="."),result)
}
#Calculate fpkm again for each group, with normalized myDGE and new set of genes
fpkm_dom<-rpkm(myDGE[,which(myDGE$samples$group==c(1,2))],gene.length=as.numeric(domLengths[,1]))
fpkm_mus<-rpkm(myDGE[,which(myDGE$samples$group==c(3,4))],gene.length=as.numeric(musLengths[,1]))
fpkm_spr<-rpkm(myDGE[,which(myDGE$samples$group==c(5,6))],gene.length=as.numeric(sprLengths[,1]))
fpkm_pah<-rpkm(myDGE[,which(myDGE$samples$group==c(7,8))],gene.length=as.numeric(pahLengths[,1]))
fpkm<-cbind(fpkm_dom,fpkm_mus,fpkm_spr,fpkm_pah)
write.table(fpkm,paste0("fpkm_filtered_table.rpkm",min_rpkm,".txt"),quote=FALSE,sep="\t")

print("Determining which genes are expressed in each species or cell type...")
include<-c()
for(i in 1:nrow(fpkm)){
      	if(FALSE %in% (samples.dom.LZ %in% colnames(fpkm)[which(fpkm[i,]>min_rpkm)])){
       	}
       	else{
       		include<-c(include,rownames(fpkm)[i])
       	}
}
DomLZ_expressed<-include
        
include<-c()
for(i in 1:nrow(fpkm)){
       	if(FALSE %in% (samples.dom.RS %in% colnames(fpkm)[which(fpkm[i,]>min_rpkm)])){
       	}
       	else{
       		include<-c(include,rownames(fpkm)[i])
       	}
}
DomRS_expressed<-include

include<-c()
for(i in 1:nrow(fpkm)){
       	if(FALSE %in% (samples.mus.LZ %in% colnames(fpkm)[which(fpkm[i,]>min_rpkm)])){
       	}
       	else{
       		include<-c(include,rownames(fpkm)[i])
       	}
}
MusLZ_expressed<-include

include<-c()
for(i in 1:nrow(fpkm)){
      	if(FALSE %in% (samples.mus.RS %in% colnames(fpkm)[which(fpkm[i,]>min_rpkm)])){
       	}
       	else{
       		include<-c(include,rownames(fpkm)[i])
       	}
}
MusRS_expressed<-include

include<-c()
for(i in 1:nrow(fpkm)){
	if(FALSE %in% (samples.spr.LZ %in% colnames(fpkm)[which(fpkm[i,]>min_rpkm)])){
        }
        else{
        	include<-c(include,rownames(fpkm)[i])
        }
}
SprLZ_expressed<-include
        
include<-c()
for(i in 1:nrow(fpkm)){
	if(FALSE %in% (samples.spr.RS %in% colnames(fpkm)[which(fpkm[i,]>min_rpkm)])){
        }
        else{
        	include<-c(include,rownames(fpkm)[i])
        }
}
SprRS_expressed<-include

include<-c()
for(i in 1:nrow(fpkm)){
        if(FALSE %in% (samples.pah.LZ %in% colnames(fpkm)[which(fpkm[i,]>min_rpkm)])){
        }
        else{
	        include<-c(include,rownames(fpkm)[i])
        }
}
PahLZ_expressed<-include
        
include<-c()
for(i in 1:nrow(fpkm)){
	if(FALSE %in% (samples.pah.RS %in% colnames(fpkm)[which(fpkm[i,]>min_rpkm)])){
        }
        else{
        	include<-c(include,rownames(fpkm)[i])
        }
}
PahRS_expressed<-include

print("Here are the numbers of genes produced in each combination of species and cell type:")
length(DomLZ_expressed)
length(DomRS_expressed)
length(MusLZ_expressed)
length(MusRS_expressed)
length(SprLZ_expressed)
length(SprRS_expressed)
length(PahLZ_expressed)
length(PahRS_expressed)

FPKM<-fpkm
#PHYLO DATA - two cell types
LZ<-c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67)
RS<-c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68)
#PHYLO DATA - NO PAHARI, two cell types, 
#LZ<-c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,42,44,46,48,50,52,54,56,58)
#RS<-c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,41,43,45,47,49,51,53,55,57,59)
FPKM_LZ<-FPKM[,LZ]
FPKM_RS<-FPKM[,RS]
#Generate file of FPKM data for all LZ_expressed and RS_expressed genes
#formated for EVEmodel input
LZ_expressed<-union(DomLZ_expressed,MusLZ_expressed)
LZ_expressed<-union(LZ_expressed,SprLZ_expressed)
LZ_expressed<-union(LZ_expressed,PahLZ_expressed)
if(length(LZ_expressed)!=length(unique(LZ_expressed))){
    print("ERROR: LZ_expressed contains repeated gene names")
}
LZ_fpkm_output<-FPKM_LZ[which(rownames(FPKM_LZ) %in% LZ_expressed),]
print(dim(LZ_fpkm_output))
error<-FALSE
for(i in LZ_expressed){
    if(!(i %in% rownames(LZ_fpkm_output))){
        print(i)
        error<-TRUE
    }
}
for(i in rownames(LZ_fpkm_output)){
    if(!(i %in% LZ_expressed)){
        print(i)
        error<-TRUE
    }
}
if(error){
    print("ERROR: LZ_fpkm_output doesn't match LZ_expressed; look at what's printed to screen for problem genes")
}
nonzero<-apply(LZ_fpkm_output, 1, function(row) all(row !=0 ))
nonzero_LZ_fpkm_output<-LZ_fpkm_output[nonzero,]
print("Here are the first few lines of the data that will go into the EVE input file:")
print(head(nonzero_LZ_fpkm_output))
print("Here's the dimensions of the table after zeros were removed:")
print(dim(nonzero_LZ_fpkm_output))
print("And before removing zeros (should be more rows):")
print(dim(LZ_fpkm_output))
write(dim(nonzero_LZ_fpkm_output)[1], file=paste0("eve_expression_input_LZ_edgeR_wholeGenome.ensemblOrthos.rpkm",min_rpkm,".txt"), append=FALSE, sep="\t")
write.table(nonzero_LZ_fpkm_output, file=paste0("eve_expression_input_LZ_edgeR_wholeGenome.ensemblOrthos.rpkm",min_rpkm,".txt"), quote=FALSE, col.names=FALSE, append=TRUE, sep="\t")#
write(colnames(nonzero_LZ_fpkm_output),file=paste0("eve_sample_order_LZ_edgeR.ensemblOrthos.rpkm",min_rpkm,".txt"),sep="\n")
write(rownames(nonzero_LZ_fpkm_output),file=paste0("gene_list_eve_LZ_edgeR_wholeGenome.ensemblOrthos.rpkm",min_rpkm,".txt"),sep="\n")

RS_expressed<-union(DomRS_expressed,MusRS_expressed)
RS_expressed<-union(RS_expressed,SprRS_expressed)
RS_expressed<-union(RS_expressed,PahRS_expressed)
if(length(RS_expressed)!=length(unique(RS_expressed))){
    print("ERROR: RS_expressed contains repeated gene names")
}
RS_fpkm_output<-FPKM_RS[which(rownames(FPKM_RS) %in% RS_expressed),]
print(dim(RS_fpkm_output))
error<-FALSE
for(i in RS_expressed){
    if(!(i %in% rownames(RS_fpkm_output))){
        print(i)
        error<-TRUE
    }
}
for(i in rownames(RS_fpkm_output)){
    if(!(i %in% RS_expressed)){
        print(i)
        error<-TRUE
    }
}
if(error){
    print("ERROR: RS_fpkm_output doesn't match RS_expressed; look at what's printed to screen for problem genes")
}
nonzero<-apply(RS_fpkm_output, 1, function(row) all(row !=0 ))
nonzero_RS_fpkm_output<-RS_fpkm_output[nonzero,]
print("Here are the data that will go into the EVE input file for RS:")
print(head(nonzero_RS_fpkm_output))
print("Here's RS dimensions after removing zeros:")
print(dim(nonzero_RS_fpkm_output))
print("And before removing zeros (should be more rows than above):")
print(dim(RS_fpkm_output))
write(dim(nonzero_RS_fpkm_output)[1], file=paste0("eve_expression_input_RS_edgeR_wholeGenome.ensemblOrthos.rpkm",min_rpkm,".txt"), append=FALSE, sep="\t")
write.table(nonzero_RS_fpkm_output, file=paste0("eve_expression_input_RS_edgeR_wholeGenome.ensemblOrthos.rpkm",min_rpkm,".txt"), quote=FALSE, col.names=FALSE, append=TRUE, sep="\t")#
write(colnames(nonzero_RS_fpkm_output),file=paste0("eve_sample_order_RS_edgeR_wholeGenome.ensemblOrthos.rpkm",min_rpkm,".txt"),sep="\n")
write(rownames(nonzero_RS_fpkm_output),file=paste0("gene_list_eve_RS_edgeR_wholeGenome.ensemblOrthos.rpkm",min_rpkm,".txt"),sep="\n")

#determine genes induced (not just expressed) in each cell type
#induced: median FPKM in one cell type > induced_cutoff*median in other cell type
#induced_cutoff set at beginning of script
LZ_meds<-apply(FPKM_LZ,1,median)
if(length(LZ_meds) != nrow(FPKM)){
    print("ERROR: number of LZ medians not equal to number of genes in analysis")
}
LZ_meds<-as.table(LZ_meds)#
RS_meds<-apply(FPKM_RS,1,median)
if(length(RS_meds) != nrow(FPKM)){
    print("ERROR: number of RS medians not equal to number of genes in analysis")
}
RS_meds<-as.table(RS_meds)#
#Two cell types
LZ_induced<-c()
for(i in rownames(LZ_meds)){
    if(i %in% rownames(RS_meds)){
        if((LZ_meds[i]>induced_cutoff*RS_meds[i]) & (LZ_meds[i]>1) & (min(FPKM_LZ[i,])>1)){
            LZ_induced<-c(LZ_induced,i)
        }
    }
    else if((!(i %in% rownames(RS_meds)))  & (LZ_meds[i]>1) & (min(FPKM_LZ[i,])>1)){
        LZ_induced<-c(LZ_induced,i)
    }
}
RS_induced<-c()
for(i in rownames(RS_meds)){
    if((i %in% rownames(LZ_meds)) & (RS_meds[i]>1) & (min(FPKM_RS[i,])>1)){
        if(RS_meds[i]>induced_cutoff*LZ_meds[i]){
            RS_induced<-c(RS_induced,i)
        }
    }
    else if((!(i %in% rownames(LZ_meds))) & (RS_meds[i]>1) & (min(FPKM_RS[i,])>1)){
        RS_induced<-c(RS_induced,i)
    }
}
print("Here is length and head for LZ and RS:")
print(length(LZ_induced))
print(length(RS_induced))
print(head(LZ_induced))
print(head(RS_induced))
for(i in LZ_induced){
    if(i %in% RS_induced){
        print(paste("ERROR: ", i, "cannot be induced in both LZ and RS"))
    }
}
for(i in RS_induced){
    if(i %in% LZ_induced){
        print(paste("ERROR: ", i, "cannot be induced in both LZ and RS"))
    }
}
write(LZ_induced, file=paste0("gene_list_LZinduced_edgeR_wholeGenome.ensemblOrthos.rpkm",min_rpkm,".txt"),append=FALSE)
write(RS_induced, file=paste0("gene_list_RSinduced_edgeR_wholeGenome.ensemblOrthos.rpkm",min_rpkm,".txt"),append=FALSE)
#FPKM_LZ_induced<-FPKM_LZ[LZ_induced,]
#FPKM_RS_induced<-FPKM_RS[RS_induced,]
FPKM_LZ_induced<-nonzero_LZ_fpkm_output[which(rownames(nonzero_LZ_fpkm_output) %in% LZ_induced),]
FPKM_RS_induced<-nonzero_RS_fpkm_output[which(rownames(nonzero_RS_fpkm_output) %in% RS_induced),]
print("Here are the final dimensions for LZ and RS induced data tables:")
print(dim(FPKM_LZ_induced))
print(dim(FPKM_RS_induced))
write(length(LZ_induced), file=paste0("eve_expression_input_LZinduced_edgeR_wholeGenome.ensemblOrthos.rpkm",min_rpkm,".txt"), append=FALSE, sep="\t")
write.table(FPKM_LZ_induced, file=paste0("eve_expression_input_LZinduced_edgeR_wholeGenome.ensemblOrthos.rpkm",min_rpkm,".txt"), quote=FALSE, col.names=FALSE, append=TRUE, sep="\t")
write(length(RS_induced), file=paste0("eve_expression_input_RSinduced_edgeR_wholeGenome.ensemblOrthos.rpkm",min_rpkm,".txt"), append=FALSE, sep="\t")
write.table(FPKM_RS_induced, file=paste0("eve_expression_input_RSinduced_edgeR_wholeGenome.ensemblOrthos.rpkm",min_rpkm,".txt"), quote=FALSE, col.names=FALSE, append=TRUE, sep="\t")

print("Done!")
