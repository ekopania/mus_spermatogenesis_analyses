#PURPOSE: Figure out if genes identified as lineage specific by zigzag and the logFC>2 method overlap

library(VennDiagram)

#Read in arguments:
#cell type
args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line arguments.\nArgument 1: Cell type ('LZ' or 'RS)")
}
ct<-args[1]

print("Reading in data...")
#Read in logFC lineage specific genes
logfc_dom<-scan(paste0("../lineage_specific_gene_lists.logFC_enrichment_2.dom",ct,"_only.txt"),what=character())
logfc_mus<-scan(paste0("../lineage_specific_gene_lists.logFC_enrichment_2.mus",ct,"_only.txt"),what=character())
logfc_spr<-scan(paste0("../lineage_specific_gene_lists.logFC_enrichment_2.spr",ct,"_only.txt"),what=character())
logfc_dom_mus<-scan(paste0("../lineage_specific_gene_lists.logFC_enrichment_2.dom",ct,"_DM.txt"),what=character())

logfc_dom_x<-scan(paste0("../lineage_specific_gene_lists.logFC_enrichment_2.dom",ct,"_only_X.txt"),what=character())
logfc_mus_x<-scan(paste0("../lineage_specific_gene_lists.logFC_enrichment_2.mus",ct,"_only_X.txt"),what=character())
logfc_spr_x<-scan(paste0("../lineage_specific_gene_lists.logFC_enrichment_2.spr",ct,"_only_X.txt"),what=character())
logfc_dom_mus_x<-scan(paste0("../lineage_specific_gene_lists.logFC_enrichment_2.dom",ct,"_DM_X.txt"),what=character())

#Get autos only for logFC lineage specific genes
library(EnsDb.Mmusculus.v75)
edb<-EnsDb.Mmusculus.v75
edb_x<-addFilter(edb, SeqNameFilter("X"))
x_genes<-genes(edb_x)
edb_auto<-addFilter(edb, SeqNameFilter(c(1:19)))
auto_genes<-genes(edb_auto)

logfc_dom_auto<-intersect(logfc_dom,auto_genes$gene_id)
logfc_mus_auto<-intersect(logfc_mus,auto_genes$gene_id)
logfc_spr_auto<-intersect(logfc_spr,auto_genes$gene_id)
logfc_dom_mus_auto<-intersect(logfc_dom_mus,auto_genes$gene_id)

print("Sanity checks; these should all be TRUE")
all.equal(logfc_dom_x,intersect(logfc_dom,x_genes$gene_id))
all.equal(logfc_mus_x,intersect(logfc_mus,x_genes$gene_id))
all.equal(logfc_spr_x,intersect(logfc_spr,x_genes$gene_id))
all.equal(logfc_dom_mus_x,intersect(logfc_dom_mus,x_genes$gene_id))

#Read in zigzag lineage specific genes
zigzag_dom_auto<-scan(paste0("lineage_specific_genes.spec.dom_only_auto.",ct,".txt"),what=character())
zigzag_mus_auto<-scan(paste0("lineage_specific_genes.spec.mus_only_auto.",ct,".txt"),what=character())
zigzag_spr_auto<-scan(paste0("lineage_specific_genes.spec.spr_only_auto.",ct,".txt"),what=character())
zigzag_dom_mus_auto<-scan(paste0("lineage_specific_genes.spec.mus_dom_only_auto.",ct,".txt"),what=character())
zigzag_dom_x<-scan(paste0("lineage_specific_genes.spec.dom_only_x.",ct,".txt"),what=character())
zigzag_mus_x<-scan(paste0("lineage_specific_genes.spec.mus_only_x.",ct,".txt"),what=character())
zigzag_spr_x<-scan(paste0("lineage_specific_genes.spec.spr_only_x.",ct,".txt"),what=character())
zigzag_dom_mus_x<-scan(paste0("lineage_specific_genes.spec.mus_dom_only_x.",ct,".txt"),what=character())

#Make venn diagrams
print("Generating venn diagrams...")
vd_dom_auto<-venn.diagram(x=list(logfc_dom_auto,zigzag_dom_auto),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.dom_auto.",ct,".pdf"),output=TRUE)
vd_mus_auto<-venn.diagram(x=list(logfc_mus_auto,zigzag_mus_auto),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.mus_auto.",ct,".pdf"),output=TRUE)
vd_spr_auto<-venn.diagram(x=list(logfc_spr_auto,zigzag_spr_auto),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.spr_auto.",ct,".pdf"),output=TRUE)
vd_dom_mus_auto<-venn.diagram(x=list(logfc_dom_mus_auto,zigzag_dom_mus_auto),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.dom_mus_auto.",ct,".pdf"),output=TRUE)
vd_dom_x<-venn.diagram(x=list(logfc_dom_x,zigzag_dom_x),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.dom_x.",ct,".pdf"),output=TRUE)
vd_mus_x<-venn.diagram(x=list(logfc_mus_x,zigzag_mus_x),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.mus_x.",ct,".pdf"),output=TRUE)
vd_spr_x<-venn.diagram(x=list(logfc_spr_x,zigzag_spr_x),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.spr_x.",ct,".pdf"),output=TRUE)
vd_dom_mus_x<-venn.diagram(x=list(logfc_dom_mus_x,zigzag_dom_mus_x),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.dom_mus_x.",ct,".pdf"),output=TRUE)

#Get overlaps
print("Printing number of genes for logFC>2, number of genes for zigzag, and number of genes overlapping both sets:")
print("dom auto")
print(paste(length(logfc_dom_auto),length(zigzag_dom_auto),length(intersect(logfc_dom_auto,zigzag_dom_auto))))
print("mus auto")
print(paste(length(logfc_mus_auto),length(zigzag_mus_auto),length(intersect(logfc_mus_auto,zigzag_mus_auto))))
print("spr auto")
print(paste(length(logfc_spr_auto),length(zigzag_spr_auto),length(intersect(logfc_spr_auto,zigzag_spr_auto))))
print("dom_mus auto")
print(paste(length(logfc_dom_mus_auto),length(zigzag_dom_mus_auto),length(intersect(logfc_dom_mus_auto,zigzag_dom_mus_auto))))
print("dom x")
print(paste(length(logfc_dom_x),length(zigzag_dom_x),length(intersect(logfc_dom_x,zigzag_dom_x))))
print("mus x")
print(paste(length(logfc_mus_x),length(zigzag_mus_x),length(intersect(logfc_mus_x,zigzag_mus_x))))
print("spr x")
print(paste(length(logfc_spr_x),length(zigzag_spr_x),length(intersect(logfc_spr_x,zigzag_spr_x))))
print("dom_mus x")
print(paste(length(logfc_dom_mus_x),length(zigzag_dom_mus_x),length(intersect(logfc_dom_mus_x,zigzag_dom_mus_x))))

#Repeat for only genes considered "active" in that species by zigzag
print("Repeating for 'active' genes only...")
dom_active<-scan(paste0("active_genes.spec.dom.",ct,".txt"),what=character())
mus_active<-scan(paste0("active_genes.spec.mus.",ct,".txt"),what=character())
spr_active<-scan(paste0("active_genes.spec.spr.",ct,".txt"),what=character())

logfc_dom_auto_active<-intersect(logfc_dom_auto,dom_active)
logfc_mus_auto_active<-intersect(logfc_mus_auto,mus_active)
logfc_spr_auto_active<-intersect(logfc_spr_auto,spr_active)
logfc_dom_mus_auto_active<-Reduce(intersect,list(logfc_dom_mus_auto,dom_active,mus_active))
logfc_dom_x_active<-intersect(logfc_dom_x,dom_active)
logfc_mus_x_active<-intersect(logfc_mus_x,mus_active)
logfc_spr_x_active<-intersect(logfc_spr_x,spr_active)
logfc_dom_mus_x_active<-Reduce(intersect,list(logfc_dom_mus_x,dom_active,mus_active))

print("dom auto")
print(paste(length(logfc_dom_auto_active),length(zigzag_dom_auto),length(intersect(logfc_dom_auto_active,zigzag_dom_auto))))
print("mus auto")
print(paste(length(logfc_mus_auto_active),length(zigzag_mus_auto),length(intersect(logfc_mus_auto_active,zigzag_mus_auto))))
print("spr auto")
print(paste(length(logfc_spr_auto_active),length(zigzag_spr_auto),length(intersect(logfc_spr_auto_active,zigzag_spr_auto))))
print("dom_mus auto")
print(paste(length(logfc_dom_mus_auto_active),length(zigzag_dom_mus_auto),length(intersect(logfc_dom_mus_auto_active,zigzag_dom_mus_auto))))
print("dom x")
print(paste(length(logfc_dom_x_active),length(zigzag_dom_x),length(intersect(logfc_dom_x_active,zigzag_dom_x))))
print("mus x")
print(paste(length(logfc_mus_x_active),length(zigzag_mus_x),length(intersect(logfc_mus_x_active,zigzag_mus_x))))
print("spr x")
print(paste(length(logfc_spr_x_active),length(zigzag_spr_x),length(intersect(logfc_spr_x_active,zigzag_spr_x))))
print("dom_mus x")
print(paste(length(logfc_dom_mus_x_active),length(zigzag_dom_mus_x),length(intersect(logfc_dom_mus_x_active,zigzag_dom_mus_x))))

vd_dom_auto_active<-venn.diagram(x=list(logfc_dom_auto_active,zigzag_dom_auto),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.dom_auto.activeOnly.",ct,".pdf"),output=TRUE)
vd_mus_auto_active<-venn.diagram(x=list(logfc_mus_auto_active,zigzag_mus_auto),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.mus_auto.activeOnly.",ct,".pdf"),output=TRUE)
vd_spr_auto_active<-venn.diagram(x=list(logfc_spr_auto_active,zigzag_spr_auto),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.spr_auto.activeOnly.",ct,".pdf"),output=TRUE)
vd_dom_mus_auto_active<-venn.diagram(x=list(logfc_dom_mus_auto_active,zigzag_dom_mus_auto),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.dom_mus_auto.activeOnly.",ct,".pdf"),output=TRUE)
vd_dom_x_active<-venn.diagram(x=list(logfc_dom_x_active,zigzag_dom_x),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.dom_x.activeOnly.",ct,".pdf"),output=TRUE)
vd_mus_x_active<-venn.diagram(x=list(logfc_mus_x_active,zigzag_mus_x),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.mus_x.activeOnly.",ct,".pdf"),output=TRUE)
vd_spr_x_active<-venn.diagram(x=list(logfc_spr_x_active,zigzag_spr_x),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.spr_x.activeOnly.",ct,".pdf"),output=TRUE)
vd_dom_mus_x_active<-venn.diagram(x=list(logfc_dom_mus_x_active,zigzag_dom_mus_x),category.names=c("logfc","zigzag"),filename=paste0("venn_diagram.dom_mus_x.activeOnly.",ct,".pdf"),output=TRUE)

print("Done!")
