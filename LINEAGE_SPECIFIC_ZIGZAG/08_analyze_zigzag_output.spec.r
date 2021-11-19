#PURPOSE: Analyze zigzag output. Zigzag output a table of probability of being active for each gene in dataset. This script parses this table to compare between two zigzag runs for a genotype and compare genotypes.

#Threschold cutoff to consider probabilities between the two runs "different"
diff_thresh<-0.01 #Running zigzag at least twice; if same gene has a difference in prob_active > diff_thresh between the two runs it may have issues and we may want to remove it from downstream analyses; not sure what this cutoff should be yet, I've tried 0.05, 0.01, and 0.005
active_thresh<-0.5 #Genes with prob_active > active_thresh for both zigzag runs are considered active; in the paper they use 0.5
inactive_thresh<-0.5 #Genes with prob_active < inactive_thresh for both zigzag runs are considered inactive; in the paper they use 0.5
ct<-"RS" #Cell type
#outfile<-"lineage_specific_gene_numbers.LZ.txt"

#Loop through all genotypes
myfiles<-list.files(pattern=paste0(ct,".spec.combined.tpm.txt$"))
myfiles<-c(myfiles,paste0("PAHNew_",ct,".combined.tpm.txt"))
print(myfiles)
mysamples<-unlist(lapply(myfiles,function(x) gsub(".combined.tpm.txt","",x)))
#mysamples[1]<-"DOM_LZ.spec" #Only when using DOM_LZ.rmProbSample
print(mysamples)

for(s in mysamples){
	print(s)
	#Read in data
	if(ct=="LZ"){
		mydata1<-read.table(paste0(s,"_zigzag_output1/",s,"_mcmc_output/",s,"_probability_active.tab"),header=TRUE)
		mydata2<-read.table(paste0(s,"_zigzag_output2/",s,"_mcmc_output/",s,"_probability_active.tab"),header=TRUE)
	} else{ #RS looked better w/ just one active component
		mydata1<-read.table(paste0(s,"_zigzag_output1_oneActiveComp/",s,"_mcmc_output/",s,"_probability_active.tab"),header=TRUE)
		mydata2<-read.table(paste0(s,"_zigzag_output2_oneActiveComp/",s,"_mcmc_output/",s,"_probability_active.tab"),header=TRUE)
	}
	print(dim(mydata1))
	print(dim(mydata2))
	#head(mydata1)
	#head(mydata2)
	#length(which(as.numeric(as.character(mydata1$prob_active))>0.9))
	
	#Check for genes with different results between the two runs (greater than <thresh> difference in prob)
	assign(paste0(s,"_diff_genes"),c())
	assign(paste0(s,"_active_genes"),c())
	assign(paste0(s,"_inactive_genes"),c())
	for(i in 1:nrow(mydata1)){
		this_gene<-as.character(mydata1$gene[i])
		prob1<-as.numeric(as.character(mydata1$prob_active[i]))
		prob2<-as.numeric(as.character(mydata2$prob_active[i]))
		#Check for genes with different results between the two runs (greater than <thresh> difference in prob)
		if( (prob1 < prob2) && ((prob1 + diff_thresh) < prob2) ){
				assign(paste0(s,"_diff_genes"),c(get(paste0(s,"_diff_genes")),this_gene))
		}
		if( (prob1 > prob2) && ((prob1 - diff_thresh) > prob2) ){
				assign(paste0(s,"_diff_genes"),c(get(paste0(s,"_diff_genes")),this_gene))
		}

	        #Zigzag paper considers >0.5 active and <0.5 active
		#Genes need to fall into same category in each run to be considered active or inactive here
		#In other words, throw out ambiguous things that have prob active 0.51 in one run and 0.49 in other
		if( (prob1 > active_thresh) && (prob2 > active_thresh) ){
                        assign(paste0(s,"_active_genes"),c(get(paste0(s,"_active_genes")),this_gene))
		}
                if( (prob1 < inactive_thresh) && (prob2 < inactive_thresh) ){
                        assign(paste0(s,"_inactive_genes"),c(get(paste0(s,"_inactive_genes")),this_gene))
                }
	}
	
	#Sanity checks and print number of active and inactive genes
	temp_diff<-get(paste0(s,"_diff_genes"))
	temp_active<-get(paste0(s,"_active_genes"))
	temp_inactive<-get(paste0(s,"_inactive_genes"))
	print("Diff:")
	#print(head(temp_diff))
	print(length(temp_diff))
	if(length(intersect(temp_active,temp_inactive)) > 0){
		stop("ERROR: Overlap between active and inactive genes")
	}
	print("Active:")
	print(length(temp_active))
	#print(head(temp_active))
	print("Inactive:")
	print(length(temp_inactive))
	#print(head(temp_inactive))
}

#Write active genes to files
if(ct=="LZ"){
	write(DOM_LZ.spec_active_genes,"active_genes.spec.dom.LZ.txt",sep="\t")
	write(MUS_LZ.spec_active_genes,"active_genes.spec.mus.LZ.txt",sep="\t")
	write(SPR_LZ.spec_active_genes,"active_genes.spec.spr.LZ.txt",sep="\t")
	write(PAHNew_LZ_active_genes,"active_genes.spec.pah.LZ.txt",sep="\t")
} else if(ct=="RS"){
	write(DOM_RS.spec_active_genes,"active_genes.spec.dom.RS.txt",sep="\t")
        write(MUS_RS.spec_active_genes,"active_genes.spec.mus.RS.txt",sep="\t")
        write(SPR_RS.spec_active_genes,"active_genes.spec.spr.RS.txt",sep="\t")
        write(PAHNew_RS_active_genes,"active_genes.spec.pah.RS.txt",sep="\t")
} else{
	stop("Invalid cell type")
}

#Genes that are active across all genotypes in a (sub)species and are uniquely active in that (sub)species
DOM_active<-get(paste0("DOM_",ct,".spec_active_genes"))
MUS_active<-get(paste0("MUS_",ct,".spec_active_genes"))
SPR_active<-get(paste0("SPR_",ct,".spec_active_genes"))
PAH_active<-get(paste0("PAHNew_",ct,"_active_genes"))

print("Domesticus only:")
not_dom<-Reduce(union, list(MUS_active,SPR_active,PAH_active))
dom_only<-DOM_active[which(!(DOM_active %in% not_dom))]
print(length(dom_only))
print(paste("Proportion of dom overlap active genes that are dom only:",length(dom_only)/length(DOM_active)))
print("Musculus only:")
not_mus<-Reduce(union, list(DOM_active,SPR_active,PAH_active))
mus_only<-MUS_active[which(!(MUS_active %in% not_mus))]
print(length(mus_only))
print(paste("Proportion of mus overlap active genes that are mus only:",length(mus_only)/length(MUS_active)))
print("Spretus only:")
not_spr<-Reduce(union, list(DOM_active,MUS_active,PAH_active))
spr_only<-SPR_active[which(!(SPR_active %in% not_spr))]
print(length(spr_only))
print(paste("Proportion of spr overlap active genes that are spr only:",length(spr_only)/length(SPR_active)))
print("Pahari only:")
not_pah<-Reduce(union, list(DOM_active,MUS_active,SPR_active))
pah_only<-PAH_active[which(!(PAH_active %in% not_pah))]
print(length(pah_only))
print(paste("Proportion of pah overlap active genes that are pah only:",length(pah_only)/length(PAH_active)))

#Genes that are active in all mus and dom strains but not spr or pah
print("Musculus and domesticus only:")
mus_and_dom<-Reduce(intersect, list(DOM_active,MUS_active))
not_mus_dom<-Reduce(union, list(SPR_active,PAH_active))
mus_dom_only<-mus_and_dom[which(!(mus_and_dom %in% not_mus_dom))]
print(length(mus_dom_only))
print(paste("Proportion of mus and dom overlap active genes that are in mus and/or dom only:", length(mus_dom_only)/length(mus_and_dom)))

print("Lengths of 'not' vectors:")
print(length(not_dom))
print(length(not_mus))
print(length(not_spr))
print(length(not_pah))
print(length(not_mus_dom))

#Separate X and auto
print("Separating X and auto genes...")
#library(EnsDb.Mmusculus.v79)
#edb<-EnsDb.Mmusculus.v79
library(EnsDb.Mmusculus.v75)
edb<-EnsDb.Mmusculus.v75
edb_y<-addFilter(edb, SeqNameFilter("Y"))
y_genes<-genes(edb_y)
edb_x<-addFilter(edb, SeqNameFilter("X"))
x_genes<-genes(edb_x)
edb_auto<-addFilter(edb, SeqNameFilter(c(1:19)))
auto_genes<-genes(edb_auto)

print("Domesticus only, auto:")
dom_overlap_auto<-intersect(DOM_active,auto_genes$gene_id)
dom_only_auto<-intersect(dom_only,auto_genes$gene_id)
print(length(dom_only_auto))
print(paste("Proportion of dom overlap active genes that are dom only - auto:",length(dom_only_auto)/length(dom_overlap_auto)))
write(dom_only_auto,paste0("lineage_specific_genes.spec.dom_only_auto.",ct,".txt"),sep="\t")
print("Musculus only, auto:")
mus_overlap_auto<-intersect(MUS_active,auto_genes$gene_id)
mus_only_auto<-intersect(mus_only,auto_genes$gene_id)
print(length(mus_only_auto))
print(paste("Proportion of mus overlap active genes that are mus only - auto:",length(mus_only_auto)/length(mus_overlap_auto)))
write(mus_only_auto,paste0("lineage_specific_genes.spec.mus_only_auto.",ct,".txt"),sep="\t")
print("Spretus only, auto:")
spr_overlap_auto<-intersect(SPR_active,auto_genes$gene_id)
spr_only_auto<-intersect(spr_only,auto_genes$gene_id)
print(length(spr_only_auto))
print(paste("Proportion of spr overlap active genes that are spr only - auto:",length(spr_only_auto)/length(spr_overlap_auto)))
write(spr_only_auto,paste0("lineage_specific_genes.spec.spr_only_auto.",ct,".txt"),sep="\t")
print("Musculus and domesticus only, auto:")
mus_and_dom_auto<-intersect(mus_and_dom,auto_genes$gene_id)
mus_dom_only_auto<-intersect(mus_dom_only,auto_genes$gene_id)
print(length(mus_dom_only_auto))
print(paste("Proportion of mus and dom overlap active genes that are in mus and/or dom only - auto:", length(mus_dom_only_auto)/length(mus_and_dom_auto)))
write(mus_dom_only_auto,paste0("lineage_specific_genes.spec.mus_dom_only_auto.",ct,".txt"),sep="\t")

print("Domesticus only, x:")
dom_overlap_x<-intersect(DOM_active,x_genes$gene_id)
dom_only_x<-intersect(dom_only,x_genes$gene_id)
print(length(dom_only_x))
print(paste("Proportion of dom overlap active genes that are dom only - x:",length(dom_only_x)/length(dom_overlap_x)))
write(dom_only_x,paste0("lineage_specific_genes.spec.dom_only_x.",ct,".txt"),sep="\t")
print("Musculus only, x:")
mus_overlap_x<-intersect(MUS_active,x_genes$gene_id)
mus_only_x<-intersect(mus_only,x_genes$gene_id)
print(length(mus_only_x))
print(paste("Proportion of mus overlap active genes that are mus only - x:",length(mus_only_x)/length(mus_overlap_x)))
write(mus_only_x,paste0("lineage_specific_genes.spec.mus_only_x.",ct,".txt"),sep="\t")
print("Spretus only, x:")
spr_overlap_x<-intersect(SPR_active,x_genes$gene_id)
spr_only_x<-intersect(spr_only,x_genes$gene_id)
print(length(spr_only_x))
print(paste("Proportion of spr overlap active genes that are spr only - x:",length(spr_only_x)/length(spr_overlap_x)))
write(spr_only_x,paste0("lineage_specific_genes.spec.spr_only_x.",ct,".txt"),sep="\t")
print("Musculus and domesticus only, x:")
mus_and_dom_x<-intersect(mus_and_dom,x_genes$gene_id)
mus_dom_only_x<-intersect(mus_dom_only,x_genes$gene_id)
print(length(mus_dom_only_x))
print(paste("Proportion of mus and dom overlap active genes that are in mus and/or dom only - x:", length(mus_dom_only_x)/length(mus_and_dom_x)))
write(mus_dom_only_x,paste0("lineage_specific_genes.spec.mus_dom_only_x.",ct,".txt"),sep="\t")

dom_col<-c(length(dom_only_auto),length(dom_overlap_auto),length(dom_only_x),length(dom_overlap_x))
mus_col<-c(length(mus_only_auto),length(mus_overlap_auto),length(mus_only_x),length(mus_overlap_x))
spr_col<-c(length(spr_only_auto),length(spr_overlap_auto),length(spr_only_x),length(spr_overlap_x))
musdom_col<-c(length(mus_dom_only_auto),length(mus_and_dom_auto),length(mus_dom_only_x),length(mus_and_dom_x))
out_table<-cbind(dom_col,mus_col,spr_col,musdom_col)
colnames(out_table)<-c("dom","mus","spr","musdom")
rownames(out_table)<-c("auto_lineage_specific","auto_all_active","x_lineage_specific","x_all_active")
write.table(out_table, file=paste0("lineage_specific_gene_numbers.",ct,".txt"), quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
q()
#Get genes with high expression but not active - DOM
print("Genes that are highly expressed but considered in active in dom:")
for(i in DOM_RS.spec_inactive_genes){
	expression_data<-read.table(myfiles[1], header = T, row.names = 1)
	logexp<-log(expression_data[,1])
	names(logexp)<-rownames(expression_data)
	if(i %in% names(logexp)){
		if(logexp[i] > 4){
			print(i)
		}
	}
}

print("Done!")
