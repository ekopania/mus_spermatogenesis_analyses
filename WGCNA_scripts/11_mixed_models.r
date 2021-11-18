#PURPOSE: Run mixed effects models to determine the relationships between module eigenvalues and lineages

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line arguments.\nArgument 1: dataset")
}

dataset<-args[1]

print(paste0("Running gene linear models for ", dataset))

library(lme4)
library(lmerTest)
library(multcomp)
library(multcompView)

#Load WGCNA data for dataset
load(paste("05-WGCNA",dataset,"RData", sep="."))

#Get column names as module colors for MEs
colname_nums<-sapply(colnames(MEs), function(x) gsub("ME","",x))
colname_cols<-sapply(colname_nums, function(x) label_to_color$col[which(label_to_color$lab==x)])
colnames(MEs)<-colname_cols
print(dim(MEs))
print(MEs[1:5,1:5])

#Merge MEs and sample data
#Get just the strain name for each sample
strain.type<-unlist(lapply(rownames(MEs), function(x) gsub("New\\.","",gsub("[0-9].*","",x))))
#print(strain.type)
#Get cell type for each sample; the first gsub deals with weird spretus sample names
cell.type<-unlist(lapply(rownames(MEs), function(x) gsub("^S.*\\.","",gsub("^.*\\.","", gsub("_counts.*","",x)))))
#print(cell.type)

#Set up to have lineage as factor
lin.type<-c()
for(i in strain.type){
        if(i %in% c("BIK","DGA","LLLL","WWWW")){
                lin.type<-c(lin.type, "dom")
        } else if(i %in% c("CCCC","MBS","PPPP")){
                lin.type<-c(lin.type, "mus")
        } else if(i %in% c("SEG","SFM","STF")){
                lin.type<-c(lin.type, "spr")
        } else if(i=="PAH"){
                lin.type<-c(lin.type, "pah")
        }
}
final_data<-as.data.frame(cbind(strain.type, lin.type, cell.type, MEs))
rownames(final_data)<-rownames(MEs)
final_data$strain.type<-as.factor(final_data$strain.type)
final_data$lin.type<-as.factor(final_data$lin.type)
final_data$cell.type<-as.factor(final_data$cell.type)

#Loop through each module
#MEs not KMEs!!! (by individual, not gene)
#pdf(paste("lmer_plots","exp1","pdf", sep="."), onefile=TRUE)
tukey_lineage<-c()
tukey_strain<-c()
tukey_ct<-c()
tukey_lineageXct<-c()
tukey_strainXct<-c()
for(m in colnames(MEs)){
        print(paste("Working on module", m))
        #Eigenvalue for each module with lineage type as fixed effect; perform Tukey test for pairwise comparisons among lineage types
        lm_linType<-lm(get(m) ~ lin.type, data=final_data)
        tukey_result<-summary(glht(lm_linType, mcp(lin.type="Tukey")))
        print(tukey_result)
        tukey_df<-as.data.frame(cbind(rep(m, length(tukey_result$test$coefficients)), tukey_result$test$coefficients, tukey_result$test$sigma, tukey_result$test$tstat, as.vector(tukey_result$test$pvalues)))
        tukey_lineage<-rbind(tukey_lineage, tukey_df)
	
	#Eigenvalue for each module with strain as fixed effect; perform Tukey test for pairwise comparisons among strains
	lm_strainType<-lm(get(m) ~ strain.type, data=final_data)
        tukey_result<-summary(glht(lm_strainType, mcp(strain.type="Tukey")))
        print(tukey_result)
        tukey_df<-as.data.frame(cbind(rep(m, length(tukey_result$test$coefficients)), tukey_result$test$coefficients, tukey_result$test$sigma, tukey_result$test$tstat, as.vector(tukey_result$test$pvalues)))
        tukey_strain<-rbind(tukey_strain, tukey_df)

	if(dataset=="all"){
		#Test for cell type effects
		lm_cellType<-lm(get(m)	~ cell.type, data=final_data)
		tukey_result<-summary(glht(lm_cellType, mcp(cell.type="Tukey")))
		print(tukey_result)
		tukey_df<-as.data.frame(cbind(rep(m, length(tukey_result$test$coefficients)), tukey_result$test$coefficients, tukey_result$test$sigma, tukey_result$test$tstat, as.vector(tukey_result$test$pvalues)))
		tukey_ct<-rbind(tukey_ct, tukey_df)

		#Test for lineage effects with cell type as random effect
		lmer_linCt<-lmer(get(m) ~ lin.type + (1|cell.type), data=final_data)
		tukey_result<-summary(glht(lmer_linCt, mcp(lin.type="Tukey")))
		print(tukey_result)
		tukey_df<-as.data.frame(cbind(rep(m, length(tukey_result$test$coefficients)), tukey_result$test$coefficients, tukey_result$test$sigma, tukey_result$test$tstat, as.vector(tukey_result$test$pvalues)))
		tukey_lineageXct<-rbind(tukey_lineageXct, tukey_df)

		#Test for strain effects with cell type as random effect
		lmer_strainCt<-lmer(get(m) ~ strain.type + (1|cell.type), data=final_data)
                tukey_result<-summary(glht(lmer_strainCt, mcp(strain.type="Tukey")))
                print(tukey_result)
                tukey_df<-as.data.frame(cbind(rep(m, length(tukey_result$test$coefficients)), tukey_result$test$coefficients, tukey_result$test$sigma, tukey_result$test$tstat, as.vector(tukey_result$test$pvalues)))
                tukey_strainXct<-rbind(tukey_strainXct, tukey_df)

	}
}
warnings()

colnames(tukey_lineage)<-c("module","estimate","st.err","z.val","p.val")
p.adj<-p.adjust(tukey_lineage$p.val, method="fdr")
tukey_final<-as.data.frame(cbind(tukey_lineage, p.adj))
write.csv(tukey_final, file=paste0("module_lmers.lineageWithTukey.",dataset,".csv"))

colnames(tukey_strain)<-c("module","estimate","st.err","z.val","p.val")
p.adj<-p.adjust(tukey_strain$p.val, method="fdr")
tukey_final<-as.data.frame(cbind(tukey_strain, p.adj))
write.csv(tukey_final, file=paste0("module_lmers.strainWithTukey.",dataset,".csv"))

if(dataset=="all"){
	colnames(tukey_ct)<-c("module","estimate","st.err","z.val","p.val")
	p.adj<-p.adjust(tukey_ct$p.val, method="fdr")
	tukey_final<-as.data.frame(cbind(tukey_ct, p.adj))
	write.csv(tukey_final, file=paste0("module_lmers.cellTypeWithTukey.",dataset,".csv"))
	
	colnames(tukey_lineageXct)<-c("module","estimate","st.err","z.val","p.val")
	p.adj<-p.adjust(tukey_lineageXct$p.val, method="fdr")
	tukey_final<-as.data.frame(cbind(tukey_lineageXct, p.adj))
	write.csv(tukey_final, file=paste0("module_lmers.lineageXcellTypeWithTukey.",dataset,".csv"))

	colnames(tukey_strainXct)<-c("module","estimate","st.err","z.val","p.val")
        p.adj<-p.adjust(tukey_strainXct$p.val, method="fdr")
        tukey_final<-as.data.frame(cbind(tukey_strainXct, p.adj))
        write.csv(tukey_final, file=paste0("module_lmers.strainXcellTypeWithTukey.",dataset,".csv"))

}
