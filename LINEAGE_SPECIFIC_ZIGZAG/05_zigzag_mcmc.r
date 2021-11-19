#PURPOSE: Run burnin and MCMC steps of zigzag
#Tutorial here: https://ammonthompson.github.io/zigzag/
#Preprint here (As of 25 MARCH 2020): https://www.biorxiv.org/content/10.1101/711630v2

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line argument. Enter the sample name for tpm and length files (output from 00_run_calculate_tpm.sh; example: WWWW_LZ)")
}

#Required packages/libraries
library(coda)
library(zigzag)

#Read in files
tpm_file<-paste0(args[1],".combined.tpm.txt")
len_file<-paste0(args[1],".length.txt")
print(tpm_file)
print(len_file)
expression_data<-read.table(tpm_file, header = T, row.names = 1)
head(round(expression_data, digits = 3))
gene_lengths<-read.table(len_file, header = F, row.names = 1)
head(round(gene_lengths, digits = 1))

#Set up zigzag object
#FOR LZ
#my_zigzag1<-zigzag$new(data=expression_data, gene_length=gene_lengths, output_directory=paste0(args[1],"_zigzag_output1"), num_active_components=2, threshold_a = c(1, 6)) #Num active componenets and threshold_a should be changed based on distribution of expression levels!
#my_zigzag2<-zigzag$new(data=expression_data, gene_length=gene_lengths, output_directory=paste0(args[1],"_zigzag_output2"), num_active_components=2, threshold_a = c(1, 6)) #Zigzag tutorial recommends doing everything at least twice and making sure parameters converge on same values for both runs
#FOR RS
#my_zigzag1<-zigzag$new(data=expression_data, gene_length=gene_lengths, output_directory=paste0(args[1],"_zigzag_output1"), num_active_components=2, threshold_a = c(-1, 4))
#my_zigzag2<-zigzag$new(data=expression_data, gene_length=gene_lengths, output_directory=paste0(args[1],"_zigzag_output2"), num_active_components=2, threshold_a = c(-1, 4))
my_zigzag1<-zigzag$new(data=expression_data, gene_length=gene_lengths, output_directory=paste0(args[1],"_zigzag_output1_oneActiveComp"), num_active_components=1)
my_zigzag2<-zigzag$new(data=expression_data, gene_length=gene_lengths, output_directory=paste0(args[1],"_zigzag_output2_oneActiveComp"), num_active_components=1)

#Run and plot burnin
#my_zigzag1$burnin(sample_frequency=50, ngen=100000, write_to_files=TRUE)
burn1<-read.table(paste0(args[1],"_zigzag_output1_oneActiveComp/output_burnin/output_model_parameters.log"), header=T, row.names=1)
pdf(paste0(args[1],"_burnin_plots1_oneActiveComp.pdf"))
for(i in seq(ncol(burn1))) plot(burn1[,i], type="l", main=colnames(burn1)[i])
dev.off()
my_zigzag2$burnin(sample_frequency=50, ngen=100000, write_to_files=TRUE)
burn2<-read.table(paste0(args[1],"_zigzag_output2_oneActiveComp/output_burnin/output_model_parameters.log"), header=T, row.names=1)
pdf(paste0(args[1],"_burnin_plots2_oneActiveComp.pdf"))
for(i in seq(ncol(burn2))) plot(burn2[,i], type="l", main=colnames(burn2)[i])
dev.off()



#Run and plot MCMC
my_zigzag1$mcmc(sample_frequency=50, ngen=100000, run_posterior_predictive=TRUE, mcmcprefix=args[1])
post_mcmc1<-read.table(paste0(args[1],"_zigzag_output1_oneActiveComp/",args[1],"_mcmc_output/",args[1],"_model_parameters.log"), header=T, row.names=1)
pdf(paste0(args[1],"_mcmc_plots1_oneActiveComp.pdf"))
for(i in seq(ncol(post_mcmc1))) plot(post_mcmc1[,i], type="l", main=colnames(post_mcmc1)[i])
dev.off()
write.table(effectiveSize(post_mcmc1),file=paste0(args[1],"_effectiveSize1_oneActiveComp.txt"),append=FALSE,quote=FALSE,row.names=TRUE,col.names=FALSE,sep="\t")
pdf(paste0(args[1],"_inactive_variance_hist1_oneActiveComp.pdf"))
hist(post_mcmc1$inactive_variance, breaks=100, xlim=c(0.01,5), main="")
dev.off()
my_zigzag2$mcmc(sample_frequency=50, ngen=100000, run_posterior_predictive=TRUE, mcmcprefix=args[1])
post_mcmc2<-read.table(paste0(args[1],"_zigzag_output2_oneActiveComp/",args[1],"_mcmc_output/",args[1],"_model_parameters.log"), header=T, row.names=1)
pdf(paste0(args[1],"_mcmc_plots2_oneActiveComp.pdf"))
for(i in seq(ncol(post_mcmc2))) plot(post_mcmc2[,i], type="l", main=colnames(post_mcmc2)[i])
dev.off()
write.table(effectiveSize(post_mcmc2),file=paste0(args[1],"_effectiveSize2.txt"),append=FALSE,quote=FALSE,row.names=TRUE,col.names=FALSE,sep="\t")
pdf(paste0(args[1],"_inactive_variance_hist2_oneActiveComp.pdf"))
hist(post_mcmc2$inactive_variance, breaks=100, xlim=c(0.01,5), main="")
dev.off()

print("Done with 05_zigzag_mcmc")
