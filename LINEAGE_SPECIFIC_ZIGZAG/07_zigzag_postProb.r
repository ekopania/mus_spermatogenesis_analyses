#PURPOSE: Use zigzag to calculate posterior probabilities that each gene is active
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

#Get probability genes are active
#prob_active1<-read.table(paste0(args[1],"_zigzag_output1/",args[1],"_mcmc_output/",args[1],"_probability_active.tab"), header=T, row.names=1)
prob_active1<-read.table(paste0(args[1],"_zigzag_output1_oneActiveComp/",args[1],"_mcmc_output/",args[1],"_probability_active.tab"), header=T, row.names=1)
head(round(cbind(expression_data, prob_active1), digits=3))
#prob_active2<-read.table(paste0(args[1],"_zigzag_output2/",args[1],"_mcmc_output/",args[1],"_probability_active.tab"), header=T, row.names=1)
prob_active2<-read.table(paste0(args[1],"_zigzag_output2_oneActiveComp/",args[1],"_mcmc_output/",args[1],"_probability_active.tab"), header=T, row.names=1)
head(round(cbind(expression_data, prob_active2), digits=3))

#Plot prob active over expression density plot
#pdf(paste0(args[1],"_prob_active1.pdf"))
pdf(paste0(args[1],"_prob_active1.oneActiveComp.pdf"))
plot(density(log(expression_data[,1])), main=paste0(args[1],": Probability Gene Is Active"), xlab="log Expression", ylim=c(0,0.25), axes=FALSE)
axis(1)
points(log(rowMeans(expression_data)), 0.25*prob_active1$prob_active, col="red", cex=0.5)
for(i in seq(ncol(expression_data))) lines(density(log(expression_data[,i])))
dev.off()
#pdf(paste0(args[1],"_prob_active2.pdf"))
pdf(paste0(args[1],"_prob_active2.oneActiveComp.pdf"))
plot(density(log(expression_data[,1])), main=paste0(args[1],": Probability Gene Is Active"), xlab="log Expression", ylim=c(0,0.25), axes=FALSE)
axis(1)
points(log(rowMeans(expression_data)), 0.25*prob_active2$prob_active, col="red", cex=0.5)
for(i in seq(ncol(expression_data))) lines(density(log(expression_data[,i])))
dev.off()


#Post prediction (make sure model is adequate and consistent between two MCMC runs)
#postpred1<-read.table(paste0(args[1],"_zigzag_output1/",args[1],"_mcmc_output/",args[1],".post_predictive_output.log"), header=T, row.names=1)
postpred1<-read.table(paste0(args[1],"_zigzag_output1_oneActiveComp/",args[1],"_mcmc_output/",args[1],".post_predictive_output.log"), header=T, row.names=1)
#pdf(paste0(args[1],"_postpred_boxplot1.pdf"))
pdf(paste0(args[1],"_postpred_boxplot1.oneActiveComp.pdf"))
boxplot(postpred1, outline=FALSE)
abline(h=0, col="Red")
dev.off()
#postpred2<-read.table(paste0(args[1],"_zigzag_output2/",args[1],"_mcmc_output/",args[1],".post_predictive_output.log"), header=T, row.names=1)
postpred2<-read.table(paste0(args[1],"_zigzag_output2_oneActiveComp/",args[1],"_mcmc_output/",args[1],".post_predictive_output.log"), header=T, row.names=1)
#pdf(paste0(args[1],"_postpred_boxplot2.pdf"))
pdf(paste0(args[1],"_postpred_boxplot2.oneActiveComp.pdf"))
boxplot(postpred2, outline=FALSE)
abline(h=0, col="Red")
dev.off()

print("Done with 07_zigzag_postProb")
