#PURPOSE: Run R package ZigZag (Kopp Lab) to determine if genes are active or inactive using a Bayesian approach
#Tutorial here: https://ammonthompson.github.io/zigzag/
#Preprint here (As of 25 MARCH 2020): https://www.biorxiv.org/content/10.1101/711630v2

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line argument. Enter the sample name for tpm and length files (output from 00_run_calculate_tpm.sh; example: WWWW_LZ)")
}
#print("DO NOT run this script with Rscript - plot won't work. Instead, turn X11 forwarding ON and copy and paste into R terminal")

#sample_group<-"WWWW_LZ" #Change this according to which genotype/cell type you want to look at!

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

#View density of log expression
#Make sure plot is bimodal (one large peak w/ left shoulder)
#Use to determine how many subgroups of active genes there should be and where cutoffs for active and inactive should be (this seems a bit arbitrary, see tutorial for visual)
pdf(paste0(args[1],".log_exp_distribution.pdf"))
plot(density(log(expression_data[,1])), main ="", xlab = "log Expression", ylim = c(0, 0.25))
for(i in seq(ncol(expression_data))) lines(density(log(expression_data[,i])))
#Use these if you want to color by sample:
#for(i in seq(ncol(expression_data))) lines(density(log(expression_data[,i])),col=i)
#legend("topright",legend=colnames(expression_data),lty=1,col=c(1:4))
dev.off()

print("Done with 03_zigzag")
