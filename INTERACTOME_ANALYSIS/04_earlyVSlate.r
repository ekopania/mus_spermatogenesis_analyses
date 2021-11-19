#PURPOSE: Test if number of interactions is different between genes induced early and genes induced late

suppressPackageStartupMessages(library(ggplot2))

print("Reading in data...")

#Get interaction counts
LZinteractions<-read.table("interaction_counts.LZinduced.txt", header=TRUE)
RSinteractions<-read.table("interaction_counts.RSinduced.txt", header=TRUE)

LZ_wLabel<-as.data.frame(cbind(ct=rep("LZ", nrow(LZinteractions)), LZinteractions))
RS_wLabel<-as.data.frame(cbind(ct=rep("RS", nrow(RSinteractions)), RSinteractions))
final_df<-as.data.frame(rbind(LZ_wLabel, RS_wLabel))
print(head(final_df))

#Print medians
print("Median number of interactions:")
print(paste("Early, all:", median(LZ_wLabel$num_interactions, na.rm=TRUE)))
print(paste("Early, mid:", median(LZ_wLabel$num_midScore_interactions, na.rm=TRUE)))
print(paste("Early, high:", median(LZ_wLabel$num_highScore_interactions, na.rm=TRUE)))
print(paste("Late, all:", median(RS_wLabel$num_interactions, na.rm=TRUE)))
print(paste("Late, mid:", median(RS_wLabel$num_midScore_interactions, na.rm=TRUE)))
print(paste("Late, high:", median(RS_wLabel$num_highScore_interactions, na.rm=TRUE)))

#Wilcoxon rank sum tests to compare early vs late
print("Wilcoxon rank sum tests (all, mid score, high score):")
wilcox_all<-wilcox.test(LZ_wLabel$num_interactions, RS_wLabel$num_interactions)
wilcox_mid<-wilcox.test(LZ_wLabel$num_midScore_interactions, RS_wLabel$num_midScore_interactions)
wilcox_high<-wilcox.test(LZ_wLabel$num_highScore_interactions, RS_wLabel$num_highScore_interactions)
print(wilcox_all)
print(wilcox_mid)
print(wilcox_high)
p.values<-c(wilcox_all$p.value, wilcox_mid$p.value, wilcox_high$p.value)
print("FDR-corrected p-values")
p.adjust(p.values, method="fdr")
print("Bonferroni-corrected p-values")
p.adjust(p.values, method="bonferroni")

#Plot
print("Plotting...")
pdf("number_interactions.earlyVSlate.boxplots.pdf", onefile=TRUE)
#All interactions
p<-ggplot(final_df, aes(x=ct, y=log(as.numeric(as.character(num_interactions))), fill=ct)) + geom_boxplot()
p<-p + labs(title="Number of interactions in genes induced early vs induced late", x="Cell type", y="log(Number of interactions)")
p<-p + theme(axis.text.y = element_text(size=20))
p<-p + theme_minimal()
p<-p + scale_fill_manual(values=c("chocolate1","lightsteelblue1"))
print(p)
#Mid-score interactions
p<-ggplot(final_df, aes(x=ct, y=log(as.numeric(as.character(num_midScore_interactions))), fill=ct)) + geom_boxplot()
p<-p + labs(title="Number of mid-score interactions in genes induced early vs induced late", x="Cell type", y="log(Number of mid-score interactions)")
p<-p + theme(axis.text.y = element_text(size=20))
p<-p + theme_minimal()
p<-p + scale_fill_manual(values=c("chocolate1","lightsteelblue1"))
print(p)
#Mid-score interactions
p<-ggplot(final_df, aes(x=ct, y=log(as.numeric(as.character(num_highScore_interactions))), fill=ct)) + geom_boxplot()
p<-p + labs(title="Number of high-score interactions in genes induced early vs induced late", x="Cell type", y="log(Number of high-score interactions)")
p<-p + theme(axis.text.y = element_text(size=20))
p<-p + theme_minimal()
p<-p + scale_fill_manual(values=c("chocolate1","lightsteelblue1"))
print(p)

dev.off()

print("Done with early vs late")
