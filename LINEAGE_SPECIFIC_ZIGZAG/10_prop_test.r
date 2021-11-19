#PURPOSE: Proportion test to check for significant difference in proportion lineage specific genes between early and alte

#Read in data
print("Reading in data...")
LZ<-read.table("lineage_specific_gene_numbers.LZ.txt",header=TRUE)
RS<-read.table("lineage_specific_gene_numbers.RS.txt",header=TRUE)
print(LZ)
print(RS)

for(i in 1:ncol(LZ)){
	print(colnames(LZ[i]))
	auto<-prop.test(c(LZ["auto_lineage_specific",i],RS["auto_lineage_specific",i]), c(LZ["auto_all_active",i],RS["auto_all_active",i]))
	print(auto)
	x<-prop.test(c(LZ["x_lineage_specific",i],RS["x_lineage_specific",i]), c(LZ["x_all_active",i],RS["x_all_active",i]))
	print(x)
}
q()
#Get total active genes for each group
print("Getting active genes...")
dom_active_LZ<-LZ[27]
dom_active_RS<-RS[27]
mus_active_LZ<-LZ[47]
mus_active_RS<-RS[47]
spr_active_LZ<-LZ[67]
spr_active_RS<-RS[67]
#Didn't actually print this for the mus-dom ancestral branch so need to back calculate using (# of lineage spec gene/proportion lineage specific)
musdom_total_lineageSpec_LZ<-as.numeric(LZ[119])
musdom_line_LZ<-unlist(strsplit(LZ[121]," "))
musdom_prop_LZ<-as.numeric(sub("\"","",musdom_line_LZ[length(musdom_line_LZ)]))
musdom_active_LZ<-round(musdom_total_lineageSpec_LZ / musdom_prop_LZ,0)
musdom_total_lineageSpec_RS<-as.numeric(RS[119])
musdom_line_RS<-unlist(strsplit(RS[121]," "))
musdom_prop_RS<-as.numeric(sub("\"","",musdom_line_RS[length(musdom_line_RS)]))
musdom_active_RS<-round(musdom_total_lineageSpec_RS / musdom_prop_RS,0)

print(paste("dom LZ active:",dom_active_LZ))
print(paste("dom RS active:",dom_active_RS))
print(paste("mus LZ active:",mus_active_LZ))
print(paste("mus RS active:",mus_active_RS))
print(paste("spr LZ active:",spr_active_LZ))
print(paste("spr RS active:",spr_active_RS))
print(paste("dom-mus ancestral branch LZ active:",musdom_active_LZ))
print(paste("dom-mus ancestral branch RS active:",musdom_active_RS))

#Get auto lineage specific genes
print("Getting autosomal lineage specific genes...")
dom_auto_LZ<-LZ[319]
dom_auto_RS<-RS[319]
mus_auto_LZ<-LZ[325]
mus_auto_RS<-RS[325]
spr_auto_LZ<-LZ[331]
spr_auto_RS<-RS[331]
musdom_auto_LZ<-LZ[337]
musdom_auto_RS<-RS[337]

print(paste("dom lineage specific LZ auto:",dom_auto_LZ))
print(paste("dom lineage specific RS auto:",dom_auto_RS))
print(paste("mus lineage specific LZ auto:",mus_auto_LZ))
print(paste("mus lineage specific RS auto:",mus_auto_RS))
print(paste("spr lineage specific LZ auto:",spr_auto_LZ))
print(paste("spr lineage specific RS auto:",spr_auto_RS))
print(paste("dom-mus ancestral branch lineage specific LZ auto:",musdom_auto_LZ))
print(paste("dom-mus ancestral branch lineage specific RS auto:",musdom_auto_RS))

#Get x lineage specific genes
print("Getting X chromosome lineage specific genes...")
dom_x_LZ<-LZ[343]
dom_x_RS<-RS[343]
mus_x_LZ<-LZ[349]
mus_x_RS<-RS[349]
spr_x_LZ<-LZ[355]
spr_x_RS<-RS[355]
musdom_x_LZ<-LZ[361]
musdom_x_RS<-RS[361]

print(paste("dom lineage specific LZ x:",dom_x_LZ))
print(paste("dom lineage specific RS x:",dom_x_RS))
print(paste("mus lineage specific LZ x:",mus_x_LZ))
print(paste("mus lineage specific RS x:",mus_x_RS))
print(paste("spr lineage specific LZ x:",spr_x_LZ))
print(paste("spr lineage specific RS x:",spr_x_RS))
print(paste("dom-mus ancestral branch lineage specific LZ x:",musdom_x_LZ))
print(paste("dom-mus ancestral branch lineage specific RS x:",musdom_x_RS))

#Proportion tests
print("Performing proportion tests...")
mygroups<-c("dom","mus","spr","musdom")
for(i in mygroups){
	print(i)
	auto_LZ<-as.numeric(get(paste0(i,"_auto_LZ")))
	#x_LZ<-as.numeric(get(paste0(i,"_x_LZ")))
	active_LZ<-as.numeric(get(paste0(i,"_active_LZ")))
	auto_RS<-as.numeric(get(paste0(i,"_auto_RS")))
        #x_RS<-as.numeric(get(paste0(i,"_x_RS")))
	active_RS<-as.numeric(get(paste0(i,"_active_RS")))
	print(prop.test(c(auto_LZ,auto_RS), c(active_LZ,active_RS)))
	#print(prop.test(c(x_LZ,x_RS), c(active_LZ,active_RS)))
	#auto<-prop.test(c(auto_LZ,auto_RS), c(active_LZ,active_RS))$p.value
	#x<-prop.test(c(x_LZ,x_RS), c(active_LZ,active_RS))$p.value
	#print(paste(i, "autos prop test p-value:", auto))
	#print(paste(i, "x prop test p-value:", x))
}

print("Done!")
