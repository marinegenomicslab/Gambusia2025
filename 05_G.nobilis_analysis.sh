####### G. nobilis last stages of filtering and beginning of analysis #######
####### Oct 30, 2024 #######

#Loading Libraries
{```{R}```
library('devtools')
library('pegas')
library('adegenet')
library('vcfR') 
library('dartR')
library('zvau')
library('geosphere')
library('stringr')
library('ggmap')
library('ggcompoplot')
library('vegan')
library('spdep')
library('adespatial')
library('igraph')
library('poppr') 
library('smatr')
library('radiator')
library('related')
library('ggcompoplot')
library('scales')
library('hierfstat')
library('VennDiagram')
library('sdmpredictors') 
library('ggord')
library('coin')
library('car')
library('reshape2')
}

#Add functions
{```{R}```
#Pulls the unique loci names from alleles selected from genind object names (e.g. dDocent_Contig_23960.001)
Loci_names <- function(NAMES, SEP="[.]", REMOVE=1){
COL <- length(strsplit(head(NAMES,n=1), SEP)[[1]])
TMP_U <- unique(matrix(unlist(strsplit(NAMES,SEP)),ncol=COL,byrow=T)[,1:(COL-REMOVE)])
if(is.matrix(TMP_U)){TMP_DF <- data.frame(TMP_U)
} else {TMP_DF <- data.frame(matrix(TMP_U, ncol=(COL-REMOVE), byrow=T))}
return(tidyr::unite(TMP_DF, "loci", 1:ncol(TMP_DF), sep=SEP))
}

} #notepad cleanup

#Removing loci for genind object from a vcf; VCF locNames need to have a "dDocent_Contig_#_Pos" format
VCF_remove <- function(genind.vcf, loci.list){
tmp.loci <- matrix(unlist(strsplit(locNames(genind.vcf),"_")), ncol=4, byrow=T)
if(sum(is.na(tail(tmp.loci, n=1)))){print("locNames(genind) in incorrect format"); break}
tmp.loci <- data.frame(Locus=paste(tmp.loci[,1], tmp.loci[,2], tmp.loci[,3],sep="_"), Pos=tmp.loci[,4])

tmp <- vector()
for(i in loci.list){
ROW <- which(tmp.loci$Locus == i)
tmp <- append(tmp,paste(tmp.loci$Locus, tmp.loci$Pos, sep="_")[ROW])
}

set.loc <- subset(locNames(genind.vcf), !locNames(genind.vcf) %in% tmp)
return(genind.vcf[ ,loc=set.loc])
}

}  #notepad cleanup

#Function for assessing the Kmeans of a group
Kmeans_heatmap <- function(gen.ind, n.clust, iter=100, plot=TRUE, out_file, verbose=TRUE, diagonal=TRUE){
#gen.ind: A genind object
#n.clust: The K-means value to test (number)
#iter: number of iterations to use (number)
#plot: should the function provide you with a plot? (T|F)
#out_file: The name of the heatmap plot (string)
#verbose: Should the function provide progress updates (T|F)
#diaganol: Should the matrix include diagonal values (T|F)

#Making multiple versions of the K value
if(verbose==TRUE){print("Starting K selections")}
for(i in 1:iter){
if(verbose==TRUE){if(i/10 == round(i/10,0)){print(paste("Starting Round",i))}}
tmp.grp <- find.clusters(gen.ind, max.n.clust=40, n.pca=400, n.clust=n.clust, method="kmeans")
if(!exists("Kmean.df")){Kmean.df <- data.frame(T1=tmp.grp$grp)
} else {Kmean.df[,paste("T",i,sep="")] <- tmp.grp$grp[match(rownames(Kmean.df),names(tmp.grp$grp))]}
}

Kmean.map <- data.frame(matrix(nrow=nrow(Kmean.df), ncol=nrow(Kmean.df)))
colnames(Kmean.map) <- rownames(Kmean.map) <- rownames(Kmean.df)

Kmean.pairs <- t(combn(colnames(Kmean.map),2))

print("Starting pairwise comparisons")
for(i in 1:dim(Kmean.pairs)[1]){
if(verbose==TRUE){if(i/1000 == round(i/1000,0)){print(paste("Starting Round",i,"of",dim(Kmean.pairs)[1]))}}
INDV1 <- Kmean.pairs[i,1] 
INDV2 <- Kmean.pairs[i,2]
tmp.df <- Kmean.df[Kmean.pairs[i,],]
Kmean.map[INDV1,INDV2] <- Kmean.map[INDV2,INDV1] <- sum(apply(Kmean.df[c(INDV1,INDV2),], 2, function(x) sum(x[1] == x[2])))
}

#Filling in diaganol
if(diagonal==TRUE){for(i in colnames(Kmean.map)){Kmean.map[i,i] <- iter}}

#Making percent
Kmean.map <- Kmean.map/iter

#Getting generic view of the heatmap
if(plot==TRUE){pheatmap::pheatmap(Kmean.map, treeheight_row = 0, treeheight_col = 0)}

png(paste(out_file),res=200, width=2000, height=2000)
heatmap(as.matrix(Kmean.map))
dev.off()

return(Kmean.map)
}

}}}}}}}} #notepad cleanup
} #notepad cleanup

#Importing Data
{```{R}```
strata <- read.csv(file = "../../filtering/Gnobilis_May2024/SNP.TRS.F06.2024-06-01/indv.csv", header = TRUE)
#Adding site to strata
tmp.v <- matrix(unlist(strsplit(as.character(as.matrix(strata$POP)), "-")), byrow=T, ncol=2)[,2]
strata$Site <- as.factor(tmp.v)

#Adjusting sites
tmp.v <- as.character(as.matrix(strata$Site))
tmp.v[tmp.v %in% c("DY", "FLS")] <- "KS"
tmp.v[tmp.v == "CHS"] <- "BAL"
strata$Site2 <- as.factor(tmp.v)

#Adding nobilis grouping inforrmation
tmp.v <- as.character(as.matrix(strata$Site))
tmp.v[tmp.v %in% c("CHS", "ES", "PHL", "BAL")] <- "Gan-WT"
tmp.v[tmp.v %in% c("EU", "FLS", "HEAD", "DY")] <- "Gan-DY"
tmp.v[tmp.v %in% c("Sink27", "Sink31", "Sink37", "Sink7", "BC", "BLNWR")] <- "Gan-NM"
strata$POP2 <- as.factor(tmp.v)
strata$POP <- factor(strata$POP, levels=c("DY","EU","FLS","HEAD","BAL","CHS","ES","PHL","BC","BLNWR","Sink27","Sink31","Sink37","Sink7", "LCCSP"))
head(strata)

vcf<-read.vcfR(file="SNP.TRS.F09.vcf")
gen.vcf<-vcfR2genind(vcf)
strata(gen.vcf) <- strata[match(indNames(gen.vcf),strata$INDV),]
head(gen.vcf@strata)
rm(vcf)

gen <- read.genepop(file = "SNP.TRS.F09.gen", ncode=3L, quiet = FALSE)
strata(gen) <- strata[match(indNames(gen),strata$INDV),]
head(gen@strata)

#Adding related ID to strata
gen.indv <- substr(gen@strata$INDV, 1, 20)
gen.indv <- data.frame(cbind(as.character(gen@strata$INDV), gen.indv))
names(gen.indv) <- c("gen.ID", "relate.ID")

gen@strata$relate.ID <- gen.indv[match(indNames(gen),gen.indv$gen.ID),"relate.ID"]
head(gen@strata)

save(gen, file="gen.gz", compress=T)
save(gen.vcf, file="gen.vcf.gz", compress=T)
#load("gen.gz")
#load("gen.vcf.gz")
} #notepad cleanup
#Setting color schemes
{```{R}```
#         "Gan-DY" "Gan-NM"  "Gan-WT"
col.All<-c("red4","grey50","dodgerblue")

c2 <- c("mediumblue", "red4")
c3 <- c("dodgerblue", "red4", "mediumblue")
c4 <- c("lightslateblue", "mediumblue", "red4", "dodgerblue")
c5 <- c("dodgerblue", "mediumblue", "royalblue1", "lightslateblue", "red4")
c6 <- c("mediumblue", "lightslateblue", "royalblue1", "red4","cornflowerblue", "dodgerblue")
} #notepad cleanup
#Removing monomorphics
{```{R}```
mono.loc <- names(which(gen@loc.n.all == 1))
} #notepad cleanup

#Removing duplicates
{```{R}```
dup.list <- names(which(table(gen.vcf@strata$Sample) > 1))

rm.list <- list()
for(i in 1:length(dup.list)){
j <- as.character(dup.list[i])

tmp <- gen.vcf@strata[gen.vcf@strata$Sample %in% j, ]
tmp$Diff <- abs(as.numeric(as.character(tmp[ , "E.HOM."])) - as.numeric(as.character(tmp[ , "O.HOM."])))

Miss.tab <- table(as.numeric(as.character(tmp$N_MISS)))
Depth.tab <- table(as.numeric(as.character(tmp$MEAN_DEPTH)))

MIN.Miss <- min(as.numeric(names(Miss.tab)))
MAX.Depth <- max(as.numeric(names(Depth.tab)))

if(Miss.tab[names(Miss.tab)==MIN.Miss]==1 && Depth.tab[names(Depth.tab)==MAX.Depth]==1){keep <- as.character(tmp[tmp$N_MISS==MIN.Miss,1])
} else {keep <- as.character(tmp[tmp$Diff == min(tmp$Diff),1])}
print(keep)

rm.list[[i]] <- as.character(tmp[grep(keep, tmp$INDV, invert=T),1])
rm(keep, tmp)
} #notepad cleanup

#Gno-LCCSP_19759-Gan_Lib3_I2 is a misidentification and therefore removed
dup.rm <- c(unlist(rm.list),"Gaa-LCCSP_MGL-19774_Lib1_I11")
set.indv <- indNames(gen)[!indNames(gen) %in% dup.rm]
set.loc <- locNames(gen)[which(!locNames(gen) %in% mono.loc)]
gen2 <- gen[set.indv, loc=set.loc, drop=T]
gen2.vcf <- VCF_remove(gen.vcf[set.indv, , drop=T], mono.loc)

gen2@strata$Site2 <- factor(gen2@strata$Site2, levels=c("HEAD","KS","EU","BAL","ES","PHL","BC","BLNWR","Sink7","Sink27","Sink31","Sink37"))
gen2.vcf@strata$Site2 <- factor(gen2.vcf@strata$Site2, levels=c("HEAD","KS","EU","BAL","ES","PHL","BC","BLNWR","Sink7","Sink27","Sink31","Sink37"))

save(gen2, file="gen2.gz", compress=T)
save(gen2.vcf, file="gen2.vcf.gz", compress=T)
#load("gen2.gz")
#load("gen2.vcf.gz")
}} #notepad cleanup
#Library Bias (Outflank)
{```{R}```
setPop(gen2.vcf) <- ~Lib
tmp <- gl.outflank(gen2.vcf, qthreshold = 0.1)
length(which(tmp$outflank$results[15]=="TRUE"))
} #notepad cleanup

#Library Bias (Bayescan)
{```{R}```
#Exporting data for Bayescan
setPop(gen2) <- ~Lib
writeGenPop(gen2, "SNP.TRS.Lib_all.gen", "All Spinner data without dups by Library")
} #notepad cleanup

#Converting to BS format
{```{bash}```
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile SNP.TRS.Lib_all.gen -inputformat GENEPOP -outputfile SNP.TRS.Lib_all.BS -outputformat BAYESCAN -spid /home/afields/bin/genepop_to_BS.spid
} #notepad cleanup

#Running Bayescan
{```{bash}```
mkdir Bayescan
bayescan SNP.TRS.Lib_all.BS -od ./Bayescan -all-trace -threads 10 -thin 100 -nbp 30 -pr_odds 100

#Analyzing results
head -n2 ../SNP.TRS.Lib_all.gen | tail -n 1 | sed 's/, /\n/g' > contigs
echo "Contigs" | cat - contigs > tmp; mv tmp contigs
paste contigs SNP.TRS.Lib_al_fst.txt | awk 'NR==1{print $0}NR!=1{$2=""; print}' > fst_lib.txt

sort -k4,4n fst_lib.txt | less
awk 'NR==1{next;} $4 < 0.05{print $0}' fst_lib.txt | wc -l
awk 'NR==1{next;} $4 < 0.05{print $0}' fst_lib.txt | cut -f1 -d" "> Bayescan_Outliers_Lib.list
} #notepad cleanup

#Checking convergance
#http://evomics.org/wp-content/uploads/2016/01/BayeScan_BayeScEnv_exercises.pdf
{```{R}```
library(coda)
#Load data
chain <- read.table("SNP.TRS.Lib_al.sel",header=TRUE)
chain <- chain[-c(1)]
chain <- mcmc(chain,thin=10)
#Plot
plot(chain)
#Summary
summary(chain)

effectiveSize(chain)

#Checking Convergance
#Geweke’s convergence diagnostic (values should be between -1.96 and 1.96 for alpha = 0.05)
geweke.diag(chain, frac1=0.1, frac2=0.5)
#Heidelberg and Welch’s convergence diagnostic
heidel.diag(chain, eps=0.1, pvalue=0.05)
} #notepad cleanup

#Removing Library effects
{```{R}```
tmp.out <- as.character(as.matrix(read.csv("Bayescan/Bayescan_Outliers_Lib.list", head=F)))
OF.Gaa <- read.table("../May_2024/Gaa_Outflank_Lib_outlier.list", head=F)[,1] 
OF.Gag <- read.table("../May_2024/Gag_Outflank_Lib_outlier.list", head=F)[,1] 
OF.Gan <- read.table("../May_2024/Gan-DY_Outflank_Lib_outlier.list", head=F)[,1] 
BS.Gan <- read.table("../May_2024/Bayescan/Bayescan_Library_DY.list", head=F)[,1] 
gen.out <- c(tmp.out, OF.Gaa, OF.Gag, OF.Gan, BS.Gan)

set.loc <- locNames(gen2)[!locNames(gen2) %in% gen.out]
gen3 <- gen2[, loc=set.loc, drop=T]
gen3.vcf <- VCF_remove(gen2.vcf, gen.out)

save(gen3, file="gen3.gz", compress=T)
save(gen3.vcf, file="gen3.vcf.gz", compress=T)
#load("gen3.gz")
#load("gen3.vcf.gz")
} #notepad cleanup
#Relatedness within site
{```{R}```
setPop(gen3) <- ~Site2
gen3@pop <- factor(gen3@pop, levels=c("HEAD","KS","EU","BAL","ES","PHL","BC","BLNWR","Sink27","Sink31","Sink37","Sink7"))
gen.list <- seppop(gen3)
par(mfrow=c(3,4))
for(i in names(gen.list)){
tmp.gen <- gen.list[[i]]
if(nInd(tmp.gen) < 10){next}
gl <- gi2gl(tmp.gen)
test.cov <- gl2related(gl, save=F)
tmp.rel <- coancestry(test.cov, wang =1)$relatedness
hist(tmp.rel$wang, breaks=50, col="red4", main=i)
assign(paste(i,"rel",sep="_"), tmp.rel)
}

MIN <- -0.6
MAX <- 0.6

tmp.break <- seq(MIN,MAX,(MAX-MIN)/100)

png("Site_specific_relatedness.png", res=200, width=8000, height=6000)
par(mfrow=c(3,3))
for(i in names(gen.list)){
tmp.gen <- gen.list[[i]]
if(nInd(tmp.gen) < 10){next}
tmp.rel <- get(paste(i,"rel",sep="_"))
hist(tmp.rel$wang, breaks=tmp.break, col="red4", main=i, xlim=c(MIN,MAX))
}
dev.off()
}
}
} #notepad cleanup
#Exploring KS
{```{R}```
par(mfrow=c(1,1))
boxplot(as.numeric(as.matrix(gen2@strata$Fis[gen2@strata$Site2 == "KS"])))

rm.list <- gen2@strata$INDV[as.numeric(as.matrix(gen2@strata$Fis)) < 0.4 & gen2@strata$Site2 == "KS"]
set.indv <- indNames(gen.list[["KS"]])[!indNames(gen.list[["KS"]]) %in% rm.list]

tmp.gen <- gen.list[["KS"]][set.indv, ,drop=T]
gl <- gi2gl(tmp.gen)
test.cov <- gl2related(gl, save=F)
tmp.rel <- coancestry(test.cov, wang =1)$relatedness
hist(tmp.rel$wang, breaks=50, col="red4", main="DY_edit")

tmp.rel[which(tmp.rel$wang > 0.3),1:6]
gen2@strata[gen2@strata$relate.ID %in% c("Gno-DY_DY11-24_Lib3_", "Gno-DY_DY27-24_Lib3_"),]
rel.rm.list <- c(as.character(as.matrix(rm.list)),"Gno-DY_DY11-24_Lib3_I7")
} #notepad cleanup
#Exploring BAL
{```{R}```
boxplot(as.numeric(as.matrix(gen2@strata$Fis[gen2@strata$Site2 == "BAL"])))

rm.list <- gen2@strata$INDV[as.numeric(as.matrix(gen2@strata$Fis)) < 0.2 & gen2@strata$Site2 == "BAL"]
set.indv <- indNames(gen.list[["BAL"]])[!indNames(gen.list[["BAL"]]) %in% rm.list]
tmp.gen <- gen.list[["BAL"]][set.indv, ,drop=T]
gl <- gi2gl(tmp.gen)
test.cov <- gl2related(gl, save=F)
tmp.rel <- coancestry(test.cov, wang =1)$relatedness
par(mfrow=c(2,1))
hist(BAL_rel$wang, breaks=50, col="red4", main="BAL")
hist(tmp.rel$wang, breaks=50, col="red4", main="BAL_edit")

par(mfrow=c(3,2))
plot(density(tmp.rel$wang))
for(i in rm.list){
set.indv <- indNames(gen.list[["BAL"]])[!indNames(gen.list[["BAL"]]) %in% rm.list]
set.indv <- c(set.indv, i)
tmp.gen <- gen.list[["BAL"]][set.indv, ,drop=T]
gl <- gi2gl(tmp.gen)
test.cov <- gl2related(gl, save=F)
tmp.rel <- coancestry(test.cov, wang =1)$relatedness
plot(density(tmp.rel$wang), main=i)
}

rel.rm.list <- c(rel.rm.list, as.character(as.matrix(rm.list)))
} #notepad cleanup
#Exploring PHL
{```{R}```
boxplot(as.numeric(as.matrix(gen2@strata$Fis[gen2@strata$Site2 == "PHL"])))

rm.list <- gen2@strata$INDV[as.numeric(as.matrix(gen2@strata$Fis)) < 0 & gen2@strata$Site2 == "PHL"]
set.indv <- indNames(gen.list[["PHL"]])[!indNames(gen.list[["PHL"]]) %in% rm.list]

tmp.gen <- gen.list[["PHL"]][set.indv, ,drop=T]
gl <- gi2gl(tmp.gen)
test.cov <- gl2related(gl, save=F)
tmp.rel <- coancestry(test.cov, wang =1)$relatedness
par(mfrow=c(2,1))
hist(PHL_rel$wang, breaks=50, col="red4", main="PHL")
hist(tmp.rel$wang, breaks=50, col="red4", main="PHL_edit")

tmp.rel[tmp.rel$wang > 0.3,1:6]
gen3@strata[gen3@strata$relate.ID %in% c("Gan-PHL_PSGn_11_Lib2","Gan-PHL_PSGn_1_Lib2_"),]
gen3@strata[gen3@strata$relate.ID %in% c("Gan-PHL_PSGn_12_Lib2","Gan-PHL_PSGn_7_Lib2_"),]

rel.rm.list <- c(rel.rm.list, as.character(as.matrix(rm.list)))
} #notepad cleanup
#Exploring Sink31
{```{R}```
par(mfrow=c(1,1))
boxplot(as.numeric(as.matrix(gen2@strata$Fis[gen2@strata$Site2 == "Sink31"])))

tmp.gen <- gen.list[["Sink31"]][set.indv, ,drop=T]
gl <- gi2gl(tmp.gen)
test.cov <- gl2related(gl, save=F)
tmp.rel <- coancestry(test.cov, wang =1)$relatedness
hist(Sink31_rel$wang, breaks=50, col="red4", main="Sink31")

Sink31_rel[Sink31_rel$wang > 0.4,1:6]
gen3@strata[gen3@strata$relate.ID %in% c("Gan-Sink31_94695-A_L","Gan-Sink31_94695-B_L", "Gan-Sink31_94695-C_L"),]
} #notepad cleanup
#Exploring EU
{```{R}```
boxplot(as.numeric(as.matrix(gen2@strata$Fis[gen2@strata$Site2 == "EU"])))

rm.list <- gen2@strata$INDV[as.numeric(as.matrix(gen2@strata$Fis)) < 0.65 & gen2@strata$Site2 == "EU"]
set.indv <- indNames(gen.list[["EU"]])[!indNames(gen.list[["EU"]]) %in% rm.list]

tmp.gen <- gen.list[["EU"]][set.indv, ,drop=T]
gl <- gi2gl(tmp.gen)
test.cov <- gl2related(gl, save=F)
tmp.rel <- coancestry(test.cov, wang =1)$relatedness
par(mfrow=c(2,1))
hist(EU_rel$wang, breaks=50, col="red4", main="EU", xlim = c(-0.25, 0.2))
abline(v=mean(EU_rel$wang), col="black")
abline(v=median(EU_rel$wang), col="black", lty=2)
hist(tmp.rel$wang, breaks=50, col="red4", main="EU_edit", xlim = c(-0.25, 0.2))
abline(v=mean(tmp.rel$wang), col="black")
abline(v=median(tmp.rel$wang), col="black", lty=2)

rel.rm.list <- c(rel.rm.list, as.character(as.matrix(rm.list)))
} #notepad cleanup
#Exploring BC
{```{R}```
boxplot(as.numeric(as.matrix(gen2@strata$Fis[gen2@strata$Site2 == "BC"])))

rm.list <- gen2@strata$INDV[as.numeric(as.matrix(gen2@strata$Fis)) < 0.35 & gen2@strata$Site2 == "BC"]
set.indv <- indNames(gen.list[["BC"]])[!indNames(gen.list[["BC"]]) %in% rm.list]

tmp.gen <- gen.list[["BC"]][set.indv, ,drop=T]
gl <- gi2gl(tmp.gen)
test.cov <- gl2related(gl, save=F)
tmp.rel <- coancestry(test.cov, wang =1)$relatedness
par(mfrow=c(2,1))
hist(BC_rel$wang, breaks=50, col="red4", main="BC", xlim=c(-0.2, 0.15))
abline(v=mean(BC_rel$wang), col="black")
abline(v=median(BC_rel$wang), col="black", lty=2)
hist(tmp.rel$wang, breaks=50, col="red4", main="BC_edit", xlim=c(-0.2, 0.15))
abline(v=mean(tmp.rel$wang), col="black")
abline(v=median(tmp.rel$wang), col="black", lty=2)

tmp.tab <- table(c(head(BC_rel[order(BC_rel$wang), 2],n=25), head(BC_rel[order(BC_rel$wang), 3],n=25)))
tail(tmp.tab[order(tmp.tab)])

tmp.tab <- table(c(head(tmp.rel[order(tmp.rel$wang), 2],n=25), head(tmp.rel[order(tmp.rel$wang), 3],n=25)))
tail(tmp.tab[order(tmp.tab)])

set.indv <- indNames(gen.list[["BC"]])[!indNames(gen.list[["BC"]]) %in% c("Gno-BC_BL49-24_Lib3_I10")]
tmp.gen2 <- gen.list[["BC"]][set.indv, ,drop=T]
gl <- gi2gl(tmp.gen2)
test.cov <- gl2related(gl, save=F)
tmp.rel2 <- coancestry(test.cov, wang =1)$relatedness
par(mfrow=c(3,1))
hist(BC_rel$wang, breaks=50, col="red4", main="BC", xlim=c(-0.2, 0.15))
abline(v=mean(BC_rel$wang), col="black")
abline(v=median(BC_rel$wang), col="black", lty=2)
hist(tmp.rel$wang, breaks=50, col="red4", main="BC_edit", xlim=c(-0.2, 0.15))
abline(v=mean(tmp.rel$wang), col="black")
abline(v=median(tmp.rel$wang), col="black", lty=2)
hist(tmp.rel2$wang, breaks=50, col="red4", main="BC_edit2", xlim=c(-0.2, 0.15))
abline(v=mean(tmp.rel2$wang), col="black")
abline(v=median(tmp.rel2$wang), col="black", lty=2)

tmp.tab <- table(c(head(tmp.rel2[order(tmp.rel2$wang), 2],n=25), head(tmp.rel2[order(tmp.rel2$wang), 3],n=25)))

tmp.tab <- table(c(tail(BC_rel[order(BC_rel$wang), 2],n=25), tail(BC_rel[order(BC_rel$wang), 3],n=25)))

G2.map <- Kmeans_heatmap(gen.list[["BC"]], 2, iter=100, plot=TRUE, "BC_G2_Kmeans_heatmap.png", verbose=TRUE, diagonal=FALSE)
G3.map <- Kmeans_heatmap(gen.list[["BC"]], 3, iter=100, plot=TRUE, "BC_G3_Kmeans_heatmap.png", verbose=TRUE, diagonal=FALSE)

rel.rm.list <- c(rel.rm.list, as.character(as.matrix(rm.list)))
} #notepad cleanup
#Exploring BLNWR
{```{R}```
par(mfrow=c(1,1))
boxplot(as.numeric(as.matrix(gen2@strata$Fis[gen2@strata$Site2 == "BLNWR"])))

rm.list <- gen2@strata$INDV[as.numeric(as.matrix(gen2@strata$Fis)) < 0.5 & gen2@strata$Site2 == "BLNWR"]
set.indv <- indNames(gen.list[["BLNWR"]])[!indNames(gen.list[["BLNWR"]]) %in% rm.list]

tmp.gen <- gen.list[["BLNWR"]][set.indv, ,drop=T]
gl <- gi2gl(tmp.gen)
test.cov <- gl2related(gl, save=F)
tmp.rel <- coancestry(test.cov, wang =1)$relatedness
par(mfrow=c(2,1))
hist(BLNWR_rel$wang, breaks=50, col="red4", main="BLNWR")
hist(tmp.rel$wang, breaks=50, col="red4", main="BLNWR_edit")

rel.rm.list <- c(rel.rm.list, as.character(as.matrix(rm.list)))
} #notepad cleanup

#Listing all related individuals for removal
{```{R}```
related.indv <- c("Gno-DY_DY11-24_Lib3_I7", "Gan-PHL_PSGn_1_Lib2_I10", "Gan-Sink31_94695-A_Lib2_I8", "Gan-Sink31_94695-B_Lib2_I10")
rel.rm.list <- c("Gno-DY_DY12-24_Lib3_I4", "Gno-DY_DY16-24_Lib3_I10", "Gno-DY_DY24-24_Lib3_I2", "Gno-DY_DY11-24_Lib3_I7", "Gno-BAL_B09-24_Lib3_I2","Gno-BAL_B11-24_Lib3_I10", "Gno-BAL_B18-24_Lib3_I4","Gan-CHS_MGL-16377_Lib1_I2", 
"Gan-CHS_MGL-16387_Lib1_I11", "Gan-PHL_PSGn_9_Lib2_I8", "Gan-EU_MGL-16445_Lib1_I11","Gno-BC_BL28-24_Lib3_I10", "Gno-BLNWR_BL43-24_Lib3_I2")
} #notepad cleanup

#Removing loci with low Fis
{```{R}```
set.ind <- indNames(gen3)[!indNames(gen3) %in% rel.rm.list]
gen4 <- gen3[set.ind, , drop=T]
gen4.vcf <- gen3.vcf[set.ind, , drop=T]

save(gen4, file="gen4.gz", compress=T)
save(gen4.vcf, file="gen4.vcf.gz", compress=T)
#load("gen4.gz")
#load("gen4.vcf.gz")
} #notepad cleanup
#Fst Outliers (Outflank)
{```{R}```
setPop(gen4.vcf) <- ~Site2
tmp <- gl.outflank(gen4.vcf, qthreshold = 0.1)
length(which(tmp$outflank$results[15]=="TRUE" & tmp$outflank$results[11]=="goodH"))

setPop(gen4.vcf) <- ~POP2
gen.list <- seppop(gen4.vcf, drop=T)
setPop(gen.list[["Gan-NM"]]) <- ~Site2
tmp <- gl.outflank(gen.list[["Gan-NM"]], qthreshold = 0.1)
length(which(tmp$outflank$results[15]=="TRUE" & tmp$outflank$results[11]=="goodH"))

setPop(gen.list[["Gan-WT"]]) <- ~Site2
tmp <- gl.outflank(gen.list[["Gan-WT"]], qthreshold = 0.1)
length(which(tmp$outflank$results[15]=="TRUE" & tmp$outflank$results[11]=="goodH"))

setPop(gen.list[["Gan-DY"]]) <- ~Site2
tmp <- gl.outflank(gen.list[["Gan-DY"]], qthreshold = 0.1)
length(which(tmp$outflank$results[15]=="TRUE" & tmp$outflank$results[11]=="goodH"))

outliers <- tmp$outflank$results[which(tmp$outflank$results[15]=="TRUE" & tmp$outflank$results[11]=="goodH"),1]
out.list <- matrix(unlist(strsplit(as.vector(outliers),"[.]")),ncol=2,byrow=T)[,1]
tmp <- matrix(unlist(strsplit(as.vector(outliers),"[_]")),ncol=4,byrow=T)
gen.out.list<-unique(paste(tmp[,1],tmp[,2],tmp[,3],sep="_"))
length(gen.out.list) 

write.table(gen.out.list, "Gan-DY_Outflank_outlier.list", quote=F, col.names=F, row.names=F)
} #notepad cleanup

#Fst Outliers (Bayescan)
{```{R}```
#Exporting data for Bayescan
setPop(gen4) <- ~Site2
writeGenPop(gen4, "SNP.TRS.Site2_all.gen", "All G. nobilis data without dups by Site2")

setPop(gen4) <- ~POP2
gen.list <- seppop(gen4, drop=T)

setPop(gen.list[["Gan-NM"]]) <- ~Site2
writeGenPop(gen.list[["Gan-NM"]], "Gan-NM.Site2_all.gen", "New Mexico G. nobilis data without dups by Site2")

setPop(gen.list[["Gan-WT"]]) <- ~Site2
writeGenPop(gen.list[["Gan-WT"]], "Gan-WT.Site2_all.gen", "West Texas G. nobilis data without dups by Site2")

setPop(gen.list[["Gan-DY"]]) <- ~Site2
writeGenPop(gen.list[["Gan-DY"]], "Gan-DY.Site2_all.gen", "Diamond Y G. nobilis data without dups by Site2")
} #notepad cleanup

#Converting to BS format
{```{bash}```
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile SNP.TRS.Site2_all.gen -inputformat GENEPOP -outputfile SNP.TRS.Site2_all.BS -outputformat BAYESCAN -spid /home/afields/bin/genepop_to_BS.spid
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gan-NM.Site2_all.gen -inputformat GENEPOP -outputfile Gan-NM.Site2_all.BS -outputformat BAYESCAN -spid /home/afields/bin/genepop_to_BS.spid
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gan-WT.Site2_all.gen -inputformat GENEPOP -outputfile Gan-WT.Site2_all.BS -outputformat BAYESCAN -spid /home/afields/bin/genepop_to_BS.spid
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gan-DY.Site2_all.gen -inputformat GENEPOP -outputfile Gan-DY.Site2_all.BS -outputformat BAYESCAN -spid /home/afields/bin/genepop_to_BS.spid
} #notepad cleanup

#Running Bayescan
{```{bash}```
mkdir Bayescan
bayescan SNP.TRS.Site2_all.BS -od ./Bayescan -all-trace -threads 15 -thin 100 -nbp 30 -pr_odds 100
bayescan Gan-NM.Site2_all.BS -od ./Bayescan -all-trace -threads 15 -thin 100 -nbp 30 -pr_odds 100
bayescan Gan-WT.Site2_all.BS -od ./Bayescan -all-trace -threads 15 -thin 100 -nbp 30 -pr_odds 100
bayescan Gan-DY.Site2_all.BS -od ./Bayescan -all-trace -threads 15 -thin 100 -nbp 30 -pr_odds 100
} #notepad cleanup

#Analyzing results
{```{bash}```
#All sites
head -n2 ../SNP.TRS.Site2_all.gen | tail -n 1 | sed 's/, /\n/g' > contigs
echo "Contigs" | cat - contigs > tmp; mv tmp contigs
paste contigs SNP.TRS.Site2_al_fst.txt | awk 'NR==1{print $0}NR!=1{$2=""; print}' > fst_out.txt

awk 'NR==1{next;} $4 < 0.05{print $0}' fst_out.txt | wc -l
awk 'NR==1{next;} $4 < 0.05 && $5 < 0{print $0}' fst_out.txt | cut -f1 -d" "> Bayescan_Outliers_Site2_neg.list
awk 'NR==1{next;} $4 < 0.05 && $5 > 0{print $0}' fst_out.txt | cut -f1 -d" "> Bayescan_Outliers_Site2_pos.list

#NM Bayescan
head -n2 ../Gan-NM.Site2_all.gen | tail -n 1 | sed 's/, /\n/g' > contigs
echo "Contigs" | cat - contigs > tmp; mv tmp contigs
paste contigs Gan-NM.Site2_al_fst.txt | awk 'NR==1{print $0}NR!=1{$2=""; print}' > fst_out.txt

awk 'NR==1{next;} $4 < 0.05{print $0}' fst_out.txt | wc -l
awk 'NR==1{next;} $4 < 0.05{print $0}' fst_out.txt | less
awk 'NR==1{next;} $4 < 0.05 && $5 < 0{print $0}' fst_out.txt | cut -f1 -d" "> Bayescan_Outliers_NM_neg.list
awk 'NR==1{next;} $4 < 0.05 && $5 > 0{print $0}' fst_out.txt | cut -f1 -d" "> Bayescan_Outliers_NM_pos.list

#WT Bayescan
head -n2 ../Gan-WT.Site2_all.gen | tail -n 1 | sed 's/, /\n/g' > contigs
echo "Contigs" | cat - contigs > tmp; mv tmp contigs
paste contigs Gan-WT.Site2_al_fst.txt | awk 'NR==1{print $0}NR!=1{$2=""; print}' > fst_out.txt

awk 'NR==1{next;} $4 < 0.05{print $0}' fst_out.txt | wc -l

#DY Bayescan
head -n2 ../Gan-DY.Site2_all.gen | tail -n 1 | sed 's/, /\n/g' > contigs
echo "Contigs" | cat - contigs > tmp; mv tmp contigs
paste contigs Gan-DY.Site2_al_fst.txt | awk 'NR==1{print $0}NR!=1{$2=""; print}' > fst_out.txt

awk 'NR==1{next;} $4 < 0.05{print $0}' fst_out.txt | wc -l
awk 'NR==1{next;} $4 < 0.05{print $0}' fst_out.txt | less
awk 'NR==1{next;} $4 < 0.05 && $5 < 0{print $0}' fst_out.txt | cut -f1 -d" "> Bayescan_Outliers_DY_neg.list
} #notepad cleanup

#Gathering Outlier loci
{```{R}```
setPop(gen4) <- ~POP2
gen.list <- seppop(gen4)

#Outlier Loci
OF.DY <- read.table("Gan-DY_Outflank_outlier.list", head=F)[,1]						#Directional
BS.S2.pos <- read.table("Bayescan/Bayescan_Outliers_Site2_pos.list", head=F)[,1]	#Directional
BS.NM.pos <- read.table("Bayescan/Bayescan_Outliers_NM_pos.list", head=F)[,1]		#Directional

BS.S2.neg <- read.table("Bayescan/Bayescan_Outliers_Site2_neg.list", head=F)[,1]	#Balancing
BS.NM.neg <- read.table("Bayescan/Bayescan_Outliers_NM_neg.list", head=F)[,1]		#Balancing
BS.DY.neg <- read.table("Bayescan/Bayescan_Outliers_DY_neg.list", head=F)[,1]		#Balancing

out.pos <- unique(c(OF.DY, BS.S2.pos, BS.NM.pos))
out.neg <- unique(c(BS.S2.neg, BS.NM.neg, BS.DY.neg))

length(out.pos)

length(out.neg)
} #notepad cleanup

#Looking at PHL shift in Directional loci by library
{```{R}```
set.indv <- gen4@strata$INDV[gen4@strata$Site == "PHL"]
gen.tmp <- gen4[set.indv, ]
set.loc <- locNames(gen.tmp)[which(!locNames(gen.tmp) %in% c(out.pos, out.neg))]
gen.tmp.net <- gen.tmp[, loc=set.loc]
set.loc <- locNames(gen.tmp)[which(locNames(gen.tmp) %in% out.pos)]
gen.tmp.out.pos <- gen.tmp[, loc=set.loc]
set.loc <- locNames(gen.tmp)[which(locNames(gen.tmp) %in% out.neg)]
gen.tmp.out.neg <- gen.tmp[, loc=set.loc]

X.tmp <- scaleGen(gen.tmp, NA.method="mean", scale=F)
pca.tmp <- dudi.pca(X.tmp,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.tmp.net <- scaleGen(gen.tmp.net, NA.method="mean", scale=F)
pca.tmp.net <- dudi.pca(X.tmp.net,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.tmp.out.pos <- scaleGen(gen.tmp.out.pos, NA.method="mean", scale=F)
pca.tmp.out.pos <- dudi.pca(X.tmp.out.pos,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.tmp.out.neg <- scaleGen(gen.tmp.out.neg, NA.method="mean", scale=F)
pca.tmp.out.neg <- dudi.pca(X.tmp.out.neg,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

png("PCA_split_data_WT_Lib.png", res=200, width=4000, height=4000)
setPop(gen.tmp) <- setPop(gen.tmp.net) <- setPop(gen.tmp.out.pos) <- setPop(gen.tmp.out.neg) <- ~Lib/Site
par(mfrow=c(2,2))
ade4::s.class(pca.tmp$li, pop(gen.tmp), col=funky(5), cstar=0, axesell=F)
mtext(paste("All data (n=",length(locNames(gen.tmp))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.tmp.net$li, pop(gen.tmp.net), col=funky(5), cstar=0, axesell=F)
mtext(paste("Neutral data (n=",length(locNames(gen.tmp.net))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.tmp.out.pos$li, pop(gen.tmp.out.pos), col=funky(5), cstar=0, axesell=F)
mtext(paste("Directional Outlier data (n=",length(locNames(gen.tmp.out.pos))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.tmp.out.neg$li, pop(gen.tmp.out.neg), col=funky(5), cstar=0, axesell=F)
mtext(paste("Balancing Outlier data (n=",length(locNames(gen.tmp.out.neg))," loci)", sep=""), 3, 2, adj = 0.95)
dev.off()

#DAPC
set.indv <- gen4@strata$INDV[gen4@strata$Site == "PHL"]
gen.tmp <- gen4[set.indv, ]
set.loc <- locNames(gen.tmp)[which(locNames(gen.tmp) %in% out.pos)]
gen.tmp.out.pos <- gen.tmp[, loc=set.loc]

setPop(gen.tmp.out.pos) <- ~Lib/Site
X.tmp.out.pos <- scaleGen(gen.tmp.out.pos, NA.method="mean", scale=F)
xval.PHL <- xvalDapc(X.tmp.out.pos, pop(gen.tmp.out.pos), n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.PHL[2:6]

scatter(xval.PHL$DAPC, cstar=0, col=c2,posi.da="bottomright")

tail(xval.PHL$DAPC$var.contr[order(xval.Site$DAPC$var.contr),]

gen.tmp <- gen.list[["Gan-WT"]]
par(mfrow=c(2,2))

set.loc <- locNames(gen.tmp)[which(locNames(gen.tmp) %in% out.pos)]
gen.tmp.out.pos <- gen.tmp[, loc=set.loc]
X.tmp.out.pos <- scaleGen(gen.tmp.out.pos, NA.method="mean", scale=F)
pca.tmp.out.pos <- dudi.pca(X.tmp.out.pos,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
setPop(gen.tmp.out.pos) <- ~Lib/Site
ade4::s.class(pca.tmp.out.pos$li, pop(gen.tmp.out.pos), col=funky(5), cstar=0, axesell=F)

set.loc <- locNames(gen.tmp)[which(locNames(gen.tmp) %in% out.pos)]
set.loc <- set.loc[!set.loc %in% c("dDocent_Contig_28522")]
gen.tmp.out.pos <- gen.tmp[, loc=set.loc]
X.tmp.out.pos <- scaleGen(gen.tmp.out.pos, NA.method="mean", scale=F)
pca.tmp.out.pos <- dudi.pca(X.tmp.out.pos,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
setPop(gen.tmp.out.pos) <- ~Lib/Site
ade4::s.class(pca.tmp.out.pos$li, pop(gen.tmp.out.pos), col=funky(5), cstar=0, axesell=F)

set.loc <- locNames(gen.tmp)[which(locNames(gen.tmp) %in% out.pos)]
set.loc <- set.loc[!set.loc %in% c("dDocent_Contig_28522", "dDocent_Contig_136448")]
gen.tmp.out.pos <- gen.tmp[, loc=set.loc]
X.tmp.out.pos <- scaleGen(gen.tmp.out.pos, NA.method="mean", scale=F)
pca.tmp.out.pos <- dudi.pca(X.tmp.out.pos,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
setPop(gen.tmp.out.pos) <- ~Lib/Site
ade4::s.class(pca.tmp.out.pos$li, pop(gen.tmp.out.pos), col=funky(5), cstar=0, axesell=F)

set.loc <- locNames(gen.tmp)[which(locNames(gen.tmp) %in% out.pos)]
set.loc <- set.loc[!set.loc %in% c("dDocent_Contig_28522", "dDocent_Contig_136448", "dDocent_Contig_111655")]
gen.tmp.out.pos <- gen.tmp[, loc=set.loc]
X.tmp.out.pos <- scaleGen(gen.tmp.out.pos, NA.method="mean", scale=F)
pca.tmp.out.pos <- dudi.pca(X.tmp.out.pos,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
setPop(gen.tmp.out.pos) <- ~Lib/Site
ade4::s.class(pca.tmp.out.pos$li, pop(gen.tmp.out.pos), col=funky(5), cstar=0, axesell=F)

#Removing these 3 loci removes signal between the PHL samples
PHL_Lib_loci <- c("dDocent_Contig_28522", "dDocent_Contig_136448", "dDocent_Contig_111655")
} #notepad cleanup

#Looking at BAL shift in all loci nby library
{```{R}```
set.indv <- gen4@strata$INDV[gen4@strata$Site2 == "BAL"]
gen.tmp <- gen4[set.indv, ]

setPop(gen.tmp) <- ~Lib/Site
X.tmp <- scaleGen(gen.tmp, NA.method="mean", scale=F)
xval.BAL <- xvalDapc(X.tmp, pop(gen.tmp), n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.BAL[2:6]

scatter(xval.BAL$DAPC, cstar=0, col=c2,posi.da="bottomright")

tail(xval.BAL$DAPC$var.contr[order(xval.BAL$DAPC$var.contr),], n=10)

gen.tmp <- gen.list[["Gan-WT"]]

set.loc <- locNames(gen.tmp)[which(locNames(gen.tmp) %in% locNames(BAL_hunt[[1]]))]
gen.tmp.net <- gen.tmp[, loc=set.loc]
set.loc <- locNames(gen.tmp)[which(locNames(gen.tmp) %in% locNames(BAL_hunt[[2]]))]
gen.tmp.out <- gen.tmp[, loc=set.loc]

X.tmp <- scaleGen(gen.tmp, NA.method="mean", scale=F)
pca.tmp <- dudi.pca(X.tmp,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.tmp.net <- scaleGen(gen.tmp.net, NA.method="mean", scale=F)
pca.tmp.net <- dudi.pca(X.tmp.net,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.tmp.out <- scaleGen(gen.tmp.out, NA.method="mean", scale=F)
pca.tmp.out <- dudi.pca(X.tmp.out,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

indv.plot <- gen4@strata$INDV[gen4@strata$Site2 == "BAL"]
setPop(gen.tmp) <- setPop(gen.tmp.net) <- setPop(gen.tmp.out) <- ~Lib/Site
par(mfrow=c(2,2))
ade4::s.class(pca.tmp$li[which(rownames(pca.tmp$li) %in% indv.plot), 1:2], pop(gen.tmp)[indNames(gen.tmp) %in% indv.plot], col=funky(5), cstar=0, axesell=F)
mtext(paste("All data (n=",length(locNames(gen.tmp))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.tmp.net$li[which(rownames(pca.tmp.net$li) %in% indv.plot), 1:2], pop(gen.tmp.net)[indNames(gen.tmp.net) %in% indv.plot], col=funky(5), cstar=0, axesell=F)
mtext(paste("Neutral data (n=",length(locNames(gen.tmp.net))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.tmp.out$li[which(rownames(pca.tmp.out$li) %in% indv.plot), 1:2], pop(gen.tmp.out)[indNames(gen.tmp.out) %in% indv.plot], col=funky(5), cstar=0, axesell=F)
mtext(paste("Outlier data (n=",length(locNames(gen.tmp.out.pos))," loci)", sep=""), 3, 2, adj = 0.95)

set.indv <- gen4@strata$INDV[gen4@strata$Site2 == "BAL"]
gen.tmp <- gen4[set.indv, ]

BAL_hunt <- Locus_hunt(gen.tmp, xval.BAL$DAPC, THRES=0.99, GROUP="~Site", AXIS=1, ALPHA=0.05, MIN=0.85, ABS=0.00001, MODEL="normal")

locNames(BAL_hunt[[2]])

gen.tmp <- gen.list[["Gan-WT"]]

set.loc <- locNames(gen.tmp)[which(locNames(gen.tmp) %in% locNames(BAL_hunt[[1]]))]
gen.tmp.net <- gen.tmp[, loc=set.loc]
set.loc <- locNames(gen.tmp)[which(locNames(gen.tmp) %in% locNames(BAL_hunt[[2]]))]
gen.tmp.out <- gen.tmp[, loc=set.loc]

X.tmp <- scaleGen(gen.tmp, NA.method="mean", scale=F)
pca.tmp <- dudi.pca(X.tmp,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.tmp.net <- scaleGen(gen.tmp.net, NA.method="mean", scale=F)
pca.tmp.net <- dudi.pca(X.tmp.net,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.tmp.out <- scaleGen(gen.tmp.out, NA.method="mean", scale=F)
pca.tmp.out <- dudi.pca(X.tmp.out,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

indv.plot <- gen4@strata$INDV[gen4@strata$Site2 == "BAL"]
setPop(gen.tmp) <- setPop(gen.tmp.net) <- setPop(gen.tmp.out) <- ~Lib/Site
par(mfrow=c(2,2))
ade4::s.class(pca.tmp$li[which(rownames(pca.tmp$li) %in% indv.plot), 1:2], pop(gen.tmp)[indNames(gen.tmp) %in% indv.plot], col=funky(5), cstar=0, axesell=F)
mtext(paste("All data (n=",length(locNames(gen.tmp))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.tmp.net$li[which(rownames(pca.tmp.net$li) %in% indv.plot), 1:2], pop(gen.tmp.net)[indNames(gen.tmp.net) %in% indv.plot], col=funky(5), cstar=0, axesell=F)
mtext(paste("Neutral data (n=",length(locNames(gen.tmp.net))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.tmp.out$li[which(rownames(pca.tmp.out$li) %in% indv.plot), 1:2], pop(gen.tmp.out)[indNames(gen.tmp.out) %in% indv.plot], col=funky(5), cstar=0, axesell=F)
mtext(paste("Outlier data (n=",length(locNames(gen.tmp.out))," loci)", sep=""), 3, 2, adj = 0.95)

BAL_Lib_loci <- locNames(BAL_hunt[[2]])

#BAL_Lib_loci <- c("dDocent_Contig_179925","dDocent_Contig_89850","dDocent_Contig_52759","dDocent_Contig_183784","dDocent_Contig_99392","dDocent_Contig_51964","dDocent_Contig_194257","dDocent_Contig_200418","dDocent_Contig_72629",
#"dDocent_Contig_1956","dDocent_Contig_118313","dDocent_Contig_152201","dDocent_Contig_110732","dDocent_Contig_77511","dDocent_Contig_12561","dDocent_Contig_115705","dDocent_Contig_39185","dDocent_Contig_142132","dDocent_Contig_155603",
#"dDocent_Contig_116642")
} #notepad cleanup

#Flagging Monomorphic loci
{```{R}```
mono.loc <- names(which(gen4@loc.n.all == 1))
length(mono.loc)
} #notepad cleanup

#Finalize datasets (gen5)
{```{R}```
set.loc <- locNames(gen4)[!locNames(gen4) %in% c(PHL_Lib_loci, BAL_Lib_loci, mono.loc)]
gen5 <- gen4[, loc=set.loc, drop=T]
gen5.vcf <- VCF_remove(gen4.vcf, c(PHL_Lib_loci, BAL_Lib_loci, mono.loc))

setPop(gen5) <- ~POP2
gen.list <- seppop(gen5)

set.loc <- locNames(gen5)[which(!locNames(gen5) %in% c(out.pos, out.neg))]
gen.net <- gen5[, loc=set.loc]
set.loc <- locNames(gen5)[which(locNames(gen5) %in% out.pos)]
gen.out.pos <- gen5[, loc=set.loc]
set.loc <- locNames(gen5)[which(locNames(gen5) %in% out.neg)]
gen.out.neg <- gen5[, loc=set.loc]

gen.net.vcf <- VCF_remove(gen5.vcf, c(out.pos, out.neg))
set.loc <- locNames(gen5)[which(!locNames(gen5) %in% out.pos)]
gen.out.pos.vcf <- VCF_remove(gen5.vcf, set.loc)
set.loc <- locNames(gen5)[which(!locNames(gen5) %in% out.neg)]
gen.out.neg.vcf <- VCF_remove(gen5.vcf, set.loc)

save(gen5, file="gen5.gz", compress=T)
save(gen5.vcf, file="gen5.vcf.gz", compress=T)
#load("gen5.gz")
#load("gen5.vcf.gz")
} #notepad cleanup

#Plotting PCA
{```{R}```
#Reassigning thel location codes
tmp.v <- as.character(as.matrix(gen5@strata$Site2))
tmp.v[tmp.v == "KS"] <- "KGS"
tmp.v[tmp.v == "BAL"] <- "SS"
tmp.v[tmp.v == "PHL"] <- "PL"
gen5@strata$Site3 <- factor(tmp.v, levels=c("HEAD", "KGS", "EU", "SS", "PL", "ES", "BLNWR", "BC", "Sink7", "Sink27", "Sink31", "Sink37"))

set.loc <- locNames(gen5)[which(!locNames(gen5) %in% c(out.pos, out.neg))]
gen.net <- gen5[, loc=set.loc]
set.loc <- locNames(gen5)[which(locNames(gen5) %in% out.pos)]
gen.out.pos <- gen5[, loc=set.loc]
set.loc <- locNames(gen5)[which(locNames(gen5) %in% out.neg)]
gen.out.neg <- gen5[, loc=set.loc]

#Neutral data
gen.net <- gen.net[,loc=locNames(gen.net)[which(gen.net@loc.n.all > 1)]]
X.net <- scaleGen(gen.net, NA.method="mean", scale=F)
pca.net <- dudi.pca(X.net,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

#Getting colors
LC_pal <- colorRampPalette(c("#000000", "#B0D8A2", "#FFFFFF"))
LC_col <- LC_pal(11)

TC_pal <- colorRampPalette(c("#000000", "#FC9003", "#FFFFFF"))
TC_col <- TC_pal(11)

NM_pal <- colorRampPalette(c("#000000", "#87BBE5", "#FFFFFF"))
NM_col <- NM_pal(11)

col.v <- as.character(as.matrix(levels(gen.net@strata$Site3)))
col.v[col.v == "HEAD"] <- LC_col[8]
col.v[col.v == "KGS"] <- LC_col[6]
col.v[col.v == "EU"] <- LC_col[4]

col.v[col.v == "SS"] <- TC_col[8]
col.v[col.v == "ES"] <- TC_col[4]
col.v[col.v == "PL"] <- TC_col[6]

col.v[col.v == "BC"] <- NM_col[8]
col.v[col.v == "BLNWR"] <- NM_col[9]
col.v[col.v == "Sink7"] <- NM_col[7]
col.v[col.v == "Sink27"] <- NM_col[6]
col.v[col.v == "Sink31"] <- NM_col[5]
col.v[col.v == "Sink37"] <- NM_col[3]

par(mfrow=c(1,1))
tmp.dat <- cbind(gen.net@strata, pca.net$li[match(gen.net@strata$INDV, rownames(pca.net$li)),])
tmp.dat$POP2 <- factor(tmp.dat$POP2, levels=c("Gan-DY","Gan-WT","Gan-NM"))

p.net_Site2_v2 <- ggplot() + geom_point(data=tmp.dat, aes(x=Axis1, y=Axis2, fill=Site3, shape=POP2), size=3) + scale_fill_manual(values=col.v) + 
scale_shape_manual(name = "Region", labels = c("Leon Creek", "Toyah Creek", "New Mexico"), values=c(21,24,22)) + theme_classic() + 
labs(x=paste("PC1: ",sprintf("%.2f",(pca.net$eig[1]/sum(pca.net$eig))*100)," %", sep=""), y=paste("PC2: ",sprintf("%.2f",(pca.net$eig[2]/sum(pca.net$eig))*100)," %", sep="")) +
guides(fill=guide_legend(title="Sites", override.aes = list(shape = 23))) + 
ggtitle(paste("Neutral data (n=",length(locNames(gen.net))," loci)", sep=""))
#scale_shape_discrete(name = "Region", labels = c("Diamond Y", "West Texas", "New Mexico"), values=)
p.net_Site2_v2

ggsave(p.net_Site2_v2, file="PCA_gennet_Site2_v2.png", device="png")
ggsave(p.net_Site2_v2, file="PCA_gennet_Site2_v2.svg", device="svg")

#Directional selection
X.out.pos <- scaleGen(gen.out.pos, NA.method="mean", scale=F)
pca.out.pos <- dudi.pca(X.out.pos,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

par(mfrow=c(1,1))
tmp.dat <- cbind(gen.out.pos@strata, pca.out.pos$li[match(gen.out.pos@strata$INDV, rownames(pca.out.pos$li)),])
tmp.dat$POP2 <- factor(tmp.dat$POP2, levels=c("Gan-DY","Gan-WT","Gan-NM"))

p.out.pos_Site2_v2 <- ggplot() + geom_point(data=tmp.dat, aes(x=Axis1, y=Axis2, fill=Site3, shape=POP2), size=3) + scale_fill_manual(values=col.v) + 
scale_shape_manual(name = "Region", labels = c("Leon Creek", "Toyah Creek", "New Mexico"),values=c(21,24,22)) + theme_classic() + 
labs(x=paste("PC1: ",sprintf("%.2f",(pca.out.pos$eig[1]/sum(pca.out.pos$eig))*100)," %", sep=""), y=paste("PC2: ",sprintf("%.2f",(pca.out.pos$eig[2]/sum(pca.out.pos$eig))*100)," %", sep="")) +
guides(fill=guide_legend(title="Sites", override.aes = list(shape = 23)), shape=guide_legend(title="Regions")) + ggtitle(paste("Directional Outlier data (n=",length(locNames(gen.out.pos))," loci)", sep=""))
p.out.pos_Site2_v2

ggsave(p.out.pos_Site2_v2, file="PCA_genDir_Site2_v2.png", device="png")
ggsave(p.out.pos_Site2_v2, file="PCA_genDir_Site2_v2.svg", device="svg")

#Balancing selection
X.out.neg <- scaleGen(gen.out.neg, NA.method="mean", scale=F)
pca.out.neg <- dudi.pca(X.out.neg,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

par(mfrow=c(1,1))
tmp.dat <- cbind(gen.out.neg@strata, pca.out.neg$li[match(gen.out.neg@strata$INDV, rownames(pca.out.neg$li)),])
tmp.dat$POP2 <- factor(tmp.dat$POP2, levels=c("Gan-DY","Gan-WT","Gan-NM"))

p.out.neg_Site2_v2 <- ggplot() + geom_point(data=tmp.dat, aes(x=Axis1, y=Axis2, fill=Site3, shape=POP2), size=3) + scale_fill_manual(values=col.v) + 
scale_shape_manual(name = "Region", labels = c("Leon Creek", "Toyah Creek", "New Mexico"),values=c(21,24,22)) + theme_classic() + 
labs(x=paste("PC1: ",sprintf("%.2f",(pca.out.neg$eig[1]/sum(pca.out.neg$eig))*100)," %", sep=""), y=paste("PC2: ",sprintf("%.2f",(pca.out.neg$eig[2]/sum(pca.out.neg$eig))*100)," %", sep="")) +
guides(fill=guide_legend(title="Sites", override.aes = list(shape = 23)), shape=guide_legend(title="Regions")) + ggtitle(paste("Balancing Outlier data (n=",length(locNames(gen.out.neg))," loci)", sep=""))
p.out.neg_Site2_v2

ggsave(p.out.neg_Site2_v2, file="PCA_genBal_Site2_v2.png", device="png")
ggsave(p.out.neg_Site2_v2, file="PCA_genBal_Site2_v2.svg", device="svg")
}

#Looking at monomorphic loci by region
{```{R}```
setPop(gen.net) <- ~POP2
gen.list <- seppop(gen.net, drop=T)
lapply(gen.list, function(x) length(which(x@loc.n.all == 1)))

lapply(gen.list, function(x) length(which(x@loc.n.all == 1))/nLoc(x))
} #notepad cleanup

#Neutral K-means and DAPC analysis
{```{R}```
#Looking for groupings in neutral data
{```{R}```
grp <- find.clusters(gen.net, max.n.clust=20, n.pca=400, method = "kmeans")

grp$Kstat

png("Kmeans_net_BIC.png", res=200, width=2000, height=2000)
plot(grp$Kstat, pch=21, col="mediumblue", bg="white", xlab="Number of clusters", ylab="BIC")
lines(grp$Kstat, col="mediumblue")
dev.off()

G3.map <- Kmeans_heatmap(gen.net, 3, iter=100, plot=TRUE, "K3_net_all.png", verbose=TRUE, diagonal=TRUE)
G4.map <- Kmeans_heatmap(gen.net, 4, iter=100, plot=TRUE, "K4_net_all.png", verbose=TRUE, diagonal=TRUE)
G5.map <- Kmeans_heatmap(gen.net, 5, iter=100, plot=TRUE, "K5_net_all.png", verbose=TRUE, diagonal=TRUE)
G6.map <- Kmeans_heatmap(gen.net, 6, iter=100, plot=TRUE, "K6_net_all.png", verbose=TRUE, diagonal=TRUE)
G7.map <- Kmeans_heatmap(gen.net, 7, iter=100, plot=TRUE, "K7_net_all.png", verbose=TRUE, diagonal=TRUE)
G8.map <- Kmeans_heatmap(gen.net, 8, iter=100, plot=TRUE, "K8_net_all.png", verbose=TRUE, diagonal=TRUE)
G9.map <- Kmeans_heatmap(gen.net, 8, iter=100, plot=TRUE, "K9_net_all.png", verbose=TRUE, diagonal=TRUE)
G10.map <- Kmeans_heatmap(gen.net, 8, iter=100, plot=TRUE, "K10_net_all.png", verbose=TRUE, diagonal=TRUE)

#Adding groups to genind
tmp.hcl <- hclust(as.dist(1-G3.map))
tmp.clust <- cutree(tmp.hcl, k = 3)
gen.net@strata$G3 <- as.factor(tmp.clust[indNames(gen.net)])

tmp.hcl <- hclust(as.dist(1-G4.map))
tmp.clust <- cutree(tmp.hcl, k = 4)
gen.net@strata$G4 <- as.factor(tmp.clust[indNames(gen.net)])

tmp.hcl <- hclust(as.dist(1-G5.map))
tmp.clust <- cutree(tmp.hcl, k = 5)
gen.net@strata$G5 <- as.factor(tmp.clust[indNames(gen.net)])

tmp.hcl <- hclust(as.dist(1-G6.map))
tmp.clust <- cutree(tmp.hcl, k = 6)
gen.net@strata$G6 <- as.factor(tmp.clust[indNames(gen.net)])

tmp.hcl <- hclust(as.dist(1-G7.map))
tmp.clust <- cutree(tmp.hcl, k = 7)
gen.net@strata$G7 <- as.factor(tmp.clust[indNames(gen.net)])

tmp.hcl <- hclust(as.dist(1-G8.map))
tmp.clust <- cutree(tmp.hcl, k = 8)
gen.net@strata$G8 <- as.factor(tmp.clust[indNames(gen.net)])
} #notepad cleanup

#DAPC crossvalidation
{```{R}```
X.net <- scaleGen(gen.net, NA.method="mean", scale=F)
tiff("Xval_grps_net.tif", res=300, height =10000, width=2000)
par(mfrow=c(6,1))
xval.3 <- xvalDapc(X.net, gen.net@strata$G3, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.4 <- xvalDapc(X.net, gen.net@strata$G4, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.5 <- xvalDapc(X.net, gen.net@strata$G5, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.6 <- xvalDapc(X.net, gen.net@strata$G6, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.7 <- xvalDapc(X.net, gen.net@strata$G7, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.8 <- xvalDapc(X.net, gen.net@strata$G8, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
dev.off()

xval.Region <- xvalDapc(X.net, gen.net@strata$POP2, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.Region[2:6]
xval.Site2 <- xvalDapc(X.net, gen.net@strata$Site2, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.Site2[2:6]

#Outputting the crossvalidation data
results.df <- data.frame(matrix(ncol=5))
colnames(results.df) <- c("K", "Success", "Success_PCs", "RMSE", "RMSE_PCs")
for(i in 3:8){
xval.tmp <- get(paste("xval",i,sep="."))
tmp.row <- c(i, max((xval.tmp[[3]])), xval.tmp[[4]], min((xval.tmp[[5]])), xval.tmp[[6]])
results.df <- rbind(results.df, tmp.row)
}

results.df <- results.df[-1,]
results.df
} #notepad cleanup

#Composition plots
{```{R}```
#Directional Selection
setPop(gen.net) <- ~Site2
gen.net@pop <- factor(gen.net@pop, levels=c("HEAD","KS","EU","BAL","ES","PHL","BC","BLNWR","Sink7","Sink27","Sink31","Sink37"))
p3.comp <- ggcompoplot(xval.3$DAPC, gen.net, pal = funky(3), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 3 Groups")
p4.comp <- ggcompoplot(xval.4$DAPC, gen.net, pal = funky(4), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 4 Groups")
p5.comp <- ggcompoplot(xval.5$DAPC, gen.net, pal = funky(5), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 5 Groups")
p6.comp <- ggcompoplot(xval.6$DAPC, gen.net, pal = funky(6), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 6 Groups")
p7.comp <- ggcompoplot(xval.7$DAPC, gen.net, pal = funky(7), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 7 Groups")
p8.comp <- ggcompoplot(xval.8$DAPC, gen.net, pal = funky(8), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 8 Groups")
pSite2.comp <- ggcompoplot(xval.Site2$DAPC, gen.net, pal = funky(12), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("Site2 Groups")

p3.comp
p4.comp
p5.comp
p6.comp
p7.comp
p8.comp
pSite2.comp

ggsave("Gnobilis_K3_Site2_net.tif", p3.comp, device="tiff")
ggsave("Gnobilis_K4_Site2_net.tif", p4.comp, device="tiff")
ggsave("Gnobilis_K5_Site2_net.tif", p5.comp, device="tiff")
ggsave("Gnobilis_K6_Site2_net.tif", p6.comp, device="tiff")
ggsave("Gnobilis_K7_Site2_net.tif", p7.comp, device="tiff")
ggsave("Gnobilis_K8_Site2_net.tif", p8.comp, device="tiff")
ggsave("Gnobilis_Site2_Site2_net.tif", p8.comp, device="tiff")
} #notepad cleanup

#Exporting data for Arlequin
{```{R}```
setPop(gen.net) <- ~Site
writeGenPop(gen.net, "Gnobilis_Site_net.gen", "Neutral loci in G. nobilis data by Site")

setPop(gen.net) <- ~Site2
writeGenPop(gen.net, "Gnobilis_Site2_net.gen", "Neutral loci in G. nobilis data by Site2")

setPop(gen.net) <- ~G5
writeGenPop(gen.net, "Gnobilis_G5_net.gen", "Neutral loci in G. nobilis data by k-means 5")

setPop(gen.net) <- ~G8
writeGenPop(gen.net, "Gnobilis_G8_net.gen", "Neutral loci in G. nobilis data by k-means 8")
} #notepad cleanup
{```{bash}```
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gnobilis_Site_net.gen -inputformat GENEPOP -outputfile Gnobilis_Site_net.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gnobilis_Site2_net.gen -inputformat GENEPOP -outputfile Gnobilis_Site2_net.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid

java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gnobilis_G5_net.gen -inputformat GENEPOP -outputfile Gnobilis_G5_net.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gnobilis_G8_net.gen -inputformat GENEPOP -outputfile Gnobilis_G8_net.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid
} #notepad cleanup

#Fst MAX; modified KRD's code
{```{R}```
#Separate into groups
setPop(gen.net) <- ~POP2
gen.list <- seppop(gen.net)

#Function to modify the numbers in each vector
modify_numbers <- function(vec, num) {
  modified_vec <- sub("^\\d", num, vec)
  return(modified_vec)
}

#Change alleles
for(i in 2:length(gen.list)){
modified_list <- lapply(gen.list[[i]]@all.names, function(x) modify_numbers(x,(i-1)))
gen.list[[i]]@all.names <- modified_list
}

#repool the geninds back together
pooled_gens <- repool(gen.list)

#Export for Arlequin
setPop(gen.net) <- ~POP2
writeGenPop(gen.net, file.name = "Gnobilis_POP2_net.gen", comment="Gnobilis neutral loci separated by Region")

setPop(pooled_gens) <- ~POP2
writeGenPop(pooled_gens, file.name = "Gnobilis_POP2_for_Fmax_net.gen", comment="Gnobilis neutral loci separated by Region for Fmax calculation")
} #notepad cleanup
{```{bash}```
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gnobilis_POP2_net.gen -inputformat GENEPOP -outputfile Gnobilis_POP2_net.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gnobilis_POP2_for_Fmax_net.gen -inputformat GENEPOP -outputfile Gnobilis_POP2_for_Fmax_net.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid
} #notepad cleanup

#Exporting data for NeEstimator with no missing data
{```{R}```
#Creating a dataset with no missing data
loc.miss <- apply(gen.net@tab,2,function(x) sum(is.na(x)))
loc.list <- names(loc.miss)[which(loc.miss == 0)]
set.loc <- unlist(Loci_names(loc.list, SEP="[.]", REMOVE=1))
gen.tmp <- gen.net[,loc=set.loc]

setPop(gen.tmp) <- ~Site
writeGenPop(gen.tmp, "Gnobilis_Site_net_nomiss.gen", "Neutral loci in G. nobilis data by Site with no missing data")

setPop(gen.tmp) <- ~Site2
writeGenPop(gen.tmp, "Gnobilis_Site2_net_nomiss.gen", "Neutral loci in G. nobilis data by Site2 with no missing data")
} #notepad cleanup
#Shortening the names of the samples for Ne estimation
{```{bash}```
awk 'length($1) > 20 && NR > 2{$1=substr($1,1,20); $1=$1","; print $0; next}{print $0}' Gnobilis_Site_net_nomiss.gen > tmp.file; mv tmp.file Gnobilis_Site_net_nomiss.gen
awk 'length($1) > 20 && NR > 2{$1=substr($1,1,20); $1=$1","; print $0; next}{print $0}' Gnobilis_Site2_net_nomiss.gen > tmp.file; mv tmp.file Gnobilis_Site2_net_nomiss.gen
} #notepad cleanup
} #notepad cleanup for neutral loci

#Directional Outlier K-means and DAPC analysis
{```{R}```
#Looking for groupings in neutral data
{```{R}```
grp <- find.clusters(gen.out.pos, max.n.clust=20, n.pca=400, method = "kmeans")

grp$Kstat

png("Kmeans_out_pos_BIC.png", res=200, width=2000, height=2000)
plot(grp$Kstat, pch=21, col="mediumblue", bg="white", xlab="Number of clusters", ylab="BIC")
lines(grp$Kstat, col="mediumblue")
dev.off()

G3.map <- Kmeans_heatmap(gen.out.pos, 3, iter=100, plot=TRUE, "K3_out_pos_all.png", verbose=TRUE, diagonal=TRUE)
G4.map <- Kmeans_heatmap(gen.out.pos, 4, iter=100, plot=TRUE, "K4_out_pos_all.png", verbose=TRUE, diagonal=TRUE)
G5.map <- Kmeans_heatmap(gen.out.pos, 5, iter=100, plot=TRUE, "K5_out_pos_all.png", verbose=TRUE, diagonal=TRUE)
G6.map <- Kmeans_heatmap(gen.out.pos, 6, iter=100, plot=TRUE, "K6_out_pos_all.png", verbose=TRUE, diagonal=TRUE)
G7.map <- Kmeans_heatmap(gen.out.pos, 7, iter=100, plot=TRUE, "K7_out_pos_all.png", verbose=TRUE, diagonal=TRUE)
G8.map <- Kmeans_heatmap(gen.out.pos, 8, iter=100, plot=TRUE, "K8_out_pos_all.png", verbose=TRUE, diagonal=TRUE)

#Adding groups to genind
tmp.hcl <- hclust(as.dist(1-G3.map))
tmp.clust <- cutree(tmp.hcl, k = 3)
gen.out.pos@strata$G3 <- as.factor(tmp.clust[indNames(gen.out.pos)])

tmp.hcl <- hclust(as.dist(1-G4.map))
tmp.clust <- cutree(tmp.hcl, k = 4)
gen.out.pos@strata$G4 <- as.factor(tmp.clust[indNames(gen.out.pos)])

tmp.hcl <- hclust(as.dist(1-G5.map))
tmp.clust <- cutree(tmp.hcl, k = 5)
gen.out.pos@strata$G5 <- as.factor(tmp.clust[indNames(gen.out.pos)])

tmp.hcl <- hclust(as.dist(1-G6.map))
tmp.clust <- cutree(tmp.hcl, k = 6)
gen.out.pos@strata$G6 <- as.factor(tmp.clust[indNames(gen.out.pos)])

tmp.hcl <- hclust(as.dist(1-G7.map))
tmp.clust <- cutree(tmp.hcl, k = 7)
gen.out.pos@strata$G7 <- as.factor(tmp.clust[indNames(gen.out.pos)])

tmp.hcl <- hclust(as.dist(1-G8.map))
tmp.clust <- cutree(tmp.hcl, k = 8)
gen.out.pos@strata$G8 <- as.factor(tmp.clust[indNames(gen.out.pos)])
} #notepad cleanup

#DAPC crossvalidation
{```{R}```
X.out.pos <- scaleGen(gen.out.pos, NA.method="mean", scale=F)
tiff("Xval_grps_out_pos.tif", res=300, height =10000, width=2000)
par(mfrow=c(6,1))
xval.3 <- xvalDapc(X.out.pos, gen.out.pos@strata$G3, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.4 <- xvalDapc(X.out.pos, gen.out.pos@strata$G4, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.5 <- xvalDapc(X.out.pos, gen.out.pos@strata$G5, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.6 <- xvalDapc(X.out.pos, gen.out.pos@strata$G6, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.7 <- xvalDapc(X.out.pos, gen.out.pos@strata$G7, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.8 <- xvalDapc(X.out.pos, gen.out.pos@strata$G8, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
dev.off()

xval.Region <- xvalDapc(X.out.pos, gen.out.pos@strata$POP2, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.Region[2:6]
xval.Site2 <- xvalDapc(X.out.pos, gen.out.pos@strata$Site2, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.Site2[2:6]

#Outputting the crossvalidation data
results.df <- data.frame(matrix(ncol=5))
colnames(results.df) <- c("K", "Success", "Success_PCs", "RMSE", "RMSE_PCs")
for(i in 3:8){
xval.tmp <- get(paste("xval",i,sep="."))
tmp.row <- c(i, max((xval.tmp[[3]])), xval.tmp[[4]], min((xval.tmp[[5]])), xval.tmp[[6]])
results.df <- rbind(results.df, tmp.row)
}

results.df <- results.df[-1,]
results.df
} #notepad cleanup

#Composition plots
{```{R}```
#Directional Selection
setPop(gen.out.pos) <- ~Site2
gen.out.pos@pop <- factor(gen.out.pos@pop, levels=c("HEAD","KS","EU","BAL","ES","PHL","BC","BLNWR","Sink7","Sink27","Sink31","Sink37"))
p3.comp <- ggcompoplot(xval.3$DAPC, gen.out.pos, pal = funky(3), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 3 Groups")
p4.comp <- ggcompoplot(xval.4$DAPC, gen.out.pos, pal = funky(4), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 4 Groups")
p5.comp <- ggcompoplot(xval.5$DAPC, gen.out.pos, pal = funky(5), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 5 Groups")
p6.comp <- ggcompoplot(xval.6$DAPC, gen.out.pos, pal = funky(6), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 6 Groups")
p7.comp <- ggcompoplot(xval.7$DAPC, gen.out.pos, pal = funky(7), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 7 Groups")
p8.comp <- ggcompoplot(xval.8$DAPC, gen.out.pos, pal = funky(8), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 8 Groups")
pSite2.comp <- ggcompoplot(xval.Site2$DAPC, gen.out.pos, pal = funky(12), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("Site2 Groups")

p3.comp
p4.comp
p5.comp
p6.comp
p7.comp
p8.comp
pSite2.comp

ggsave("Gnobilis_K3_Site2_out_pos.tif", p3.comp, device="tiff")
ggsave("Gnobilis_K4_Site2_out_pos.tif", p4.comp, device="tiff")
ggsave("Gnobilis_K5_Site2_out_pos.tif", p5.comp, device="tiff")
ggsave("Gnobilis_K6_Site2_out_pos.tif", p6.comp, device="tiff")
ggsave("Gnobilis_K7_Site2_out_pos.tif", p7.comp, device="tiff")
ggsave("Gnobilis_K8_Site2_out_pos.tif", p8.comp, device="tiff")
ggsave("Gnobilis_Site2_Site2_out_pos.tif", p8.comp, device="tiff")
} #notepad cleanup

#Exporting data for Arlequin
{```{R}```
setPop(gen.out.pos) <- ~Site
writeGenPop(gen.out.pos, "Gnobilis_Site_out_pos.gen", "Neutral loci in G. nobilis data by Site")

setPop(gen.out.pos) <- ~Site2
writeGenPop(gen.out.pos, "Gnobilis_Site2_out_pos.gen", "Neutral loci in G. nobilis data by Site2")

setPop(gen.out.pos) <- ~G3
writeGenPop(gen.out.pos, "Gnobilis_G3_out_pos.gen", "Neutral loci in G. nobilis data by k-means 3")

setPop(gen.out.pos) <- ~G4
writeGenPop(gen.out.pos, "Gnobilis_G8_out_pos.gen", "Neutral loci in G. nobilis data by k-means 4")
} #notepad cleanup
{```{bash}```
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gnobilis_Site_out_pos.gen -inputformat GENEPOP -outputfile Gnobilis_Site_out_pos.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gnobilis_Site2_out_pos.gen -inputformat GENEPOP -outputfile Gnobilis_Site2_out_pos.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid

java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gnobilis_G5_out_pos.gen -inputformat GENEPOP -outputfile Gnobilis_G5_out_pos.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gnobilis_G8_out_pos.gen -inputformat GENEPOP -outputfile Gnobilis_G8_out_pos.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid
} #notepad cleanup
} #directional selection notepad cleanup

#Balancing Outlier K-means and DAPC analysis
{```{R}```
#Looking for groupings in neutral data
{```{R}```
grp <- find.clusters(gen.out.neg, max.n.clust=20, n.pca=400, method = "kmeans")

grp$Kstat

png("Kmeans_out_neg_BIC.png", res=200, width=2000, height=2000)
plot(grp$Kstat, pch=21, col="mediumblue", bg="white", xlab="Number of clusters", ylab="BIC")
lines(grp$Kstat, col="mediumblue")
dev.off()

G3.map <- Kmeans_heatmap(gen.out.neg, 3, iter=100, plot=TRUE, "K3_out_neg_all.png", verbose=TRUE, diagonal=TRUE)
G4.map <- Kmeans_heatmap(gen.out.neg, 4, iter=100, plot=TRUE, "K4_out_neg_all.png", verbose=TRUE, diagonal=TRUE)
G5.map <- Kmeans_heatmap(gen.out.neg, 5, iter=100, plot=TRUE, "K5_out_neg_all.png", verbose=TRUE, diagonal=TRUE)
G6.map <- Kmeans_heatmap(gen.out.neg, 6, iter=100, plot=TRUE, "K6_out_neg_all.png", verbose=TRUE, diagonal=TRUE)
G7.map <- Kmeans_heatmap(gen.out.neg, 7, iter=100, plot=TRUE, "K7_out_neg_all.png", verbose=TRUE, diagonal=TRUE)
G8.map <- Kmeans_heatmap(gen.out.neg, 8, iter=100, plot=TRUE, "K8_out_neg_all.png", verbose=TRUE, diagonal=TRUE)

#Adding groups to genind
tmp.hcl <- hclust(as.dist(1-G3.map))
tmp.clust <- cutree(tmp.hcl, k = 3)
gen.out.neg@strata$G3 <- as.factor(tmp.clust[indNames(gen.out.neg)])

tmp.hcl <- hclust(as.dist(1-G4.map))
tmp.clust <- cutree(tmp.hcl, k = 4)
gen.out.neg@strata$G4 <- as.factor(tmp.clust[indNames(gen.out.neg)])

tmp.hcl <- hclust(as.dist(1-G5.map))
tmp.clust <- cutree(tmp.hcl, k = 5)
gen.out.neg@strata$G5 <- as.factor(tmp.clust[indNames(gen.out.neg)])

tmp.hcl <- hclust(as.dist(1-G6.map))
tmp.clust <- cutree(tmp.hcl, k = 6)
gen.out.neg@strata$G6 <- as.factor(tmp.clust[indNames(gen.out.neg)])

tmp.hcl <- hclust(as.dist(1-G7.map))
tmp.clust <- cutree(tmp.hcl, k = 7)
gen.out.neg@strata$G7 <- as.factor(tmp.clust[indNames(gen.out.neg)])

tmp.hcl <- hclust(as.dist(1-G8.map))
tmp.clust <- cutree(tmp.hcl, k = 8)
gen.out.neg@strata$G8 <- as.factor(tmp.clust[indNames(gen.out.neg)])
} #notepad cleanup

#DAPC crossvalidation
{```{R}```
X.out.neg <- scaleGen(gen.out.neg, NA.method="mean", scale=F)
tiff("Xval_grps_out_neg.tif", res=300, height =10000, width=2000)
par(mfrow=c(6,1))
xval.3 <- xvalDapc(X.out.neg, gen.out.neg@strata$G3, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.4 <- xvalDapc(X.out.neg, gen.out.neg@strata$G4, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.5 <- xvalDapc(X.out.neg, gen.out.neg@strata$G5, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.6 <- xvalDapc(X.out.neg, gen.out.neg@strata$G6, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.7 <- xvalDapc(X.out.neg, gen.out.neg@strata$G7, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.8 <- xvalDapc(X.out.neg, gen.out.neg@strata$G8, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
dev.off()

xval.Region <- xvalDapc(X.out.neg, gen.out.neg@strata$POP2, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.Region[2:6]
xval.Site2 <- xvalDapc(X.out.neg, gen.out.neg@strata$Site2, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.Site2[2:6]

#Outputting the crossvalidation data
results.df <- data.frame(matrix(ncol=5))
colnames(results.df) <- c("K", "Success", "Success_PCs", "RMSE", "RMSE_PCs")
for(i in 3:8){
xval.tmp <- get(paste("xval",i,sep="."))
tmp.row <- c(i, max((xval.tmp[[3]])), xval.tmp[[4]], min((xval.tmp[[5]])), xval.tmp[[6]])
results.df <- rbind(results.df, tmp.row)
}
results.df <- results.df[-1,]
results.df
} #notepad cleanup

#Composition plots
{```{R}```
#Directional Selection
setPop(gen.out.neg) <- ~Site2
gen.out.neg@pop <- factor(gen.out.neg@pop, levels=c("HEAD","KS","EU","BAL","ES","PHL","BC","BLNWR","Sink7","Sink27","Sink31","Sink37"))
p3.comp <- ggcompoplot(xval.3$DAPC, gen.out.neg, pal = funky(3), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 3 Groups")
p4.comp <- ggcompoplot(xval.4$DAPC, gen.out.neg, pal = funky(4), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 4 Groups")
p5.comp <- ggcompoplot(xval.5$DAPC, gen.out.neg, pal = funky(5), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 5 Groups")
p6.comp <- ggcompoplot(xval.6$DAPC, gen.out.neg, pal = funky(6), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 6 Groups")
p7.comp <- ggcompoplot(xval.7$DAPC, gen.out.neg, pal = funky(7), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 7 Groups")
p8.comp <- ggcompoplot(xval.8$DAPC, gen.out.neg, pal = funky(8), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("K-means 8 Groups")
pSite2.comp <- ggcompoplot(xval.Site2$DAPC, gen.out.neg, pal = funky(12), cols=3) + theme(axis.text.x = element_blank()) + ggtitle("Site2 Groups")

p3.comp
p4.comp
p5.comp
p6.comp
p7.comp
p8.comp
pSite2.comp

ggsave("Gnobilis_K3_Site2_out_neg.tif", p3.comp, device="tiff")
ggsave("Gnobilis_K4_Site2_out_neg.tif", p4.comp, device="tiff")
ggsave("Gnobilis_K5_Site2_out_neg.tif", p5.comp, device="tiff")
ggsave("Gnobilis_K6_Site2_out_neg.tif", p6.comp, device="tiff")
ggsave("Gnobilis_K7_Site2_out_neg.tif", p7.comp, device="tiff")
ggsave("Gnobilis_K8_Site2_out_neg.tif", p8.comp, device="tiff")
ggsave("Gnobilis_Site2_Site2_out_neg.tif", p8.comp, device="tiff")
} #notepad cleanup

#Exporting data for Arlequin
{```{R}```
setPop(gen.out.neg) <- ~Site2
writeGenPop(gen.out.neg, "Gnobilis_Site2_out_neg.gen", "Neutral loci in G. nobilis data by Site2")
} #notepad cleanup
{```{bash}```
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gnobilis_Site2_out_neg.gen -inputformat GENEPOP -outputfile Gnobilis_Site2_out_neg.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid
} #notepad cleanup
} #balancing selection notepad cleanup

### Diversity Comparisons ###
#Making a list of samples by Site
{```{bash}```
mkdir pi
} #notepad cleanup
#Exporting data for determining pi by Site2
{```{R}```
#Exporting samples by Site2
for(i in levels(gen5@strata$Site2)){if(length(gen5@strata$INDV[gen5@strata$Site2 == i]) > 0){length(write.table(gen5@strata$INDV[gen5@strata$Site2 == i], file=paste("pi/",i,"_indv.txt",sep=""), col.names=F, row.names=F, quote=F))}}

#Exporting neutral loci
tmp.m <- matrix(unlist(strsplit(names(gen5.vcf@loc.n.all),"_")), byrow=T, ncol=4)
tmp.df <- data.frame(CHROM=apply(tmp.m[,1:3], 1, function(x) paste(x,collapse="_")),POS=tmp.m[,4])
write.table(tmp.df, file="neutral_loci.txt", col.names=F, row.names=F, quote=F, sep="\t")
}}}} #notebook cleanup
#Getting pi for neutral loci
{```{bash}```
#Splitting populations
ls pi/*_indv.txt | while read i; do
POP=$(echo $i | cut -f1,2 -d"_")
vcftools --vcf SNP.TRS.F09.vcf --out $POP --recode --recode-INFO-all --keep $i --positions neutral_loci.txt
done
#Getting nucleotide diversity for each locus
ls pi/*_indv.txt | while read i; do
POP=$(echo $i | cut -f1,2 -d"_")
vcftools --vcf $POP.recode.vcf --out $POP --window-pi 10000
done
#Getting nucleotide diversity for the entire dataset
vcftools --vcf SNP.TRS.F09.vcf --out SNP.TRS.F09.net --recode --recode-INFO-all --positions neutral_loci.txt
vcftools --vcf SNP.TRS.F09.net.recode.vcf --out SNP.TRS.F09.net --window-pi 10000
mv SNP.TRS.F09.net.windowed.pi Site2.pi

cp ../Gnobilis_Jun2024/dDocent_reference.fasta .
python ~/bin/seq_length.py dDocent_reference.fasta > ref_length.txt
} #notepad cleanup

##Sorting pi files into a single file##
#Loading data
{```{R}```
ref_len <- read.table("ref_length.txt", head=F)

FILES <- list.files(path = "pi", pattern="*.windowed.pi")
LOC <- matrix(unlist(strsplit(FILES, "[.]")), ncol=4, byrow=T)[,1]
LOCI <- unique(read.table("neutral_loci.txt", head=F)[,1])
POP2.pi <- read.table("Site2.pi", head=T)
} #notepad cleanup

#If there are loci not present in the pi, checking to see if they are missing for the data or have 0 variation
{```{R}```
for(i in LOC){
print(paste("Processing", i, "data"))
assign(paste(i, "pi", sep="_"), read.table(paste("pi/",i,".txt.windowed.pi", sep=""), head=T))
tmp.df <- get(paste(i, "pi", sep="_"))
tmp.df[,1] <- as.character(as.matrix(tmp.df[,1]))
vcf <- read.vcfR(file=paste("pi/", i, ".txt.recode.vcf", sep=""), verbose=F)
gen.vcf <- vcfR2genind(vcf)
vcf.loci <- Loci_names(locNames(gen.vcf), SEP="_")
vcf.loci <- vcf.loci[!vcf.loci[,1] %in% as.character(as.matrix(tmp.df$CHROM)),]
for(j in vcf.loci){

if(is.integer(gen.vcf@tab[,grep(paste(j,"_",sep=""), colnames(gen.vcf@tab))])){tmp.m <- table(gen.vcf@tab[,grep(paste(j,"_",sep=""), colnames(gen.vcf@tab))], useNA="always")
}else{tmp.m <- apply(gen.vcf@tab[,grep(paste(j,"_",sep=""), colnames(gen.vcf@tab))], 2, function(x) table(x,useNA="always")) }

if(nrow(tmp.m) > 1){tmp.df <- rbind(tmp.df, c(j, 1, 10000, ncol(tmp.m), 0))}
}
assign(paste(i, "pi", sep="_"), tmp.df)
}

for(i in LOC){
tmp_pi <- get(paste(i, "pi", sep="_"))
tmp_pi$PI <- as.numeric(as.matrix(tmp_pi$PI))
cat(i,"\t", round(mean(tmp_pi$PI),7), "\n")
}

rm(vcf)

#Testing the number of loci per 
for(i in ls()[grep("_pi", ls())]){
print(paste(i, nrow(get(i)), sep=":"))
}

#Testing for NA values
for(i in ls()[grep("_pi", ls())]){
print(paste(i, sum(is.na(get(i)[,5])), sep=":"))
}
}}} #notepad cleanup

#Getting POP2 pi data into a data.frame
{```{R}```
POP_pi <- POP2.pi[,c(1,4)]
POP_pi <- cbind(POP_pi,matrix(ncol=length(levels(gen5@strata$Site2)), nrow=nrow(POP_pi)))
colnames(POP_pi) <- c(colnames(POP_pi)[1:2], levels(gen5@strata$Site2))
for(i in levels(gen5@strata$Site2)){POP_pi[,i] <- get(paste(i, "indv_pi", sep="_"))[match(POP_pi$CHROM, get(paste(i, "indv_pi", sep="_"))[,1]),5]}

#Scaling by the number of bases in the locus
for(i in 3:ncol(POP_pi)){POP_pi[,i] <- as.numeric(as.matrix(POP_pi[,i]))*10000}
POP_pi[,1] <- as.character(as.matrix(POP_pi[,1]))

for(i in 1:nrow(POP_pi)){POP_pi[i,3:ncol(POP_pi)] <- POP_pi[i,3:ncol(POP_pi)]/ref_len[ref_len[,1] == POP_pi[i,1],2]}

write.table(POP_pi, file="Site2_pi.txt", col.names=T, row.names=F, quote=F, sep="\t")
} #notepad cleanup

#Making variables
{```{R}```
dna <- ape::read.dna("dDocent_reference.fasta", format="fasta")
fa <- seqinr::read.fasta("dDocent_reference.fasta")
vcf<-read.vcfR(file="SNP.TRS.F09.vcf")

LOCI_tmp <- matrix(unlist(str_split(locNames(gen5), "_")),ncol=3,byrow=T)
LOCI_COUNT <- table(paste(LOCI_tmp[,1], LOCI_tmp[,2], LOCI_tmp[,3], sep="_"))
LOCI <- unique(paste(LOCI_tmp[,1], LOCI_tmp[,2], LOCI_tmp[,3], sep="_"))
LOCI_LEN <- sapply(fa, function(x) length(x[x != "n"]))
gff.df <- data.frame(seqid=LOCI, source=".", type="gene", start=1, end=as.numeric(LOCI_LEN[match(LOCI,names(LOCI_LEN))]), score=".", strand="+", phase=".", attributes=paste("ID=gene_",LOCI,sep=""))
} #notepad cleanup

#Preparing the final data frame
{```{R}```
sum.data <- data.frame(matrix(ncol=5))
colnames(sum.data) <- c("Pop","n","Parameter", "mean","sd")
} #notepad cleanup
#Loop gathering summary data
{```{R}```
STRATA <- "Site2"
for(j in levels(gen5@strata[,STRATA])){
if(length(which(gen5@strata[,STRATA] == j)) == 0){next}
data.df <- data.frame(matrix(ncol=10))
colnames(data.df) <- c("pop","locus","n","Nuc_div","Hap_div","Seg_sites","theta.s","tajima_D","tajima_pval_nor","tajima_pval_beta")
print(paste("Processing population", j))
INDV <- as.character(as.matrix(gen5@strata$INDV[gen5@strata[, STRATA] == j]))
for(i in 1:nrow(gff.df)){
tmp.vcf <- vcf[, LOCI[i], INDV]
if(sum(is.na(extract.gt(tmp.vcf[tmp.vcf@fix[,1] == gff.df[i,1]]))) > 0){
INDV.NA <- apply(extract.gt(tmp.vcf[tmp.vcf@fix[,1] == gff.df[i,1]]),2,function(x) sum(is.na(x)))
tmp.INDV <- names(INDV.NA)[which(INDV.NA == 0)]
if(length(tmp.INDV) < 2 ){
	data.df <- rbind(data.df, c(j, as.character(gff.df[i,1]), length(INDV), "NA", "NA", "NA", "NA", "NA", "NA", "NA"))
next
}
tmp.vcf <- vcf[, LOCI[i], tmp.INDV]
}
my_dnabin1 <- vcfR2DNAbin(tmp.vcf[tmp.vcf@fix[,1] == gff.df[i,1]], consensus = F, extract.haps = T, ref.seq=dna[paste(gff.df[i,1])], start.pos=gff.df[i,4], verbose=F, unphased_as_NA =F)
tmp.n <- nrow(my_dnabin1)
tmp.nuc <- nuc.div(my_dnabin1)
tmp.hap <- hap.div(my_dnabin1)
tmp.seg <- length(seg.sites(my_dnabin1))
tmp.theta <- theta.s(tmp.seg,nrow(my_dnabin1))
tmp.TD <- tajima.test(my_dnabin1)
data.df <- rbind(data.df, c(j, as.character(gff.df[i,1]), tmp.n/2, tmp.nuc, tmp.hap, tmp.seg, tmp.theta, tmp.TD$D, tmp.TD$Pval.normal, tmp.TD$Pval.beta))
if(i/500 == round(i/500,0)){print(paste("Finished locus", i, "out of", nrow(gff.df)))}
}
data.df <- data.df[-1,]
write.table(data.df, file=paste(j,"theta_plus.txt", sep="_"), col.names=T, row.names=F, quote=F, sep="\t")

for(i in c("Nuc_div","Hap_div","Seg_sites","theta.s","tajima_D")){
sum.data <- rbind(sum.data, c(j, data.df[1,3], i, mean(as.numeric(data.df[,i]), na.rm=T), sd(as.numeric(data.df[,i]), na.rm=T)))
}
}

sum.data <- sum.data[-1,]
head(sum.data)
}}}}}  #notepad cleanup

#Preparing the data tables
{```{R}```
HEAD_data.df <- read.table(paste("HEAD","theta_plus.txt", sep="_"), head=T, sep="\t")
KS_data.df <- read.table(paste("KS","theta_plus.txt", sep="_"), head=T, sep="\t")
EU_data.df <- read.table(paste("EU","theta_plus.txt", sep="_"), head=T, sep="\t")
BAL_data.df <- read.table(paste("BAL","theta_plus.txt", sep="_"), head=T, sep="\t")
ES_data.df <- read.table(paste("ES","theta_plus.txt", sep="_"), head=T, sep="\t")
PHL_data.df <- read.table(paste("PHL","theta_plus.txt", sep="_"), head=T, sep="\t")
BC_data.df <- read.table(paste("BC","theta_plus.txt", sep="_"), head=T, sep="\t")
BLNWR_data.df <- read.table(paste("BLNWR","theta_plus.txt", sep="_"), head=T, sep="\t")
Sink27_data.df <- read.table(paste("Sink27","theta_plus.txt", sep="_"), head=T, sep="\t")
Sink31_data.df <- read.table(paste("Sink31","theta_plus.txt", sep="_"), head=T, sep="\t")
Sink37_data.df <- read.table(paste("Sink37","theta_plus.txt", sep="_"), head=T, sep="\t")
Sink7_data.df <- read.table(paste("Sink7","theta_plus.txt", sep="_"), head=T, sep="\t")
} #notepad cleanup

#Allelic Richness
{```{R}```
setPop(gen5) <- ~Site2
Ar_data.df <- allelic.richness(gen5)$Ar
} #notepad cleanup

#Nucleotide, Haplotype, Segregation sites and theta s (Watterson's theta)
{```{R}```
STRATA <- "Site2"
for(i in 2:4){
VAR <- c("Nuc_div","Hap_div","Seg_sites","theta.s")[i]
data.df <- NULL
for(j in levels(gen5@strata[,STRATA])){
if(j == "CHS"){next}
if(length(which(gen5@strata[,STRATA] == j)) == 0){next}
PREFIX <- matrix(unlist(strsplit(j,"-")),ncol=2)[1,2]
data.df <- cbind(data.df, get(paste(PREFIX,"_data.df",sep=""))[,VAR])
colnames(data.df)[ncol(data.df)] <- j
}
assign(paste(VAR, "_data.df", sep=""), data.df)
}

Nuc_div_data.df <- read.table("Site2_pi.txt", head=T)
Nuc_div_data.df <- Nuc_div_data.df[3:ncol(Nuc_div_data.df)]
colnames(Nuc_div_data.df) <- gsub("[.]","-",colnames(Nuc_div_data.df))
}}}}}  #notepad cleanup

#Expected Heterozygosity from Arlequin
{```{R}```
He_data.df <- read.csv("Site2_He.csv", head=T)
He_data.df <- He_data.df[,2:13]
} #notepad cleanup

##Comparing Site groups ##
#Running the Friedmans
{```{R}```
for(i in 1:6){
VAR <- c("Nuc_div","Hap_div","Seg_sites","theta.s", "Ar", "He")[i]
data.df <- get(paste(VAR, "_data.df", sep=""))

#Removing rows with NA
na.test <- apply(data.df, 1, function(x) sum(is.na(x)))
data.df <- data.df[which(na.test == 0),]

#"Gathering" data
data.m <- matrix(ncol=3)
for(j in 1:length(locNames(gen5)[which(na.test == 0)])){
tmp.m <- matrix(c(rep(locNames(gen5)[j],ncol(data.df)), colnames(data.df), as.numeric(data.df[j,])),ncol=3,byrow=F)
data.m <- rbind(data.m, tmp.m)
}
test.df <- as.data.frame(data.m[-1, ])
colnames(test.df) <- c("Locus","POP","VAR")

#Friedman's test
tmp.dat <- friedman.test(VAR ~ POP | Locus, data=test.df)
print(paste("Friedman's test on", VAR, "( n=", nrow(data.df),")"))
print(tmp.dat)
print(apply(get(paste(VAR, "_data.df", sep="")), 2, function(x) mean(x, na.rm=T)))
}
} #notepad cleanup

#Saving image
{```{R}```
save.image("Fried_Site2.RData.gz", compress=T)
#load("Fried_Site2.RData.gz")
} #notepad cleanup

#Running the Wilcox Test
{```{R}```
error.list <- NULL

for(i in 1:6){
VAR <- c("Nuc_div","Hap_div","Seg_sites","theta.s", "Ar", "He")[i]
data.df <- get(paste(VAR, "_data.df", sep=""))

data.m <- matrix(ncol=3)
for(j in 1:length(locNames(gen5))){
tmp.m <- matrix(c(rep(locNames(gen5)[j],ncol(data.df)), colnames(data.df), as.numeric(data.df[j,])),ncol=3,byrow=F)
data.m <- rbind(data.m, tmp.m)
}
test.df <- as.data.frame(data.m[-1, ])
colnames(test.df) <- c("Locus","POP","VAR")
test.df$Locus <- factor(test.df$Locus)
test.df$POP <- factor(test.df$POP)
test.df$VAR <- as.numeric(as.matrix(test.df$VAR))
result_m <- t(combn(levels(test.df$POP), 2))
result_m <- cbind(result_m, matrix(ncol=3, nrow=nrow(result_m)))

for(k in 1:nrow(result_m)){
tmp.df <- test.df[test.df$POP %in% result_m[k,1:2],]
tmp.df$POP <- gdata::drop.levels(tmp.df$POP)
if(length(table(tmp.df$POP)) != 2){next}
rm.loci <- as.character(as.matrix(unique(tmp.df$Locus[is.na(tmp.df$VAR)])))
tmp.df <- tmp.df[!tmp.df$Locus %in% rm.loci,]
tmp.df$Locus <- gdata::drop.levels(tmp.df$Locus)
possibleError <- tryCatch(tmp.out <- wilcoxsign_test(VAR ~ POP | Locus, tmp.df), error=function(e) e)
if(inherits(possibleError, "error")){error.list <- c(error.list,c(VAR,result_m[k,1:2])); next}

tmp.df$POP <- factor(tmp.df$POP, levels=c(result_m[k,1], result_m[k,2]))
tmp.out <- wilcoxsign_test(VAR ~ POP | Locus,tmp.df)
tmp.stat <- statistic(tmp.out)
tmp.p <- pvalue(tmp.out)
result_m[k,3:4] <- c(tmp.stat, tmp.p)

tmp.td <- tidyr::spread(tmp.df, POP, VAR)
tmp.td$diff <- tmp.td[,2] - tmp.td[,3]
tmp.td$rank <- rank(abs(tmp.td$diff))
tmp.td$sgn <- tmp.td$diff
tmp.td$sgn[tmp.td$sgn > 0] <- 1
tmp.td$sgn[tmp.td$sgn < 0] <- -1

result_m[k,5] <- sum(tmp.td$sgn*tmp.td$rank)
}

result_df <- as.data.frame(result_m)
colnames(result_df) <- c("Loc1", "Loc2", "stat", "pvalue", "T_stat")
for(j in 3:5){result_df[,j] <- as.numeric(as.matrix(result_df[,j]))}
result_df$pvalue_adj <- p.adjust(result_df$pvalue, method="fdr")
print(paste("Wilcox test on", VAR))
print(paste("Number of significant results:", length(which(result_df$pvalue < 0.05))))
print(paste("Number of significant results adj:", length(which(result_df$pvalue_adj < 0.05))))
#print(result_df[result_df$pvalue < 0.05])
assign(paste(VAR,"_wilcox.df",sep=""), result_df)
}

#Putting results into a table
if(exists("all_df")){rm(all_df)}

for(i in 1:6){
VAR <- c("Nuc_div","Hap_div","Seg_sites","theta.s", "Ar", "He")[i]
result_df <- get(paste(VAR, "wilcox.df", sep="_"))
colnames(result_df)[3:5] <- paste(VAR, colnames(result_df)[3:5], sep="_")
if(!exists("all_df")){all_df <- result_df
}else{all_df <- cbind(all_df, result_df[,3:5])}
}

#Viewing results
all_df

#Saving image
save.image("Wilcox_Loc2_net.RData.gz", compress=T)
#load("Wilcox_Loc2_net.RData.gz")
}}}}}  #notepad cleanup
