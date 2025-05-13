####### Gambusia last stages of filtering and beginning of analysis #######

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
} #notepad cleanup

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

} #notepad cleanup
} #notepad cleanup

#Importing Data
{```{R}```
strata <- read.csv(file = "../../filtering/filtering_v2_May2024/SNP.TRS.F06.2024-05-17/indv.csv", header = TRUE)
#Adding nobilis grouping inforrmation
tmp.v <- as.character(as.matrix(strata$POP))
tmp.v[tmp.v %in% c("Gaa-GSW1", "Gaa-HEAD", "Gaa-LCCSP", "Gaf-LCCSP")] <- "Gaa"
tmp.v[tmp.v %in% c("Gag-CHS", "Gag-ES","Gag-SMR")] <- "Gag"
tmp.v[tmp.v %in% c("Gan-CHS", "Gan-ES", "Gan-PHL", "Gan-BAL", "Gno-BAL")] <- "Gan-WT"
tmp.v[tmp.v %in% c("Gan-EU", "Gan-FLS", "Gan-HEAD", "Gan-DY", "Gno-DY")] <- "Gan-DY"
tmp.v[tmp.v %in% c("Gan-Sink27", "Gan-Sink31", "Gan-Sink37", "Gan-Sink7", "Gan-BC", "Gan-BLNWR", "Gno-BLNWR", "Gno-BC")] <- "Gan-NM"
strata$POP2 <- as.factor(tmp.v)

vcf <- read.vcfR(file="SNP.TRS.F09.vcf")
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
}

#Setting color schemes
{```{R}```
#           "Gaa"   "Gag"    "Gan-DY"     "Gan-NM"      "Gan-WT"
col.All<-c("red4","grey50","dodgerblue","mediumblue","mediumorchid")

c2 <- c("mediumblue", "red4")
c3 <- c("red4", "mediumblue", "chocolate2")
c4 <- c("lightslateblue", "mediumblue", "red4", "dodgerblue")
c5 <- c("dodgerblue", "mediumblue", "royalblue1", "lightslateblue", "red4")
c6 <- c("mediumblue", "lightslateblue", "royalblue1", "red4","cornflowerblue", "dodgerblue")
}

#Library Effects (Outflank)
{```{R}```
setPop(gen.vcf) <- ~Lib
tmp <- gl.outflank(gen.vcf, qthreshold = 0.1)
length(which(tmp$outflank$results[15]=="TRUE" & tmp$outflank$results[11]=="goodH"))

tmp.v <- as.character(as.matrix(gen.vcf@strata$POP2))
tmp.v[tmp.v == "Gaa-GSW"] <- "Gaa"
gen.vcf@strata$POP2 <- as.factor(tmp.v)

setPop(gen.vcf) <- ~POP2
gen.list <- seppop(gen.vcf, drop=T)

#Gambusia affinis
{```{R}```
gen.Gaa <- gen.list[["Gaa"]]
mono.loc <- names(which(gen.Gaa@loc.n.all == 1))
length(mono.loc)

tmp.df <- data.frame(matrix(ncol=3, nrow=ncol(gen.Gaa@tab)), row.names=colnames(gen.Gaa@tab))
for(i in 1:3){
set.indv <- gen.Gaa@strata$INDV[gen.Gaa@strata$Lib == levels(gen.Gaa@strata$Lib)[i]]
tmp.df[,i] <- apply(gen.Gaa@tab[set.indv,], 2, function(x) sum(is.na(x))/length(set.indv))
}
tmp.v <- apply(tmp.df, 1, function(x) length(which(x == 1)))
loc.na <- Loci_names(names(which(tmp.v > 0)))[,1]
length(loc.na)

set.loc <- locNames(gen.Gaa)[!locNames(gen.Gaa) %in% c(mono.loc,loc.na)]
gen.Gaa <- gen.Gaa[,loc=set.loc]

setPop(gen.Gaa) <- ~Lib
tmp.Gaa <- gl.outflank(gen.Gaa, qthreshold = 0.1)
length(which(tmp.Gaa$outflank$results[15]=="TRUE" & tmp.Gaa$outflank$results[11]=="goodH"))

outliers <- tmp.Gaa$outflank$results[which(tmp.Gaa$outflank$results[15]=="TRUE" & tmp.Gaa$outflank$results[11]=="goodH"),1]
out.list <- matrix(unlist(strsplit(as.vector(outliers),"[.]")),ncol=2,byrow=T)[,1]
tmp.Gaa <- matrix(unlist(strsplit(as.vector(outliers),"[_]")),ncol=4,byrow=T)
gen.out.list<-unique(paste(tmp.Gaa[,1],tmp.Gaa[,2],tmp.Gaa[,3],sep="_"))
length(gen.out.list) 

write.table(gen.out.list, "Gaa_Outflank_Lib_outlier.list", quote=F, col.names=F, row.names=F)
} #notepad cleanup

#Gambusia geiseri
{```{R}```
gen.Gag <- gen.list[["Gag"]]
mono.loc <- names(which(gen.Gag@loc.n.all == 1))
length(mono.loc)

gen.Gag@strata$Lib <- droplevels(gen.Gag@strata$Lib)
tmp.df <- data.frame(matrix(ncol=length(levels(gen.Gag@strata$Lib)), nrow=ncol(gen.Gag@tab)), row.names=colnames(gen.Gag@tab))
for(i in 1:length(levels(gen.Gag@strata$Lib))){
set.indv <- as.character(as.matrix(gen.Gag@strata$INDV[gen.Gag@strata$Lib == levels(gen.Gag@strata$Lib)[i]]))
tmp.df[,i] <- apply(gen.Gag@tab[set.indv,], 2, function(x) sum(is.na(x))/length(set.indv))
}
tmp.v <- apply(tmp.df, 1, function(x) length(which(x == 1)))
loc.na <- Loci_names(names(which(tmp.v > 0)))[,1]
length(loc.na)

set.loc <- locNames(gen.Gag)[!locNames(gen.Gag) %in% c(mono.loc,loc.na)]
gen.Gag <- gen.Gag[,loc=set.loc]

setPop(gen.Gag) <- ~Lib
tmp.Gag <- gl.outflank(gen.Gag, qthreshold = 0.1)
length(which(tmp.Gag$outflank$results[15]=="TRUE" & tmp.Gag$outflank$results[11]=="goodH"))

outliers <- tmp.Gag$outflank$results[which(tmp.Gag$outflank$results[15]=="TRUE" & tmp.Gag$outflank$results[11]=="goodH"),1]
out.list <- matrix(unlist(strsplit(as.vector(outliers),"[.]")),ncol=2,byrow=T)[,1]
tmp.Gag <- matrix(unlist(strsplit(as.vector(outliers),"[_]")),ncol=4,byrow=T)
gen.out.list<-unique(paste(tmp.Gag[,1],tmp.Gag[,2],tmp.Gag[,3],sep="_"))
length(gen.out.list) 

write.table(gen.out.list, "Gag_Outflank_Lib_outlier.list", quote=F, col.names=F, row.names=F)
} #notepad cleanup

#Gambusia nobilis
{```{R}```
gen.Gan <- repool(gen.list[["Gan-WT"]], gen.list[["Gan-DY"]], gen.list[["Gan-NM"]])
mono.loc <- names(which(gen.Gan@loc.n.all == 1))
length(mono.loc)

gen.Gan@strata$Lib <- droplevels(gen.Gan@strata$Lib)
tmp.df <- data.frame(matrix(ncol=length(levels(gen.Gan@strata$Lib)), nrow=ncol(gen.Gan@tab)), row.names=colnames(gen.Gan@tab))
for(i in 1:length(levels(gen.Gan@strata$Lib))){
set.indv <- as.character(as.matrix(gen.Gan@strata$INDV[gen.Gan@strata$Lib == levels(gen.Gan@strata$Lib)[i]]))
tmp.df[,i] <- apply(gen.Gan@tab[set.indv,], 2, function(x) sum(is.na(x))/length(set.indv))
}
tmp.v <- apply(tmp.df, 1, function(x) length(which(x == 1)))
length(which(tmp.v > 0))

set.loc <- locNames(gen.Gan)[!locNames(gen.Gan) %in% mono.loc]
gen.Gan <- gen.Gan[,loc=set.loc]

setPop(gen.Gan) <- ~Lib
tmp.Gan <- gl.outflank(gen.Gan, qthreshold = 0.1)
length(which(tmp.Gan$outflank$results[15]=="TRUE" & tmp.Gan$outflank$results[11]=="goodH"))
} #notepad cleanup

#Gambusia nobilis - West Texas
{```{R}```
gen.WT <- gen.list[["Gan-WT"]]
mono.loc <- names(which(gen.WT@loc.n.all == 1))
length(mono.loc)

gen.WT@strata$Lib <- droplevels(gen.WT@strata$Lib)
tmp.df <- data.frame(matrix(ncol=length(levels(gen.WT@strata$Lib)), nrow=ncol(gen.WT@tab)), row.names=colnames(gen.WT@tab))
for(i in 1:length(levels(gen.WT@strata$Lib))){
set.indv <- as.character(as.matrix(gen.WT@strata$INDV[gen.WT@strata$Lib == levels(gen.WT@strata$Lib)[i]]))
tmp.df[,i] <- apply(gen.WT@tab[set.indv,], 2, function(x) sum(is.na(x))/length(set.indv))
}
tmp.v <- apply(tmp.df, 1, function(x) length(which(x == 1)))
length(which(tmp.v > 0))

set.loc <- locNames(gen.WT)[!locNames(gen.WT) %in% mono.loc]
gen.WT <- gen.WT[,loc=set.loc]

setPop(gen.WT) <- ~Lib
tmp.WT <- gl.outflank(gen.WT, qthreshold = 0.1)
length(which(tmp.WT$outflank$results[15]=="TRUE" & tmp.WT$outflank$results[11]=="goodH"))
} #notepad cleanup

#Gambusia nobilis - New Mexico
{```{R}```
gen.NM <- gen.list[["Gan-NM"]]
mono.loc <- names(which(gen.NM@loc.n.all == 1))
length(mono.loc)

gen.NM@strata$Lib <- droplevels(gen.NM@strata$Lib)
tmp.df <- data.frame(matrix(ncol=length(levels(gen.NM@strata$Lib)), nrow=ncol(gen.NM@tab)), row.names=colnames(gen.NM@tab))
for(i in 1:length(levels(gen.NM@strata$Lib))){
set.indv <- as.character(as.matrix(gen.NM@strata$INDV[gen.NM@strata$Lib == levels(gen.NM@strata$Lib)[i]]))
tmp.df[,i] <- apply(gen.NM@tab[set.indv,], 2, function(x) sum(is.na(x))/length(set.indv))
}
tmp.v <- apply(tmp.df, 1, function(x) length(which(x == 1)))
length(which(tmp.v > 0))

set.loc <- locNames(gen.NM)[!locNames(gen.NM) %in% mono.loc]
gen.NM <- gen.NM[,loc=set.loc]

setPop(gen.NM) <- ~Lib
tmp.NM <- gl.outflank(gen.NM, qthreshold = 0.1)
length(which(tmp.NM$outflank$results[15]=="TRUE" & tmp.NM$outflank$results[11]=="goodH"))
} #notepad cleanup

#Gambusia nobilis - Diamond Y
{```{R}```
gen.DY <- gen.list[["Gan-DY"]]
mono.loc <- names(which(gen.DY@loc.n.all == 1))
length(mono.loc)
#[1] 12963

gen.DY@strata$Lib <- droplevels(gen.DY@strata$Lib)
tmp.df <- data.frame(matrix(ncol=length(levels(gen.DY@strata$Lib)), nrow=ncol(gen.DY@tab)), row.names=colnames(gen.DY@tab))
for(i in 1:length(levels(gen.DY@strata$Lib))){
set.indv <- as.character(as.matrix(gen.DY@strata$INDV[gen.DY@strata$Lib == levels(gen.DY@strata$Lib)[i]]))
tmp.df[,i] <- apply(gen.DY@tab[set.indv,], 2, function(x) sum(is.na(x))/length(set.indv))
}
tmp.v <- apply(tmp.df, 1, function(x) length(which(x == 1)))
length(which(tmp.v > 0))

loc.na <- Loci_names(names(which(tmp.v > 0)))[,1]
length(loc.na)

set.loc <- locNames(gen.DY)[!locNames(gen.DY) %in% c(mono.loc, loc.na)]
gen.DY <- gen.DY[,loc=set.loc]

setPop(gen.DY) <- ~Lib
tmp.DY <- gl.outflank(gen.DY, qthreshold = 0.1)
length(which(tmp.DY$outflank$results[15]=="TRUE" & tmp.DY$outflank$results[11]=="goodH"))

outliers <- tmp.DY$outflank$results[which(tmp.DY$outflank$results[15]=="TRUE" & tmp.DY$outflank$results[11]=="goodH"),1]
out.list <- matrix(unlist(strsplit(as.vector(outliers),"[.]")),ncol=2,byrow=T)[,1]
tmp.DY <- matrix(unlist(strsplit(as.vector(outliers),"[_]")),ncol=4,byrow=T)
gen.out.list<-unique(paste(tmp.DY[,1],tmp.DY[,2],tmp.DY[,3],sep="_"))
length(gen.out.list)

write.table(gen.out.list, "Gan-DY_Outflank_Lib_outlier.list", quote=F, col.names=F, row.names=F)
} #notepad cleanup
} #notepad cleanup

#Outlier Loci (BAYESCAN)
{```{R}```
#Exporting data for Bayescan
setPop(gen.Gaa) <- ~Lib
writeGenPop(gen.Gaa, "Gaa_Lib.gen", "G. affinis data by Library")

setPop(gen.Gag) <- ~Lib
writeGenPop(gen.Gag, "Gag_Lib.gen", "G. geiseri data by Library")

setPop(gen.WT) <- ~Lib
writeGenPop(gen.WT, "Gan-WT_Lib.gen", "G. nobilis in West TX data by Library")

setPop(gen.NM) <- ~Lib
writeGenPop(gen.NM, "Gan-NM_Lib.gen", "G. nobilis in NM data by Library")

setPop(gen.DY) <- ~Lib
writeGenPop(gen.DY, "Gan-DY_Lib.gen", "G. nobilis in DY data by Library")
} #notepad cleanup

#Converting to BS format
{```{bash}```
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gaa_Lib.gen -inputformat GENEPOP -outputfile Gaa_Lib.BS -outputformat BAYESCAN -spid /home/afields/bin/genepop_to_BS.spid
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gag_Lib.gen -inputformat GENEPOP -outputfile Gag_Lib.BS -outputformat BAYESCAN -spid /home/afields/bin/genepop_to_BS.spid
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gan-WT_Lib.gen -inputformat GENEPOP -outputfile Gan-WT_Lib.BS -outputformat BAYESCAN -spid /home/afields/bin/genepop_to_BS.spid
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gan-NM_Lib.gen -inputformat GENEPOP -outputfile Gan-NM_Lib.BS -outputformat BAYESCAN -spid /home/afields/bin/genepop_to_BS.spid
java8 -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile Gan-DY_Lib.gen -inputformat GENEPOP -outputfile Gan-DY_Lib.BS -outputformat BAYESCAN -spid /home/afields/bin/genepop_to_BS.spid
} #notepad cleanup

#Running BAYESCAN
{```{bash}```
mkdir Bayescan
bayescan Gaa_Lib.BS -od ./Bayescan -all-trace -threads 15 -thin 100 -nbp 30 -pr_odds 100
bayescan Gag_Lib.BS -od ./Bayescan -all-trace -threads 15 -thin 100 -nbp 30 -pr_odds 100
bayescan Gan-WT_Lib.BS -od ./Bayescan -all-trace -threads 30 -thin 100 -nbp 30 -pr_odds 100
bayescan Gan-NM_Lib.BS -od ./Bayescan -all-trace -threads 30 -thin 100 -nbp 30 -pr_odds 100
bayescan Gan-DY_Lib.BS -od ./Bayescan -all-trace -threads 30 -thin 100 -nbp 30 -pr_odds 100
} #notepad cleanup

#Analyzing results
{```{bash}```
head -n2 ../Gaa_Lib.gen | tail -n 1 | sed 's/, /\n/g' > contigs
echo "Contigs" | cat - contigs > tmp; mv tmp contigs
paste contigs Gaa_Li_fst.txt | awk 'NR==1{print $0}NR!=1{$2=""; print}' > fst_Gaa.txt
sort -k4,4n fst_Gaa.txt | less
awk 'NR==1{next;} $4 < 0.05{print $0}' fst_Gaa.txt | wc -l

head -n2 ../Gag_Lib.gen | tail -n 1 | sed 's/, /\n/g' > contigs
echo "Contigs" | cat - contigs > tmp; mv tmp contigs
paste contigs Gag_Li_fst.txt | awk 'NR==1{print $0}NR!=1{$2=""; print}' > fst_Gag.txt
sort -k4,4n fst_Gag.txt | less
awk 'NR==1{next;} $4 < 0.05{print $0}' fst_Gag.txt | wc -l

head -n2 ../Gan-NM_Lib.gen | tail -n 1 | sed 's/, /\n/g' > contigs
echo "Contigs" | cat - contigs > tmp; mv tmp contigs
paste contigs Gan-NM_Li_fst.txt | awk 'NR==1{print $0}NR!=1{$2=""; print}' > fst_Gan-NM.txt
sort -k4,4n fst_Gan-NM.txt | less
awk 'NR==1{next;} $4 < 0.05{print $0}' fst_Gan-NM.txt | wc -l

head -n2 ../Gan-WT_Lib.gen | tail -n 1 | sed 's/, /\n/g' > contigs
echo "Contigs" | cat - contigs > tmp; mv tmp contigs
paste contigs Gan-WT_Li_fst.txt | awk 'NR==1{print $0}NR!=1{$2=""; print}' > fst_Gan-WT.txt
sort -k4,4n fst_Gan-WT.txt | less
awk 'NR==1{next;} $4 < 0.05{print $0}' fst_Gan-WT.txt | wc -l

head -n2 ../Gan-DY_Lib.gen | tail -n 1 | sed 's/, /\n/g' > contigs
echo "Contigs" | cat - contigs > tmp; mv tmp contigs
paste contigs Gan-DY_Li_fst.txt | awk 'NR==1{print $0}NR!=1{$2=""; print}' > fst_Gan-DY.txt
sort -k4,4n fst_Gan-DY.txt | less
awk 'NR==1{next;} $4 < 0.05{print $0}' fst_Gan-DY.txt | wc -l

awk 'NR==1{next;} $4 < 0.05{print $0}' fst_Gan-DY.txt | cut -f1 -d" "> Bayescan_Library_DY.list
} #notepad cleanup

#Checking convergance
#http://evomics.org/wp-content/uploads/2016/01/BayeScan_BayeScEnv_exercises.pdf
{```{R}```
library(coda)

##G. affinis
{```{R}```
chain <- read.table("Gaa_Li.sel",header=TRUE)
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

#G. geiseri
{```{R}```
chain <- read.table("Gag_Li.sel",header=TRUE)
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

#G. nobilis NM
{```{R}```
chain <- read.table("Gan-NM_Li.sel",header=TRUE)
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

#G. nobilis WTX
{```{R}```
chain <- read.table("Gan-WT_Li.sel",header=TRUE)
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
}} #notepad cleanup

#Loading BAYESCAN outliers
{```{R}```
out.Gaa <- read.table("Gaa_Outflank_Lib_outlier.list", head=F)[,1]
out.Gag <- read.table("Gag_Outflank_Lib_outlier.list", head=F)[,1]
out.Gan <- read.table("Gan-DY_Outflank_Lib_outlier.list", head=F)[,1]
gen.out.loc <- unique(c(out.Gaa,out.Gag,out.Gan))
length(gen.out.loc)
} #notepad cleanup

#Removing monomorphics
{```{R}```
mono.loc <- names(which(gen@loc.n.all == 1))
length(mono.loc)
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
}

#Removing duplicates not picked up by strata due to naming difference
dup.list <- list(
Sample1=c("Gaa-LCCSP_MGL-19747_Lib1_I2", "Gaf-LCCSP_19747-Gaa_Lib3_I7"),
Sample2=c("Gaa-LCCSP_MGL-19758_Lib1_I9", "Gaf-LCCSP_19758-Gaa_Lib3_I2"),
Sample3=c("Gaa-LCCSP_MGL-19765_Lib1_I2", "Gaf-LCCSP_19765-Gaa_Lib3_I10"),
Sample4=c("Gaf-LCCSP_19746-Gaa_Lib3_I2", "Gaa-LCCSP_MGL-19746_Lib1_I9"),
Sample5=c("Gaf-LCCSP_19755-Gaa_Lib3_I10", "Gaa-LCCSP_MGL-19755_Lib1_I9"),
Sample6=c("Gaf-LCCSP_19766-Gaa_Lib3_I4", "Gaa-LCCSP_MGL-19766_Lib1_I7"),
Sample7=c("Gaa-LCCSP_MGL-19757_Lib1_I8", "Gaf-LCCSP_19757-Gaa_Lib3_I4"))

rel.list <- list()
for(i in 1:length(dup.list)){
j <- dup.list[[i]]

tmp <- gen.vcf@strata[gen.vcf@strata$INDV %in% j, ]
tmp$Diff <- abs(as.numeric(as.character(tmp[ , "E.HOM."])) - as.numeric(as.character(tmp[ , "O.HOM."])))

Miss.tab <- table(as.numeric(as.character(tmp$N_MISS)))
Depth.tab <- table(as.numeric(as.character(tmp$MEAN_DEPTH)))

MIN.Miss <- min(as.numeric(names(Miss.tab)))
MAX.Depth <- max(as.numeric(names(Depth.tab)))

if(Miss.tab[names(Miss.tab)==MIN.Miss]==1 && Depth.tab[names(Depth.tab)==MAX.Depth]==1){keep <- as.character(tmp[tmp$N_MISS==MIN.Miss,1])
} else {keep <- as.character(tmp[tmp$Diff == min(tmp$Diff),1])}
print(keep)

rel.list[[i]] <- as.character(tmp[grep(keep, tmp$INDV, invert=T),1])
rm(keep, tmp)
}

#Gno-LCCSP_19759-Gan_Lib3_I2 and Gaa-LCCSP_MGL-19774_Lib1_I11 were confused during sampling and therefore removed
OF.Gaa <- read.table("Gaa_Outflank_Lib_outlier.list", head=F)[,1] 
OF.Gag <- read.table("Gag_Outflank_Lib_outlier.list", head=F)[,1] 
OF.Gan <- read.table("Gan-DY_Outflank_Lib_outlier.list", head=F)[,1] 
BS.Gan <- read.table("Bayescan/Bayescan_Library_DY.list", head=F)[,1] 

set.indv <- indNames(gen)[!indNames(gen) %in% c(unlist(rm.list),"Gno-LCCSP_19759-Gan_Lib3_I2",unlist(rel.list),"Gaa-LCCSP_MGL-19774_Lib1_I11")]
set.loc <- locNames(gen)[which(!locNames(gen) %in% c(mono.loc, OF.Gaa, OF.Gag, OF.Gan, BS.Gan))]

gen2 <- gen[set.indv, loc=set.loc, drop=T]
gen2.vcf <- VCF_remove(gen.vcf[set.indv, ,drop=T], c(mono.loc, OF.Gaa, OF.Gag, OF.Gan, BS.Gan))

save(gen2, file="gen2.gz", compress=T)
save(gen2.vcf, file="gen2.vcf.gz", compress=T)
#load("gen2.gz")
#load("gen2.vcf.gz")
}}} #notepad cleanup

#Relatedness within species/region
{```{R}```
setPop(gen2) <- ~ POP2
tmp.list <- seppop(gen2, drop=T)

#Looking for relatedness between of the full dataset
related.out <- list()
related.review <- NULL

for(i in names(tmp.list)){
if(nInd(tmp.list[[i]]) < 2){next}
print(paste("Processing",i))
mono.loc <- names(which(tmp.list[[i]]@loc.n.all == 1))
print(paste("Monomorphic loci:",length(mono.loc)))
set.loc <- locNames(tmp.list[[i]])[!locNames(tmp.list[[i]]) %in% mono.loc]
tmp.gl <- gi2gl(tmp.list[[i]][,loc=set.loc,drop=T])
tmp.rel <- gl2related(tmp.gl, save=F)
tmp.coan <- coancestry(tmp.rel, wang =1)
related.out[[i]] <- tmp.coan
related.review <- rbind(related.review,tmp.coan$relatedness[which(tmp.coan$relatedness$wang > 0.5),c(2:3,6)])
}
}} #notepad cleanup

#PCA of raw data
{```{R}```
#Standard
X <- scaleGen(gen2, NA.method="mean", scale=F)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
setPop(gen2) <- ~POP2
ade4::s.class(pca1$li, pop(gen2), col=funky(5), cstar=0)

tiff("PCA_POP_raw_all_loci.tif", res=300, height =2000, width=2000)
ade4::s.class(pca1$li, pop(gen2), col=col.All, cstar=0)
dev.off()

tmp.v <- as.character(as.matrix(gen2@strata$POP2))
tmp.v[tmp.v %in% c("Gan-WT", "Gan-DY", "Gan-NM")] <- "Gan"
gen2@strata$POP3 <- as.factor(tmp.v)
head(gen2@strata)
}

######### Hybrid analysis #########
### Adegenet ###
#G. geiseri and nobilis-WT
{```{R}```
#Define Strata
gen.tmp <- gen2[indNames(gen2)[gen2@strata$POP2 %in% c("Gag","Gan-WT")], , drop=T]
setPop(gen.tmp) <- ~ POP2
temp <- seppop(gen.tmp)
geis <- temp$`Gag`
nob <- temp$`Gan-WT`

#Test for missing loci
geis.miss <- apply(geis@tab, 2, function(x) sum(is.na(x))/nInd(geis))
nob.miss <- apply(nob@tab, 2, function(x) sum(is.na(x))/nInd(nob))

rm.loci <- NULL
if(length(which(geis.miss  == 1))){rm.loci <- c(rm.loci, names(which(geis.miss == 1)))}
if(length(which(nob.miss  == 1))){rm.loci <- c(rm.loci, names(which(nob.miss == 1)))}

if(length(rm.loci) > 0){
rm.loci <- unique(matrix(unlist(strsplit(rm.loci,"[.]")), ncol=2, byrow=T)[,1])
set.loc <- locNames(geis)[!locNames(geis) %in% rm.loci]
gen.vcf <- gen.vcf[, loc=set.loc]
geis <- geis[, loc=set.loc]
nob <- nob[, loc=set.loc]
}

nLoc(nob)

#Simulate samples
F1 <- hybridize(geis, nob, n = 100, pop = "geisxnob")
geis_bx <- hybridize(geis, F1, n = 100, pop = "geis_bx")
nob_bx <- hybridize(nob, F1, n = 100, pop = "nob_bx")
F2 <- hybridize(F1, F1, n= 100, pop = "F1XF1")
pooled_gens <- repool(gen.tmp, F1, geis_bx, nob_bx, F2)

save(pooled_gens, file="G_x_N_WT_sims.gz", compress=T)
#load("G_x_N_WT_sims.gz")
#PCA
table(pop(pooled_gens))
x <- scaleGen(pooled_gens, NA.method = "mean")

#PCA with all loci
pca <- dudi.pca(x,center = FALSE, scale = FALSE, scannf = FALSE, nf = 4)

#plot2
col = c("#68a691", "#ddd392", "#ed6a5a", "#e59a2e", "#bdbdbd", "#636363", '#000000', '#799c41')
s.class(pca$li, pop(pooled_gens), xax=1, yax=2, cellipse=3, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0)
MIN.x <- par("usr")[1]
MAX.x <- par("usr")[2]
MIN.y <- par("usr")[3]
MAX.y <- par("usr")[4]

png("ade_hybrid/GgxGn_WT_sim.png", res=200, width=2000, height=2000)
s.class(pca$li[pop(pooled_gens) %in% c("geisxnob", "geis_bx", "nob_bx", "F1XF1"), 1:2], pop(pooled_gens)[pop(pooled_gens) %in% c("geisxnob", "geis_bx", "nob_bx", "F1XF1")], 
xax=1, yax=2, cellipse=10, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0, xlim=c(MIN.x*1.25,MAX.x), ylim=c(MIN.y,MAX.y))
points(pca$li[1:nInd(gen.tmp),1:2], pch=19, col=col[as.numeric(pop(pooled_gens)[1:nInd(gen.tmp)])])
points(pca$li[14:15,1:2], pch=19, col=col[as.numeric(pop(pooled_gens)[14:15])])
mtext("G. geiseri and G. nobilis (West Texas)", 3, adj=0.05, cex=1.5)
dev.off()

#Finding samples around the sims
pop_center <- vector()
tmp.range <- vector()
for(i in c("geisxnob", "geis_bx", "nob_bx", "F1XF1")){
pop_center <- c(pop_center, mean(pca$li[which(pop(pooled_gens) == i),1]))
tmp.low <- min(pca$li[which(pop(pooled_gens) == i),1])
tmp.high <- max(pca$li[which(pop(pooled_gens) == i),1])
tmp.range <- c(min(tmp.range,tmp.low), max(tmp.range,tmp.high))
}
names(pop_center) <- c("geisxnob", "geis_bx", "nob_bx", "F1XF1")

tmp.pca <- pca$li[1:nInd(gen.tmp),1:2]
for(i in which(tmp.pca[,1] > tmp.range[1]*1.5 & tmp.pca[,1] < tmp.range[2]*1.5)){
tmp.err <- (tmp.pca[i,1] - pop_center)^2
tmp.sim <- names(tmp.err)[which(tmp.err == min(tmp.err))]
print(paste(rownames(pca$li)[i],": ",tmp.sim,"    ","pos=",round(pca$li[i,1],3),", ",round(pca$li[i,2],3),"    ",tmp.sim,"=",round(pop_center[tmp.sim],3),sep=""))
}

#Field MisIdentified or library mix ups
rownames(pca$li)[pca$li[,1] > pop_center["nob_bx"] & rownames(pca$li) %in% indNames(geis)]
rownames(pca$li)[pca$li[,1] < pop_center["geis_bx"] & rownames(pca$li) %in% indNames(nob)]
}}}} #notepad cleanup

#G. geiseri and nobilis-DY
{```{R}```
#Define Strata
gen.tmp <- gen2[indNames(gen2)[gen2@strata$POP2 %in% c("Gag","Gan-DY")], , drop=T]
setPop(gen.tmp) <- ~ POP2
temp <- seppop(gen.tmp)
geis <- temp$`Gag`
nob <- temp$`Gan-DY`

#Test for missing loci
geis.miss <- apply(geis@tab, 2, function(x) sum(is.na(x))/nInd(geis))
nob.miss <- apply(nob@tab, 2, function(x) sum(is.na(x))/nInd(nob))

rm.loci <- NULL
if(length(which(geis.miss  == 1))){rm.loci <- c(rm.loci, names(which(geis.miss == 1)))}
if(length(which(nob.miss  == 1))){rm.loci <- c(rm.loci, names(which(nob.miss == 1)))}

if(length(rm.loci) > 0){
rm.loci <- unique(matrix(unlist(strsplit(rm.loci,"[.]")), ncol=2, byrow=T)[,1])
set.loc <- locNames(geis)[!locNames(geis) %in% rm.loci]
gen.vcf <- gen.vcf[, loc=set.loc]
geis <- geis[, loc=set.loc]
nob <- nob[, loc=set.loc]
}

nLoc(nob)

#Simulate samples
F1 <- hybridize(geis, nob, n = 100, pop = "geisxnob")
geis_bx <- hybridize(geis, F1, n = 100, pop = "geis_bx")
nob_bx <- hybridize(nob, F1, n = 100, pop = "nob_bx")
F2 <- hybridize(F1, F1, n= 100, pop = "F1XF1")
pooled_gens <- repool(gen.tmp, F1, geis_bx, nob_bx, F2)

save(pooled_gens, file="G_x_N_DY_sims.gz", compress=T)
#load("G_x_N_DY_sims.gz")

#PCA
x <- scaleGen(pooled_gens, NA.method = "mean")

#PCA with all loci
pca <- dudi.pca(x,center = FALSE, scale = FALSE, scannf = FALSE, nf = 4)

#plot2
col = c("#68a691", "#ddd392", "#ed6a5a", "#e59a2e", "#bdbdbd", "#636363", '#000000', '#799c41')
s.class(pca$li, pop(pooled_gens), xax=1, yax=2, cellipse=3, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0)
MIN.x <- par("usr")[1]
MAX.x <- par("usr")[2]
MIN.y <- par("usr")[3]
MAX.y <- par("usr")[4]

png("ade_hybrid/GgxGn_DY_sim.png", res=200, width=2000, height=2000)
s.class(pca$li[pop(pooled_gens) %in% c("geisxnob", "geis_bx", "nob_bx", "F1XF1"), 1:2], pop(pooled_gens)[pop(pooled_gens) %in% c("geisxnob", "geis_bx", "nob_bx", "F1XF1")], 
xax=1, yax=2, cellipse=10, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0, xlim=c(MIN.x,MAX.x), ylim=c(MIN.y*2,MAX.y*1.25))
points(pca$li[1:nInd(gen.tmp),1:2], pch=19, col=col[as.numeric(pop(pooled_gens)[1:nInd(gen.tmp)])])
legend(150, 180, legend=c("G. geiseri","G. nobilis"), col=col[1:2], pch=19, bty="n", cex=1, title=as.expression(bquote(bold("Field ID"))))
legend(150, 120, legend=c("F1 cross","G. geiseri BX","G. nobilis BX","F1 x F1"), col=col[3:7], pch=21, bty="n", cex=1, title = as.expression(bquote(bold("Simulation"))))
mtext("G. geiseri and G. nobilis (Diamond Y)", 3, adj=0.05, cex=1.5)
dev.off()

#Finding samples around the sims
pop_center <- vector()
tmp.range <- vector()
for(i in c("geisxnob", "geis_bx", "nob_bx", "F1XF1")){
pop_center <- c(pop_center, mean(pca$li[which(pop(pooled_gens) == i),1]))
tmp.low <- min(pca$li[which(pop(pooled_gens) == i),1])
tmp.high <- max(pca$li[which(pop(pooled_gens) == i),1])
tmp.range <- c(min(tmp.range,tmp.low), max(tmp.range,tmp.high))
}
names(pop_center) <- c("geisxnob", "geis_bx", "nob_bx", "F1XF1")

tmp.pca <- pca$li[1:nInd(gen.tmp),1:2]
for(i in which(tmp.pca[,1] > tmp.range[1]*1.5 & tmp.pca[,1] < tmp.range[2]*1.5)){
tmp.err <- (tmp.pca[i,1] - pop_center)^2
tmp.sim <- names(tmp.err)[which(tmp.err == min(tmp.err))]
print(paste(rownames(pca$li)[i],": ",tmp.sim,"    ","pos=",round(pca$li[i,1],3),", ",round(pca$li[i,2],3),"    ",tmp.sim,"=",round(pop_center[tmp.sim],3),sep=""))
}

#Samples matching the wrong side of the plot
rownames(pca$li)[pca$li[,1] > pop_center["nob_bx"] & rownames(pca$li) %in% indNames(geis)]
rownames(pca$li)[pca$li[,1] < pop_center["geis_bx"] & rownames(pca$li) %in% indNames(nob)]
}}}} #notepad cleanup

#G. geiseri and nobilis-NM
{```{R}```
#Define Strata
gen.tmp <- gen2[indNames(gen2)[gen2@strata$POP2 %in% c("Gag","Gan-NM")], , drop=T]
setPop(gen.tmp) <- ~ POP2
temp <- seppop(gen.tmp)
geis <- temp$`Gag`
nob <- temp$`Gan-NM`

#Test for missing loci
geis.miss <- apply(geis@tab, 2, function(x) sum(is.na(x))/nInd(geis))
nob.miss <- apply(nob@tab, 2, function(x) sum(is.na(x))/nInd(nob))

rm.loci <- NULL
if(length(which(geis.miss  == 1))){rm.loci <- c(rm.loci, names(which(geis.miss == 1)))}
if(length(which(nob.miss  == 1))){rm.loci <- c(rm.loci, names(which(nob.miss == 1)))}

if(length(rm.loci) > 0){
rm.loci <- unique(matrix(unlist(strsplit(rm.loci,"[.]")), ncol=2, byrow=T)[,1])
set.loc <- locNames(geis)[!locNames(geis) %in% rm.loci]
gen.vcf <- gen.vcf[, loc=set.loc]
geis <- geis[, loc=set.loc]
nob <- nob[, loc=set.loc]
}

nLoc(nob)

#Simulate samples
F1 <- hybridize(geis, nob, n = 100, pop = "geisxnob")
geis_bx <- hybridize(geis, F1, n = 100, pop = "geis_bx")
nob_bx <- hybridize(nob, F1, n = 100, pop = "nob_bx")
F2 <- hybridize(F1, F1, n= 100, pop = "F1XF1")
pooled_gens <- repool(gen.tmp, F1, geis_bx, nob_bx, F2)

save(pooled_gens, file="G_x_N_NM_sims.gz", compress=T)
#load("G_x_N_NM_sims.gz")

#PCA
x <- scaleGen(pooled_gens, NA.method = "mean")

#PCA with all loci
pca <- dudi.pca(x,center = FALSE, scale = FALSE, scannf = FALSE, nf = 4)

#plot2
col = c("#68a691", "#ddd392", "#ed6a5a", "#e59a2e", "#bdbdbd", "#636363", '#000000', '#799c41')
s.class(pca$li, pop(pooled_gens), xax=1, yax=2, cellipse=3, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0)
MIN.x <- par("usr")[1]
MAX.x <- par("usr")[2]
MIN.y <- par("usr")[3]
MAX.y <- par("usr")[4]

png("ade_hybrid/GgxGn_NM_sim.png", res=200, width=2000, height=2000)
s.class(pca$li[pop(pooled_gens) %in% c("geisxnob", "geis_bx", "nob_bx", "F1XF1"), 1:2], pop(pooled_gens)[pop(pooled_gens) %in% c("geisxnob", "geis_bx", "nob_bx", "F1XF1")], 
xax=1, yax=2, cellipse=10, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0, xlim=c(MIN.x,MAX.x), ylim=c(MIN.y*1.1,MAX.y*1.2))
points(pca$li[1:nInd(gen.tmp),1:2], pch=19, col=col[as.numeric(pop(pooled_gens)[1:nInd(gen.tmp)])])
dev.off()

#Finding samples around the sims
pop_center <- vector()
tmp.range <- vector()
for(i in c("geisxnob", "geis_bx", "nob_bx", "F1XF1")){
pop_center <- c(pop_center, mean(pca$li[which(pop(pooled_gens) == i),1]))
tmp.low <- min(pca$li[which(pop(pooled_gens) == i),1])
tmp.high <- max(pca$li[which(pop(pooled_gens) == i),1])
tmp.range <- c(min(tmp.range,tmp.low), max(tmp.range,tmp.high))
}
names(pop_center) <- c("geisxnob", "geis_bx", "nob_bx", "F1XF1")

tmp.pca <- pca$li[1:nInd(gen.tmp),1:2]
for(i in which(tmp.pca[,1] > tmp.range[1]*1.5 & tmp.pca[,1] < tmp.range[2]*1.5)){
tmp.err <- (tmp.pca[i,1] - pop_center)^2
tmp.sim <- names(tmp.err)[which(tmp.err == min(tmp.err))]
print(paste(rownames(pca$li)[i],": ",tmp.sim,"    ","pos=",round(pca$li[i,1],3),", ",round(pca$li[i,2],3),"    ",tmp.sim,"=",round(pop_center[tmp.sim],3),sep=""))
}

#Samples matching the wrong side of the plot
rownames(pca$li)[pca$li[,1] > pop_center["nob_bx"] & rownames(pca$li) %in% indNames(geis)]
rownames(pca$li)[pca$li[,1] < pop_center["geis_bx"] & rownames(pca$li) %in% indNames(nob)]
}}}} #notepad cleanup

#G. affinis and nobilis-WT
{```{R}```
#Define Strata
gen.tmp <- gen2[indNames(gen2)[gen2@strata$POP2 %in% c("Gaa","Gan-WT")], , drop=T]
setPop(gen.tmp) <- ~ POP2
temp <- seppop(gen.tmp)
affin <- temp$`Gaa`
nob <- temp$`Gan-WT`

#Test for missing loci
affin.miss <- apply(affin@tab, 2, function(x) sum(is.na(x))/nInd(affin))
nob.miss <- apply(nob@tab, 2, function(x) sum(is.na(x))/nInd(nob))

rm.loci <- NULL
if(length(which(affin.miss  == 1))){rm.loci <- c(rm.loci, names(which(affin.miss == 1)))}
if(length(which(nob.miss  == 1))){rm.loci <- c(rm.loci, names(which(nob.miss == 1)))}

if(length(rm.loci) > 0){
rm.loci <- unique(matrix(unlist(strsplit(rm.loci,"[.]")), ncol=2, byrow=T)[,1])
set.loc <- locNames(affin)[!locNames(affin) %in% rm.loci]
gen.vcf <- gen.vcf[, loc=set.loc]
affin <- affin[, loc=set.loc]
nob <- nob[, loc=set.loc]
}

nLoc(nob)

#Simulate samples
F1 <- hybridize(affin, nob, n = 100, pop = "affinxnob")
affin_bx <- hybridize(affin, F1, n = 100, pop = "affin_bx")
nob_bx <- hybridize(nob, F1, n = 100, pop = "nob_bx")
F2 <- hybridize(F1, F1, n= 100, pop = "F1XF1")
pooled_gens <- repool(gen.tmp, F1, affin_bx, nob_bx, F2)

save(pooled_gens, file="A_x_N_WT_sims.gz", compress=T)
#load("A_x_N_WT_sims.gz")

#PCA
x <- scaleGen(pooled_gens, NA.method = "mean")

#PCA with all loci
pca <- dudi.pca(x,center = FALSE, scale = FALSE, scannf = FALSE, nf = 4)

#plot2
col = c("#68a691", "#ddd392", "#ed6a5a", "#e59a2e", "#bdbdbd", "#636363", '#000000', '#799c41')
s.class(pca$li, pop(pooled_gens), xax=1, yax=2, cellipse=3, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0)
MIN.x <- par("usr")[1]
MAX.x <- par("usr")[2]
MIN.y <- par("usr")[3]
MAX.y <- par("usr")[4]

png("ade_hybrid/GaxGn_WT_sim.png", res=200, width=2000, height=2000)
s.class(pca$li[pop(pooled_gens) %in% c("affinxnob", "affin_bx", "nob_bx", "F1XF1"), 1:2], pop(pooled_gens)[pop(pooled_gens) %in% c("affinxnob", "affin_bx", "nob_bx", "F1XF1")], 
xax=1, yax=2, cellipse=10, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0, xlim=c(MIN.x,MAX.x), ylim=c(MIN.y*2.1,MAX.y*1.25))
points(pca$li[1:nInd(gen.tmp),1:2], pch=19, col=col[as.numeric(pop(pooled_gens)[1:nInd(gen.tmp)])])
dev.off()

#Finding samples around the sims
pop_center <- vector()
tmp.range <- vector()
for(i in c("affinxnob", "affin_bx", "nob_bx", "F1XF1")){
pop_center <- c(pop_center, mean(pca$li[which(pop(pooled_gens) == i),1]))
tmp.low <- min(pca$li[which(pop(pooled_gens) == i),1])
tmp.high <- max(pca$li[which(pop(pooled_gens) == i),1])
tmp.range <- c(min(tmp.range,tmp.low), max(tmp.range,tmp.high))
}
names(pop_center) <- c("affinxnob", "affin_bx", "nob_bx", "F1XF1")

tmp.pca <- pca$li[1:nInd(gen.tmp),1:2]
for(i in which(tmp.pca[,1] > tmp.range[1]*1.5 & tmp.pca[,1] < tmp.range[2]*1.5)){
tmp.err <- (tmp.pca[i,1] - pop_center)^2
tmp.sim <- names(tmp.err)[which(tmp.err == min(tmp.err))]
print(paste(rownames(pca$li)[i],": ",tmp.sim,"    ","pos=",round(pca$li[i,1],3),", ",round(pca$li[i,2],3),"    ",tmp.sim,"=",round(pop_center[tmp.sim],3),sep=""))
}

#Samples matching the wrong side of the plot
rownames(pca$li)[pca$li[,1] > pop_center["nob_bx"] & rownames(pca$li) %in% indNames(affin)]
rownames(pca$li)[pca$li[,1] < pop_center["affin_bx"] & rownames(pca$li) %in% indNames(nob)]
}}}} #notepad cleanup

#G. affinis and nobilis-DY
{```{R}```
#Define Strata
gen.tmp <- gen2[indNames(gen2)[gen2@strata$POP2 %in% c("Gaa","Gan-DY")], , drop=T]
setPop(gen.tmp) <- ~ POP2
temp <- seppop(gen.tmp)
affin <- temp$`Gaa`
nob <- temp$`Gan-DY`

#Test for missing loci
affin.miss <- apply(affin@tab, 2, function(x) sum(is.na(x))/nInd(affin))
nob.miss <- apply(nob@tab, 2, function(x) sum(is.na(x))/nInd(nob))

rm.loci <- NULL
if(length(which(affin.miss  == 1))){rm.loci <- c(rm.loci, names(which(affin.miss == 1)))}
if(length(which(nob.miss  == 1))){rm.loci <- c(rm.loci, names(which(nob.miss == 1)))}

if(length(rm.loci) > 0){
rm.loci <- unique(matrix(unlist(strsplit(rm.loci,"[.]")), ncol=2, byrow=T)[,1])
set.loc <- locNames(affin)[!locNames(affin) %in% rm.loci]
gen.vcf <- gen.vcf[, loc=set.loc]
affin <- affin[, loc=set.loc]
nob <- nob[, loc=set.loc]
}

nLoc(nob)

#Simulate samples
F1 <- hybridize(affin, nob, n = 100, pop = "affinxnob")
affin_bx <- hybridize(affin, F1, n = 100, pop = "affin_bx")
nob_bx <- hybridize(nob, F1, n = 100, pop = "nob_bx")
F2 <- hybridize(F1, F1, n= 100, pop = "F1XF1")
pooled_gens <- repool(gen.tmp, F1, affin_bx, nob_bx, F2)

save(pooled_gens, file="A_x_N_DY_sims.gz", compress=T)
#load("A_x_N_DY_sims.gz")

#PCA
x <- scaleGen(pooled_gens, NA.method = "mean")

#PCA with all loci
pca <- dudi.pca(x,center = FALSE, scale = FALSE, scannf = FALSE, nf = 4)

#plot2
col = c("#68a691", "#ddd392", "#ed6a5a", "#e59a2e", "#bdbdbd", "#636363", '#000000', '#799c41')
s.class(pca$li, pop(pooled_gens), xax=1, yax=2, cellipse=3, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0)
MIN.x <- par("usr")[1]
MAX.x <- par("usr")[2]
MIN.y <- par("usr")[3]
MAX.y <- par("usr")[4]

png("ade_hybrid/GaxGn_DY_sim.png", res=200, width=2000, height=2000)
s.class(pca$li[pop(pooled_gens) %in% c("affinxnob", "affin_bx", "nob_bx", "F1XF1"), 1:2], pop(pooled_gens)[pop(pooled_gens) %in% c("affinxnob", "affin_bx", "nob_bx", "F1XF1")], 
xax=1, yax=2, cellipse=10, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0, xlim=c(MIN.x,MAX.x), ylim=c(MIN.y*1.25,MAX.y*1.35))
points(pca$li[1:nInd(gen.tmp),1:2], pch=19, col=col[as.numeric(pop(pooled_gens)[1:nInd(gen.tmp)])])
points(pca$li[35,1:2], pch=19, col=col[as.numeric(pop(pooled_gens)[35])])
legend(150, 150, legend=c("G. affinis","G. nobilis"), col=col[1:2], pch=19, bty="n", cex=1, title=as.expression(bquote(bold("Field ID"))))
legend(150, 120, legend=c("F1 cross","G. affineri BX","G. nobilis BX","F1 x F1"), col=col[3:7], pch=21, bty="n", cex=1, title = as.expression(bquote(bold("Simulation"))))
#mtext("G. affinis and G. nobilis (Diamond Y)", 3, adj=0.05, cex=1.5)
dev.off()

#Finding samples around the sims
pop_center <- vector()
tmp.range <- vector()
for(i in c("affinxnob", "affin_bx", "nob_bx", "F1XF1")){
pop_center <- c(pop_center, mean(pca$li[which(pop(pooled_gens) == i),1]))
tmp.low <- min(pca$li[which(pop(pooled_gens) == i),1])
tmp.high <- max(pca$li[which(pop(pooled_gens) == i),1])
tmp.range <- c(min(tmp.range,tmp.low), max(tmp.range,tmp.high))
}
names(pop_center) <- c("affinxnob", "affin_bx", "nob_bx", "F1XF1")

tmp.pca <- pca$li[1:nInd(gen.tmp),1:2]
for(i in which(tmp.pca[,1] > tmp.range[1]*1.5 & tmp.pca[,1] < tmp.range[2]*1.5)){
tmp.err <- (tmp.pca[i,1] - pop_center)^2
tmp.sim <- names(tmp.err)[which(tmp.err == min(tmp.err))]
print(paste(rownames(pca$li)[i],": ",tmp.sim,"    ","pos=",round(pca$li[i,1],3),", ",round(pca$li[i,2],3),"    ",tmp.sim,"=",round(pop_center[tmp.sim],3),sep=""))
}

#Samples matching the wrong side of the plot
rownames(pca$li)[pca$li[,1] > pop_center["nob_bx"] & rownames(pca$li) %in% indNames(affin)]
rownames(pca$li)[pca$li[,1] < pop_center["affin_bx"] & rownames(pca$li) %in% indNames(nob)]
}}}} #notepad cleanup

#G.  affinis and nobilis-NM
{```{R}```
#Define Strata
gen.tmp <- gen2[indNames(gen2)[gen2@strata$POP2 %in% c("Gaa","Gan-NM")], , drop=T]
setPop(gen.tmp) <- ~ POP2
temp <- seppop(gen.tmp)
affin <- temp$`Gaa`
nob <- temp$`Gan-NM`

#Test for missing loci
affin.miss <- apply(affin@tab, 2, function(x) sum(is.na(x))/nInd(affin))
nob.miss <- apply(nob@tab, 2, function(x) sum(is.na(x))/nInd(nob))

rm.loci <- NULL
if(length(which(affin.miss  == 1))){rm.loci <- c(rm.loci, names(which(affin.miss == 1)))}
if(length(which(nob.miss  == 1))){rm.loci <- c(rm.loci, names(which(nob.miss == 1)))}

if(length(rm.loci) > 0){
rm.loci <- unique(matrix(unlist(strsplit(rm.loci,"[.]")), ncol=2, byrow=T)[,1])
set.loc <- locNames(affin)[!locNames(affin) %in% rm.loci]
gen.vcf <- gen.vcf[, loc=set.loc]
affin <- affin[, loc=set.loc]
nob <- nob[, loc=set.loc]
}

nLoc(nob)

#Simulate samples
F1 <- hybridize(affin, nob, n = 100, pop = "affinxnob")
affin_bx <- hybridize(affin, F1, n = 100, pop = "affin_bx")
nob_bx <- hybridize(nob, F1, n = 100, pop = "nob_bx")
F2 <- hybridize(F1, F1, n= 100, pop = "F1XF1")
pooled_gens <- repool(gen.tmp, F1, affin_bx, nob_bx, F2)

save(pooled_gens, file="A_x_N_NM_sims.gz", compress=T)
#load("A_x_N_NM_sims.gz")

#PCA
x <- scaleGen(pooled_gens, NA.method = "mean")

#PCA with all loci
pca <- dudi.pca(x,center = FALSE, scale = FALSE, scannf = FALSE, nf = 4)

#plot2
col = c("#68a691", "#ddd392", "#ed6a5a", "#e59a2e", "#bdbdbd", "#636363", '#000000', '#799c41')
s.class(pca$li, pop(pooled_gens), xax=1, yax=2, cellipse=3, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0)
MIN.x <- par("usr")[1]
MAX.x <- par("usr")[2]
MIN.y <- par("usr")[3]
MAX.y <- par("usr")[4]

png("ade_hybrid/GaxGn_NM_sim.png", res=200, width=2000, height=2000)
s.class(pca$li[pop(pooled_gens) %in% c("affinxnob", "affin_bx", "nob_bx", "F1XF1"), 1:2], pop(pooled_gens)[pop(pooled_gens) %in% c("affinxnob", "affin_bx", "nob_bx", "F1XF1")], 
xax=1, yax=2, cellipse=10, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0, xlim=c(MIN.x,MAX.x), ylim=c(MIN.y*1.6,MAX.y*1.25))
points(pca$li[1:nInd(gen.tmp),1:2], pch=19, col=col[as.numeric(pop(pooled_gens)[1:nInd(gen.tmp)])])
legend(150, 150, legend=c("G. affinis","G. nobilis"), col=col[1:2], pch=19, bty="n", cex=1, title=as.expression(bquote(bold("Field ID"))))
legend(150, 120, legend=c("F1 cross","G. affineri BX","G. nobilis BX","F1 x F1"), col=col[3:7], pch=21, bty="n", cex=1, title = as.expression(bquote(bold("Simulation"))))
mtext("G. affinis and G. nobilis (New Mexico)", 3, adj=0.05, cex=1.5)
dev.off()

#Finding samples around the sims
pop_center <- vector()
tmp.range <- vector()
for(i in c("affinxnob", "affin_bx", "nob_bx", "F1XF1")){
pop_center <- c(pop_center, mean(pca$li[which(pop(pooled_gens) == i),1]))
tmp.low <- min(pca$li[which(pop(pooled_gens) == i),1])
tmp.high <- max(pca$li[which(pop(pooled_gens) == i),1])
tmp.range <- c(min(tmp.range,tmp.low), max(tmp.range,tmp.high))
}
names(pop_center) <- c("affinxnob", "affin_bx", "nob_bx", "F1XF1")

tmp.pca <- pca$li[1:nInd(gen.tmp),1:2]
for(i in which(tmp.pca[,1] > tmp.range[1]*1.5 & tmp.pca[,1] < tmp.range[2]*1.5)){
tmp.err <- (tmp.pca[i,1] - pop_center)^2
tmp.sim <- names(tmp.err)[which(tmp.err == min(tmp.err))]
print(paste(rownames(pca$li)[i],": ",tmp.sim,"    ","pos=",round(pca$li[i,1],3),", ",round(pca$li[i,2],3),"    ",tmp.sim,"=",round(pop_center[tmp.sim],3),sep=""))
}

#Samples matching the wrong side of the plot
rownames(pca$li)[pca$li[,1] > pop_center["nob_bx"] & rownames(pca$li) %in% indNames(affin)]
rownames(pca$li)[pca$li[,1] < pop_center["affin_bx"] & rownames(pca$li) %in% indNames(nob)]
}}}} #notepad cleanup

### New Hybrids ###
# geiseri and nobilis-WT
{```{R}```
write.table(indNames(gen2)[gen2@strata$POP2 == "Gaa"], file="affinis_indv.txt", quote=F, col.names=F, row.names=F)
write.table(indNames(gen2)[gen2@strata$POP2 == "Gag"], file="geiseri_indv.txt", quote=F, col.names=F, row.names=F)
write.table(indNames(gen2)[gen2@strata$POP2 == "Gan-DY"], file="DY_indv.txt", quote=F, col.names=F, row.names=F)
write.table(indNames(gen2)[gen2@strata$POP2 == "Gan-NM"], file="NM_indv.txt", quote=F, col.names=F, row.names=F)
write.table(indNames(gen2)[gen2@strata$POP2 == "Gan-WT"], file="WTX_indv.txt", quote=F, col.names=F, row.names=F)

names.m <- matrix(unlist(strsplit(locNames(gen2.vcf),"_")), ncol=4, byrow=T)
names.df <- data.frame(CHROM=apply(names.m[,1:3],1,function(x) paste(x,collapse="_")), POS=names.m[,4])
write.table(names.df, file="gen2_loci", quote=F, col.names=F, row.names=F, sep="\t")
} #notepad cleanup
{```{bash}```
mkdir new_hybrids
mv *_indv.txt new_hybrids
cd new_hybrids

#Make affinis and nobilis pairs
cat affinis_indv.txt DY_indv.txt > Ga_DY.txt
cat affinis_indv.txt WTX_indv.txt > Ga_WTX.txt
cat affinis_indv.txt NM_indv.txt > Ga_NM.txt

#Make geiseri and nobilis pairs
cat geiseri_indv.txt DY_indv.txt > Gg_DY.txt
cat geiseri_indv.txt WTX_indv.txt > Gg_WTX.txt
cat geiseri_indv.txt NM_indv.txt > Gg_NM.txt

#Make affinis and nobilis pairs
cat affinis_indv.txt geiseri_indv.txt > Ga_Gg.txt

#Make vcf files
vcftools --vcf ../SNP.TRS.F09.vcf --out Ga_DY --recode --recode-INFO-all --keep Ga_DY.txt --positions ../gen2_loci
vcftools --vcf ../SNP.TRS.F09.vcf --out Ga_WTX --recode --recode-INFO-all --keep Ga_WTX.txt --positions ../gen2_loci
vcftools --vcf ../SNP.TRS.F09.vcf --out Ga_NM --recode --recode-INFO-all --keep Ga_NM.txt --positions ../gen2_loci

vcftools --vcf ../SNP.TRS.F09.vcf --out Gg_DY --recode --recode-INFO-all --keep Gg_DY.txt --positions ../gen2_loci
vcftools --vcf ../SNP.TRS.F09.vcf --out Gg_WTX --recode --recode-INFO-all --keep Gg_WTX.txt --positions ../gen2_loci
vcftools --vcf ../SNP.TRS.F09.vcf --out Gg_NM --recode --recode-INFO-all --keep Gg_NM.txt --positions ../gen2_loci

vcftools --vcf ../SNP.TRS.F09.vcf --out Ga_Gg --recode --recode-INFO-all --keep Ga_Gg.txt --positions ../gen2_loci

#Adding the characterizations for the analysis
nano ~/bin/genotype_frequencies.txt
#Pasted info without the #
#   5
# 1_Bx      0.00000    0.250000    0.250000   0.50000
# 0_Bx      0.50000    0.250000    0.250000   0.00000
# F1        0.00000     .5          .5        0.00000
# Pure_1    0.00000    0.00000     0.00000    1.00000
# Pure_0    1.00000    0.00000     0.00000    0.00000
} #notepad cleanup

## G. geseri vs New Mexico G. nobilis ##
{```{bash}```
#Generating data
{```{bash}```
cd /home/afields/Workspace/misc/Gambusia/analysis/May_2024/new_hybrids
#Preparing dir
mkdir -p geiseri/NM
cd geiseri/NM
for i in $(seq 1 5); do mkdir run_$i; done

for i in $(seq 1 5); do 
echo "Processing run $i"
cd run_$i
Rscript ~/bin/NewHybrids_vcf_to_txt.r /home/afields/Workspace/misc/Gambusia/analysis/May_2024/new_hybrids/Gg_NM.recode.vcf newhybrids_input_GxN_NM.txt
cd ..
done

#Adding header
for i in $(seq 1 5); do 
INDV=$(vcfsamplenames /home/afields/Workspace/misc/Gambusia/analysis/May_2024/new_hybrids/Gg_NM.recode.vcf | wc -l);
sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" run_$i/newhybrids_input_GxN_NM.txt;
done

#Running iterations
for i in $(seq 1 5); do
cd run_$i
RAND1=$(echo $((RANDOM)))
RAND2=$(echo $((RANDOM)))
echo $RAND1 $RAND2 > starting_seeds.txt
newhybrids -d newhybrids_input_GxN_NM.txt -c ~/bin/genotype_frequencies.txt -s $RAND1 $RAND2 --no-gui 2>&1 | tee run.log &
sleep 5
cd ..
done
} #notepad cleanup
#Formatting the output data
{```{bash}```
for i in $(seq 1 5); do
cd run_$i
vcfsamplenames ../../../Gg_NM.recode.vcf > tmp.indv
sed -i "1 s/^/Indv\n/" tmp.indv
paste tmp.indv aa-ScaledLikelihood.txt > aa-ScaledLikelihood_indv.txt
rm tmp.indv
Rscript ~/bin/NewHybrids_results_format.r G.geiseri ../../../geiseri_indv.txt G.nobilis ../../../NM_indv.txt
cd ..
done
} #notepad cleanup
#Putting all the data into an excel workbook
{```{R}```
library(openxlsx)
dat <- read.table("run_1/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat$Max <- apply(dat[,2:ncol(dat)],1,max)
wb <- write.xlsx(dat, "GxN_NM.xlsx", colNames = T, sheetName="run1")
dat2 <- read.table("run_2/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat2$Max <- apply(dat2[,2:ncol(dat2)],1,max)
addWorksheet(wb, sheetName = "run2")
writeData(wb, "run2", dat2)
dat3 <- read.table("run_3/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat3$Max <- apply(dat3[,2:ncol(dat3)],1,max)
addWorksheet(wb, sheetName = "run3")
writeData(wb, "run3", dat3)
dat4 <- read.table("run_4/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat4$Max <- apply(dat4[,2:ncol(dat4)],1,max)
addWorksheet(wb, sheetName = "run4")
writeData(wb, "run4", dat4)
dat5 <- read.table("run_5/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat5$Max <- apply(dat5[,2:ncol(dat5)],1,max)
addWorksheet(wb, sheetName = "run5")
writeData(wb, "run5", dat5)
addWorksheet(wb, sheetName = "Summary")
saveWorkbook(wb, "GxN_NM.xlsx", overwrite = TRUE)
q("no")
} #notepad cleanup
} #notepad cleanup

## G. geseri vs Diamond Y G. nobilis ##
{```{bash}```
#Generating data
{```{bash}```
#Preparing dir
mkdir ../DY
cd ../DY
for i in $(seq 1 5); do mkdir run_$i; done

for i in $(seq 1 5); do 
echo "Processing run $i"
cd run_$i
Rscript ~/bin/NewHybrids_vcf_to_txt.r /home/afields/Workspace/misc/Gambusia/analysis/May_2024/new_hybrids/Gg_DY.recode.vcf newhybrids_input_GxN_DY.txt
cd ..
done

#Adding header
for i in $(seq 1 5); do 
INDV=$(vcfsamplenames /home/afields/Workspace/misc/Gambusia/analysis/May_2024/new_hybrids/Gg_DY.recode.vcf | wc -l);
sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" run_$i/newhybrids_input_GxN_DY.txt;
done

#Running iterations
for i in $(seq 1 5); do
cd run_$i
RAND1=$(echo $((RANDOM)))
RAND2=$(echo $((RANDOM)))
echo $RAND1 $RAND2 > starting_seeds.txt
newhybrids -d newhybrids_input_GxN_DY.txt -c ~/bin/genotype_frequencies.txt -s $RAND1 $RAND2 --no-gui 2>&1 | tee run.log &
sleep 5
cd ..
done
} #notepad cleanup
#Formatting the output data
{```{bash}```
for i in $(seq 1 5); do
cd run_$i
vcfsamplenames ../../../Gg_DY.recode.vcf > tmp.indv
sed -i "1 s/^/Indv\n/" tmp.indv
paste tmp.indv aa-ScaledLikelihood.txt > aa-ScaledLikelihood_indv.txt
rm tmp.indv
Rscript ~/bin/NewHybrids_results_format.r G.geiseri ../../../geiseri_indv.txt G.nobilis ../../../DY_indv.txt
cd ..
done
} #notepad cleanup
#Putting all the data into an excel workbook
{```{R}```
library(openxlsx)
dat <- read.table("run_1/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat$Max <- apply(dat[,2:ncol(dat)],1,max)
wb <- write.xlsx(dat, "GxN_DY.xlsx", colNames = T, sheetName="run1")
dat2 <- read.table("run_2/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat2$Max <- apply(dat2[,2:ncol(dat2)],1,max)
addWorksheet(wb, sheetName = "run2")
writeData(wb, "run2", dat2)
dat3 <- read.table("run_3/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat3$Max <- apply(dat3[,2:ncol(dat3)],1,max)
addWorksheet(wb, sheetName = "run3")
writeData(wb, "run3", dat3)
dat4 <- read.table("run_4/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat4$Max <- apply(dat4[,2:ncol(dat4)],1,max)
addWorksheet(wb, sheetName = "run4")
writeData(wb, "run4", dat4)
dat5 <- read.table("run_5/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat5$Max <- apply(dat5[,2:ncol(dat5)],1,max)
addWorksheet(wb, sheetName = "run5")
writeData(wb, "run5", dat5)
addWorksheet(wb, sheetName = "Summary")
saveWorkbook(wb, "GxN_DY.xlsx", overwrite = TRUE)
q("no")
} #notepad cleanup
} #notepad cleanup

## G. geseri vs West Texas G. nobilis ##
{```{bash}```
#Generating data
{```{bash}```
#Preparing dir
mkdir ../WTX
cd ../WTX
for i in $(seq 1 5); do mkdir run_$i; done

for i in $(seq 1 5); do 
echo "Processing run $i"
cd run_$i
Rscript ~/bin/NewHybrids_vcf_to_txt.r /home/afields/Workspace/misc/Gambusia/analysis/May_2024/new_hybrids/Gg_WTX.recode.vcf newhybrids_input_GxN_WTX.txt
cd ..
done

#Adding header
for i in $(seq 1 5); do 
INDV=$(vcfsamplenames /home/afields/Workspace/misc/Gambusia/analysis/May_2024/new_hybrids/Gg_WTX.recode.vcf | wc -l);
sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" run_$i/newhybrids_input_GxN_WTX.txt;
done

#Running iterations
for i in $(seq 1 5); do
cd run_$i
RAND1=$(echo $((RANDOM)))
RAND2=$(echo $((RANDOM)))
echo $RAND1 $RAND2 > starting_seeds.txt
newhybrids -d newhybrids_input_GxN_WTX.txt -c ~/bin/genotype_frequencies.txt -s $RAND1 $RAND2 --no-gui 2>&1 | tee run.log &
sleep 5
cd ..
done
} #notepad cleanup
#Formatting the output data
{```{bash}```
for i in $(seq 1 5); do
cd run_$i
vcfsamplenames ../../../Gg_WTX.recode.vcf > tmp.indv
sed -i "1 s/^/Indv\n/" tmp.indv
paste tmp.indv aa-ScaledLikelihood.txt > aa-ScaledLikelihood_indv.txt
rm tmp.indv
Rscript ~/bin/NewHybrids_results_format.r G.geiseri ../../../geiseri_indv.txt G.nobilis ../../../WTX_indv.txt
cd ..
done
} #notepad cleanup
#Putting all the data into an excel workbook
{```{R}```
library(openxlsx)
dat <- read.table("run_1/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat$Max <- apply(dat[,2:ncol(dat)],1,max)
wb <- write.xlsx(dat, "GxN_WTX.xlsx", colNames = T, sheetName="run1")
dat2 <- read.table("run_2/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat2$Max <- apply(dat2[,2:ncol(dat2)],1,max)
addWorksheet(wb, sheetName = "run2")
writeData(wb, "run2", dat2)
dat3 <- read.table("run_3/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat3$Max <- apply(dat3[,2:ncol(dat3)],1,max)
addWorksheet(wb, sheetName = "run3")
writeData(wb, "run3", dat3)
dat4 <- read.table("run_4/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat4$Max <- apply(dat4[,2:ncol(dat4)],1,max)
addWorksheet(wb, sheetName = "run4")
writeData(wb, "run4", dat4)
dat5 <- read.table("run_5/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat5$Max <- apply(dat5[,2:ncol(dat5)],1,max)
addWorksheet(wb, sheetName = "run5")
writeData(wb, "run5", dat5)
addWorksheet(wb, sheetName = "Summary")
saveWorkbook(wb, "GxN_WTX.xlsx", overwrite = TRUE)
q("no")
} #notepad cleanup
} #notepad cleanup

## G. affinis vs New Mexico G. nobilis ##
{```{bash}```
#Generating data
{```{bash}```
cd /home/afields/Workspace/misc/Gambusia/analysis/May_2024/new_hybrids
#Preparing dir
mkdir -p affinis/NM
cd affinis/NM
for i in $(seq 1 5); do mkdir run_$i; done

for i in $(seq 1 5); do 
echo "Processing run $i"
cd run_$i
Rscript ~/bin/NewHybrids_vcf_to_txt.r /home/afields/Workspace/misc/Gambusia/analysis/May_2024/new_hybrids/Ga_NM.recode.vcf newhybrids_input_AxN_NM.txt
cd ..
done

#Adding header
for i in $(seq 1 5); do 
INDV=$(vcfsamplenames /home/afields/Workspace/misc/Gambusia/analysis/May_2024/new_hybrids/Ga_NM.recode.vcf | wc -l);
sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" run_$i/newhybrids_input_AxN_NM.txt;
done

#Running iterations
for i in $(seq 1 5); do
cd run_$i
RAND1=$(echo $((RANDOM)))
RAND2=$(echo $((RANDOM)))
echo $RAND1 $RAND2 > starting_seeds.txt
newhybrids -d newhybrids_input_AxN_NM.txt -c ~/bin/genotype_frequencies.txt -s $RAND1 $RAND2 --no-gui 2>&1 | tee run.log &
sleep 5
cd ..
done
} #notepad cleanup
#Formatting the output data
{```{bash}```
for i in $(seq 1 5); do
cd run_$i
vcfsamplenames ../../../Ga_NM.recode.vcf > tmp.indv
sed -i "1 s/^/Indv\n/" tmp.indv
paste tmp.indv aa-ScaledLikelihood.txt > aa-ScaledLikelihood_indv.txt
rm tmp.indv
Rscript ~/bin/NewHybrids_results_format.r G.affinis ../../../affinis_indv.txt G.nobilis ../../../NM_indv.txt
cd ..
done
} #notepad cleanup
#Putting all the data into an excel workbook
{```{R}```
library(openxlsx)
dat <- read.table("run_1/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat$Max <- apply(dat[,2:ncol(dat)],1,max)
wb <- write.xlsx(dat, "AxN_NM.xlsx", colNames = T, sheetName="run1")
dat2 <- read.table("run_2/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat2$Max <- apply(dat2[,2:ncol(dat2)],1,max)
addWorksheet(wb, sheetName = "run2")
writeData(wb, "run2", dat2)
dat3 <- read.table("run_3/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat3$Max <- apply(dat3[,2:ncol(dat3)],1,max)
addWorksheet(wb, sheetName = "run3")
writeData(wb, "run3", dat3)
dat4 <- read.table("run_4/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat4$Max <- apply(dat4[,2:ncol(dat4)],1,max)
addWorksheet(wb, sheetName = "run4")
writeData(wb, "run4", dat4)
dat5 <- read.table("run_5/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat5$Max <- apply(dat5[,2:ncol(dat5)],1,max)
addWorksheet(wb, sheetName = "run5")
writeData(wb, "run5", dat5)
addWorksheet(wb, sheetName = "Summary")
saveWorkbook(wb, "AxN_NM.xlsx", overwrite = TRUE)
q("no")
} #notepad cleanup
} #notepad cleanup

## G. affinis vs Diamond Y G. nobilis ##
{```{bash}```
#Generating data
{```{bash}```
#Preparing dir
mkdir -p ../DY
cd ../DY
for i in $(seq 1 5); do mkdir run_$i; done

for i in $(seq 1 5); do 
echo "Processing run $i"
cd run_$i
Rscript ~/bin/NewHybrids_vcf_to_txt.r /home/afields/Workspace/misc/Gambusia/analysis/May_2024/new_hybrids/Ga_DY.recode.vcf newhybrids_input_AxN_DY.txt
cd ..
done

#Adding header
for i in $(seq 1 5); do 
INDV=$(vcfsamplenames /home/afields/Workspace/misc/Gambusia/analysis/May_2024/new_hybrids/Ga_DY.recode.vcf | wc -l);
sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" run_$i/newhybrids_input_AxN_DY.txt;
done

#Running iterations
for i in $(seq 1 5); do
cd run_$i
RAND1=$(echo $((RANDOM)))
RAND2=$(echo $((RANDOM)))
echo $RAND1 $RAND2 > starting_seeds.txt
newhybrids -d newhybrids_input_AxN_DY.txt -c ~/bin/genotype_frequencies.txt -s $RAND1 $RAND2 --no-gui 2>&1 | tee run.log &
sleep 5
cd ..
done
} #notepad cleanup
#Formatting the output data
{```{bash}```
for i in $(seq 1 5); do
cd run_$i
vcfsamplenames ../../../Ga_DY.recode.vcf > tmp.indv
sed -i "1 s/^/Indv\n/" tmp.indv
paste tmp.indv aa-ScaledLikelihood.txt > aa-ScaledLikelihood_indv.txt
rm tmp.indv
Rscript ~/bin/NewHybrids_results_format.r G.affinis ../../../affinis_indv.txt G.nobilis ../../../DY_indv.txt
cd ..
done
} #notepad cleanup
#Putting all the data into an excel workbook
{```{R}```
library(openxlsx)
dat <- read.table("run_1/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat$Max <- apply(dat[,2:ncol(dat)],1,max)
wb <- write.xlsx(dat, "AxN_DY.xlsx", colNames = T, sheetName="run1")
dat2 <- read.table("run_2/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat2$Max <- apply(dat2[,2:ncol(dat2)],1,max)
addWorksheet(wb, sheetName = "run2")
writeData(wb, "run2", dat2)
dat3 <- read.table("run_3/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat3$Max <- apply(dat3[,2:ncol(dat3)],1,max)
addWorksheet(wb, sheetName = "run3")
writeData(wb, "run3", dat3)
dat4 <- read.table("run_4/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat4$Max <- apply(dat4[,2:ncol(dat4)],1,max)
addWorksheet(wb, sheetName = "run4")
writeData(wb, "run4", dat4)
dat5 <- read.table("run_5/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat5$Max <- apply(dat5[,2:ncol(dat5)],1,max)
addWorksheet(wb, sheetName = "run5")
writeData(wb, "run5", dat5)
addWorksheet(wb, sheetName = "Summary")
saveWorkbook(wb, "AxN_DY.xlsx", overwrite = TRUE)
q("no")
} #notepad cleanup
} #notepad cleanup

## G. affinis vs West Texas G. nobilis ##
{```{bash}```
#Generating data
{```{bash}```
#Preparing dir
mkdir -p ../WTX
cd ../WTX
for i in $(seq 1 5); do mkdir run_$i; done

for i in $(seq 1 5); do 
echo "Processing run $i"
cd run_$i
Rscript ~/bin/NewHybrids_vcf_to_txt.r /home/afields/Workspace/misc/Gambusia/analysis/May_2024/new_hybrids/Ga_WTX.recode.vcf newhybrids_input_AxN_WTX.txt
cd ..
done

#Adding header
for i in $(seq 1 5); do 
INDV=$(vcfsamplenames /home/afields/Workspace/misc/Gambusia/analysis/May_2024/new_hybrids/Ga_WTX.recode.vcf | wc -l);
sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" run_$i/newhybrids_input_AxN_WTX.txt;
done

#Running iterations
for i in $(seq 1 5); do
cd run_$i
RAND1=$(echo $((RANDOM)))
RAND2=$(echo $((RANDOM)))
echo $RAND1 $RAND2 > starting_seeds.txt
newhybrids -d newhybrids_input_AxN_WTX.txt -c ~/bin/genotype_frequencies.txt -s $RAND1 $RAND2 --no-gui 2>&1 | tee run.log &
sleep 5
cd ..
done
} #notepad cleanup
#Formatting the output data
{```{bash}```
for i in $(seq 1 5); do
cd run_$i
vcfsamplenames ../../../Ga_WTX.recode.vcf > tmp.indv
sed -i "1 s/^/Indv\n/" tmp.indv
paste tmp.indv aa-ScaledLikelihood.txt > aa-ScaledLikelihood_indv.txt
rm tmp.indv
Rscript ~/bin/NewHybrids_results_format.r G.affinis ../../../affinis_indv.txt G.nobilis ../../../WTX_indv.txt
cd ..
done
} #notepad cleanup
#Putting all the data into an excel workbook
{```{R}```
library(openxlsx)
dat <- read.table("run_1/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat$Max <- apply(dat[,2:ncol(dat)],1,max)
wb <- write.xlsx(dat, "AxN_WTX.xlsx", colNames = T, sheetName="run1")
dat2 <- read.table("run_2/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat2$Max <- apply(dat2[,2:ncol(dat2)],1,max)
addWorksheet(wb, sheetName = "run2")
writeData(wb, "run2", dat2)
dat3 <- read.table("run_3/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat3$Max <- apply(dat3[,2:ncol(dat3)],1,max)
addWorksheet(wb, sheetName = "run3")
writeData(wb, "run3", dat3)
dat4 <- read.table("run_4/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat4$Max <- apply(dat4[,2:ncol(dat4)],1,max)
addWorksheet(wb, sheetName = "run4")
writeData(wb, "run4", dat4)
dat5 <- read.table("run_5/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat5$Max <- apply(dat5[,2:ncol(dat5)],1,max)
addWorksheet(wb, sheetName = "run5")
writeData(wb, "run5", dat5)
addWorksheet(wb, sheetName = "Summary")
saveWorkbook(wb, "AxN_WTX.xlsx", overwrite = TRUE)
q("no")
} #notepad cleanup
} #notepad cleanup

#Picking samples to use in the species tree
{```{R}```
ind.v <- NULL
ind.v <- c(ind.v, indNames(gen2)[gen2@strata$POP == "Gaa-GSW1"])
ind.v <- c(ind.v, sample(indNames(gen2)[gen2@strata$POP %in% c("Gaa-LCCSP", "Gaf-LCCSP")],3))
ind.v <- c(ind.v, sample(indNames(gen2)[gen2@strata$POP == "Gag-CHS"],2))
ind.v <- c(ind.v, sample(indNames(gen2)[gen2@strata$POP == "Gag-ES"],2))
ind.v <- c(ind.v, sample(indNames(gen2)[gen2@strata$POP == "Gag-SMR"],2))

write.table(ind.v, file="Tree_samples.txt", col.names=F, row.names=F, quote=F)
} #notepad cleanup
