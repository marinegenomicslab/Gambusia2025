### Gambusia Filtering ###

#Filters allelic balance, quality vs depth, strand representation and paired read representation
{```{bash}```
dDocent_filters SNP.TRS.QC.recode.vcf SNP.TRS.dDocent
#no
#100000
}

#Filtering singleton and doubleton loci for depth (Singletons >20 reads; Doubletons > 10 reads)
{```{bash}```
vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf --out out --singletons

awk ' $3=="S" {print $1, $2}' out.singletons > sing.loci
awk ' $3=="D" {print $1, $2}' out.singletons > doub.loci

vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf --out SNP.TRS.F01a --recode --recode-INFO-all --exclude-positions sing.loci
vcftools --vcf SNP.TRS.F01a.recode.vcf --out SNP.TRS.F01b --recode --recode-INFO-all --exclude-positions doub.loci

vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf --out SNP.TRS.F01.sing --recode --recode-INFO-all --positions sing.loci
vcftools --vcf SNP.TRS.F01.sing.recode.vcf --out SNP.TRS.F02.sing --recode --recode-INFO-all --minDP 20

vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf --out SNP.TRS.F01.doub --recode --recode-INFO-all --positions doub.loci
vcftools --vcf SNP.TRS.F01.doub.recode.vcf --out SNP.TRS.F02.doub --recode --recode-INFO-all --minDP 10

#No sites were found to have limited depth, therefore no filtering was applied
}

#Filter loci that have high variation in depth across a locus with an individual
{```{bash}```
vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf --out out --geno-depth
}
{```{R}```
gdepth<-read.table(file="out.gdepth", head=T)
gdepth[gdepth==-1]<-NA

for (i in 3:dim(gdepth)[2]) {
temp<-aggregate(gdepth[,i],by=list(gdepth[,1]), sd)
if(i==3){indv.site.sd<-data.frame(temp,row.names=1)} else
{indv.site.sd[,(i-2)]<-temp[,2]}
}
colnames(indv.site.sd)<-colnames(gdepth[3:dim(gdepth)[2]])
tmp<-apply(indv.site.sd, 1, mean, na.rm=T)
tmp2<-unique(c(names(which(tmp>25))))
length(tmp)
length(tmp2)
write.table(tmp2,file="bad.loci.sd", quote=F, col.names=F, row.names=F)
}} #Notepad cleanup
{```{bash}```
grep "dDocent" SNP.TRS.F02.vcf | cut -f 1,2 | uniq | tail -n +2 > contigs.txt
grep -wf bad.loci.sd contigs.txt > bad.loci
vcftools --vcf SNP.TRS.F02.vcf  --out SNP.TRS.F03 --exclude-positions bad.loci --recode-INFO-all --recode

charts.sh SNP.TRS.F03.recode.vcf &
}

#Finding monomorphic loci
{```{R}```
#Loading Libraries
{```{R}```
library('adegenet')
library('vcfR') 
} #notepad cleanup

#Importing Data
{```{R}```
strata <- read.csv(file = "SNP.TRS.F03.2024-04-25/indv.csv", header = TRUE)
strata$Species <- matrix(unlist(strsplit(as.character(as.matrix(strata$POP)), "-")), ncol=2, byrow=T)[,1]
strata$Species[strata$Species == "Gaf"] <- "Gaa"
strata$Species[strata$Species == "Gno"] <- "Gan"
strata$Species[strata$Species == "Gge"] <- "Gag"

vcf <- read.vcfR(file="SNP.TRS.F03.recode.vcf")
gen.vcf<-vcfR2genind(vcf)
strata(gen.vcf) <- strata[match(indNames(gen.vcf),strata$INDV),]
head(strata(gen.vcf))
rm(vcf)
} #notepad cleanup

#Recording monomorphics
{```{R}```
mono.loc  <- which(gen.vcf@loc.n.all == 1)

set.loc <- locNames(gen.vcf)[which(!locNames(gen.vcf) %in% names(mono.loc))]
gen2.vcf <- gen.vcf[, loc=set.loc, drop=T]

write.table(mono.loci, file="mono_loci.txt")
} #notepad cleanup

#Recoding individuals with high missing data
{```{R}```
tmp.indv <- gen.vcf@strata$INDV[gen.vcf@strata$F_MISS > 0.99]
write.table(tmp.indv, file="high_miss_indv.txt")
}} #notepad cleanup

#Removing monomorphic loci and individuals and loci with high missing data
{```{bash}```
vcftools --vcf SNP.TRS.F03.recode.vcf  --out SNP.TRS.F03a --exclude-positions mono_loci.txt --recode-INFO-all --recode

vcftools --vcf SNP.TRS.F03a.recode.vcf  --out SNP.TRS.F03b --exclude high_miss_indv.txt --recode-INFO-all --recode

vcftools --vcf SNP.TRS.F03b.recode.vcf  --out SNP.TRS.F04 --max-missing 0.8 --recode-INFO-all --recode
}

#Removing outgroups
{```{bash}```
cp ../filtering_Apr2024/outgroup.txt .
vcftools --vcf SNP.TRS.F4.recode.vcf  --out SNP.TRS.F05 --exclude outgroup.txt --recode-INFO-all --recode
}

#Filter data using decision tree
{```{bash}```
mkdir SOL
cp SNP.TRS.F05.recode.vcf SOL/
cd SOL
SOL.filter1.no_mac.sh SNP.TRS.F05.recode.vcf

cp vcf/B.1.3.3.2.3.2.1.2.SNP.finalb.1.recode.vcf ../
cd ..
charts.sh B.1.3.3.2.3.2.1.2.SNP.finalb.1.recode.vcf
}

#Removing loci with super high depth or over 50% heterozygosity
{```{R}```
dat <- read.table("B.1.3.3.2.3.2.1.2.SNP.finalb.1.2024-05-07/loci.csv", head=T, sep=",")
hist(dat$MEAN_DEPTH, breaks=100, col="red4")
#The central portion of this graphic looks like it should be symetrical

#finding the highest point in the density plot
peak <- which.max(density(dat$MEAN_DEPTH[dat$MEAN_DEPTH > 150])$y)
density(dat$MEAN_DEPTH[dat$MEAN_DEPTH > 150])$x[peak]
abline(v = density(dat$MEAN_DEPTH[dat$MEAN_DEPTH > 150])$x[peak])

#Minimum depth
min(dat$MEAN_DEPTH)

#New cutoff
cutoff <- 2*density(dat$MEAN_DEPTH[dat$MEAN_DEPTH > 150])$x[peak] - min(dat$MEAN_DEPTH)
cutoff
abline(v = cutoff)

write.table(dat[dat$MEAN_DEPTH > cutoff, c("CHROM","POS")], file="high_depth_loci.txt", quote=F, col.names=F, row.names=F, sep="\t")

hist(dat$PER_HET, breaks=100, col="red4")
write.table(dat[which(dat$PER_HET > 0.5), c("CHROM","POS")], file="high_het_loci.txt", quote=F, col.names=F, row.names=F, sep="\t")
}
{```{bash}```
cat high_depth_loci.txt high_het_loci.txt | sort -u > high_loci.txt
vcftools --vcf B.1.3.3.2.3.2.1.2.SNP.finalb.1.recode.vcf --out SNP.TRS.F06 --recode --recode-INFO-all --exclude-positions high_loci.txt

charts.sh SNP.TRS.F06.recode.vcf
}

#Filtering individuals with > 10% missing data
{```{bash}```
awk -v FS="," '$7 > 0.1 {print $1}' SNP.TRS.F06.2024-05-07/indv.csv | sed 's/"//g' > lowDP.indv
vcftools --vcf SNP.TRS.F06.recode.vcf --out SNP.TRS.F07 --remove lowDP.indv --recode --recode-INFO-all

charts.sh SNP.TRS.F07.recode.vcf
}

#MAF filtering
{```{R}```
#Load Libraries
{```{R}```
library('vcfR')
library('adegenet')
library('dartR')
library('related')
}

#Load data
{```{R}```
#Load sample information
strata <- read.csv(file = "SNP.TRS.F06.2024-05-17/indv.csv", header = TRUE)
#Adding nobilis grouping inforrmation
tmp.v <- as.character(as.matrix(strata$POP))
tmp.v[tmp.v %in% c("Gaa-GSW1", "Gaa-HEAD", "Gaa-LCCSP", "Gaf-LCCSP")] <- "Gaa"
tmp.v[tmp.v %in% c("Gag-CHS", "Gag-ES","Gag-SMR")] <- "Gag"
tmp.v[tmp.v %in% c("Gan-CHS", "Gan-ES", "Gan-PHL", "Gan-BAL", "Gno-BAL")] <- "Gan-WT"
tmp.v[tmp.v %in% c("Gan-EU", "Gan-FLS", "Gan-HEAD", "Gan-DY", "Gno-DY")] <- "Gan-DY"
tmp.v[tmp.v %in% c("Gan-Sink27", "Gan-Sink31", "Gan-Sink37", "Gan-Sink7", "Gan-BC", "Gan-BLNWR", "Gno-BLNWR", "Gno-BC")] <- "Gan-NM"
strata$POP2 <- as.factor(tmp.v)

#Checking the sample information
head(strata)

#Loading the vcf file
vcf <- read.vcfR(file="SNP.TRS.F07.recode.vcf")
gen.vcf<-vcfR2genind(vcf)
#Assing the sample information to the genind object
strata(gen.vcf) <- strata[match(indNames(gen.vcf),strata$INDV),]
#Checking the sample information of the genind
head(strata(gen.vcf))
}

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
set.indv <- indNames(gen.vcf)[!indNames(gen.vcf) %in% unlist(rm.list)]
gen2.vcf <- gen.vcf[set.indv, , drop=T]
}} #notepad cleanup

#Looking for minor alleles in each group
{```{R}```
setPop(gen2.vcf) <- ~POP2
ESU <- seppop(gen2.vcf)

minor.all <- list()
for(i in 1:length(ESU)){
gen.tmp <- ESU[[i]]
tmp.v <- apply(gen.tmp@tab, 2, function(x) sum(x,na.rm=T)/(2*length(which(!is.na(x)))))
minor.all[[i]] <- names(tmp.v)[tmp.v < 0.01]
}
minor.all <- lapply(minor.all, function(x) x[!is.na(x)])

lapply(minor.all, length)
}

#Getting the overlap between lists of minor alleles
{```{R}```
x <- minor.all
A <- x[[1]]
B <- x[[2]]
C <- x[[3]]
D <- x[[4]]
E <- x[[5]]
F <- x[[6]]
ncd <- intersect(C, D)
ncde <- intersect(ncd, E)
ncdea <- intersect(ncde, A)
ncdeb <- intersect(ncde, B)
ncdeab <- intersect(ncdea, B)

min_all.m <- matrix(unlist(strsplit(ncdeab, "[.]")), ncol=2, byrow=T)
loc.parts <- matrix(unlist(strsplit(min_all.m[,1], "_")), ncol=4, byrow=T)
loc.df <- data.frame(CHROM=apply(loc.parts[,1:3],1,function(x) paste(x, collapse="_")), POS=loc.parts[,4], Allele=min_all.m[,2])

vcf.rec <- vcf

#Replacing low maf values with missing data indicators
for(i in 1:nrow(loc.df)){
#Record of progress
if(i/50 == round(i/50,0)){print(paste("Processing locus number", i, "out of", length(ncdeab)))}
#Pulling Locus
tmp.loc <- which(vcf@fix[,"CHROM"] == loc.df$CHROM[i] & vcf@fix[,"POS"] == loc.df$POS[i])
#Pulling genotype calls
tmp.gt <- vcf@gt[tmp.loc,]
#Seperating the vcf information for each individual
tmp.list <- strsplit(tmp.gt,":")
#Replacing minor alleles
for(j in 2:length(tmp.list)){if(length(grep(loc.df$Allele[i], tmp.list[[j]][1])) > 0){tmp.list[[j]][1] <- ".|."}}
#Reassembling vcf information
tmp.gt <- lapply(tmp.list, function(x) paste(x,collapse=":"))
#Putting it back in the vcf object
vcf@gt[tmp.loc,] <- unlist(tmp.gt)[match(colnames(vcf@gt), names(tmp.gt))]
}

write.vcf(vcf, file="SNP.TRS.F08.vcf.gz")
}}} #notepad cleanup
}

#Haplotyping data
{```{bash}```
cd $WORK/Workspace/Gambusia/haplotyping
cp -s ../reference/ref_subset.fasta ./reference.fasta
cp -s ../pear_remap/Lib*/*-RG.bam .
ls *-RG.bam | while read i; do echo "Indexing $i"; samtools index $i; done
cp -s ../pear_SNPs/popmap .
cp ../filtering/SNP.TRS.F08.vcf .

#Haplotyping
sbatch $WORK/slurm/haplotyping.slurm
}

{```{bash}```
rsync -av haplotyping afields@hactar.ad.tamucc.edu:/home/afields/Workspace/misc/Gambusia/filtering_May2024
charts.sh SNP.TRS.F08.vcf
}
