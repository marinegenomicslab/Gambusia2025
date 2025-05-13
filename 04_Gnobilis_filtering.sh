### Gambusia nobilis filtering ###

#Selecting only the G. nobilis samples
{```{bash}```
vcftools --vcf TotalRawSNPs.vcf --out SNP.TRS.F01 --keep G.nobilis_samples.txt --recode-INFO-all --recode
} #notepad cleanup

#Filters allelic balance, quality vs depth, strand representation and paired read representation
{```{bash}```
dDocent_filters SNP.TRS.F01 SNP.TRS.dDocent
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

rm out.singletons
#No sites were found to have limited depth, therefore no filtering was applied
} #notepad cleanup

#Filter loci that have high variation in depth across a locus with an individual
{```{bash}```
vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf --out out --geno-depth
} #notepad cleanup
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
q("no")
}} #notepad cleanup
{```{bash}```
grep "dDocent" SNP.TRS.dDocent.FIL.recode.vcf | cut -f 1,2 | uniq | tail -n +2 > contigs.txt
grep -wf bad.loci.sd contigs.txt > bad.loci
vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf  --out SNP.TRS.F03 --exclude-positions bad.loci --recode-INFO-all --recode

charts.sh SNP.TRS.F03.recode.vcf &
} #notepad cleanup

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
vcftools --vcf SNP.TRS.F03b.recode.vcf  --out SNP.TRS.F04 --max-missing 0.8 --recode-INFO-all --recode
} #notepad cleanup

#Filter data using decision tree
{```{bash}```
mkdir SOL
cp SNP.TRS.F03.recode.vcf SOL/
cd SOL
SOL.filter1.no_mac.sh SNP.TRS.F03.recode.vcf

cp vcf/B.1.3.3.2.3.2.1.2.SNP.finalb.1.recode.vcf ../
cd ..
charts.sh B.1.3.3.2.3.2.1.2.SNP.finalb.1.recode.vcf
} #notepad cleanup

#Removing loci with super high depth
{```{R}```
dat <- read.table("B.1.3.3.2.3.2.1.2.SNP.finalb.1.2024-06-01/loci.csv", head=T, sep=",")
hist(dat$MEAN_DEPTH, breaks=100, col="red4")
#The central portion of this graphic looks like it should be symetrical

#finding the highest point in the density plot
peak <- which.max(density(dat$MEAN_DEPTH)$y)
density(dat$MEAN_DEPTH)$x[peak]
abline(v = density(dat$MEAN_DEPTH)$x[peak])

#Minimum depth
min(dat$MEAN_DEPTH)

#New cutoff
cutoff <- 2*density(dat$MEAN_DEPTH)$x[peak] - min(dat$MEAN_DEPTH)
abline(v = cutoff)

write.table(dat[dat$MEAN_DEPTH > cutoff, c("CHROM","POS")], file="high_depth_loci.txt", quote=F, col.names=F, row.names=F, sep="\t")
} #notepad cleanup
{```{bash}```
vcftools --vcf B.1.3.3.2.3.2.1.2.SNP.finalb.1.recode.vcf --out SNP.TRS.F06 --recode --recode-INFO-all --exclude-positions high_depth_loci.txt

charts.sh SNP.TRS.F06.recode.vcf
} #notepad cleanup

#Filtering individuals with > 10% missing data
{```{bash}```
awk -v FS="," '$7 > 0.1 {print $1}' SNP.TRS.F06.2024-06-01/indv.csv | sed 's/"//g' > lowDP.indv
vcftools --vcf SNP.TRS.F06.recode.vcf --out SNP.TRS.F07 --remove lowDP.indv --recode --recode-INFO-all

charts.sh SNP.TRS.F07.recode.vcf
} #notepad cleanup

#MAF filtering
#Did not actually filter based upon this yet
{```{R}```
#Load Libraries
library('vcfR')
library('adegenet')
library('dartR')
library('related')
library('VennDiagram')

#Load data
strata <- read.csv(file = "SNP.TRS.F06.2024-06-01/indv.csv", header = TRUE)
strata$POP[is.na(strata$POP)] <- "Gan-CHS"
strata$Species <- matrix(unlist(strsplit(as.character(as.matrix(strata$POP)), "-")),ncol=2, byrow=T)[,1]
head(strata)

vcf <- read.vcfR(file="SNP.TRS.F07.recode.vcf")
gen.vcf<-vcfR2genind(vcf)
strata(gen.vcf) <- strata[match(indNames(gen.vcf),strata$INDV),]
head(strata(gen.vcf))

ESU <- seppop(gen.vcf)

minor.all <- list()
for(i in 1:length(ESU)){
gen.tmp <- ESU[[i]]
tmp.v <- apply(gen.tmp@tab, 2, function(x) sum(x,na.rm=T)/(2*length(which(!is.na(x)))))
minor.all[[i]] <- names(tmp.v)[tmp.v < 0.01]
}
minor.all <- lapply(minor.all, function(x) x[!is.na(x)])

lapply(minor.all, length)

POP.tab <- lapply(species, function(x) table(x$strata$POP))
list.names <- unlist(lapply(POP.tab, function(x) names(x)[x == max(x)]))
venn.diagram(minor.all,
category.names=c(
paste(list.names[1],"\nloci\n(n=", length(minor.all[[1]]), ")", sep = ""),
paste(list.names[2],"\nloci\n(n=", length(minor.all[[2]]), ")", sep = ""),
paste(list.names[3],"\nloci\n(n=", length(minor.all[[3]]), ")", sep = "")
), "Venn_min_all.tiff", main="Minor alleles by species", fill=c("red4", "mediumblue", "darkgreen"), alpha=rep(0.5,length(list.names)), cat.cex=rep(0.8,length(list.names))
)

x <- minor.all
A <- x[[1]]
B <- x[[2]]
C <- x[[3]]
nab <- intersect(A, B)
nabc <- intersect(nab, C)
length(nabc)

#Removing overlapping minor alleles
min_all.m <- matrix(unlist(strsplit(nabc, "[.]")), ncol=2, byrow=T)
loc.parts <- matrix(unlist(strsplit(min_all.m[,1], "_")), ncol=4, byrow=T)
loc.df <- data.frame(CHROM=apply(loc.parts[,1:3],1,function(x) paste(x, collapse="_")), POS=loc.parts[,4], Allele=min_all.m[,2])
vcf.rec <- vcf

#Replacing low maf values with missing data indicators
for(i in 1:length(nabc)){
	if(i/50 == round(i/50,0)){print(paste("Processing locus number", i, "out of", length(nabc)))}
	tmp.loc <- which(vcf@fix[,"CHROM"] == loc.df$CHROM[i] & vcf@fix[,"POS"] == loc.df$POS[i])
	tmp.gt <- vcf@gt[tmp.loc,]
	tmp.list <- strsplit(tmp.gt,":")
	for(j in 2:length(tmp.list)){if(length(grep(loc.df$Allele[i], tmp.list[[j]][1])) > 0){tmp.list[[j]][1] <- ".|."}}
	tmp.gt <- lapply(tmp.list, function(x) paste(x,collapse=":"))
	vcf@gt[tmp.loc,] <- unlist(tmp.gt)[match(colnames(vcf@gt), names(tmp.gt))]
}

#Saving the data (Note: vcfR saves it in a vcf.gz format)
write.vcf(vcf, file="SNP.TRS.F08.vcf.gz")
}}} #notepad cleanup

#Removing paralogs and indivduals with high missing data#
{```{bash}```
mkdir Gnobilis_haplotyping
cd Gnobilis_haplotyping
cp -s /home/shared/Gambusia/Gambusia/reference/Gambusia_reference.fasta ./reference.fasta
cp -s ../popmap .

rad_haplotyper.pl -v SNP.TRS.F08.recode.vcf -r reference.fasta -p popmap -x 20 -m 0.8 -o SNP.TRS.F09.vcf -g SNP.TRS.F09.gen
} #notepad cleanup
