####### Gambusia demultiplexing and SNP calling #######

#Demultiplexing
{```{bash}```
ls -d Lib* | while read i; do
cd $i
FILE=$(ls Extract*)
sbatch --export=FILE=$FILE $WORK/slurm/Demultiplex.slurm
cd ..
done
}

#Trimming
{```{bash}```
ls | while read i; do cd $i; cp -sf $WORK/slurm/trim_config.file .; cp -sf $WORK/Workspace/Gambusia/reference/reference.fasta .; cd ..; done
ls | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch $WORK/slurm/dDocent_trimming.slurm; cd ..; done
}

#Mapping
{```{bash}```
find . -type d > ../pear_mapping/dirs.txt; cd ../pear_mapping; xargs mkdir -p < dirs.txt; rm dirs.txt
ls | while read i; do cd $i; cp -sf $WORK/Workspace/Gambusia/reference/reference.fasta .; cd ..; done
ls | while read i; do cd $i; bwa index reference.fasta & cd ..; done

#Looking for overlap
ls | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch $WORK/slurm/pear/PEAR_trimmed_reads.slurm; cd ..; sleep 2; done

#Mapping reads
ls | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch $WORK/slurm/pear/PEAR_mapping.slurm; cd ..; sleep 2; done

#Filtering bam files
ls | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch $WORK/slurm/pear/PEAR_PE_filter.slurm; cd ..; sleep 2; done

#Combining SE and PE reads
ls | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch $WORK/slurm/pear/PEAR_combo.slurm; cd ..; sleep 2; done
}

#Subsetting the reference for later hapotyping
{```{bash}```
mkdir ../pear_cat && cd ../pear_cat
cp -s ../pear_mapping/Lib*/*-RG.bam .
cp -s ../reference/reference.fasta .
cp ../pear_mapping/Lib1/reference.fasta.* .
ls *.bam | sed 's/-RG.bam//g' | sed '/cat-RRG/d' > namelist

ls *-RG.bam | grep -v cat > bamlist.list
sbatch $WORK/slurm/pear/PEAR_cat.slurm

bedtools getfasta -fi reference.fasta -bed ../pear_SNPs/mapped.bed > ../reference/ref_subset.fasta
}

#Remapping the reads
{```{bash}```
mkdir ../pear_remap && cd ../pear_remap
mkdir Lib1 Lib2 Lib3
ls | while read i; do echo "Processing $i"; cd $i; cp -s ../../pear_mapping/$i/*.fastq ../../pear_mapping/$i/namelist .; cd ..; done
ls | while read i; do echo "Processing $i"; cd $i; cp -s ../../reference/ref_subset.fasta ./reference.fasta; samtools faidx reference.fasta; bwa index reference.fasta; cd ..; done

#Mapping reads
ls | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch $WORK/slurm/pear/PEAR_mapping.slurm; cd ..; sleep 2; done

#Filtering bam files
ls | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch $WORK/slurm/pear/PEAR_PE_filter.slurm; cd ..; sleep 2; done

#Combining SE and PE reads
ls | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch $WORK/slurm/pear/PEAR_combo.slurm; cd ..; sleep 2; done
}

#SNP calling
{```{bash}```
mkdir ../pear_SNPs && cd ../pear_SNPs

cp -s ../pear_remap/Lib*/*-RG.bam .
cp ../pear_remap/Lib1/reference.fasta* .
ls *.bam | sed 's/-RG.bam//g' | sed '/cat-RRG/d' > namelist
ls *-RG.bam | grep -v cat > bamlist.list
sbatch $WORK/slurm/pear/PEAR_cat.slurm
cp ../SNP_calling/popmap .

#Splitting cat.bam for SNP calling
sbatch --export=NODES=8 $WORK/slurm/dDocent_split.slurm

#SNP calling on normal nodes, scripted to limit to 8 nodes
ulimit -s 81920
ls -d *.node | while read i; do cd $i; ulimit -s 81920; sbatch $WORK/slurm/dDocent_freebayes.slurm; cd ..; sleep 2; done

#Checking to see if nodes have the correct number of vcf files with data in them
ls -d *.node | while read i; do echo "checking $i"; cd $i; BEDS=$(ls mapped.*.bed | wc -l);	VCF=$(find . -name "*.vcf" -size +62k | wc -l); if [ $VCF -lt $BEDS ]; then echo $i "did not complete all the vcf files properly"; fi; if [ $(find . -name "*.vcf" | wc -l) -gt $BEDS ]; then echo $i "has too many vcf files present"; fi; cd ..; done

#Checking to see if vcf file have the same number of contigs in the header as the reference.fasta has
ls -d *.node | while read i; do echo "checking $i"; cd $i; ls raw.*.vcf | while read j; do VCF=$(head -n1 $j| grep "##" | wc -l); if [ $VCF -eq 0 ]; then echo $i $j "missing complete header"; echo ${i}/${j} >> ../bad.vcf; fi; done; cd ..; done

#Applying a basic filter to each 
ls -d *.node | while read i; do cd $i; sbatch $WORK/slurm/dDocent_basic_filters.slurm; cd ..; sleep 2; done

#Combine all vcfs in each node
ls -d *.node | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch $WORK/slurm/dDocent_combine_filter_node.slurm; cd ..; sleep 2; done

#Preparing to combine all the vcf files
mkdir vcf && cd vcf

#Move data to a central directory
ls -d ../*.node | while read i; do 
NODE=$(echo $i | sed 's:../::g' | cut -f1 -d .)
cp -fs $i/cat_f.vcf ./raw.$NODE.vcf
done

#Combining all the vcf files
sbatch $WORK/slurm/dDocent_combine.slurm
}
