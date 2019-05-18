## Clean raw .fastq files

The following script uses the program Trimmomatic v. 0.33.

```
#!/bin/bash 
 
# Below is a script to loop through the files in the /filtered directory, identify matches,  
# and clean the fastq files, and direct output to ~/cleanreads 
 
cd /data/project_data/fastq/mp_process/
 
MyLeftReads=`find . -maxdepth 2 -name "*_R1*"`
  
#Making absolutely sure that reads are in the correct order: 
for myLeft in $MyLeftReads
do
    #The mate file will have a _2 instead of a _1. This changes the first occurence of _1 with _2. 
    myRight=${myLeft/_R1/_R2}
    #Is there a mate file? 
    if [-f $myRight]
    then
        echo "found a match $myLeft and $myRight ; now cleaning"
        java -classpath /data/popgen/Trimmomatic-0.33/trimmomatic-0.33.jar org.usadellab.trimmomatic.TrimmomaticPE \
                -threads 10 \
                -phred33 \
                 "$myLeft" \
                 "$myRight" \
                 /data/project_data/fastq/cleanreads/"$myLeft"_clean_paired.fq \
                 /data/project_data/fastq/cleanreads/"$myLeft"_clean_unpaired.fq \
                 /data/project_data/fastq/cleanreads/"$myRight"_clean_paired.fq \
                 /data/project_data/fastq/cleanreads/"$myRight"_clean_unpaired.fq \
                 ILLUMINACLIP:/data/popgen/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
                 LEADING:28 \
                 TRAILING:28 \
                 SLIDINGWINDOW:6:28 \
                 HEADCROP:12 \
                 MINLEN:35
# latest Illumina runs are phred 33 [[http://en.wikipedia.org/wiki/FASTQ_format]] 
# HEADCROP cuts the number of bases chosen from the head 
# SLIDINGWINDOW average the quality score across (6) bases and trims if avg is less than (28) 
# MINLEN excludes reads that are less than 40 bp long 
# 4 output files paired and unpaired for each side 
    fi
done
exit 0
```



## Map clean reads to reference transcriptome

The reference trancriptome was downloaded from the i5k website to our server using the following command.

```
# Note the need to "escape" out of the "(" characters using the forward slash "\".
wget https://i5k.nal.usda.gov/data/Arthropoda/onttau-\(Onthophagus_taurus\)/Current%20Genome%20Assembly/2.Official%20or%20Primary%20Gene%20Set/BCM_version_0.5.3/consensus_gene_set/OTAU.fna
```

Map the clean reads to the reference trancriptome using the following:

```
#!/bin/bash

cd ~/bwamaps/

# This indexing step only needs to be done once for the reference file.
bwa index -p ref -a is /data/project_data/OTAU.fna

MyLeftReads=`find . -maxdepth 1 -name "*_R1*.fastq.gz"`

#Making absolutely sure that reads are in the correct order:
for myLeft in $MyLeftReads
do
    #The mate file will have a _2 instead of a _1. This changes the first occurence of _1 with _2.
    myRight=${myLeft/_R1/_R2}
    #Is there a mate file?
    if [ -f $myRight ]
    then
        #echo $myLeft $myRight
        myShort=`echo $myLeft | cut -c3-12`
        #echo $myShort
         bwa aln /data/otau/reference/ref /data/otau/cleanreads/$myLeft > $myLeft".sai"
         bwa aln /data/otau/reference/ref /data/otau/cleanreads/$myRight > $myRight".sai"
         bwa sampe -r '@RG\tID:'"$myShort"'\tSM:'"$myShort"'\tPL:Illumina' \
                  -P /data/otau/reference/ref $myLeft".sai" $myRight".sai" \
                  /data/otau/cleanreads/$myLeft \
                  /data/otau/cleanreads/$myRight > $myShort"_bwaaln.sam"
    fi
done
```

## SNP Identification and genotyping

We used a combination of samtools and sambamba tools.  Sambamba was substantially faster for the sorting step; samtools sort would run out of memory and not complete.  The samtools remove duplicates function in version 1.2 was incompatible with our data.  For this step we used samtools version 0.19, then returned to version 1.2 for the mpileup SNP identification and genotyping step.

```
#!bin/bash/

cd /data/otau/cleanreads

for file in *.sam

do

samtools view -@ 16 -bS "$file" >"$file.bam"

done

#Using version 1.2 of samtools
samtools merge merged.bam *.bam

#Finding mate of reads
samtools fixmate merged.bam merged_fixmate.bam

# Switch to using sambamba for sorting - less memory intensive and it completes!
sambamba_v0.6.0 sort -m 8G -t 8 -p --tmpdir=/data/otau/cleanreads/tmp/ /data/otau/cleanreads/merged_fixmate.bam merged.fixmate.sorted.sambada

#Using Samtools version 0.19 (bug with newer versions on this step)
samtools_019 rmdup merged_fixmate_sorted.bam merged_fixmate_sorted_rmdup.bam

#Back to version 1.2
samtools index merged_fixmate_sorted_rmdup.bam

# SNP identification and genotyping
samtools mpileup -uf /data/otau/reference/OTAU.fna \
	merged_fixmate_sorted_rmdup.bam \ 
	| bcftools call -mv -> OTAU_merged.raw.vcf
```



Filter the raw vcf file - filter SNPs for quality, depth, biallelic, remove indels

VCFtools - 0.1.14

```markdown
vcftools --vcf OTAU_merged.raw.vcf --remove-indels --minDP 10 --minGQ 20 --min-alleles 2 --max-alleles 2 --max-missing 0.8 --recode --out OTAU-filtered
```

 ## Population genomics calculations and analyses

For each population:

```
vcftools --vcf OTAU_refiltered_IT.recode.vcf --TajimaD 100000 --out IT_TajD100000
vcftools --vcf OTAU_refiltered_IT.recode.vcf --window-pi 100000 --out IT_pi100000
vcftools --vcf OTAU_refiltered_IT.recode.vcf --SNPdensity 100 --out IT_SNPden100
```

For each population pair:

​	Calculate pairwise Weir and Cockerham's FST:

```
vcftools --vcf OTAU_refiltered_out_IT.recode.vcf --weir-fst-pop WA.txt --weir-fst-pop NC.txt --out WAvsNC_Fst
```

For functional enrichment analysis, we want a single FST value per gene.  We use the maximum FST per gene - calculated using the ddply function in R.

```R
library(plyr)

dt <- data.frame(fstWAvsNC)
# include the na.omit so that it will still calculate the average FST even if some of the SNPs in a gene have NAs; mean will still be NA if all snps in the gene are NA

maxFst_WAvsNC <- ddply(dt,~CHROM,summarise,max=max(na.omit(WEIR_AND_COCKERHAM_FST)))  

write.table(maxFst_WAvsNC, file = "maxFst_WAvsNC.txt", row.names = F, quote = F)
```

  

To output individual genotype data for plots:

```less
vcftools --vcf OTAU-filtered.recode.vcf --012
```



### BayeScan2.1

To be able to run BayeScan, first use PGDSpider (v. 2.0.9.0) to make eigensoft file from vcf. For each population pair:

```
java -Xmx6G -jar PGDSpider_2.0.9.0/PGDSpider2-cli.jar \
-inputfile OTAU_refiltered_wout_NC.recode.vcf -inputformat VCF \
-outputfile OTAU_refiltered_wout_NC.recode.bayescan \
-outputformat GESTE_BAYE_SCAN -spid beetle.noNC.spid
```

```
BayeScan2.1_linux64bits OTAU_refiltered_wout_NC.recode.bayescan -threads 5 -pr_odds 1000 -out_pilot beetallbayes -out_freq beetallbayes
```
