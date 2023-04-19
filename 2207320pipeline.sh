#make the directory in which we’ll be building the pipeline. We use the mkdir command to make a directory. 
mkdir advanced.assessment.NGS

#change directory
cd ~/advanced.assessment.NGS

#We can again use the mkdir command to create a sub-directory within the advanced.assessment.NGS directory. We’ll call this subdirectory DNAseq
mkdir ~/advanced.assessment.NGS/DNAseq

#We can now change our directory to the DNAseq directory. We do this using the cd command 
cd ~/advanced.assessment.NGS/DNAseq

#we make these different subdirectories to keep everything organised down the line
mkdir data metadata results log.command 

#we can use these commands to verify that we have created the four sub-directories
ls -lF 

#change directory to the data directory
cd ~/advanced.assessment.NGS/DNAseq/data

#this creates a directory for the untrimmed files
mkdir fastq_untrimmed 

#this creates a directory for the trimmed fastq files
mkdir fastq_trimmed 

#we can use this command to confirm that we have made the two directories for the trimmed and untrimmed fastq files
ls -a 

#we now need to download our two fastq.qz files and the annotation.bed file. To do this, we use the wget function.

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

# We can now move all the files ending with ‘fastq.qz’ to the ‘fastq_untrimmed’ folder
mv *fastq.qz ~/advanced.assessment.NGS/DNAseq/data/fastq_untrimmed

#We now also need to download the reference genome. Again, we use the wget command
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

#We now need to change our directory back to the home directory so that we can install anaconda and all other tools we might need during the pipeline
cd ~/

#We can now install Anaconda
wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh


chmod +x ./Anaconda3-2022.10-Linux-x86_64.sh

bash ./Anaconda3-2022.10-Linux-x86_64.sh

source ~/.bashrc

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

#We can now begin to install the tools that we will need throughout the pipeline
conda install Samtools
conda install bwa
conda install Freebayes
conda install picard
conda install bedtools
conda install Trimmomatic
conda install fastqc
conda install vcflib

#change directory to the fastq_untrimmed directory
cd ~/advanced.assessment.NGS/DNAseq/data/fastq_untrimmed

#double check that both fastq.qz files are in that directory
ls -a

#we now need to rename the fastq.qz files to fastq.gz files in order for fastqc to work. We can do this using the mv command
mv NGS0001.R1.fastq.qz NGS0001.R1.fastq.gz

mv NGS0001.R2.fastq.qz NGS0001.R2.fastq.gz

#we can run the fastqc on all the fastq.gz files using this command
fastqc *fastq.gz
#We can make a sub-directory within the results directory for the fastqc results on the untrimmed reads
mkdir ~/advanced.assessment.NGS/DNAseq/results/fastqc_untrimmed_reads 

#we can now move all of the fastqc results to this new fastqc_untrimmed_reads file
mv *fastqc* ~/advanced.assessment.NGS/DNAseq/results/fastqc_untrimmed_reads

#we can confirm that it has moved all of the fastqc files 
ls -lh ~/advanced.assessment.NGS/DNAseq/results/fastqc_untrimmed_reads

#now we change our directory to this fastqc_untrimmed_reads directory
cd ~/advanced.assessment.NGS/DNAseq/results/fastqc_untrimmed_reads

#the command unzip could not be found and needs to be installed via this command
sudo apt install unzip

#use a for loop to unzip the zip files in the fastqc results
for zip in *.zip
> do
> unzip $zip
> done

#use the below command to see all files in the fastqc_untrimmed_reads directory
ls -a

#we can use ls -lh to see what is stored inside of the fastqc files and then use head to see the first couple of lines within the summary.txt files of each of the fastqc files
ls -lh NGS0001.R1.fastq.gz_fastqc
ls -lh NGS0001.R2.fastq.gz_fastqc

head NGS0001.R1.fastq.gz_fastqc/summary.txt
head NGS0001.R2.fastq.gz_fastqc/summary.txt

#we can then use the cat command to create a concatenated summary.txt file of both the fastqc/summary.txt files and move them to the metadata directory
cat */summary.txt > ~/advanced.assessment.NGS/DNAseq/metadata/fastqc_untrimmed_summaries.txt

#we can now begin with Trimmomatic on the untrimmed reads

#change directory to the fastq_untrimmed directory
cd ~/advanced.assessment.NGS/DNAseq/data/fastq_untrimmed

#use ls -a to confirm that both fastq.gz untrimmed reads are in this directory
ls -a

#we can type the command Trimmomatic to show all of its instructions
trimmomatic

#we then need to use the gzip -d command to convert the fastq.gz files into .fastq files 
gzip -d NGS0001.R1.fastq.gz

gzip -d NGS0001.R2.fastq.gz

#we can now run Trimmomatic using the following command
trimmomatic PE -threads 4 -phred33 ~/advanced.assessment.NGS/DNAseq/data/fastq_untrimmed/NGS0001.R1.fastq ~/advanced.assessment.NGS/DNAseq/data/fastq_untrimmed/NGS0001.R2.fastq -baseout ~/advanced.assessment.NGS/DNAseq/data/fastq_trimmed/NGS0001_trimmed_R ILLUMINACLIP:~/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 TRAILING:25 MINLEN:50

#we can now change directories to the fastq_trimmed directory seeing as the Trimmomatic resulting files were sent there
cd ~/advanced.assessment.NGS/DNAseq/data/fastq_trimmed

#we then use ls -a again to confirm that the Trimmomatic files were sent to this directory
ls -a

#we can now perform fastqc on the paired end reads that were produced from the Trimmomatic output
fastqc NGS0001_trimmed_R_1P NGS0001_trimmed_R_2P

#we can now make a sub-directory in the ‘results’ directory for the fastqc results of the trimmed reads
mkdir ~/advanced.assessment.NGS/DNAseq/results/fastqc.trimmed.reads

#we can move the fastqc results of the trimmed reads to this fastqc.trimmed.reads folder
mv *fastqc* ~/advanced.assessment.NGS/DNAseq/results/fastqc.trimmed.reads

#we can now change our directory to the fastqc.trimmed.reads folder
cd ~/advanced.assessment.NGS/DNAseq/results/fastqc.trimmed.reads

#use ls -a to view all files in this fastqc.trimmed.reads folder
ls -a


#again use a for loop to unzip all of the .zip files in the fastqc results
for zip in *.zip
> do
> unzip $zip
> done

#we again use ls -lh to view the contents of the fastqc files
ls -lh NGS0001_trimmed_R_1P_fastqc NGS0001_trimmed_R_2P_fastqc

#we can then use the head command to view the first lines of the fastqc summary.txt files 
head NGS0001_trimmed_R_1P_fastqc/summary.txt

head NGS0001_trimmed_R_2P_fastqc/summary.txt

#again use the cat command to create a concatenated file of all the summary.txt files in the two fastqc files and move them to the metadata directory
cat */summary.txt > ~/advanced.assessment.NGS/DNAseq/metadata/fastqc_trimmed_summaries.txt

#we now need to change our directory to the data directory
cd ~/advanced.assessment.NGS/DNAseq/data

#use ls -a to see contents of the directory
ls -a

#we now need to make a sub-directory within the data directory for where we want to send the reference genome 
mkdir ref

#we now need to move the hg.19.fa.gz from the data directory into the ref sub-directory that we just created
mv ~/advanced.assessment.NGS/DNAseq/data/hg19.fa.gz ~/advanced.assessment.NGS/DNAseq/data/ref/

#we can now begin with the indexing process of bwa using the hg19 reference genome
bwa index ~/advanced.assessment.NGS/DNAseq/data/ref/hg19.fa.gz 

#we can use the ls command to see the contents of the ref directory
ls ~/advanced.assessment.NGS/DNAseq/data/ref

#we now need to make a sub-directory in the data directory for the aligned data 
mkdir ~/advanced.assessment.NGS/DNAseq/data/aligned_data


#in order to obtain some information to later on edit the bwa mem step to include read group information, we need to use the head command to see the first lines of the untrimmed fastq files
head ~/advanced.assessment.NGS/DNAseq/data/fastq_untrimmed/NGS0001.R1.fastq

head ~/advanced.assessment.NGS/DNAseq/data/fastq_untrimmed/NGS0001.R2.fastq


#we can now perform the bwa mem step, and alter the output sam file to include read group information 
bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1.111.D1375ACXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001\tDT:2017-02-23\tPU:11V6WR1' -I 250,50 ~/advanced.assessment.NGS/DNAseq/data/ref/hg19.fa.gz ~/advanced.assessment.NGS/DNAseq/data/fastq_trimmed/NGS0001_trimmed_R_1P ~/advanced.assessment.NGS/DNAseq/data/fastq_trimmed/NGS0001_trimmed_R_2P > ~/advanced.assessment.NGS/DNAseq/data/aligned_data/NGS0001_trimmed.sam

#we can now change our directory to the aligned_data directory
cd ~/advanced.assessment.NGS/DNAseq/data/aligned_data

#we can use Samtools view -h -b to convert the .sam file to a .bam file
samtools view -h -b NGS0001_trimmed.sam > NGS0001_trimmed.bam

#we can use Samtools sort to sort the trimmed .bam file
samtools sort NGS0001_trimmed.bam > NGS0001_trimmed_sorted.bam

#this step will create a .bai index file
samtools index NGS0001_trimmed_sorted.bam

#use ls to see the contents of the align_data directory
ls

#perform duplicate marking on the sorted bam file 
picard MarkDuplicates I=NGS0001_trimmed_sorted.bam O=NGS0001_trimmed_sorted_marked.bam M=marked_dup_metrics.txt

#index the sorted and marked bam file
samtools index NGS0001_trimmed_sorted_marked.bam

#quality filter the sorted and marked bam file using samtools view
samtools view -F 1796  -q 20 -o NGS0001_trimmed_sorted_filtered.bam NGS0001_trimmed_sorted_marked.bam

#index the sorted and filtered bam file
samtools index NGS0001_trimmed_sorted_filtered.bam

#perform flagstat statistical analysis on the sorted and filtered bam file and move the resulting file to a flagstat.results.bam file
samtools flagstat NGS0001_trimmed_sorted_filtered.bam > NGS0001.flagstat.results.bam

#perform idxstats statistical analysis on the sorted and filtered bam file and move the resulting file to a idxstats.results.bam file 
samtools idxstats NGS0001_trimmed_sorted_filtered.bam > NGS0001.idxstats.results.bam

#perform depth of coverage analysis of the sorted and filtered bam file using the annotation.bed that we have and move the resulting file to a coverage_depth_results.txt file
bedtools coverage -a NGS0001_trimmed_sorted_filtered.bam -b ~/advanced.assessment.NGS/DNAseq/data/annotation.bed > coverage_depth_results.txt

#we can also collect insert size metric using picard CollectInsertSizeMetrics
picard CollectInsertSizeMetrics I=NGS0001_trimmed_sorted_filtered.bam O=NGS0001_trimmed_sorted_filtered.bam.hist H=NGS0001_trimmed_sorted_filtered.bam.pdf

#use the zcat command to change the ref genome from hg19.fa.gz to hg19.fa
zcat ~/advanced.assessment.NGS/DNAseq/data/ref/hg19.fa.gz > ~/advanced.assessment.NGS/DNAseq/data/ref/hg19.fa

#create an index of the ref genome
samtools faidx ~/advanced.assessment.NGS/DNAseq/data/ref/hg19.fa

#call variants using freebayes
freebayes --bam ~/advanced.assessment.NGS/DNAseq/data/aligned_data/NGS0001_trimmed_sorted_filtered.bam --fasta-reference ~/advanced.assessment.NGS/DNAseq/data/ref/hg19.fa --vcf ~/advanced.assessment.NGS/DNAseq/results/NGS0001.vcf

#use bgzip to change the vcf file from .vcf to a .vcf.gz
bgzip ~/advanced.assessment.NGS/DNAseq/results/NGS0001.vcf

#use the tabix command to the turn the .vcf.gz file into a tbi file
tabix -p vcf ~/advanced.assessment.NGS/DNAseq/results/NGS0001.vcf.gz

#perform quality filtering using vcffilter
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
	> ~/advanced.assessment.NGS/DNAseq/results/NGS0001.vcf.gz > ~/advanced.assessment.NGS/DNAseq/results/NGS0001_filtered.vcf

#restrict quality filtering analysis to regions included in the bed file
bedtools intersect -header -wa -a ~/advanced.assessment.NGS/DNAseq/results/NGS0001_filtered.vcf -b ~/advanced.assessment.NGS/DNAseq/data/annotation.bed \
	> > ~/advanced.assessment.NGS/DNAseq/results/NGS0001_filtered_annotated.vcf

#use bgzip to change the vcf file into a .vcf.gz file 
bgzip ~/advanced.assessment.NGS/DNAseq/results/NGS0001_filtered_annotated.vcf


#use tabix to generate a .tbi file 
tabix -p vcf ~/advanced.assessment.NGS/DNAseq/results/NGS0001_filtered_annotated.vcf.gz

#Download annovar 
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz


#Once annovar is downloaded, unpack it and set up annovar using the following command
tar -zxvf annovar.latest.tar.gz

#Before using Annovar you need to download the databases it uses for the annotation. The following commands will download some nasic ones:
cd annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

#convert the .vcf.gz file into a .avinput file
./convert2annovar.pl -format vcf4 ~/advanced.assessment.NGS/DNAseq/results/NGS0001_filtered_annotated.vcf.gz > ~/advanced.assessment.NGS/DNAseq/results/NGS0001_filtered_annotated.avinput

#perform variant prioritisation using annovar
./table_annovar.pl ~/advanced.assessment.NGS/DNAseq/results/NGS0001_filtered_annotated.avinput humandb/ -buildver hg19 -out ~/advanced.assessment.NGS/DNAseq/results/NGS0001_filtered_annotated -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout


#download snpEFF
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

#unzip the snpEFF zip file using unzip
unzip snpEff_latest_core.zip

#change the directory to the snpEFF directory
cd snpEff  

#run snpEff variant annotation on filtered VCF file 
java -Xmx8g -jar snpEff.jar GRCh37.75 ~/advanced.assessment.NGS/DNAseq/results/NGS0001_filtered_annotated.vcf.gz > ~/advanced.assessment.NGS/DNAseq/results/NGS0001_filtered_snp_annotated.vcf.gz

#filter for exonic variants 
java -jar SnpSift.jar filter “ANN[0].EFFECT has ‘exon_variant’” ~/advanced.assessment.NGS/DNAseq/results/NGS0001_filtered_snp_annotated.vcf.gz > ~/advanced.assessment.NGS/DNAseq/results/NGS0001_filtered_exon_annotated.vcf.gz

