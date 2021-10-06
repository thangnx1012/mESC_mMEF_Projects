cd /media/hkh/8TB/XUANTHANG/INO80_ESC_MEF/RNAseq

mkdir -p Data/sra Data/Fastq Data/FastQC Data/Trimming
mkdir -p Results/STAR Results/feartureCounts Results/Rplot

# Loading all the modules
export PATH=/media/hkh/8TB/XUANTHANG/TOOLS/sratoolkit.2.10.9-ubuntu64/bin:$PATH

# ChIPseq ESC-SINGLE-END; MEF-PAIREND
prefetch --option-file Data/SRA_RNAseq_INO80_ESC.txt --output-directory Data/sra/ESC
prefetch --option-file Data/SRA_RNAseq_INO80_MEF.txt --output-directory Data/sra/MEF
fastq-dump Data/sra/ESC/SRR*/SRR* --outdir Data/Fastq/ESC --gzip
fastq-dump --split-files Data/sra/MEF/SRR*/SRR* --outdir Data/Fastq/MEF --gzip
fastqc -t 6 -o Data/FastQC --noextract -f fastq Data/Fastq/*/SRR*
multiqc Data/FastQC/. 

#2 Alignment (USING UCSC MOUSE mm10)
cd /media/hkh/8TB/XUANTHANG/INO80_ESC_MEF/RNAseq
** REFERENCE aligment **
# STAR --runThreadN 6 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles chr1.fa --sjdbGTFfile chr1-hg19_genes.gtf --sjdbOverhang 99


# ESC RNAseq (Single-end Using STAR)
for i in $(ls Data/Fastq/ESC/*.fastq.gz | cut -c 16- | rev | cut -c 10- | rev | uniq); do \
STAR 	--runThreadN 6 \
	--genomeDir /media/hkh/8TB/XUANTHANG/References/Reference_Mouse/Mus_musculus_UCSC_mm10/Sequence/STAR_Index \
	--runMode alignReads \
	--readFilesIn Data/Fastq/ESC/${i}.fastq.gz \
	--outSAMtype BAM SortedByCoordinate \
	--readFilesCommand zcat \
	--outSAMunmapped Within \
	--outSAMattributes Standard \
	--outFileNamePrefix Results/STAR/${i} ; done

# MEF RNA-seq (Pair-end, Using STAR)      
for i in $(ls Data/Fastq/MEF/*.fastq.gz | cut -c 16- | rev | cut -c 12- | rev | uniq); do \
STAR 	--runThreadN 6 \
	--genomeDir /media/hkh/8TB/XUANTHANG/References/Reference_Mouse/Mus_musculus_UCSC_mm10/Sequence/STAR_Index \
	--runMode alignReads \
	--readFilesIn Data/Fastq/MEF/${i}_1.fastq.gz Data/Fastq/MEF/${i}_2.fastq.gz  \
	--outSAMtype BAM SortedByCoordinate \
	--readFilesCommand zcat \
	--outSAMunmapped Within \
	--outSAMattributes Standard \
	--outFileNamePrefix Results/STAR/${i} ; done
	
#3 Quality check of aligned reads and Index file
cd /media/hkh/8TB/XUANTHANG/INO80_ESC_MEF/RNAseq
#Quality check
multiqc -o Results/STAR/ --pdf Results/STAR/.
#Index
for file in Results/STAR/*.bam; do \
samtools index $file ; done

#5  Downstream Analysis 
#5.1 RNAseq - Counting reads with SubRead (can detect and sort paired-end/single-end)

featureCounts \
      -t exon \
      -g gene_id \
      --primary \
      -a /media/hkh/8TB/XUANTHANG/References/Reference_Mouse/Mus_musculus_UCSC_mm10/Annotation/Genes/genes.gtf  \
      -o  Results/featureCounts/INO80_KD_ESC_MEF_counts.tsv Results/STAR/*.bam 

#Then move to Rstudio
cat Results/featureCounts/featureCounts.tsv | column | head


