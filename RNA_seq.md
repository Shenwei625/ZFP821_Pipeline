# RNA-seq pipeline
## Trim adapter and Fastqc
```bash
SAMPLE=

mkdir -p test_data/${SAMPLE}/fastqc
fastqc -t 20 -o test_data/${SAMPLE}/fastqc test_data/${SAMPLE}/${SAMPLE}_1.fastq.gz
fastqc -t 20 -o test_data/${SAMPLE}/fastqc test_data/${SAMPLE}/${SAMPLE}_2.fastq.gz

mkdir -p test_data/${SAMPLE}/trim
fastp --thread 20 --detect_adapter_for_pe --trim_poly_x \
    -i test_data/${SAMPLE}/${SAMPLE}_1.fastq.gz \
    -I test_data/${SAMPLE}/${SAMPLE}_2.fastq.gz \
    -o test_data/${SAMPLE}/trim/${SAMPLE}_1.fastq.gz \
    -O test_data/${SAMPLE}/trim/${SAMPLE}_2.fastq.gz \
    -h test_data/${SAMPLE}/trim/report.html \
    -j test_data/${SAMPLE}/trim/report.json

mkdir -p test_data/${SAMPLE}/trim/fastqc
fastqc -t 20 -o test_data/${SAMPLE}/trim/fastqc test_data/${SAMPLE}/trim/${SAMPLE}_1.fastq.gz
fastqc -t 20 -o test_data/${SAMPLE}/trim/fastqc test_data/${SAMPLE}/trim/${SAMPLE}_2.fastq.gz
```
## Align
```bash
mkdir -p reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz -O reference/mm10.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz -O reference/mm10.ncbiRefSeq.gtf.gz
pigz -dcf reference/mm10.fa.gz > reference/mm10.fa
pigz -dcf reference/mm10.ncbiRefSeq.gtf.gz > reference/mm10.ncbiRefSeq.gtf
rm reference/mm10.fa.gz reference/mm10.ncbiRefSeq.gtf.gz

mkdir -p STAR_index/mm10
STAR --runThreadN 20 \
     --runMode genomeGenerate \
     --genomeDir STAR_index/mm10 \
     --genomeFastaFiles reference/mm10.fa \
     --sjdbGTFfile reference/mm10.ncbiRefSeq.gtf \
     --sjdbOverhang 149

mkdir -p map_output_STAR
STAR --runThreadN 20 \
     --genomeDir STAR_index/mm10 \
     --readFilesIn test_data/${SAMPLE}/trim/${SAMPLE}_1.fastq.gz test_data/${SAMPLE}/trim/${SAMPLE}_2.fastq.gz \
     --readFilesCommand zcat \
     --outFileNameSAMPLE map_output_STAR/${SAMPLE}_ \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMstrandField intronMotif \
     --alignEndsType EndToEnd \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outMultimapperOrder Random \
     --outSAMmultNmax 1 \
     --outFilterType BySJout \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000
```
## Count
```bash
mkdir -p tecount
TEcount --format BAM -b map_output_STAR/${SAMPLE}_Aligned.sortedByCoord.out.bam --GTF reference/mm10.ncbiRefSeq.gtf --TE GRCm38_GENCODE_rmsk_TE.gtf --sortByPos --mode multi --project tecount/${SAMPLE}
```
