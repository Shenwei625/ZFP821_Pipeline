# RNA-seq pipeline
## Trim adapter and Fastqc
```bash
SAMPLE=

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
# For gene count
mkdir -p featurecount
featureCounts -p -Q 10 -s 2 --countReadPairs -T 40 \
-a reference/mm10.ncbiRefSeq.gtf -g gene_id -t exon -o featurecount/${SAMPLE}_gene_count.txt \
map_LINE/${BAM_FILE}_Aligned.sortedByCoord.out.bam

# For TE count
mkdir -p tecount
TEcount --format BAM -b map_output_STAR/${SAMPLE}_Aligned.sortedByCoord.out.bam --GTF reference/mm10.ncbiRefSeq.gtf --TE GRCm38_GENCODE_rmsk_TE.gtf --sortByPos --mode multi --project tecount/${SAMPLE}
## The file "GRCm38_GENCODE_rmsk_TE.gtf" should be downloaded from the TEtranscripts official website.
```
## Extraction of differentially expressed genes
```R
library(DESeq2) 
data <- read.table("data/RNA_seq/RS_count.tsv", header = TRUE, sep = "\t", row.names = 1)
colData <- data.frame(
  condition = factor(c("WT", "WT", "KO", "KO")),
  batch = factor(c("batch1","batch2","batch1","batch2"))
)
countdata <- data[rowMeans(data) > 0, ]
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = colData,
  design = ~ batch + condition
)
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized=TRUE)
norm_df <- as.data.frame(normalized_counts)
norm_df$Gene <- rownames(norm_df)
norm_df <- norm_df[, c("Gene", setdiff(colnames(norm_df), "Gene"))]

res <- results(dds, contrast = c("condition", "KO", "WT"))
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)
res_df <- res_df[, c("Gene", setdiff(colnames(res_df), "Gene"))]
write.table(res_df, "Output/RNA_seq/RS_deseq2.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
```
