# CHIP-seq pipeline
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
mkdir -p index
bowtie2-build --threads 20 reference/mhaESC_v1.1_with_mT2T-Y_v1.0.250617.fasta index/mouse_t2t

mkdir -p map_output
bowtie2 \
    -x index/mouse_t2t \
    -1 test_data/${SAMPLE}/trim/${SAMPLE}_1.fastq.gz \
    -2 test_data/${SAMPLE}/trim/${SAMPLE}_2.fastq.gz \
    -p 20 \
    --end-to-end \
    --very-sensitive \
    --phred33 \
    --no-unal \
    -S map_output/${SAMPLE}.sam
samtools view -@ 20 -h -b map_output/${SAMPLE}.sam > map_output/${SAMPLE}.bam
rm map_output/${SAMPLE}.sam
```
## Remove PCR duplicates
```bash
mkdir -p dedup
samtools collate -@ 20 -O -u map_output/${SAMPLE}.bam | samtools fixmate -@ 20 -m -u - - | samtools sort -@ 20 -u - | samtools markdup -@ 20 -r - dedup/${SAMPLE}_dedup.bam
```
## Bam2bw
```bash
mkdir -p bw_file
samtools index -@ 20 dedup/${SAMPLE}_dedup.bam
bamCoverage \
    -b dedup/${SAMPLE}_dedup.bam \
    -o bw_file/${SAMPLE}.bigwig \
    --binSize 5 \
    --normalizeUsing RPKM \
    --numberOfProcessors 20
```
## Callpeak
```bash
mkdir -p peak
# Narrow peak
macs2 callpeak \
      -t dedup/CHIP_IP_dedup.bam \
      -c dedup/CHIP_INPUT_dedup.bam \
      -g mm \
      --nomodel \
      -p 0.00001 \
      -n SAMPLE_NAME \
      --outdir peak -f BAMPE 

# Broad peak
macs2 callpeak \
      -t dedup/CHIP_IP_dedup.bam \
      -c dedup/CHIP_INPUT_dedup.bam \
      -g mm \
      --nomodel \
      -p 0.00001 \
      -n SAMPLE_NAME \
      --broad \
      --broad-cutoff 0.05 \
      --outdir peak -f BAMPE
```