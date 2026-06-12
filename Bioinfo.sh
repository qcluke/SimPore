### ILLUMINA ----------------------------------------------------------------------------------------------------------------

## Truncating and Merging
module load VSEARCH/2.18
# 16S V4
for R1 in /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/raw_data/16sv4/*_1.fastq; do
    R2=${R1/_1.fastq/_2.fastq}
    BASE=$(basename "$R1" _1.fastq)
    echo "======== Merging reads for: ${R1/_1.fastq/}"
    vsearch --fastq_filter "$R1" \
        --reverse "$R2" \
        --fastq_trunclen 220 \
        --fastqout /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/merged/16sv4/${BASE}_1_truncated.fastq \
        --fastqout_rev /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/merged/16sv4/${BASE}_2_truncated.fastq \
    vsearch --fastq_mergepairs /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/merged/16sv4/${BASE}_1_truncated.fastq \
        --reverse /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/merged/16sv4/${BASE}_2_truncated.fastq \
        --fastqout /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/merged/16sv4/${BASE}.fastq 
done
# 18S V9
for R1 in /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/raw_data/18sv9/*_1.fastq; do
    R2=${R1/_1.fastq/_2.fastq}
    BASE=$(basename "$R1" _1.fastq)
    echo "======== Merging reads for: ${R1/_1.fastq/}"
    vsearch --fastq_filter "$R1" \
        --reverse "$R2" \
        --fastq_trunclen 110 \
        --fastqout /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/merged/18sv9/${BASE}_1_truncated.fastq \
        --fastqout_rev /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/merged/18sv9/${BASE}_2_truncated.fastq \
    vsearch --fastq_mergepairs /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/merged/18sv9/${BASE}_1_truncated.fastq \
        --reverse /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/merged/18sv9/${BASE}_2_truncated.fastq \
        --fastqout /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/merged/18sv9/${BASE}.fastq 
done

## Primer Trimming
module load cutadapt/3.5
# 16S V4
for i in *.fastq; do
  INPUT_FILE="/scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/merged/16sv4/${i}"
  OUTPUT_FILE="/scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/removed/16sv4/${i}"
  LOG_FILE="/scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/removed/log_cutadapt_16sv4.txt"
  echo "Removing primer for: ${i}"
  cutadapt ${INPUT_FILE} \
    --cores 8 \
    --adapter GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC \
    --revcomp \
    --error-rate 0.1 \
    --match-read-wildcards \
    --action trim \
    --discard-untrimmed \
    --minimum-length 230 \
    --maximum-length 280 \
    --report=full \
    --output ${OUTPUT_FILE} >> ${LOG_FILE}
done
# 18S V9
for i in *.fastq; do
  INPUT_FILE="/scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/merged/18sv9/${i}"
  OUTPUT_FILE="/scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/removed/18sv9/${i}"
  LOG_FILE="/scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/removed/log_cutadapt_18sv9.txt"
  echo "Removing primer for: ${i}"
  cutadapt ${INPUT_FILE} \
    --cores 8 \
    --adapter GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC \
    --revcomp \
    --error-rate 0.1 \
    --match-read-wildcards \
    --action trim \
    --discard-untrimmed \
    --minimum-length 110 \
    --maximum-length 150 \
    --report=full \
    --output ${OUTPUT_FILE} >> ${LOG_FILE}
done

## Length and Quality Filtering
module load VSEARCH/2.18
# 16S V4
for i in /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/removed/16sv4/*.fastq; do
  BASE=$(basename "$i" .fastq)
  echo "Removing primer for: ${BASE}"
  vsearch --fastq_filter /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/removed/16sv4/${BASE}.fastq \
    --fastq_maxee 1.0 \
    --fastq_maxlen 280 \
    --fastq_minlen 230 \
    --fastq_maxns 0 \
    --fastaout /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/filtered/16sv4/${BASE}.fasta \
    --relabel ${BASE}.
done
# 18S V9
for i in /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/removed/18sv9/*.fastq; do
  BASE=$(basename "$i" .fastq)
  echo "Removing primer for: ${BASE}"
  vsearch --fastq_filter /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/removed/18sv9/${BASE}.fastq \
    --fastq_maxee 1.0 \
    --fastq_maxlen 150 \
    --fastq_minlen 110 \
    --fastq_maxns 0 \
    --fastaout /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/filtered/18sv9/${BASE}.fasta \
    --relabel ${BASE}.
done

## Conversion and Dereplication
module load VSEARCH/2.18
# 16S V4
cat /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/filtered/16sv4/*.fasta \
    > /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/dereplicated/16sv4_combined.fasta
vsearch --derep_fulllength /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/dereplicated/16sv4_combined.fasta \
    --output /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/dereplicated/16sv4_uniques.fasta \
    --relabel uniq. \
    --sizeout \
    --uc /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/dereplicated/16sv4_derep.uc
# 18S V9
cat /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/filtered/18sv9/*.fasta \
    > /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/dereplicated/18sv9_combined.fasta
vsearch --derep_fulllength /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/dereplicated/18sv9_combined.fasta \
    --output /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/dereplicated/18sv9_uniques.fasta \
    --relabel uniq. \
    --sizeout \
    --uc /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/dereplicated/18sv9_derep.uc

## Clustering
module load VSEARCH/2.18
# 16S V4
vsearch --cluster_size /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/dereplicated/16sv4_uniques.fasta \
    --id 0.97 \
    --sizein \
    --relabel clustered. \
    --sizeout \
    --consout /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/clustered/16sv4_clustered.fasta \
    --uc /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/clustered/16sv4_cluster.uc
# 18S V9
vsearch --cluster_size /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/dereplicated/18sv9_uniques.fasta \
    --id 0.97 \
    --sizein \
    --relabel clustered. \
    --sizeout \
    --consout /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/clustered/18sv9_clustered.fasta \
    --uc /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/clustered/18sv9_cluster.uc

## Chimera Removal
module load VSEARCH/2.18
# 16S V4
vsearch --uchime_denovo /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/clustered/16sv4_clustered.fasta \
    --sizein \
    --nonchimeras /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/chimera_removed/16sv4_otu.fasta \
    --relabel otu.
# 18S V9
vsearch --uchime_denovo /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/clustered/18sv9_clustered.fasta \
    --sizein \
    --nonchimeras /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/chimera_removed/18sv9_otu.fasta \
    --relabel otu.

## Count Table
module load VSEARCH/2.18
# 16S V4
vsearch --usearch_global /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/dereplicated/16sv4_combined.fasta \
    --db /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/chimera_removed/16sv4_otu.fasta \
    --strand plus \
    --id 0.97 \
    --otutabout /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/frequency_table/illumina_16sv4_otu_table.tsv
# 18S V9
vsearch --usearch_global /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/dereplicated/18sv9_combined.fasta \
    --db /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/chimera_removed/18sv9_otu.fasta \
    --strand plus \
    --id 0.97 \
    --otutabout /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/frequency_table/illumina_18sv9_otu_table.tsv

## Sequence Import
conda activate qiime2-amplicon-2025.7
# 16S V4
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/tax_assigned/16sv4_otu.fasta \
  --output-path /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/tax_assigned/16sv4_otu.qza
# 18S V9
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/tax_assigned/18sv9_otu.fasta \
  --output-path /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/tax_assigned/18sv9_otu.qza

## Taxonomic Assignment
conda activate qiime2-amplicon-2025.7
# 16S V4
qiime feature-classifier classify-sklearn \
  --i-classifier /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/tax_assigned/16sv4_classifier.qza \
  --i-reads /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/tax_assigned/16sv4_otu.qza \
  --p-n-jobs 8 \
  --o-classification /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/tax_assigned/illumina_16sv4_taxonomy.qza
# 18S V9
qiime feature-classifier classify-sklearn \
  --i-classifier /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/tax_assigned/18sv9_classifier.qza \
  --i-reads /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/tax_assigned/18sv9_otu.qza \
  --p-n-jobs 8 \
  --o-classification /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/tax_assigned/illumina_18sv9_taxonomy.qza

## Taxonomic Table
conda activate qiime2-amplicon-2025.7
# 16S V4
qiime metadata tabulate \
  --m-input-file illumina_16sv4_taxonomy.qza \
  --o-visualization illumina_16sv4_taxonomy.qzv
# 18S V9
qiime metadata tabulate \
  --m-input-file illumina_18sv9_taxonomy.qza \
  --o-visualization illumina_18sv9_taxonomy.qzv



### ONT SHORT ----------------------------------------------------------------------------------------------------------------

## Basecalling and Demultiplexing
# 16S V4 and 18S V9
dorado basecaller \
    /scr/u/qcluke/dorado-1.0.0-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    /lustre1/g/sbs_sey/Luke/raw_data/nanopore-chapter2-raw-sequences/ONT-16SV4_18SV9/20241205_1726_MN47038_FAZ20309_2d4e7130/pod5/ \
    --device cuda:all \
    --kit-name SQK-NBD114-96 \
    --emit-fastq \
    --no-trim \
    --output-dir /scr/u/qcluke/basecalled/chapter2_16sv4_18sv9/demultiplexed/ \

## Primer Trimming
module load cutadapt/3.5
# 16S V4
for i in {01..48}; do
  INPUT_FILE="/scr/u/qcluke/Chapter2/16sv4_18sv9/demultiplexed/16sv4/2d4e7130-76e2-4f24-b39b-268aba7fc398_SQK-NBD114-96_barcode${i}.fastq"
  OUTPUT_FILE="/scr/u/qcluke/Chapter2/16sv4_18sv9/removed/16sv4/barcode${i}_removed.fastq"
  LOG_FILE="/scr/u/qcluke/Chapter2/16sv4_18sv9/removed/log_cutadapt_16sv4.txt"
  echo "Trimming barcode${i}"
  cutadapt ${INPUT_FILE} \
    --cores 8 \
    --adapter GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC \
    --revcomp \
    --error-rate 0.3 \
    --match-read-wildcards \
    --action trim \
    --discard-untrimmed \
    --minimum-length 230 \
    --maximum-length 280 \
    --report=full \
    --output ${OUTPUT_FILE} >> ${LOG_FILE}
done
# 18S V9
for i in {49..96}; do
  INPUT_FILE="/scr/u/qcluke/Chapter2/16sv4_18sv9/demultiplexed/18sv9/2d4e7130-76e2-4f24-b39b-268aba7fc398_SQK-NBD114-96_barcode${i}.fastq"
  OUTPUT_FILE="/scr/u/qcluke/Chapter2/16sv4_18sv9/removed/18sv9/barcode${i}_removed.fastq"
  LOG_FILE="/scr/u/qcluke/Chapter2/16sv4_18sv9/removed/log_cutadapt_18sv9.txt"
  echo "Removing primer for: barcode${i}"
  cutadapt ${INPUT_FILE} \
    --cores 8 \
    --adapter CCCTGCCHTTTGTACACAC...GTAGGTGAACCTGCRGAAGG \
    --revcomp \
    --error-rate 0.3 \
    --match-read-wildcards \
    --action trim \
    --discard-untrimmed \
    --minimum-length 110 \
    --maximum-length 150 \
    --report=full \
    --output ${OUTPUT_FILE} >> ${LOG_FILE}
done

## Length and Quality Filtering
conda activate chopper
# 16S V4
for i in {01..48}; do
    INPUT_FILE="/scr/u/qcluke/Chapter2/16sv4_18sv9/removed/16sv4/barcode${i}_removed.fastq"
    OUTPUT_FILE="/scr/u/qcluke/Chapter2/16sv4_18sv9/filtered/16sv4/barcode${i}_filtered.fastq"
    echo "filtering for: barcode${i}"
    chopper --quality 20 \
    --maxqual 50 \
    --minlength 230 \
    --maxlength 280 \
    --threads 4 \
    --input ${INPUT_FILE} > ${OUTPUT_FILE}
done
# 18S V9
for i in {49..96}; do
    INPUT_FILE="/scr/u/qcluke/Chapter2/16sv4_18sv9/removed/18sv9/barcode${i}_removed.fastq"
    OUTPUT_FILE="/scr/u/qcluke/Chapter2/16sv4_18sv9/filtered/18sv9/barcode${i}_filtered.fastq"
    echo "filtering for: barcode${i}"
    chopper --quality 20 \
    --maxqual 50 \
    --minlength 110 \
    --maxlength 150 \
    --threads 4 \
    --input ${INPUT_FILE} > ${OUTPUT_FILE}
done

## Conversion and Dereplication
module load VSEARCH/2.18
# 16S V4
for i in /scr/u/qcluke/Chapter2/16sv4_18sv9/filtered/16sv4/*_filtered.fastq; do
  BASE=$(basename "$i" _filtered.fastq)
  echo "========== Converting ${BASE} =========="
  vsearch --fastq_filter /scr/u/qcluke/Chapter2/16sv4_18sv9/filtered/16sv4/${BASE}_filtered.fastq \
    --fastq_qmax 50 \
    --fastaout /scr/u/qcluke/Chapter2/16sv4_18sv9/generated_fasta/16sv4/${BASE}.fasta \
    --relabel ${BASE}.
done
cat /scr/u/qcluke/Chapter2/16sv4_18sv9/generated_fasta/16sv4/*.fasta \
    > /scr/u/qcluke/Chapter2/16sv4_18sv9/dereplicated/16sv4_combined.fasta
vsearch --derep_fulllength /scr/u/qcluke/Chapter2/16sv4_18sv9/dereplicated/16sv4_combined.fasta \
    --output /scr/u/qcluke/Chapter2/16sv4_18sv9/dereplicated/16sv4_uniques.fasta \
    --relabel uniq. \
    --sizeout \
    --uc /scr/u/qcluke/Chapter2/16sv4_18sv9/dereplicated/16sv4_derep.uc
# 18S V9
for i in /scr/u/qcluke/Chapter2/16sv4_18sv9/filtered/18sv9/*_filtered.fastq; do
  BASE=$(basename "$i" _filtered.fastq)
  echo "========== Converting ${BASE} =========="
  vsearch --fastq_filter /scr/u/qcluke/Chapter2/16sv4_18sv9/filtered/18sv9/${BASE}_filtered.fastq \
    --fastq_qmax 50 \
    --fastaout /scr/u/qcluke/Chapter2/16sv4_18sv9/generated_fasta/18sv9/${BASE}.fasta \
    --relabel ${BASE}.
done
cat /scr/u/qcluke/Chapter2/16sv4_18sv9/generated_fasta/18sv9/*.fasta \
    > /scr/u/qcluke/Chapter2/16sv4_18sv9/dereplicated/18sv9_combined.fasta
vsearch --derep_fulllength /scr/u/qcluke/Chapter2/16sv4_18sv9/dereplicated/18sv9_combined.fasta \
    --output /scr/u/qcluke/Chapter2/16sv4_18sv9/dereplicated/18sv9_uniques.fasta \
    --relabel uniq. \
    --sizeout \
    --uc /scr/u/qcluke/Chapter2/16sv4_18sv9/dereplicated/18sv9_derep.uc

## Clustering
module load VSEARCH/2.18
# 16S V4
vsearch --cluster_size /scr/u/qcluke/Chapter2/16sv4_18sv9/dereplicated/16sv4_uniques.fasta \
    --id 0.97 \
    --sizein \
    --relabel clustered. \
    --sizeout \
    --consout /scr/u/qcluke/Chapter2/16sv4_18sv9/clustered/16sv4_clustered.fasta \
    --uc /scr/u/qcluke/Chapter2/16sv4_18sv9/clustered/16sv4_cluster.uc
# 18S V9
vsearch --cluster_size /scr/u/qcluke/Chapter2/16sv4_18sv9/dereplicated/18sv9_uniques.fasta \
    --id 0.97 \
    --sizein \
    --relabel clustered. \
    --sizeout \
    --consout /scr/u/qcluke/Chapter2/16sv4_18sv9/clustered/18sv9_clustered.fasta \
    --uc /scr/u/qcluke/Chapter2/16sv4_18sv9/clustered/18sv9_cluster.uc

## Chimera Removal
module load VSEARCH/2.18
# 16S V4
vsearch --uchime_denovo /scr/u/qcluke/Chapter2/16sv4_18sv9/clustered/16sv4_clustered.fasta \
    --sizein \
    --relabel otu.\
    --nonchimeras /scr/u/qcluke/Chapter2/16sv4_18sv9/chimera_removed/16sv4_otu.fasta \
# 18S V9
vsearch --uchime_denovo /scr/u/qcluke/Chapter2/16sv4_18sv9/clustered/18sv9_clustered.fasta \
    --sizein \
    --relabel otu. \
    --nonchimeras /scr/u/qcluke/Chapter2/16sv4_18sv9/chimera_removed/18sv9_otu.fasta \

## Count Table
module load VSEARCH/2.18
# 16S V4
vsearch --usearch_global /scr/u/qcluke/Chapter2/16sv4_18sv9/dereplicated/16sv4_combined.fasta \
    --db /scr/u/qcluke/Chapter2/16sv4_18sv9/chimera_removed/16sv4_otu.fasta \
    --id 0.97 \
    --otutabout /scr/u/qcluke/Chapter2/16sv4_18sv9/frequency_table/ont_16sv4_otu_table.tsv
# 18S V9
vsearch --usearch_global /scr/u/qcluke/Chapter2/16sv4_18sv9/dereplicated/18sv9_combined.fasta \
    --db /scr/u/qcluke/Chapter2/16sv4_18sv9/chimera_removed/18sv9_otu.fasta \
    --id 0.97 \
    --otutabout /scr/u/qcluke/Chapter2/16sv4_18sv9/frequency_table/ont_18sv9_otu_table.tsv

## Sequence Import
conda activate qiime2-amplicon-2025.7
# 16S V4
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /scr/u/qcluke/Chapter2/16sv4_18sv9/tax_assigned/16sv4_otu.fasta \
  --output-path /scr/u/qcluke/Chapter2/16sv4_18sv9/tax_assigned/16sv4_otu.qza
# 18S V9
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /scr/u/qcluke/Chapter2/16sv4_18sv9/tax_assigned/18sv9_otu.fasta \
  --output-path /scr/u/qcluke/Chapter2/16sv4_18sv9/tax_assigned/18sv9_otu.qza

## Taxonomic Assignment
conda activate qiime2-amplicon-2025.7
# 16S V4
qiime feature-classifier classify-sklearn \
  --i-classifier /scr/u/qcluke/Chapter2/16sv4_18sv9/tax_assigned/16sv4_classifier.qza \
  --i-reads /scr/u/qcluke/Chapter2/16sv4_18sv9/tax_assigned/16sv4_otu.qza \
  --p-n-jobs 8 \
  --o-classification /scr/u/qcluke/Chapter2/16sv4_18sv9/tax_assigned/ont_16sv4_taxonomy.qza
# 18S V9
qiime feature-classifier classify-sklearn \
  --i-classifier /scr/u/qcluke/Chapter2/16sv4_18sv9/tax_assigned/18sv9_classifier.qza \
  --i-reads /scr/u/qcluke/Chapter2/16sv4_18sv9/tax_assigned/18sv9_otu.qza \
  --p-n-jobs 8 \
  --o-classification /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/tax_assigned/ont_18sv9_taxonomy.qza

## Taxonomic Table
conda activate qiime2-amplicon-2025.7
# 16S V4
qiime metadata tabulate \
  --m-input-file ont_16sv4_taxonomy.qza \
  --o-visualization ont_16sv4_taxonomy.qzv
# 18S V9
qiime metadata tabulate \
  --m-input-file ont_18sv9_taxonomy.qza \
  --o-visualization ont_18sv9_taxonomy.qzv



### ONT Long ----------------------------------------------------------------------------------------------------------------

## Basecalling and Demultiplexing
# FL 16S and FL 18S
dorado basecaller \
    /scr/u/qcluke/dorado-1.0.0-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    /lustre1/g/sbs_sey/Luke/raw_data/nanopore-chapter2-raw-sequences/ONT-FL16S_FL18S/20241213_1439_MN47038_FAZ20237_a3800d91/pod5/ \
    --device cuda:all \
    --kit-name SQK-NBD114-96 \
    --emit-fastq \
    --no-trim \
    --output-dir /scr/u/qcluke/basecalled/chapter2_fl16s_fl18s/demultiplexed/ \

## Primer Trimming
module load cutadapt/3.5
# FL 16S
for i in {01..48}; do
  INPUT_FILE="/scr/u/qcluke/Chapter2/fl16s_fl18s/demultiplexed/fl16s/a3800d91-d794-452e-b1f0-fb55e70220c3_SQK-NBD114-96_barcode${i}.fastq"
  OUTPUT_FILE="/scr/u/qcluke/Chapter2/fl16s_fl18s/removed/fl16s/barcode${i}_removed.fastq"
  LOG_FILE="/scr/u/qcluke/Chapter2/fl16s_fl18s/removed/log_cutadapt_fl16s.txt"
  echo "=========== Trimming barcode${i} ==========="
  cutadapt ${INPUT_FILE} \
    --cores 8 \
    --adapter AGRGTTYGATYMTGGCTCAG...AAGTCGTAACAAGGTARCY \
    --revcomp \
    --error-rate 0.3 \
    --match-read-wildcards \
    --action trim \
    --discard-untrimmed \
    --minimum-length 1350 \
    --maximum-length 1550 \
    --report=full \
    --output ${OUTPUT_FILE} >> ${LOG_FILE}
done
# FL 18S
for i in {49..96}; do
  INPUT_FILE="/scr/u/qcluke/Chapter2/fl16s_fl18s/demultiplexed/fl18s/a3800d91-d794-452e-b1f0-fb55e70220c3_SQK-NBD114-96_barcode${i}.fastq"
  OUTPUT_FILE="/scr/u/qcluke/Chapter2/fl16s_fl18s/removed/fl18s/barcode${i}_removed.fastq"
  LOG_FILE="/scr/u/qcluke/Chapter2/fl16s_fl18s/removed/log_cutadapt_fl18s.txt"
  echo "=========== Trimming barcode${i} ==========="
  cutadapt ${INPUT_FILE} \
    --cores 8 \
    --adapter AACCTGGTTGATCCTGCCAGT...GTAGGTGAACCTGCAGAAGGATCA \
    --revcomp \
    --error-rate 0.3 \
    --match-read-wildcards \
    --action trim \
    --discard-untrimmed \
    --minimum-length 1660 \
    --maximum-length 1850 \
    --report=full \
    --output ${OUTPUT_FILE} >> ${LOG_FILE}
done

## Length and Quality Filtering
conda activate chopper
# FL 16S
for i in {01..48}; do
    INPUT_FILE="/scr/u/qcluke/Chapter2/fl16s_fl18s/removed/fl16s/barcode${i}_removed.fastq"
    OUTPUT_FILE="/scr/u/qcluke/Chapter2/fl16s_fl18s/filtered/fl16s/barcode${i}_filtered.fastq"
    echo "filtering for: barcode${i}"
    chopper --quality 20 \
    --maxqual 50 \
    --minlength 1350 \
    --maxlength 1550 \
    --threads 4 \
    --input ${INPUT_FILE} > ${OUTPUT_FILE}
done
# FL 18S
for i in {49..96}; do
    INPUT_FILE="/scr/u/qcluke/Chapter2/fl16s_fl18s/removed/fl18s/barcode${i}_removed.fastq"
    OUTPUT_FILE="/scr/u/qcluke/Chapter2/fl16s_fl18s/filtered/fl18s/barcode${i}_filtered.fastq"
    echo "filtering for: barcode${i}"
    chopper --quality 20 \
    --maxqual 50 \
    --minlength 1650 \
    --maxlength 1850 \
    --threads 4 \
    --input ${INPUT_FILE} > ${OUTPUT_FILE}
done

## Conversion and Dereplication
module load VSEARCH/2.18
# FL 16S
for i in /scr/u/qcluke/Chapter2/fl16s_fl18s/filtered/fl16s/*_filtered.fastq; do
  BASE=$(basename "$i" _filtered.fastq)
  echo "========== Converting ${BASE} =========="
  vsearch --fastq_filter /scr/u/qcluke/Chapter2/fl16s_fl18s/filtered/fl16s/${BASE}_filtered.fastq \
    --fastq_ascii 33 \
    --fastq_qmax 50 \
    --fastaout /scr/u/qcluke/Chapter2/fl16s_fl18s/generated_fasta/fl16s/${BASE}.fasta \
    --relabel ${BASE}.
done
cat /scr/u/qcluke/Chapter2/fl16s_fl18s/generated_fasta/fl16s/*.fasta \
    > /scr/u/qcluke/Chapter2/fl16s_fl18s/dereplicated/fl16s_combined.fasta
vsearch --derep_fulllength /scr/u/qcluke/Chapter2/fl16s_fl18s/dereplicated/fl16s_combined.fasta \
    --output /scr/u/qcluke/Chapter2/fl16s_fl18s/dereplicated/fl16s_uniques.fasta \
    --relabel uniq. \
    --sizeout \
    --uc /scr/u/qcluke/Chapter2/fl16s_fl18s/dereplicated/fl16s_derep.uc
# FL 18S
for i in /scr/u/qcluke/Chapter2/fl16s_fl18s/filtered/fl18s/*_filtered.fastq; do
  BASE=$(basename "$i" _filtered.fastq)
  echo "========== Converting ${BASE} =========="
  vsearch --fastq_filter /scr/u/qcluke/Chapter2/fl16s_fl18s/filtered/fl18s/${BASE}_filtered.fastq \
    --fastq_ascii 33 \
    --fastq_qmax 50 \
    --fastaout /scr/u/qcluke/Chapter2/fl16s_fl18s/generated_fasta/fl18s/${BASE}.fasta \
    --relabel ${BASE}.
done
cat /scr/u/qcluke/Chapter2/fl16s_fl18s/generated_fasta/fl18s/*.fasta \
    > /scr/u/qcluke/Chapter2/fl16s_fl18s/dereplicated/fl18s_combined.fasta
vsearch --derep_fulllength /scr/u/qcluke/Chapter2/fl16s_fl18s/dereplicated/fl18s_combined.fasta \
    --output /scr/u/qcluke/Chapter2/fl16s_fl18s/dereplicated/fl18s_uniques.fasta \
    --relabel uniq. \
    --sizeout \
    --uc /scr/u/qcluke/Chapter2/fl16s_fl18s/dereplicated/fl18s_derep.uc

## Clustering
module load VSEARCH/2.18
# FL 16S
vsearch --cluster_size /scr/u/qcluke/Chapter2/fl16s_fl18s/dereplicated/fl16s_uniques.fasta \
    --id 0.97 \
    --sizein \
    --sizeout \
    --relabel clustered. \
    --consout /scr/u/qcluke/Chapter2/fl16s_fl18s/clustered/fl16s_clustered.fasta \
    --uc /scr/u/qcluke/Chapter2/fl16s_fl18s/clustered/fl16s_cluster.uc
# FL 18S
vsearch --cluster_size /scr/u/qcluke/Chapter2/fl16s_fl18s/dereplicated/fl18s_uniques.fasta \
    --id 0.97 \
    --sizein \
    --sizeout \
    --relabel clustered. \
    --consout /scr/u/qcluke/Chapter2/fl16s_fl18s/clustered/fl18s_clustered.fasta \
    --uc /scr/u/qcluke/Chapter2/fl16s_fl18s/clustered/fl18s_cluster.uc

## Chimera Removal
module load VSEARCH/2.18
# FL 16S
vsearch --uchime_denovo /scr/u/qcluke/Chapter2/fl16s_fl18s/clustered/fl16s_clustered.fasta \
    --sizein \
    --relabel otu. \
    --nonchimeras /scr/u/qcluke/Chapter2/fl16s_fl18s/chimera_removed/fl16s_otu.fasta \
# FL 18S
vsearch --uchime_denovo /scr/u/qcluke/Chapter2/fl16s_fl18s/clustered/fl18s_clustered.fasta \
    --sizein \
    --relabel otu. \
    --nonchimeras /scr/u/qcluke/Chapter2/fl16s_fl18s/chimera_removed/fl18s_otu.fasta \

## Count Table
module load VSEARCH/2.18
# FL 16S
vsearch --usearch_global /scr/u/qcluke/Chapter2/fl16s_fl18s/dereplicated/fl16s_combined.fasta \
    --db /scr/u/qcluke/Chapter2/fl16s_fl18s/chimera_removed/fl16s_otu.fasta \
    --id 0.97 \
    --otutabout /scr/u/qcluke/Chapter2/fl16s_fl18s/frequency_table/ont_fl16s_otu_table.tsv
# FL 18S
vsearch --usearch_global /scr/u/qcluke/Chapter2/fl16s_fl18s/dereplicated/fl18s_combined.fasta \
    --db /scr/u/qcluke/Chapter2/fl16s_fl18s/chimera_removed/fl18s_otu.fasta \
    --id 0.97 \
    --otutabout /scr/u/qcluke/Chapter2/fl16s_fl18s/frequency_table/ont_fl18s_otu_table.tsv

## Sequence Import
conda activate qiime2-amplicon-2025.7
# FL 16S
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /scr/u/qcluke/Chapter2/fl16s_fl18s/tax_assigned/fl16s_otu.fasta \
  --output-path /scr/u/qcluke/Chapter2/fl16s_fl18s/tax_assigned/fl16s_otu.qza
# FL 18S
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /scr/u/qcluke/Chapter2/fl16s_fl18s/tax_assigned/fl18s_otu.fasta \
  --output-path /scr/u/qcluke/Chapter2/fl16s_fl18s/tax_assigned/fl18s_otu.qza

## Reference Sequence Extraction
conda activate qiime2-amplicon-2025.7
# FL 16S
qiime feature-classifier extract-reads \
  --i-sequences /lustre1/g/sbs_sey/Luke/database/2024.09.backbone.full-length.fna.qza \
  --p-f-primer AGRGTTYGATYMTGGCTCAG \
  --p-r-primer RGYTACCTTGTTACGACTT \
  --p-n-jobs 8 \
  --o-reads fl16s-ref-seqs.qza
# FL 18S
qiime feature-classifier extract-reads \
  --i-sequences /lustre1/g/sbs_sey/Luke/database/pr2-5.0.0-seqs.qza \
  --p-f-primer AACCTGGTTGATCCTGCCAGT \
  --p-r-primer TGATCCTTCTGCAGGTTCACCTAC \
  --p-n-jobs 8 \
  --o-reads fl18s-ref-seqs.qza

## Classifier Traning
conda activate qiime2-amplicon-2025.7
# FL 16S
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads fl16s-ref-seqs.qza \
  --i-reference-taxonomy /lustre1/g/sbs_sey/Luke/database/2024.09.backbone.tax.qza \
  --o-classifier fl16s-classifier.qza
# FL 18S
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads fl18s-ref-seqs.qza \
  --i-reference-taxonomy /lustre1/g/sbs_sey/Luke/database/pr2-5.0.0-tax.qza \
  --o-classifier fl18s-classifier.qza

## Taxonomic Assignment
conda activate qiime2-amplicon-2025.7
# FL 16S
qiime feature-classifier classify-sklearn \
  --i-classifier /scr/u/qcluke/Chapter2/fl16s_fl18s/tax_assigned/fl16s_classifier.qza \
  --i-reads /scr/u/qcluke/Chapter2/fl16s_fl18s/tax_assigned/fl16s_otu.qza \
  --p-n-jobs 16 \
  --o-classification /scr/u/qcluke/Chapter2/16sv4_18sv9/tax_assigned/ont_fl16s_taxonomy.qza
# FL 18S
qiime feature-classifier classify-sklearn \
  --i-classifier /scr/u/qcluke/Chapter2/fl16s_fl18s/tax_assigned/fl18s-classifier.qza \
  --i-reads /scr/u/qcluke/Chapter2/16svfl16s_fl18s4_18sv9/tax_assigned/fl18s_otu.qza \
  --p-n-jobs 16 \
  --o-classification /scr/u/qcluke/Chapter2/illumina_16sv4_18sv9/tax_assigned/ont_fl18s_taxonomy.qza

## Taxonomic Table
conda activate qiime2-amplicon-2025.7
# FL 16S
qiime metadata tabulate \
  --m-input-file ont_fl16s_taxonomy.qza \
  --o-visualization ont_fl16s_taxonomy.qzv
# FL 18S
qiime metadata tabulate \
  --m-input-file ont_fl18s_taxonomy.qza \
  --o-visualization ont_fl18s_taxonomy.qzv