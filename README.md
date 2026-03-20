# SimPore
A concise bioinformatic workflow for ONT amplicon processing, both short and long reads.

Find SimPore in our paper



---------------

## Flowchart
<img src="https://raw.githubusercontent.com/qcluke/main/SimPore/Flowchart.png" width="800" />

## Step 1. Basecalling and Demultiplexing
Software: Dorado https://software-docs.nanoporetech.com/dorado/latest/

```sh
dorado basecaller \
  /PATH/TO/BASECALLING_MODEL \
  /PATH/TO/POD5_DIRECTORY/ \
  --device auto \
  --kit-name KIT \
  --emit-fastq \
  --no-trim \
  --output-dir /PATH/TO/DEMULTIPLXED_DIRECTORY/

```
* Remember to isolate unclassified file

## Step 2. Primer Trimming
Software: Cutadapt https://cutadapt.readthedocs.io/en/stable/

```sh
for item in /PATH/TO/DEMULTIPLXED_DIRECTORY/*; do
  barcode=${i##*_}
  barcode=${barcode%.fastq} 
  input_file="/PATH/TO/DEMULTIPLXED_DIRECTORY/PREFIX_${barcode}.fastq"
  output_file="/PATH/TO/TRIMMED_DIRECTORY/${barcode}.fastq"
  log_file="PATH/TO/TRIMMED_DIRECTORY/LOG.txt"
  echo "========== Trimming ${barcode} =========="
  cutadapt ${input_file} \
    --cores NUMBER \
    --adapter FORWARD_PRIMER...REVERSED_COMPLEMENTARY_REVERSE_PRIMER \
    --revcomp \
    --error-rate 0.3 \
    --action trim \
    --discard-untrimmed \
    --report=full \
    --output ${output_file} >> ${log_file}
done
```

## Step 3. Length and Quality Filtering
Software: Chopper https://github.com/wdecoster/chopper

```sh
for item in //PATH/TO/TRIMMED_DIRECTORY/*; do
  barcode=${barcode%.fastq} 
  input_file="/PATH/TO/TRIMMED_DIRECTORY/${barcode}.fastq"
  output_file="/PATH/TO/FILTERED_DIRECTORY/${barcode}.fastq"
  echo "========== Filtering ${barcode} =========="
  chopper --quality 20 \
    --maxqual 50 \
    --minlength NUMBER \
    --maxlength NUMBER \
    --threads NUMBER \
    --input ${input_file} > ${output_file}
done
```

## Step 4. Transformation and Dereplication
Software: Vsearch https://github.com/torognes/vsearch

```sh
for item in /PATH/TO/FILTERED_FILES/*; do
  barcode=${barcode%_.fastq}
  input_file="/PATH/TO/FILTERED_DIRECTORY/${barcode}.fastq"
  output_file="/PATH/TO/TRAMSFORMED_DIRECTORY/${barcode}.fasta"
  echo "========== Transforming ${BARCODE} =========="
  vsearch --fastq_filter ${input_file} \
    --fastq_ascii 33 \
    --fastq_qmax 50 \
    --fastaout ${output_file} \
    --relabel ${barcode}.
done
```

```sh
cat //PATH/TO/TRAMSFORMED_DIRECTORY/*。fasta \
    > /PATH/TO/TRAMSFORMED_DIRECTORY/POOL.fasta
```

```sh
vsearch --derep_fulllength /PATH/TO/TRAMSFORMED_DIRECTORY/POOL.fasta \
  --output /PATH/TO/DEREPLICATED_DIRECTORY/UNIQUE.fasta \
  --relabel unique. \
  --sizeout \
  --uc /PATH/TO/DEREPLICATED_DIRECTORY/DEREPLICATION.uc
```

## Step 5. Clustering
Software: Vsearch https://github.com/torognes/vsearch

```sh
vsearch --cluster_size /PATH/TO/DEREPLICATED_DIRECTORY/UNIQUE.fasta \
  --id 0.97 \
  --sizein \
  --sizeout \
  --relabel cluster. \
  --consout /PATH/TO/CLUSTERED_DIRECTORY/CLUSTER.fasta \
  --uc /PATH/TO/CLUSTERED_DIRECTORY/CLUSTERING.uc
```

## Step 6. Chimera Removal
Software: Vsearch https://github.com/torognes/vsearch

```sh
vsearch --uchime_denovo /PATH/TO/CLUSTERED_DIRECTORY/CLUSTER.fasta \
  --sizein \
  --nonchimeras /PATH/TO/OUTPUT_DIRECTORY/OTU.fasta \
  --relabel otu.
```

## Step 7. Count Table Generation
Software: Vsearch https://github.com/torognes/vsearch

```sh
vsearch --usearch_global /PATH/TO/TRAMSFORMED_DIRECTORY/POOL.fasta \
  --db /PATH/TO/OUTPUT_DIRECTORY/OTU.fasta \
  --id 0.97 \
  --otutabout /PATH/TO/OUTPUT_DIRECTORY/COUNT_TABLE.tsv
```

## Step 8. Consensus Sequence Import
Software: Qiime 2 https://amplicon-docs.qiime2.org/en/stable/

```sh
Qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /PATH/TO/OUTPUT_DIRECTORY/OTU.fasta \
  --output-path /PATH/TO/TAXONOMY_DIRECTORY/OTU.qza
```

## Step 9. Reference Datebase Extraction
Software: Qiime 2 https://amplicon-docs.qiime2.org/en/stable/

```sh
qiime feature-classifier extract-reads \
  --i-sequences /PATH/TO/DATABASE_DIRECTORY/DATABSE_SEQUENCE.qza \
  --p-f-primer FORWARD_PRIMER \
  --p-r-primer REVERSE_PRIMER \
  --p-n-jobs NUMBER \
  --o-reads /PATH/TO/TAXONOMY_DIRECTORY/REFERENCE_SEQUENCE.qza 
```

## Step 10. Classifier Traning
Software: Qiime 2 https://amplicon-docs.qiime2.org/en/stable/

```sh
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /PATH/TO/TAXONOMY_DIRECTORY/REFERENCE_SEQUENCE.qza  \
  --i-reference-taxonomy /PATH/TO/DATABASE_DIRECTORY/DATABSE_TAXONOMY.qza  \
  --o-classifier /PATH/TO/TAXONOMY_DIRECTORY/CLASSIFIER.qza
```

## Step 11. Taxonomic Assignment
Software: Qiime 2 https://amplicon-docs.qiime2.org/en/stable/

```sh
qiime feature-classifier classify-sklearn \
  --i-classifier /PATH/TO/TAXONOMY_DIRECTORY/CLASSIFIER.qza \
  --i-reads  /PATH/TO/TAXONOMY_DIRECTORY/OTU.qza \
  --p-n-jobs NUMBER \
  --o-classification /PATH/TO/TAXONOMY_DIRECTORY/TAXONOMY.qza
```

## Step 12. Taxonomic Table Generation
Software: Qiime 2 https://amplicon-docs.qiime2.org/en/stable/

```sh
qiime metadata tabulate \
  --m-input-file /PATH/TO/TAXONOMY_DIRECTORY/TAXONOMY.qza \
  --o-visualization /PATH/TO/OUTPUT_DIRECTORY/TAXONOMIC_TABLE.qzv
```
* Import TAXONOMIC_TABLE.qzv to QIIME 2 View (https://view.qiime2.org/) for tsv format



