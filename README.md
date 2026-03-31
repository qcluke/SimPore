# SimPore
A concise bioinformatic workflow for ONT amplicon processing, both short and long reads.

Find SimPore in our paper



---------------

## Flowchart
<img src="https://raw.githubusercontent.com/qcluke/SimPore/main/.png" width="1000" />

## Step 1.1 Basecalling and Demultiplexing
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

## Step 1.2 Primer Trimming
Software: Cutadapt https://cutadapt.readthedocs.io/en/stable/

```sh
for item in /PATH/TO/DEMULTIPLXED_DIRECTORY/*; do
  barcode=$(basename "$item")
  input="/PATH/TO/3.1DEMULTIPLXED_DIRECTORY/${barcode}.fastq"
  output="/PATH/TO/TRIMMED_DIRECTORY/${barcode}.fastq"
  log="PATH/TO/TRIMMED_DIRECTORY/LOG.txt"
  echo "========== Trimming ${barcode} =========="
  cutadapt ${input} \
    --cores NUMBER \
    --adapter FORWARD_PRIMER...REVERSED_COMPLEMENTARY_REVERSE_PRIMER \
    --revcomp \
    --error-rate 0.3 \
    --match-read-wildcards \
    --action trim \
    --discard-untrimmed \
    --minimum-length NUMBER \
    --maximum-length NUMBER \
    --report=full \
    --output ${output} >> ${log}
done
```

## Step 1.3 Length and Quality Filtering
Software: Chopper https://github.com/wdecoster/chopper

```sh
for item in //PATH/TO/TRIMMED_DIRECTORY/*; do
  barcode=$(basename "$item")
  barcode=${barcode%.fastq}
  input="/PATH/TO/TRIMMED_DIRECTORY/${barcode}.fastq"
  output="/PATH/TO/FILTERED_DIRECTORY/${barcode}.fastq"
  echo "========== Filtering ${barcode} =========="
  chopper --quality 20 \
    --maxqual 50 \
    --minlength NUMBER \
    --maxlength NUMBER \
    --threads NUMBER \
    --input ${input} > ${output}
done
```

## Step 2.1 Conversion and Dereplication
Software: Vsearch https://github.com/torognes/vsearch

```sh
for item in /PATH/TO/FILTERED_FILES/*; do
  barcode=$(basename "$item")
  barcode=${barcode%.fastq}
  input="/PATH/TO/FILTERED_DIRECTORY/${barcode}.fastq"
  output="/PATH/TO/CONVERTED_DIRECTORY/${barcode}.fasta"
  echo "========== Converting ${barcode} =========="
  vsearch --fastq_filter ${input} \
    --fastq_ascii 33 \
    --fastq_qmax 50 \
    --fastaout ${output} \
    --relabel ${barcode}.
done
```

```sh
cat //PATH/TO/CONVERTED_DIRECTORY/* \
    > /PATH/TO/CONVERTED_DIRECTORY/POOL.fasta
```

```sh
vsearch --derep_fulllength /PATH/TO/CONVERTED_DIRECTORY/POOL.fasta \
  --output /PATH/TO/DEREPLICATED_DIRECTORY/UNIQUE.fasta \
  --relabel unique. \
  --sizeout \
  --uc /PATH/TO/DEREPLICATED_DIRECTORY/DEREPLICATION.uc
```

## Step 2.2 Clustering
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

## Step 2.3 Chimera Removal
Software: Vsearch https://github.com/torognes/vsearch

```sh
vsearch --uchime_denovo /PATH/TO/CLUSTERED_DIRECTORY/CLUSTER.fasta \
  --sizein \
  --relabel otu. \
  --nonchimeras /PATH/TO/REMOVED_DIRECTORY/OTU.fasta 
```

## Step 2.4 Count Table Generation
Software: Vsearch https://github.com/torognes/vsearch

```sh
vsearch --usearch_global /PATH/TO/CONVERTED_DIRECTORY/POOL.fasta \
  --db /PATH/TO/REMOVED_DIRECTORY/OTU.fasta \
  --id 0.97 \
  --otutabout /PATH/TO/TABLE_DIRECTORY/COUNT_TABLE.tsv
```

## Step 3.1 Consensus Sequence Import
Software: Qiime 2 https://amplicon-docs.qiime2.org/en/stable/

```sh
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /PATH/TO/REMOVED_DIRECTORY/OTU.fasta \
  --output-path /PATH/TO/TAXONOMY_DIRECTORY/OTU.qza
```

## Step 3.2 Reference Sequence Extraction
Software: Qiime 2 https://amplicon-docs.qiime2.org/en/stable/

```sh
qiime feature-classifier extract-reads \
  --i-sequences /PATH/TO/DATABASE_DIRECTORY/DATABSE_SEQUENCE.qza \
  --p-f-primer FORWARD_PRIMER \
  --p-r-primer REVERSE_PRIMER \
  --p-n-jobs NUMBER \
  --o-reads /PATH/TO/TAXONOMY_DIRECTORY/REFERENCE_SEQUENCE.qza 
```

## Step 3.3 Classifier Traning
Software: Qiime 2 https://amplicon-docs.qiime2.org/en/stable/

```sh
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /PATH/TO/TAXONOMY_DIRECTORY/REFERENCE_SEQUENCE.qza  \
  --i-reference-taxonomy /PATH/TO/DATABASE_DIRECTORY/DATABSE_TAXONOMY.qza  \
  --o-classifier /PATH/TO/TAXONOMY_DIRECTORY/CLASSIFIER.qza
```

## Step 3.4 Taxonomic Assignment
Software: Qiime 2 https://amplicon-docs.qiime2.org/en/stable/

```sh
qiime feature-classifier classify-sklearn \
  --i-classifier /PATH/TO/TAXONOMY_DIRECTORY/CLASSIFIER.qza \
  --i-reads  /PATH/TO/TAXONOMY_DIRECTORY/OTU.qza \
  --p-n-jobs NUMBER \
  --o-classification /PATH/TO/TAXONOMY_DIRECTORY/TAXONOMY.qza
```

## Step 3.5 Taxonomic Table Generation
Software: Qiime 2 https://amplicon-docs.qiime2.org/en/stable/

```sh
qiime metadata tabulate \
  --m-input-file /PATH/TO/TAXONOMY_DIRECTORY/TAXONOMY.qza \
  --o-visualization /PATH/TO/TABLE_DIRECTORY/TAXONOMY_TABLE.qzv
```
* Import the TAXONOMIC_TABLE.qzv to QIIME 2 View (https://view.qiime2.org/) for tsv format



