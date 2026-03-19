# SimPore
A concise bioinformatic workflow for ONT amplicon processing, both short and long reads.

## Flowchart




---------------
# Step 1. Basecalling and Demultiplexing
Software: Dorado https://software-docs.nanoporetech.com/dorado/latest/

```sh
dorado basecaller \
  /PATH/BASECALLING_MODEL \
  /PATH/POD5_DIRECTORY/ \
  --device auto \
  --kit-name KIT \
  --emit-fastq \
  --no-trim \
  --output-dir /PATH/DEMULTIPLXED_DIRECTORY/

```
Remember to isolate unclassified file

# Step 2. Primer Trimming
Software: Cutadapt https://cutadapt.readthedocs.io/en/stable/

```sh
for item in /PATH/DEMULTIPLXED_DIRECTORY/*; do
  BARCODE=${i##*_}
  BARCODE=${barcode%.fastq} 
  INPUT_FILE="/PATH/DEMULTIPLXED_DIRECTORY/FILE_NAME_PREFIX_${BARCODE}.fastq"
  OUTPUT_FILE="/PATH/TO/OUTPUT_DIRECTORY/${BARCODE}_trimmed.fastq"
  LOG_FILE="PATH/LOG_FILE"
  echo "==========Trimming ${BARCODE} =========="
  cutadapt ${INPUT_FILE} \
    --cores NUMBER \
    --adapter FORWARD_PRIMER...REVERSED_COMPLEMENTARY_REVERSE_PRIMER \
    --revcomp \
    --error-rate 0.3 \
    --action trim \
    --discard-untrimmed \
    --report=full \
    --output ${OUTPUT_FILE} >> ${LOG_FILE}
done
```

# Step 3. Length and Quality Filtering
Software: Chopper https://github.com/wdecoster/chopper

```sh
for item in //PATH/TO/TRIMMED_FILES/*; do
  BARCODE=${barcode%_trimmed.fastq} 
  INPUT_FILE="/PATH/TO/TRIMMED_DIRECTORY/${BARCODE}_trimmed.fastq"
  OUTPUT_FILE="/PATH/TO/OUTPUT_DIRECTORY/${BARCODE}_filtered.fastq"
  echo "========== Filtering ${BARCODE} =========="
  chopper --quality 20 \
    --maxqual NUMBER \
    --minlength NUMBER \
    --maxlength NUMBER \
    --threads NUMBER \
    --input ${INPUT_FILE} > ${OUTPUT_FILE}
done
```

# Step 4. Transformation and Dereplication
Software: Vsearch https://github.com/torognes/vsearch

```sh
for item in /PATH/TO/FILTERED_FILES/*; do
  BARCODE=${barcode%_filtered.fastq}
  INPUT_FILE="/PATH/TO/FILTERED_DIRECTORY/${BARCODE}_filtered.fastq"
  OUTPUT_FILE="/PATH/TO/TRAMSFORMED_DIRECTORY/${BARCODE}.fasta"
  echo "========== Transforming ${BARCODE} =========="
  vsearch --fastq_filter ${INPUT_FILE} \
    --fastq_ascii 33 \
    --fastq_qmax 50 \
    --fastaout ${OUTPUT_FILE} \
    --relabel ${BARCODE}.
done
```

```sh
cat //PATH/TO/TRAMSFORMED_DIRECTORY/*。fasta \
    > /PATH/TO/TRAMSFORMED_DIRECTORY/combined.fasta
```

```sh
vsearch --derep_fulllength /PATH/TO/TRAMSFORMED_DIRECTORY/combined.fasta \
  --output /PATH/TO/DEREPLICATED_DIRECTORY/unique.fasta \
  --relabel unique. \
  --sizeout \
  --uc /PATH/TO/DEREPLICATED_DIRECTORY/dereplication.uc
```

# Step 5. Clustering
Software: Vsearch https://github.com/torognes/vsearch

```sh
vsearch --cluster_size /PATH/TO/DEREPLICATED_DIRECTORY/unique.fasta \
  --id 0.97 \
  --sizein \
  --sizeout \
  --relabel cluster. \
  --consout /PATH/TO/CLUSTERED_DIRECTORY/cluster.fasta \
  --uc /PATH/TO/CLUSTERED_DIRECTORY/clustering.uc
```

# Step 6. Chimera Removal
Software: Vsearch https://github.com/torognes/vsearch

```sh
vsearch --uchime_denovo /PATH/TO/CLUSTERED_DIRECTORY/cluster.fasta \
  --sizein \
  --nonchimeras /PATH/TO/OUTPUT_DIRECTORY/otu.fasta \
  --relabel otu.
```

# Step 7. Count Table Generation
Software: Vsearch https://github.com/torognes/vsearch

```sh
vsearch --usearch_global /PATH/TO/TRAMSFORMED_DIRECTORY/combined.fasta \
    --db /PATH/TO/OUTPUT_DIRECTORY/otu.fasta \
    --id 0.97 \
    --otutabout /PATH/TO/OUTPUT_DIRECTORY/count_table.tsv
```

# Step 8. Taxonomic Table Generation
Software: Qiime 2 https://amplicon-docs.qiime2.org/en/stable/

```sh
qiime feature-classifier extract-reads \
  --i-sequences /PATH/TO/TAXONOMY_DIRECTORY/database_sequence.qza \
  --p-f-primer FORWARD_PRIMER \
  --p-r-primer REVERSE_PRIMER\
  --p-n-jobs NUMBER \
  --o-reads /PATH/TO/TAXONOMY_DIRECTORY/reference_sequence.qza 
```

```sh
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /PATH/TO/TAXONOMY_DIRECTORY/reference_sequence.qza  \
  --i-reference-taxonomy /PATH/TO/TAXONOMY_DIRECTORY/database_taxonomy.qza  \
  --o-classifier /PATH/TO/TAXONOMY_DIRECTORY/classifier.qza
```

```sh
Qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /PATH/TO/OUTPUT_DIRECTORY/otu.fasta \
  --output-path /PATH/otu.qza
```

```sh
qiime feature-classifier classify-sklearn \
  --i-classifier /PATH/CLASSIFIER.qza \
  --i-reads  /PATH/OTU.fasta \
  --p-n-jobs NUMBER \
  --o-classification /PATH/TAXONOMY.qza
```

```sh
qiime metadata tabulate \
  --m-input-file /PATH/TAXONOMY.qza \
  --o-visualization /PATH/TAXONOMY.qzv
```
Import TAXONOMY.qzv to QIIME 2 View (https://view.qiime2.org/) for tsv format


# Citation

