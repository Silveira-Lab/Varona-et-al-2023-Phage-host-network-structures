# Varona-et-al-2023-Phage-host-network-structures
Supplementary materials and code utilised in the the publication by Varona et al 2023


# Sequence data quality control
Prior to any read quality control and filtration the reads were visualised with <i>FastQC</i>, head and tail trimming parameters were identified and quality control began with <i>bbduk</i>.

```bash
# Identify trimming and filtering parameters
fastqc sample001_R1.fastq.gz sample001_R2.fastq.gz

# QC reads with bbduk
for f in *_R1.fastq.gz; do
      name=$(basename $f R1_001.fastq.gz)
   # Adaptor trimming 
      bbduk.sh -Xmx512m -da \
      in1=${name}_R1.fastq.gz in2=${name}_R2.fastq.gz \
      out1=${name}_QC_R1.fastq.gz out2=${name}_QC_R2.fastq.gz \
      ktrim=rl k=23 mink=11 hdist=1 qtrim=rl trimq=30 \
      ref=/path/to/adapters.fa # make sure the path to the Illumina adaptors is correct!
      
    # Quality filtering: This will discard reads with average quality below 30.
      # If quality-trimming is enabled, the average quality will be calculated
      # on the trimmed read.
      bbduk.sh -Xmx512m -da \
      in1=${name}_1out1.fastq.gz in2=${name}_2out1.fastq.gz \
      out1=${name}_1out2.fastq.gz out2=${name}_2out2.fastq.gz maq=30
      # Entropy filtering
      bbduk.sh -Xmx512m -da \
      in1=${name}_1out2.fastq.gz in2=${name}_2out2.fastq.gz \
      out1=${name}_1_outfinal.fastq.gz out2=${name}_2_outfinal.fastq.gz \
      entropy=0.90
  done
  
 # After QC with bbduk, the reads were passed through FastQC a final time to confirm the QC was sucessful.
for R1 in *_QC_R1.fastq.gz; do
fastqc sample001_R1.fastq.gz sample001_R2.fastq.gz
done
```
