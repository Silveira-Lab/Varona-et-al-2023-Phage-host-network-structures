# Varona-et-al-2023-Phage-host-bioinformatic-preprocessing
Supplementary materials and code utilised in the the publication by Varona et al 2023


# Sequence data quality control
Prior to any read quality control and filtration the reads were visualised with <i>FastQC</i>, head and tail trimming parameters were identified and quality control began with <i>bbduk</i>.

```bash
# Identify trimming and filtering parameters
fastqc sample001_R1.fastq.gz sample001_R2.fastq.gz

## QC reads with bbduk
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

## Contig assembly with SPAdes
Assembling qc'd reads to contigs occured with SPAdes using --meta flag.
The Assembler was run as a loop for each sample to assemble by sample.
```bash
for f in *_1_outfinal.fastq; do 
        name=$(echo ${f} | sed 's/_1_outfinal.fastq//g') 
        echo ${name}
        mkdir ${name}contigs.fastq 
        spades.py --meta --pe1-1 /PATH/TO/QC_READS/${name}_1_outfinal.fastq  --pe1-2 /PATH/TO/QC_READS/${name}_2_outfinal.fastq --only-assembler -o /PATH/TO/CONTIGS/${name}contigs.fasta
done
```

## Mapping reads to contigs
qc'd reads from each sample were mapped to contigs using bowtie2 to obtain sequence alignment information and then sorted using samtools
```bash

#map qc'd reads to contigs
for f in CUR21-*_1.fastq ; do
        echo ${f}
        name=$(echo ${f} | sed 's/_1.fastq//g')
        echo ${name}
bowtie2 -x /PATH/TO/CONTIGS/${name}contigs \
-1 ${name}_1.fastq \
-2 ${name}_2.fastq \
> /OUTPUT/PATH/${name}.alignment.sam
samtools sort -o /OUTPUT/PATH2/${name}.alignment.sorted.bam /OUTPUT/PATH/${name}.alignment.sam;
samtools index /OUTPUT/PATH2/${name}.alignment.sorted.bam
done

```

## Binning bacterial contigs
contigs were binned using concoct and maxbin

```bash

#to bin, 2 binners were used, concoct & metabat
#note that these binners were also carried out with contigs from cross-sample-assemblies that were generated from:
megahit --presets meta-large -1 List,of,R1,qc-d,fatqs \
-2 list,of,R2,qc-d,fastqs --min-contig-len 1000 -m 0.85 \
-o megahit/ -t 16

### CONCOCT
  # split the reads from your contgs into 10KB chunks and create a bed file
  # to latter re-merge them
  cut_up_fasta.py flye_${sample}/assembly.fasta -c 10000 -o 0 --merge_last -b ${sample}_contigs_10K.bed > ${sample}_contigs_10K.fa
  # generate the coverage table
  concoct_coverage_table.py ${sample}_contigs_10K.bed ${sample}.aln.sort.bam > ${sample}_coverage_table.csv
  # run concoct
  concoct --composition_file ${sample}_ccontigs_10K.fa --coverage_file ${sample}_coverage_table.csv -b ${sample}_concoct-bins/
  # Merge subcontig clustering into original contig clustering
  merge_cutup_clustering.py ${sample}_concoct-bins/clustering_gt1000.csv > ${sample}_concoct-bins/clustering_merged.csv
  # Extract bins as individual FASTA
  extract_fasta_bins.py flye_${sample}/assembly.fasta ${sample}_concoct-bins/clustering_merged.csv --output_path ${sample}_concoct-bins/
        # CONCOCT output the bin name as a "number".fa (e.g. "1.fa")
        # which a lot of software is NOT going to like filenames starting
        # number -_- (change to "bin1.fa")
        for f in ${sample}_concoct-bins//*.fa; do
          mv $f ${sample}_concoct-bins/concoct-bin-"$(basename $f)"
        done
       
  ### MetaBAT2
 # generate the depth profiles
  jgi_summarize_bam_contig_depths --outputDepth ${sample}_depth.txt ${sample}.aln.sort.bam
  # run metabat2
  metabat2 -m 1500 -t $threads -i flye_${sample}/assembly.fasta -a ${sample}_depth.txt -o ${sample}_metabat2-bins/

### MetaWRAP bin refinement
  metawrap bin_refinement -o ${sample}_refined-bins -c 10 -x 70 -t $threads -A ${sample}_concoct-bins/ -B {sample}_maxbin2-bins/ -C ${sample}_metabat2-bins/

# CheckM
checkm lineage_wf -x fa --pplacer_threads 64 -t 64 bins/ checkm-out --tab_table --file concoct_checkM.csv
conda deactivate

```
## Dereplicating and generating bacterial clusters
Dereplicate bMAGs at 95%
```bash
# performed using anvio dereplicate genomes feature
anvi-dereplicate-genomes --ani-dir /PATH/TO/bacterial_bins/ \ 
                         -o /PATH/TO/95_perc_bacterial_bins \
                         --program fastani --similiarity-threshold 0.95

```

## Identification of viral contigs and vMAGs
Viral contigs were identified with VIBRANT, dereplicated using cd-hit-est and then assembled to bins using vRhyme
```bash
#before inputting to VIBRANT, catenate all contigs from different samples into one file
#run vibrant
cat *contigs.fasta > all.contigs.fasta
VIBRANT_run.py -i /PATH/TO/CONTIGS/all.contigs.fasta

#dereplicate the viral contigs at 95%
#!/bin/bash
cd-hit-est -i /PATH/TO/VIBRANT/CONTIGS -o /PATH/TO/95_perc_viral_contigs.fasta -M 20000 -c 0.95 -aS 0.85

#run vRhyme to bin your viral contigs
# -i: "input contigs:
# --method longest : dereplicates same scaffolds and keeps longest
# -r: path/to/your/forwardandreverse.reads (note vRhyme can also use .sam or bam files)
# -o: output/directory/
vRhyme -i /PATH/TO/95_perc_viral_contigs.fasta -o /PATH/TO/VIBRANT/OUTPUT/
-r /PATH/TO/QC_READS/*.fastq.gz \
--method longest
```

## Dereplicating and generating viral clusters
```bash
#inorder to dereplicate using virathon, bins need to be N-linked using vRhyme's N linkage script
link_bin_sequences.py -i vRhyme_best_bins_fasta \
-o vRhyme_vMAGs_N_linked/ \
-e fasta \
-n 1000 \
-c N

#quality check viruses by running checkV
cat vRhyme_vMAGs_linked/* > checkv/all_linked_vRhyme_vMAGs.fasta
checkv end_to_end all_linked_vRhyme_vMAGs.fasta checkv/ -t 16 -d /PATH/TO/CHECKV/checkv-db-v1.4

#note: to handle N-linkages programs using prodigal were modified to use -m flag 

#dereplicate using virathon
python3 /nethome/nsv19/anaconda3/envs/virathon/share/virathon/Virathon.py --genome_files all_linked_vRhyme_vMAGs.fasta --make_pops True --threads 24

```

## bMAG and vMAG abundances
Fraction of reads mapped to viral and bacterial bin were caluclated by mapping qc'd reads to bins with smalt
```bash

#make a smalt index 
smalt index YOUR_INDEX BINs.fasta

#then run smalt
for R1 in *_1.fastq.gz; do
R2=$(echo $R1 | sed 's/_1.fastq.gz/_2.fastq.gz/g')
name=$(basename $R1 _1.fastq.gz)
echo $name 
echo $R2 
        smalt map -n 40 -y 0.95 -o $name.aln.sam YOUR_INDEX $R1 $R2
        samtools sort -o $name.aln.sort.bam $name.aln.sam
        rm $name.aln.bam
        samtools index $name.aln.sort.bam
        samtools idxstats $name.aln.sort.bam > $name.idxstats
done

```


