These `BAM` files were created using Bismark (`v0.12.2`) and mapping against the hg19 reference genome. The `FASTQ` files and the __Bismark__ mapping reports are also included for completeness. Only the first 10,000 reads or read-pairs from each file were used for these examples.

## Single-end directional
* __Filename__: `se_directional.fq.gz_bismark_bt2.bam`
* __Short Read Archive accession number__: [SRR921747](http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR921747)
* __Organism:__ _Homo sapiens_
* __Number of reads (mapped/total)__: 8,120/10,000
* __Bismark command__: `bismark --bam --bowtie2 /usr/local/work/hickey/genomes/hg19_bowtie2/ se_directional.fq.gz`

## Single-end non-directional
None included.

## Paired-end directional
* __Filename__: `pe_directional_1.fq.gz_bismark_bt2_pe.bam`
* __Short Read Archive accession number__: [SRR400564](http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR400564)
* __Organism__: _Homo sapiens_
* __Number of readpairs (mapped/total)__: 6,655/10,000
* __Bismark command__: `bismark --bam --bowtie2 /usr/local/work/hickey/genomes/hg19_bowtie2/ -1 pe_directional_1.fq.gz -2 pe_directional_2.fq.gz`

## Paired-end non-directional
None included.
