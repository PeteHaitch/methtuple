These `BAM` files were created using Bismark (`v0.10.1`) by mapping against the hg19 reference genome.  

## Single-end directional
* __Filename__: `se_directional.bam`
* __SRA__: [SRR921747](http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR921747)
* __Organism:__ _Homo sapiens_
* __Number of reads__: 100,000
* __Bismark command__: `bismark --bam --bowtie2 se_directional.fq.gz /usr/local/work/hickey/genomes/hg19_bowtie2/`



## Single-end non-directional
None included.


## Paired-end directional
* __Filename__: `pe_directional.bam`
* __SRA__: [SRR400564](http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR400564)
* __Organism__: _Homo sapiens_
* __Number of reads__: 100,000
* __Bismark command__: `bismark --bam -1 pe_directional_1.fq.gz -2 pe_directional_2.fq.gz /usr/local/work/hickey/genomes/hg19_bowtie2/`

## Paired-end non-directional
None included.
