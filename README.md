[![Build Status](https://travis-ci.org/PeteHaitch/comethylation.png?branch=master)](https://travis-ci.org/PeteHaitch/comethylation)

# comethylation

`comethylation` is a methylation caller for methylation events that co-occur on the same DNA fragment from high-throughput bisulfite sequencing data, such as `methylC-seq`. A typical read from such an experiment reports a binary methylated or unmethylated measurement at multiple loci, where each read originates from a single cell. `comethylation` allows us to investigate the co-occurence of methylation marks at the level of a single cell.

# Installation and dependencies
Simply running 

```
python setup.py install
```
in the root `comethylation` directory should work for most systems.

`comethylation` is written in Python and relies upon the pysam module. Running `python setup.py install` will attempt to install pysam if it isn't found on your system. Alternatively, pysam is available from [https://code.google.com/p/pysam/](https://code.google.com/p/pysam/).

`comethylation` has only been tested on Python 2.7 with pysam v0.7.7.

# Usage
## Method
`comethylation` extracts _m-tuples_ of methylation loci that co-occur on the same read. The simplest _m-tuple_ is the 1-tuple (m = 1), which corresponds to counting the number of reads that are methylated (_M_) and unmethylated (_U_) for each cytosine in the genome; 1-tuples are the type of methylation calling performed by most methylation calling software such as Bismark's `bismark_methylation_extractor`. 

A 2-tuple (m = 2) is a pair of methylation loci; `comethylation` tabulates the number of reads that methylated at each locus in the pair (_MM_), both unmethylated (_UU_) or methylated at one locus but not the other (_MU_ or _UM_). This idea extends to 3-tuples, 4-tuples, etc. 

For a chosen value of _m_, all _m-tuples_ are extracted and tabulated across the genome. The key point is that _m-tuples_ are only constructed if all _m_ loci are within a single read, so we can investigate the co-occurence of methylation events at the level of a single cell.

## Basic usage
`comethylation` processes a single `BAM` file and works for both single-end and paired-end sequencing data. Example `BAM` files from single-end directional and paired-end directional bisulfite-sequencing experiments are available in the `data/` directory. 

Methylation measurements may be filtered by base quality or other criteria such as the mapping quality of the read or whether the read is marked as a PCR duplicate. For a full list of filtering options, please run `comethylation --help` or see the __Advanced Usage__ section below. 

Currently, the BAM file must have been created with [Bismark](http://www.bioinformatics.bbsrc.ac.uk/projects/download.html#bismark). If the data were aligned with Bismark version < 0.8.3 please use the `--oldBismark` flag. A future version of this software will support the use of BAM files created with other popular bisulfite aligners such as [BSMAP](https://code.google.com/p/bsmap/), [BSmooth](https://github.com/BenLangmead/bsmooth-align), [bwa-meth](https://github.com/brentp/bwa-meth/), [gsnap](http://research-pub.gene.com/gmap/), [last](http://last.cbrc.jp/) and [Novoalign](http://www.novocraft.com/). Support will either be provided natively in `comethylation` or via an pre-processing script `bismarkify`.

The main options to pass `comethylation` are the size of the m-tuple (`--mTuple`); the type of methylation, which is some combination of _CG_, _CHG_, _CHH_ and _CNN_ (`--methylationType`); any filters to be applied to reads or positions within reads (see below); the BAM file; and the sample name, which will be used as a prefix for all output files. 

To simultaneously study multiple methylation types the `--methylationType` parameter must be specified multiple times, e.g. to study CG and CHH methylation `--methylationType CG --methylationType CHH`.

## Output
Three output files are created and summary information is written to `STDOUT`.

The main output file is a tab-delimited file of all m-tuples (`<sampleName>.<--methylationType>.<--mTuple>.tsv`). 

Here are the 5 rows (along with the header row) from `data/se_directional.fq.gz_bismark_bt2.CG.2.tsv`, which is created by running the single-end directional example shown below:

```
chr		pos1		pos2		MM	MU	UM	UU
chr10	42383637	42383660	1	0	0	0
chr10	42383660	42383672	0	2	0	0
chr10	42383672	42383686	0	0	2	0
chr10	42385004	42385040	1	0	0	0
chr10	42385040	42385063	1	0	0	0
```
So, for the CpG 2-tuple chr10:(42,383,637, 42,383,660) we observed 1 read that was methylated at chr10:42,383,637 and methylated at chr10:42,383,660. For CpG 2-tuple chr10:(42,383,660, 42,383,672) we observed 2 reads that were both methylated at chr10:42,383,660 but unmethylated at chr10:42,383,672.

The second file (`<sampleName>.<--methylationType>_per_read.hist`) is a text histogram of the number of methylation loci per read or readpair (of the type specified by `--methylationType`) that passed the filters specified at runtime of `comethylation`.

Here is the file `data/se_directional.fq.gz_bismark_bt2.CG_per_read.hist`, which is created by running the single-end directional example shown below:

```
n       count
0       4389
1       2262
2       750
3       287
4       131
5       55
6       27
7       18
8       3
9       4
10      2
11      1
12      3
13      4
14      1
18      2
```
So, 4,389 reads aligned to a position with no CpGs while 2 reads aligned to a position with 18 CpGs.

An optional third and final file (`<sampleName.reads_that_failed_QC.txt>`) records the querynames (`QNAME`) of all reads that failed to pass quality control filters and which filter the read failed. This file may be omitted by use of the `--noFailedQCFile` flag.

Here are the first 2 rows of `data/se_directional.fq.gz_bismark_bt2.reads_that_failed_QC.txt`, which is created by running the single-end directional example shown below:

```
SRR921747.5_JONAS:2105:C0G2AACXX:4:2205:18672:56112_length=101  read has indel
SRR921747.81_JONAS:2105:C0G2AACXX:4:2205:18836:56196_length=101 read has indel
```
So, both of these reads were filtered out because they contained indels (see section "__Limitations and notes__").

## Examples
Two small example datasets are included in the `data/` directory. Included are the `FASTQ` files and the `BAM` files generated with __Bismark__ in __Bowtie2__ mode. More details of the example datasets can be found in `data/README.md`

Although the example datasets are both from directional bisulfite-sequencing protocols, `comethylation` works with data from both directional and non-directional bisulfite-sequencing protocols.


### Single-end
The following command will extract all CpG 2-tuples from the file `data/se_directional.bam`:

```
comethylation --mTuple 2 --methylationType CG data/se_directional.fq.gz_bismark_bt2.bam data/se_directional.fq.gz_bismark_bt2
```

This results in 3 files: 

* `data/se_directional.fq.gz_bismark_bt2.CG.2.tsv`
* `data/se_directional.fq.gz_bismark_bt2.CG_per_read.hist`
* `data/se_directional.fq.gz_bismark_bt2.reads_that_failed_QC.txt`



### Paired-end
Paired-end data must firstly be sorted by queryname prior to running `comethylation`. `BAM` files created by Bismark, such as `data/pe_directional.bam`, are already sorted by queryname. So, to extract all CG/CHH 3-tuples we would simply run:

```
comethylation --mTuple 3 --methylationType CG --methylationType CHH --strandSpecific data/pe_directional_1.fq.gz_bismark_bt2_pe.bam data/pe_directional_1.fq.gz_bismark_bt2_pe
```
Note that we had to specify `--strandSpecific` because __CHH__ methylation is not symmetric between the Watson and Crick strands.

This results in 3 files: 

* `data/pe_directional_1.fq.gz_bismark_bt2_pe.CG_CHH.3.tsv`
* `data/pe_directional_1.fq.gz_bismark_bt2_pe.CG_CHH_per_read.hist` 
* `data/pe_directional_1.fq.gz_bismark_bt2_pe.reads_that_failed_QC.txt`

If your paired-end BAM file is sorted by genomic coordinates, then you must first sort the `BAM` by queryname and then run `comethylation` on the queryname-sorted `BAM`:

```
# Create a coordinate-sorted BAM for the sake of argument
samtools sort data/pe_directional_1.fq.gz_bismark_bt2_pe.bam data/cs_pe_directional_1.fq.gz_bismark_bt2_pe
# Re-sort the coordinate-sorted BAM by queryname
samtools sort -n data/cs_pe_directional_1.fq.gz_bismark_bt2_pe.bam data/qs_pe_directional_1.fq.gz_bismark_bt2_pe
# Run comethylation on the queryname sorted BAM
comethylation --mTuple 3 --strandSpecific --methylationType CG --methylationType CHH data/qs_pe_directional_1.fq.gz_bismark_bt2_pe.bam data/qs_pe_directional
```

## Memory usage and running time
Memory usage is dependent upon the number of methylation loci on the chromosome (more methylation loci means increased memory usage) and the value of `--mTuple` (roughly, larger values means increased memory usage), but largely independent of the number of reads in the `BAM` file. In contrast, running time is dependent on the number of reads in the `BAM` file and largely independent of the choice of `--mTuple`.

I will include more detailed performance benchmarks in future releases. For a rough indication of performance, here are the results for processing approximately 41,000,000 100bp paired-end reads from chr1 of a 20-30x coverage whole-genome methylC-seq experiment of human data. This analysis used a single AMD Opteron 6276 CPU (2.3GHz) on a shared memory system.

### `--mTuple 2`
Memory usage peaked at 2.9GB and the running time was approximately 1.5 hours. 

### `--mTuple 5`
Memory usage peaked at 8.8GB and the running time was approximately 1.5 hours.

## Helper script
I frequently work with large, coordinate-sorted `BAM` files. To speed up the extraction of m-tuples I use a simple (naive) parallelisation strategy. The idea is to split the `BAM` file into chromosome-level `BAM` files, process each chromosome-level `BAM` separately and then recombine these chromosome-level files into a genome-level file. The script `helper_scripts/run_comethylation.sh` implements this strategy; simply edit the key variables in this script.

### Warnings

* __WARNING__: This simple strategy uses as many cores as there are chromosomes. This can result in __very__ large memory usage, depending on the value of `--mTuple`, and may cause problems if you have more chromosomes than available cores.
* __WARNING__: The script `tabulate_hist.R` must be in the same directory as `run_comethylation.sh`

## Advanced usage

A full list of options is available by running `comethylation --help`:

```
usage: comethylation [-h] [--mTuple <int>] [--methylationType <string>]
                     [--oldBismark] [--ignoreDuplicates]
                     [--ignoreStart_r1 <int>] [--ignoreStart_r2 <int>]
                     [--ignoreEnd_r1 <int>] [--ignoreEnd_r2 <int>]
                     [--minQual <int>] [--minMapQ <int>] [--phred64]
                     [--overlappingPairedEndFilter <string>]
                     [--strandSpecific] [--useImproperPairs]
                     [--noFailedQCFile] [--version]
                     BAM sampleName

Extract within-fragment co-methylation measurements at methylation loci from
the aligned reads of a bisulfite-sequencing experiment. WARNING: Requires
Bismark-style BAM files including XG-, XR- and XM-tags and corrected SAM
flags. Currently only supports the directional (aka 2-strand) bisulfite-
sequencing protocol.

positional arguments:
  BAM                   The path to the SAM/BAM file
  sampleName            The name of the sample. All output files will have
                        this prefix.

optional arguments:
  -h, --help            show this help message and exit
  --mTuple <int>        The size of the methylation-loci m-tuples (i.e. the
                        choice of m); must be an integer > 1 (default: 2).
  --methylationType <string>
                        The type of methylation loci to study: CG, CHG, CHH or
                        CNN. This option may be specified multiple times in
                        order to study multiple methylation types
                        simultaneously, e.g. --methylationType CG
                        --methylationType CHG
  --oldBismark          SAM/BAM created with Bismark version < 0.8.3. The FLAG
                        and QNAME field in SAM/BAM files created by these
                        older versions of Bismark differed from the SAM
                        specifications and need to be adjusted on the fly by
                        comethylation
  --ignoreDuplicates    Ignore reads that have been flagged as PCR duplicates
                        by Picard's MarkDuplicates function
  --ignoreStart_r1 <int>
                        Ignore first <int> bases from start (5' end) of read
                        (respectively, read_1) for single-end data
                        (respectively, paired-end data) (default: 0). WARNING:
                        Parameter value not sanity checked by program.
  --ignoreStart_r2 <int>
                        Only used if data are paired-end. Ignore first <int>
                        bases from start (5' end) of read_2 (default: 0).
                        WARNING: Parameter value not sanity checked by
                        program.
  --ignoreEnd_r1 <int>  Ignore last <int> bases from end (3' end) of read
                        (respectively, read_1) for single-end data
                        (respectively, paired-end data) (default: 0). WARNING:
                        Parameter value not sanity checked by program.
  --ignoreEnd_r2 <int>  Only used if data are paired-end. Ignore last <int>
                        bases from end (5' end) of read_2 (default: 0).
                        WARNING: Parameter value not sanity checked by
                        program.
  --minQual <int>       Minimum base-quality (default: 0). Any base with a
                        lower base-quality is ignored.
  --minMapQ <int>       Minimum mapping quality mapQ (default: 0). Any read
                        with a lower mapQ is ignored. WARNING: This option has
                        no effect with SAM/BAM files created with Bismark.
                        Bismark sets all mapQ values to 255, which indicates
                        that the mapping quality is not available.
  --phred64             Quality scores are encoded as Phred64 (default:
                        Phred33).
  --overlappingPairedEndFilter <string>
                        What filter should be applied to any overlapping
                        paired-end reads. Read-pairs that don't pass the
                        filter are not used for methylation calling (options
                        listed by most-to-least stringent): check the entire
                        overlapping sequence is identical (sequence), check
                        the XM-tag is identical for the overlapping region
                        (XM), do no check of the overlapping bases but use the
                        read with the higher quality basecalls in the
                        overlapping region (quality), do no check of the
                        overlapping bases and just use the overlapping bases
                        from read_1 a la bismark_methylation_extractor
                        (bismark) (default: XM).
  --strandSpecific      Produce strand-specific counts, i.e. don't collapse
                        methylation calls across Watson and Crick strands
  --useImproperPairs    Do not filter out improper readpairs. The definition
                        of a proper readpair is aligner-specific and the value
                        set with the 0x2 bit in the SAM flag
  --noFailedQCFile      Do not create the file listing the reads that failed
                        to pass a QC filter and which filter they failed
  --version             show program's version number and exit
```


# Limitations and notes
These are current limitations and their statuses:

## Only works with data aligned with the __Bismark__ mapping software 
`comethylation` makes use of Bismark's custom SAM tags `XM`, `XR` and `XG`. The `XM` tag is used to infer the methylation state of each sequenced cytosine while the `XR` and `XG` tags are used to infer the orientation and strand of the alignment. If the data were aligned with Bismark version < 0.8.3 please use the `--oldBismark` flag.

A future version of this software will support the use of `BAM` files created with other popular bisulfite aligners such as [BSMAP](https://code.google.com/p/bsmap/), [BSmooth](https://github.com/BenLangmead/bsmooth-align), [bwa-meth](https://github.com/brentp/bwa-meth/), [gsnap](http://research-pub.gene.com/gmap/), [last](http://last.cbrc.jp/) and [Novoalign](http://www.novocraft.com/). Support will either be provided natively in `comethylation` or via an pre-processing script `bismarkify`. See [Issue #30](https://github.com/PeteHaitch/comethylation/issues/30)

## Paired-end data must be sorted by queryname
This is required in order to avoid lookups when finding the mate of a paired-end read.

The `BAM` file created by Bismark is natively in queryname order so this is not a problem. If the file is not in queryname order then use `samtools sort` with the `-n` option to sort you `BAM` by queryname. The helper script `helper_scripts/run_comethylation.sh` works with a coordinate-sorted `BAM` file and so includes a step to sort the chromosome-level `BAM` files by queryname.

## `comethylation` will skip any read containing an indel 
It is difficult, although not impossible, to assign coordinates to a cytosine within an indel. To avoid this complication, `comethylation` currently skips any reads containing an indel. I aim to improve the handling of indels in the next release.

## The `--oldBismark` option is a bit crude
Specifically, it assumes that there are no '/' characters in the read names (`QNAME`) and that the BAM has not been processed with any other programs, e.g. Picard's MarkDuplicates, that might change the `FLAG` field. I am happy to improve this functionality if requested.

## Memory usage can be large 
For instance, 5-20Gb per chromosome for a typical 20-30x coverage whole-genome methylC-seq experiment of human data. I think this is largely due to inefficiencies in how I store the m-tuples internally in `comethylation` (which is basically as a dictionary of dictionaries). See [Issue #64](https://github.com/PeteHaitch/comethylation/issues/64)

## Other notes

* Bismark always sets the mapping quality (`mapQ`) as the value 255, which means unavailable in the SAM format specification. Thus the `--minMapQ` option will not have any effect for Bismark data.
* `comethylation` skip paired-end reads where either mate is unmapped.

# Acknowledgements
A big thank you to [Felix Krueger](http://www.bioinformatics.babraham.ac.uk/people.html) (the author of Bismark) for his help in understanding mapping of bisulfite-sequencing data and for answering my many questions along the way.

Thanks also to Tobias Sargeant ([@folded](https://github.com/folded)) for his help in turning the original `comethylation.py` script into the current Python module `comethylation` and for help in setting up a testing framework.


# Questions and comments
Please use the GitHub Issue Tracker (available at [www.github.com/PeteHaitch/comethylation](www.github.com/PeteHaitch/comethylation)) to file bug reports or request new functionality. I welcome questions and comments; you can email me at peter.hickey@gmail.com. 
