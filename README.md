# methtuple

![Python package](https://github.com/PeteHaitch/methtuple/actions/workflows/python.yaml/badge.svg)

## Overview

### What does it do?

`methtuple` allows the user to investigate the co-occurence of methylation marks at the level of individual DNA fragments. It does this by performing methylation calling at _m-tuples_ of methylation loci from high-throughput bisulfite sequencing data, such as _methylC-seq_. In short, `methtuple` extracts and tabulates the methylation states of all m-tuples from a SAM/BAM file (for a user-defined value of _m_).

### Why would I want to do that?

A typical read from a bisulfite-sequencing experiment reports a binary methylated or unmethylated measurement at multiple loci. Each read originates from a single cell. Because methylation calls are made from individual reads/read-pairs, we can investigate the co-occurence of methylation events at the level of individual DNA fragments.

I have been using `methtuple` to investigate the spatial dependence of DNA methylation at the level of individual DNA fragments by studying methylation patterns of CpG 2-tuples. `methtuple` can also be used as a drop-in replacement for `bismark_methylation_extractor` while also providing enhanced filtering options and a slightly faster runtime (10-20% faster, albeit with an increased memory usage).

### What is an m-tuple?

The simplest _m-tuple_ is the 1-tuple (_m_ = 1). `methtuple` tabulates the number of reads that are methylated (_M_) and unmethylated (_U_) for each methylation 1-tuple in the genome. 1-tuples are the type of methylation calling performed by most methylation calling software such as Bismark's `bismark_methylation_extractor`.

A 2-tuple (_m_ = 2) is a pair of methylation loci. `methtuple` tabulates the number of reads that methylated at each locus in the pair (_MM_), both unmethylated (_UU_) or methylated at one locus but not the other (_MU_ or _UM_). This idea readily extends to 3-tuples, 4-tuples, etc.

In its default settings, and with _m_ > 1, `methtuple` tries to create only m-tuples made of "neighbouring" loci. However, please see the example below for why I say this only "tries" to create m-tuples of neighbouring loci. For a DNA fragment containing _k_ methylation loci there are _k - m + 1_ m-tuples made of neighbouring loci.

Alternatively, we can create all combinations of m-tuples by using the `--all-combinations` flag. For a DNA fragment containing _k_ methylation loci there are "_k_ choose _m_" m-tuples when using `--all-combinations`, a number that grows rapidly in _k_, particularly when _m_ is close to _k/2_.

Regardless of how m-tuples are constructed, `methtuple` always takes care to only count each methylation locus once when it has been twice-sequenced by overlapping paired-end reads.

### Draw me a picture

Well, I hope ASCII art will do.

Suppose we sequence a region of the genome containing five methylation loci with three paired-end reads (`A`, `B` and `C`):

```text
ref: 1    2   3 4 5
A_1: |----->
A_2:         <------|
B_1: |----->
B_2:           <----|
C_1:    |----->
C_2:      <------|
```

If we are interested in 1-tuples, then we would obtain the following from each read by running `methtuple`:

```text
A: {1}, {2}, {3}, {4}, {5}
B: {1}, {2}, {4}, {5}
C: {2}, {3}, {4}
```

This result is true regardless of whether the `--all-combinations` flag is set.

If we are interested in 3-tuples, then we would obtain the following from each read by running `methtuple` in its default mode:

```text
A: {1, 2, 3}, {2, 3, 4}, {3, 4, 5}
B: {1, 2, 4}, {2, 4, 5}
C: {2, 3, 4}
```

Things to note:

* Read-pair `A` sequences all three (= 5 - 3 + 1) "neighbouring" 3-tuples
* Read-pair `B` sequences none of the "neighbouring" 3-tuples but does "erroneously" construct two non-neighbouring 3-tuples. This happens because m-tuples are created independently from each read-pair; effectively, read-pair `B` is "unaware" of methylation locus `3`. Depending on the downstream analysis, you may want to _post-hoc_ filter out these "non-neighbouring" m-tuples.
* The twice-sequenced methylation loci, `2` and `3`, in read-pair `C` are not double counted.

However, if we were to run `methtuple` with `--all-combinations` then we would obtain:

```text
A: {1, 2, 3}, {2, 3, 4}, {3, 4, 5}, {1, 2, 4}, {1, 2, 5}, {1, 3, 4}, {1, 3, 5}, {1, 4, 5}, {2, 3, 5}, {2, 4, 5}
B: {1, 2, 4}, {2, 4, 5}, {1, 2, 5}, {1, 4, 5}
C: {2, 3, 4}
```

## Installation and dependencies

`methtuple` is written in Python and relies upon the `pysam` module. __NOTE: `methtuple` now requires `pysam v0.8.4` or greater.__

Running `python setup.py install` will attempt to install `pysam` if it isn't found on your system. Alternatively, instructions for installing `pysam` are available from [https://github.com/pysam-developers/pysam](https://github.com/pysam-developers/pysam).

I have extensively used and tested `methtuple` with Python 2.7. It should also work on Python 3.2, 3.3, 3.4, and 3.5 with the current version of `pysam` (`v0.8.4`), as indicated by the [Travis-CI builds](https://travis-ci.org/PeteHaitch/methtuple).

### Using `pip`

The simplest way:

```sh
pip install methtuple
```

`methtuple` is written in Python and requires the `pysam` module. __NOTE: `methtuple` now requires `pysam v0.8.4` or greater.__

Alternatively, after cloning or downloading the `methtuple` git repositority, simply run:

```sh
python setup.py install
```

in the root `methtuple` directory should work for most systems.

## Usage

### Basic usage

`methtuple` processes a single SAM or BAM file and works for both single-end and paired-end sequencing data. Example BAM files from single-end directional and paired-end directional bisulfite-sequencing experiments are available in the `data/` directory.

Methylation measurements may be filtered by base quality or other criteria such as the mapping quality of the read or whether the read is marked as a PCR duplicate. For a full list of filtering options, please run `methtuple --help` or see the __Advanced Usage__ section below.

Currently, the SAM/BAM file must have been created with [Bismark](http://www.bioinformatics.bbsrc.ac.uk/projects/download.html#bismark). If the data were aligned with Bismark version < 0.8.3 please use the `--aligner Bismark_old` flag. Please file an issue if you would like to use a SAM/BAM file created with another aligner and I will do my best to support it.

The main options to pass `methtuple` are the size of the m-tuple (`-m`); the type of methylation, which is some combination of _CG_, _CHG_, _CHH_ and _CNN_ (`--methylation-type`); any filters to be applied to reads or positions within reads (see below); the SAM/BAM file; and the sample name, which will be used as a prefix for all output files. Multiple methylation types may be specified jointly, e.g., `--methylation-type CG --methylation-type CHG`

### Output

Three output files are created and summary information is written to `STDOUT`. The main output file is a tab-delimited file of all m-tuples, `<in>.<--methylation-type>.<-m>[ac].tsv`, where `<in>` is the prefix of the `<in.bam>|<in.sam>` SAM/BAM file and `ac` is added if the `--all-combinations` flag was used, e.g., `SRR949207.CG.2ac.tsv`. Output files may be gzipped (`--gzip`) or bzipped (`--bzip2`).

Here are the first 5 rows (including with the header row) from `data/se_directional.fq.gz_bismark_bt2.CG.2.tsv`, which is created by running the single-end directional example shown below:

```text
chr     strand  pos1    pos2    MM      MU      UM      UU
chr1    +       6387768 6387783 1       0       0       0
chr1    +       7104116 7104139 1       0       0       0
chr1    +       7104139 7104152 1       0       0       0
chr1    +       9256170 9256179 0       0       0       1
```

So, for example, at the CpG 2-tuple chr1:+:(6,387,768, 6,387,783) we observed 1 read that was methylated at chr1:+:6,387,768 and methylated at chr1:+:6,387,783.

The `strand` is recorded as `+` (forward strand, "OT" in Bismark), `-` (reverse strand, "OB" in Bismark) or `*`, meaning not applicable (if the `--strand-collapse` option is set). The position of all methylation loci is always with respect to the forward strand.

The second file (`<in>.<--methylation-type>_per_read.hist`) is a text histogram of the number of methylation loci per read/readpair (of the type specified by `--methylation-type`) that passed the filters specified at runtime of `methtuple`.

Here is the file `data/se_directional.fq.gz_bismark_bt2.CG_per_read.hist`, which is created by running the single-end directional example shown below:

```text
n       count
0       4561
1       2347
2       789
3       296
4       144
5       61
6       29
7       19
8       3
9       4
10      2
11      1
12      3
13      4
14      1
18      2
```

So, 4,561 reads aligned to a position containing no CpGs while 2 reads aligned to a position containing 18 CpGs.

An optional third and final file (`<in>.reads_that_failed_QC.txt>`) records the querynames (`QNAME`) of all reads that failed to pass quality control filters and which filter the read failed. This file may be omitted by use of the `--no-failed-filter-file` flag.

In this case we didn't set any quality control filters and so this file is empty.

### Examples

Two small example datasets are included in the `data/` directory. Included are the FASTQ files and the SAM/BAM files generated with __Bismark__ in __Bowtie2__ mode. More details of the example datasets can be found in `data/README.md`

Although the example datasets are both from directional bisulfite-sequencing protocols, `methtuple` also works with data from non-directional bisulfite-sequencing protocols.

#### Single-end reads

The following command will extract all CpG 2-tuples from the file `data/se_directional.bam`:

```sh
methtuple -m 2 --methylation-type CG data/se_directional.fq.gz_bismark_bt2.bam
```

This results in 3 files:

* `data/se_directional.fq.gz_bismark_bt2.CG.2.tsv`
* `data/se_directional.fq.gz_bismark_bt2.CG_per_read.hist`
* `data/se_directional.fq.gz_bismark_bt2.reads_that_failed_QC.txt`

#### Paired-end reads

Paired-end data must firstly be sorted by queryname prior to running `methtuple`. BAM files created by Bismark, such as `data/pe_directional.bam`, are already sorted by queryname. So, to extract all CG/CHH 3-tuples we would simply run:

```sh
methtuple -m 3 --methylation-type CG --methylation-type CHH data/pe_directional_1.fq.gz_bismark_bt2_pe.bam
```

This results in 3 files:

* `data/pe_directional_1.fq.gz_bismark_bt2_pe.CG_CHH.3.tsv`
* `data/pe_directional_1.fq.gz_bismark_bt2_pe.CG_CHH_per_read.hist`
* `data/pe_directional_1.fq.gz_bismark_bt2_pe.reads_that_failed_QC.txt`

##### Note on sort-order of paired-end SAM/BAM files

If your paired-end SAM/BAM file is sorted by genomic coordinates, then you must first sort the SAM/BAM by queryname and then run `methtuple` on the queryname-sorted SAM/BAM. This can be done by using `samtools sort` with the `-n` option or Picard's `SortSam` function with the `SO=queryname` option:

```sh
# Create a coordinate-sorted SAM/BAM for the sake of argument
samtools sort data/pe_directional_1.fq.gz_bismark_bt2_pe.bam data/cs_pe_directional_1.fq.gz_bismark_bt2_pe
# Re-sort the coordinate-sorted BAM by queryname
samtools sort -n data/cs_pe_directional_1.fq.gz_bismark_bt2_pe.bam data/qs_pe_directional_1.fq.gz_bismark_bt2_pe
# Run methtuple on the queryname sorted BAM
methtuple -m 3 --methylation-type CG --methylation-type CHG data/qs_pe_directional_1.fq.gz_bismark_bt2_pe.bam
```

### Memory usage and running time

For a rough indication of performance, here are the results for processing approximately 40,000,000 100bp paired-end reads from chr1 of a 20-30x coverage whole-genome methylC-seq experiment of human data. This analysis used a single AMD Opteron 6276 CPU (2.3GHz) on a shared memory system.

#### `-m 2`

Memory usage peaked at 1.9GB and the running time was approximately 5 hours.

#### `-m 2 --all-combinations`

Memory usage peaked at 7GB and the running time was approximately 5.5 hours.

Use of the `--all-combinations` flag creates all possible m-tuples, including non-neighbouring ones. This produces many more m-tuples and so increases the memory usage.

#### `-m 5`

Memory usage peaked at 1.5GB and the running time was approximately 4.3 hours.

### Helper script

I frequently work with large, coordinate-sorted SAM/BAM files. To speed up the extraction of m-tuples, I use a simple parallelisation strategy with [GNU parallel](http://www.gnu.org/software/parallel/). The idea is to split the SAM/BAM file into chromosome-level SAM/BAM files, process each chromosome-level SAM/BAM separately and then recombine these chromosome-level files into a genome-level file. The script `helper_scripts/run_methtuple.sh` implements this strategy; simply edit the key variables in this script or adapt it to your own needs. Please check the requirements listed in `helper_scripts/run_methtuple.sh`.

#### Warnings

* __WARNING__: This simple strategy uses as many cores as there are chromosomes. This can result in __very__ large memory usage, depending on the value of `-m`, and may cause problems if you have more chromosomes than available cores.
* __WARNING__: The script `tabulate_hist.R` must be in the same directory as `run_methtuple.sh`

### Advanced usage

A full list of options is available by running `methtuple --help`:

```text
usage: methtuple [options] <in.bam>|<in.sam>
Please run 'methtuple -h' for a full list of options.

Extract methylation patterns at m-tuples of methylation loci from the aligned
reads of a bisulfite-sequencing experiment. Currently only supports SAM/BAM
files created with Bismark.

Input options:
  <in.bam>|<in.sam>     Input file in BAM or SAM format. Use - to specify STDIN.
                        The header must be included and alignments must have
                        been done using Bismark.
  --aligner {Bismark,Bismark_old}
                        The aligner used to generate the SAM/BAM file.
                        Bismark_old refers to Bismark version < 0.8.3 (default:
                        Bismark)
  --Phred64             Quality scores are encoded as Phred64 rather than
                        Phred33 (default: False)

Output options:
  -o PREFIX, --output-prefix PREFIX
                        By default, all output files have the same prefix as
                        that of the input file. This will override the prefix of
                        output file names
  --sc, --strand-collapse
                        Collapse counts across across Watson and Crick strands.
                        Only possible for CG methylation type. The strand is
                        recorded as '*' if this option is selected. (default:
                        False)
  --nfff, --no-failed-filter-file
                        Do not create the file listing the reads that failed to
                        pass to pass the filters and which filter it failed
                        (default: False)
  --gzip                gzip all output files. --gzip and --bzip2 are mutually
                        exclusive (default: False)
  --bzip2               bzip2 all output files. --gzip and --bzip2 are mutually
                        exclusive (default: False)

Construction of methylation loci m-tuples:
  --mt {CG,CHG,CHH,CNN}, --methylation-type {CG,CHG,CHH,CNN}
                        The methylation type. Multiple methylation types may be
                        analysed jointly by repeated use of this argument, e.g.,
                        --methylation-type CG --methylation-type CHG. The
                        default ('None') corresponds to CG (default: None)
  -m <int>              The size of the m-tuples, i.e., the 'm' in m-tuples
                        (default: 1)
  --ac, --all-combinations
                        Create all combinations of m-tuples, including non-
                        neighbouring m-tuples. WARNING: This will greatly
                        increase the memory usage, particularly for larger
                        values of -m and when analysing non-CG methylation
                        (default: False)

Filtering of reads:
  Applied before filtering of bases

  --id, --ignore-duplicates
                        Ignore reads that have been flagged as PCR duplicates
                        by, for example, Picard's MarkDuplicates function. More
                        specifically, ignore reads with the 0x400 bit in the
                        FLAG (default: False)
  --mmq <int>, --min-mapq <int>
                        Ignore reads with a mapping quality score (mapQ) less
                        than <int> (default: 0)
  --of {sequence_strict,sequence,XM_strict,XM,XM_ol,quality,Bismark}, --overlap-filter {sequence_strict,sequence,XM_strict,XM,XM_ol,quality,Bismark}
                        The type of check to be performed (listed roughly from
                        most-to-least stringent): Ignore the read-pair if the
                        sequence in the overlap differs between mates
                        (sequence_strict); Ignore the overlapping region if the
                        sequence in the overlap differs between mates
                        (sequence); Ignore the read-pair if the XM-tag in the
                        overlap differs (XM_strict); Ignore the overlapping
                        region if the XM-tag in the overlap differs between
                        mates (XM); Ignore any positions in the overlapping
                        region where the XM-tags differ between the mates
                        (XM_ol); Use the mate with the higher average quality
                        basecalls in the overlapping region (quality); Use the
                        first mate of each read-pair, i.e., the method used by
                        bismark_methylation_extractor with the --no_overlap flag
                        (Bismark) (default: XM_ol)
  --uip, --use-improper-pairs
                        Use the improper read-pairs, i.e. don't filter them.
                        More specifically, check the 0x2 FLAG bit of each read;
                        the exact definition of an improper read-pair depends on
                        the aligner and alignment parameters (default: False)

Filtering of bases:
  Applied after filtering of reads

  --ir1p VALUES, --ignore-read1-positions VALUES
                        If single-end data, ignore these read positions from all
                        reads. If paired-end data, ignore these read positions
                        from just read_1 of each pair. Multiple values should be
                        comma-delimited, ranges can be specified by use of the
                        hyphen and all positions should use 1-based co-
                        ordinates. For example, 1-5,80,95-100 corresponds to
                        ignoring read-positions 1, 2, 3, 4, 5, 80, 98, 99, 100.
                        (default: None)
  --ir2p VALUES, --ignore-read2-positions VALUES
                        Ignore these read positions from just read_2 of each
                        pair if paired-end sequencing. Multiple values should be
                        comma-delimited, ranges can be specified by use of the
                        hyphen and all positions should use 1-based co-
                        ordinates. For example, 1-5,80,95-100 corresponds to
                        ignoring read-positions 1, 2, 3, 4, 5, 80, 98, 99, 100.
                        (default: None)
  --mbq <int>, --min-base-qual <int>
                        Ignore read positions with a base quality score less
                        than <int> (default: 0)

Other:
  -v, --version         show program's version number and exit
  -h, --help            show this help message and exit

methtuple (v1.7.0) by Peter Hickey (peter.hickey@gmail.com,
https://github.com/PeteHaitch/methtuple/)
```

## Limitations and notes

These are current limitations and their statuses:

### Only works with data aligned with the __Bismark__ mapping software

`methtuple` makes use of Bismark's custom SAM tags `XM`, `XR` and `XG`. The `XM` tag is used to infer the methylation state of each sequenced cytosine while the `XR` and `XG` tags are used to infer the orientation and strand of the alignment. If the data were aligned with Bismark version < 0.8.3 please use the `--oldBismark` flag.

Please file an issue if you would like to use a SAM/BAM file created with another aligner and I will do my best to support it; also, see [Issue #30](https://github.com/PeteHaitch/methtuple/issues/30)

### Paired-end data must be sorted by queryname and it is recommended that SAMtools is used to do this

Sorting paired-end data by queryname (`QNAME`) is required in order to avoid expensive lookups when finding the mate of a paired-end read. `methtuple` requires that read_1 always occurs before read_2 in each pair. Unfortunately, the SAM specifications provides no guarantee about the relative order of `read_1` and `read_2` when they have identical `QNAME` and the file is sorted by queryname.

The SAM/BAM file created by Bismark is already in the required order and so this is not an issue. However, if the SAM/BAM file produced by Bismark has sunsequently been re-sorted in some way (e.g., sorted by genomic co-ordinates for use with Picard's MarkDuplicates utility), then the file will need to be re-sorted. According to Heng Li, one of the developer's of SAMtools, "samtools sort starts to put `read_1` before `read_2` since 0.1.19, released on 03/19/2013." ([https://github.com/samtools/hts-specs/issues/5#issuecomment-106797588](https://github.com/samtools/hts-specs/issues/5#issuecomment-106797588)). Therefore, `samtools sort -n` is the recommended way to sort the SAM/BAM file by queryname in order to ensure that `read_1` occurs before `read_2`. The helper script `helper_scripts/run_methtuple.sh` works with a coordinate-sorted SAM/BAM file and does so by including a step to sort the chromosome-level SAM/BAM files by queryname using `samtools sort -n`.

### The `--aligner Bismark_old` option is a bit crude

Specifically, it assumes that there are no '/' characters in the read names (`QNAME`) and that the SAM/BAM has not been processed with any other programs, e.g. Picard's MarkDuplicates, that might change the `FLAG` field. Please file an issue or submit a pull request if you would like this improved.

### Construction of "non-neighbouring" m-tuples

As discussed in the above example, `methtuple` tries not to create "non-neighbouring" m-tuples, however, these do occur due to m-tuples being created independently from each read/read-pair. I do not make use of non-neighbouring m-tuples in my downstream analyses and so I _post-hoc_ filter these out.

If you would like the option to create all possible m-tuples, both "neighbouring" and "non-neighbouring", please let me know at [https://github.com/PeteHaitch/methtuple/issues/85](https://github.com/PeteHaitch/methtuple/issues/85) as there is a simple solution that just awaits motivation for me to implement it.

### Choice of `--overlap-filter`

The two mates of a paired-end read, `read_1` and `read_2`, often overlap in bisulfite-sequencing data. `methtuple` ensures that the overlapping sequence isn't double-counted and offers several different choices of how overlapping paired-end reads are processed via the `--overlap-filter` flag. Listed roughly from most-to-least stringent these are:

1. `sequence_strict`: Check that the entire overlapping sequence is identical; if not identical then do not use any methylation calls from the entire read-pair.
2. `sequence`: Check that the entire overlapping sequence is identical; if not identical then do not use any methylation calls from the overlap.
3. `XM_strict`: Check that the XM-tag is identical for the overlapping region; if not identical then do not use any methylation calls from the entire read-pair.
4. `XM`: Check that the XM-tag is identical for the overlapping region; if not identical then do not use any methylation calls from the overlap.
5. `XM_ol`: Check that the XM-tag is identical for the overlapping region; if not identical then exclude those positions of disagreement and count once the remaining positions in the overlap.
6. `quality`: No check of the overlapping bases; simply use the read with the higher average quality basecalls in the overlapping region.
7. `Bismark`: No check of the overlapping bases; simply use the overlapping bases from read_1, i.e., the method used by `bismark_methylation_extractor` with the `--no_overlap` flag.

### Other notes

* Bismark-Bowtie1 always sets the mapping quality (`mapQ`) as the value 255, which means unavailable in the SAM format specification. Thus the `--min-mapq` option will not have any effect for Bismark-Bowtie1 data.
* `methtuple` skips paired-end reads where either mate is unmapped.

## Acknowledgements

A big thank you to [Felix Krueger](http://www.bioinformatics.babraham.ac.uk/people.html) (the author of Bismark) for his help in understanding mapping of bisulfite-sequencing data and for answering my many questions along the way.

Thanks also to Tobias Sargeant ([@folded](https://github.com/folded)) for his help in turning the original `methtuple.py` script into the current Python module `methtuple` and for help in setting up a testing framework.

## Questions and comments

Please use the [GitHub Issue Tracker](www.github.com/PeteHaitch/methtuple) to file bug reports or request new functionality. I welcome questions and comments; you can email me at <peter.hickey@gmail.com>.
