Comethylation
=============

__Comethylation__ is a methylation caller for methylation events that co-occur on the same DNA fragment from high-throughput bisulfite sequening data, such as methylC-seq. A typical read from such an experiment reports a binary methylated or unmethylated measurement at multiple loci, where each read originates from a single cell. In theory this allows us to investigate the co-occurence of methylation marks at the level of a single cell, however, there is no publicly available software for such an analysis. __Comethylation__ is the first such software.

Method
=============
__Comethylation__ extracts _m-tuples_ of methylation loci that co-occur on the same read. The simplest _m-tuple_ is the 1-tuple (m = 1), which corresponds to counting the number of reads that are methylated (M) and unmethylated (U) for each cytosine in the genome; this is the basis for the type of methylation calling performed by current software such as Bismark's `bismark_methylation_extractor`. 

A 2-tuple (m = 2) is a pair of methylation loci; __Comethylation__ tabulates the number of reads that methylated at each locus in the pair (MM), both unmethylated (UU) or methylated at one locus but not the other (MU or UM). This idea extends to 3-tuples, 4-tuples, etc. 

For a chosen value of _m_, all _m-tuples_ are extracted and tabulated across the genome. The key point is that m-tuples are only constructed if all m loci are within a single read, so we can investigate the co-occurence of methylation events at the level of a single cell.

__Comethylation__ works with both single-end and paired-end sequencing data. Methylation measurements may be filtered by base quality or other criteria such as the mapping quality of the read or whether the read is marked as a PCR duplicate. For a full list of filtering options, please run `python comethylation.py --help`.


Installation and requirements
=============
__Comethylation__ is written in Python and relies upon the Pysam module. Users must install Pysam before using __Comethylation__. Pysam is available from https://code.google.com/p/pysam/.

The script `comethylation.py` extracts and tabules m-tuples from a single SAM/BAM file of bisulfite sequencing data. 

The SAM/BAM file must have been created with Bismark because `comethylation.py` uses Bismark's custom `XM`, `XR` and `XG` SAM tags. The `XM` tag is used to infer the methylation state of each sequenced cytosine while the `XR` and `XG` tags are used to infer the orientation and strand of the alignment. If the data were aligned with Bismark version < 0.8.3 please use the `--oldBismark` flag. 

__TODO__
But don't despair if your data were not aligned with Bismark. Althought other aligners are not directly supported, the included script `bismarkify.py` will add the required `XR`, `XG` and `XM` tags to a SAM/BAM file that was created using other major aligners. `bismarkify.py` will work with SAM/BAM files created with __TODO: List the supported aligners__, 

Usage
=============
__TODO__


Examples
=============
__TODO__

Limitations
=============
__TODO__
* __Comethylation__ only works with data aligned with the Bismark mapping software.
* __Comethylation__ can only process _directional_ (aka _2-strand_ or _Lister protocol_) bisulfite-sequencing data. It will not work with _non-directional_ (aka _4-strand_ or _Cokus protocol_) bisulfite-sequencing data, nor will it work with PBAT data.

Questions
=============
Please use the GitHub Issue Tracker (available at www.github.com/PeteHaitch/Comethylation) to file bug reports or request new functionality. I welcome questions and comments; you can email me at peter.hickey@gmail.com. 
