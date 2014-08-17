# Version 1.3.0

## New or improved features

* Support reads containing indels (closes [#25](https://github.com/PeteHaitch/comethylation/issues/25) and [#45](https://github.com/PeteHaitch/comethylation/issues/45)).
* Refined process for dealing with overlapping reads (closes [#15](https://github.com/PeteHaitch/comethylation/issues/15)).
* Improved documentation (closes [#82](https://github.com/PeteHaitch/comethylation/issues/82) and [#84](https://github.com/PeteHaitch/comethylation/issues/84)_.
* Python 3.4 support.

## Minor improvements

* Added NEWS file.

## Bug fixes

* Provide default `--ir1p` and `--ir2p` (closes [#79](https://github.com/PeteHaitch/comethylation/issues/79)).
* Switched to *nix style line endings (\n), rather than the Windows line endings (\r\n), for the .tsv file (closes [#83](https://github.com/PeteHaitch/comethylation/issues/83)).
* Switched back to samtools instead of Picard in helper_scipts/run_comethylation.sh (closes [#78](https://github.com/PeteHaitch/comethylation/issues/78)).
* Updated examples to use Bismark v0.12.5.

## Internal

* Simplified creation of m-tuples.
* Added `get_positions` function to get the reference positions of all bases in a read (including `None` if no reference base, such as inserted bases or soft-clipped read-positions).
* Added `process_overlap` function to process overlapping paired-end reads using the given `--overlap-check`

## :(

* Running time of `v1.3` is 2-3x slower than `v1.2`. Partly because reads with indels are no longer skipped and partly because of code simplifications. Should be able to improve performance now that some internals have been simplified.

# Version 1.2.0

## New or improved features

* Python 3 support (not yet Python 3.4 compatible due to limitations of current version of Pysam).
* Include strand of m-tuple in output.

## Minor improvements
* Updated examples to use Bismark v0.12.2 (closes [#65](https://github.com/PeteHaitch/comethylation/issues/65)).
* Update `data/run_comethylation.sh` to use command line interface introduced in `v1.1`.
* Improved filtering of specific read-positions.

## Bug fixes
* Fixed bug that meant the sum of the overlap scores was incorrect for paired-end reads aligned to the OT-strand.

## Internal

# Version 1.1.0

## New or improved features

* Complete re-design and simplification of command line interface (closes [#74](https://github.com/PeteHaitch/comethylation/issues/74))
* Reduced memory usage (addresses [#64](https://github.com/PeteHaitch/comethylation/issues/64)).
* Added support for gzip and bzip2 output files (closes [#69](https://github.com/PeteHaitch/comethylation/issues/69)).
* Added parallelisation via GNU parallel in `helper_scripts/run_comethylation.sh`(closes [#68](https://github.com/PeteHaitch/comethylation/issues/68)).
* Ignore specific read-positions rather than just start or end (closes [#71](https://github.com/PeteHaitch/comethylation/issues/71))

## Bug fixes
* `--methylationType CHG` requires `--strandSpecific` (closes [#72](https://github.com/PeteHaitch/comethylation/issues/72)).

## Internal

* Refactored `WithinFragmentComethylationMTuple` class as part of memory usage improvements.

# Version 1.0.0

Initial public release
