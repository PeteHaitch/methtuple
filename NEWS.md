# Version 1.8.0

- Fixed parsing of reads containing insertions (closes [#104](https://github.com/PeteHaitch/methtuple/issues/104))
- Fixed handling of gzip and bzip output files and filenames
- Remove support for Python 2.7. It may still work, but I'm no longer going to try to maintain compatability with both Python 2 and 3.

# Version 1.7.0

- Fixed parsing of `--methylation-type` argument (closes [#102](https://github.com/PeteHaitch/methtuple/issues/102))

# Version 1.5.4

`methtuple` now requires pysam `v0.8.4` or greater.
`methtuple` now support Python 3.5.

## Bug fixes

- Fixed failing tests due to changes between pysam v0.8.3 and v0.8.4.

## Internal

- `methtuple` now requires pysam `v0.8.4` or greater

# Version 1.5.3

`methtuple` now requires pysam `v0.8.3` or greater.

## Bug fixes

- Fixed failing tests due to changes between pysam v0.8.2 and v0.8.3.

## Internal

- `methtuple` now requires pysam `v0.8.3` or greater

# Version 1.5.0

`methtuple` now requires `pysam v0.8.1` or greater.

## Bug fixes

* Fixed Python 3 compatibility issues due to changes in the `pysam` API (closes [#94](https://github.com/PeteHaitch/methtuple/issues/94)).

## Internal

Moved to the new `pysam` API. As a result, `methtuple` now requires `pysam v0.8.1` or greater.

# Version 1.4.0

This version of the software changes the name from `comethylation` to `methtuple` to better reflect what it actually does.

I have also switched to the MIT licence.

## New or improved features

* Added `--all-combinations` option to create all possible m-tuples including non-neighbouring ones. This will greatly increase the memory usage if combined with larger values of `--m` (e.g. > 2).

# Version 1.3.0

## New or improved features

* Support reads containing indels (closes [#25](https://github.com/PeteHaitch/comethylation/issues/25) and [#45](https://github.com/PeteHaitch/comethylation/issues/45)).
* Refined process for dealing with overlapping reads (closes [#15](https://github.com/PeteHaitch/comethylation/issues/15)).
* Improved documentation (closes [#82](https://github.com/PeteHaitch/comethylation/issues/82) and [#84](https://github.com/PeteHaitch/comethylation/issues/84)).
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
