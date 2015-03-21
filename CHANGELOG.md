# Change Log

## [Unreleased](https://github.com/PeteHaitch/methtuple/tree/HEAD)

[Full Changelog](https://github.com/PeteHaitch/methtuple/compare/v1.5...HEAD)

**Implemented enhancements:**

- Support coveralls.io [\#96](https://github.com/PeteHaitch/methtuple/issues/96)

- Upgrade to pysam v0.8.2 [\#95](https://github.com/PeteHaitch/methtuple/issues/95)

## [v1.5](https://github.com/PeteHaitch/methtuple/tree/v1.5) (2015-02-05)

[Full Changelog](https://github.com/PeteHaitch/methtuple/compare/v1.4.1...v1.5)

**Closed issues:**

- Update to pysam 0.8.1 [\#94](https://github.com/PeteHaitch/methtuple/issues/94)

## [v1.4.1](https://github.com/PeteHaitch/methtuple/tree/v1.4.1) (2014-09-04)

[Full Changelog](https://github.com/PeteHaitch/methtuple/compare/v1.4...v1.4.1)

## [v1.4](https://github.com/PeteHaitch/methtuple/tree/v1.4) (2014-08-29)

[Full Changelog](https://github.com/PeteHaitch/methtuple/compare/v1.3...v1.4)

**Implemented enhancements:**

- Option to create all possible m-tuples [\#85](https://github.com/PeteHaitch/methtuple/issues/85)

- Improve README [\#82](https://github.com/PeteHaitch/methtuple/issues/82)

- Add support for reads containing INDELs [\#25](https://github.com/PeteHaitch/methtuple/issues/25)

- Refine process for dealing with overlapping read-pairs \(part 2\) [\#15](https://github.com/PeteHaitch/methtuple/issues/15)

**Fixed bugs:**

- `comethylation` may not work with soft-clipped reads or reads containing indels [\#45](https://github.com/PeteHaitch/methtuple/issues/45)

**Merged pull requests:**

- Add ac to filename of tsp files created with —all-combinations. [\#92](https://github.com/PeteHaitch/methtuple/pull/92) ([PeteHaitch](https://github.com/PeteHaitch))

- All pairs [\#91](https://github.com/PeteHaitch/methtuple/pull/91) ([PeteHaitch](https://github.com/PeteHaitch))

## [v1.3](https://github.com/PeteHaitch/methtuple/tree/v1.3) (2014-08-17)

[Full Changelog](https://github.com/PeteHaitch/methtuple/compare/v1.2...v1.3)

**Implemented enhancements:**

- Note regarding --of Bismark [\#84](https://github.com/PeteHaitch/methtuple/issues/84)

- Add native parallel support to comethylation [\#76](https://github.com/PeteHaitch/methtuple/issues/76)

- Supporting aligners besides Bismark [\#73](https://github.com/PeteHaitch/methtuple/issues/73)

- Allow `bismarkify.py` to write to STDOUT [\#32](https://github.com/PeteHaitch/methtuple/issues/32)

- bismarkify.py [\#30](https://github.com/PeteHaitch/methtuple/issues/30)

**Fixed bugs:**

- Set default ir1p or ir2p [\#79](https://github.com/PeteHaitch/methtuple/issues/79)

- Make sure `bismarify.py` can deal with a BAM that doesn't already have a `@PG` tag in the header [\#33](https://github.com/PeteHaitch/methtuple/issues/33)

- Add automatic updating of `@PG` tag when updates are made to `bismarkify.py` [\#29](https://github.com/PeteHaitch/methtuple/issues/29)

**Closed issues:**

- More options for run\_comethylation.sh [\#78](https://github.com/PeteHaitch/methtuple/issues/78)

- Improve memory usage [\#64](https://github.com/PeteHaitch/methtuple/issues/64)

- Notes on SAM/BAM files created by BSmooth  [\#19](https://github.com/PeteHaitch/methtuple/issues/19)

**Merged pull requests:**

- Indel support [\#89](https://github.com/PeteHaitch/methtuple/pull/89) ([PeteHaitch](https://github.com/PeteHaitch))

## [v1.2](https://github.com/PeteHaitch/methtuple/tree/v1.2) (2014-06-20)

[Full Changelog](https://github.com/PeteHaitch/methtuple/compare/v1.1...v1.2)

**Implemented enhancements:**

- Python3 support [\#75](https://github.com/PeteHaitch/methtuple/issues/75)

- Simplify command line arguments [\#74](https://github.com/PeteHaitch/methtuple/issues/74)

**Closed issues:**

- Update examples using newest version of bismark [\#65](https://github.com/PeteHaitch/methtuple/issues/65)

**Merged pull requests:**

- Add strand [\#81](https://github.com/PeteHaitch/methtuple/pull/81) ([PeteHaitch](https://github.com/PeteHaitch))

- Ignore by seq cycle [\#77](https://github.com/PeteHaitch/methtuple/pull/77) ([PeteHaitch](https://github.com/PeteHaitch))

## [v1.1](https://github.com/PeteHaitch/methtuple/tree/v1.1) (2014-06-11)

[Full Changelog](https://github.com/PeteHaitch/methtuple/compare/v1.0...v1.1)

**Implemented enhancements:**

- More general filtering by cycle [\#71](https://github.com/PeteHaitch/methtuple/issues/71)

- Add option to gzip output [\#69](https://github.com/PeteHaitch/methtuple/issues/69)

- Use GNU parallel in helper script [\#68](https://github.com/PeteHaitch/methtuple/issues/68)

## [v1.0](https://github.com/PeteHaitch/methtuple/tree/v1.0) (2014-03-20)

**Implemented enhancements:**

- Handle PBAT libraries [\#62](https://github.com/PeteHaitch/methtuple/issues/62)

- Handle non-directional protocols [\#61](https://github.com/PeteHaitch/methtuple/issues/61)

- Make a strand/orientation checking function [\#60](https://github.com/PeteHaitch/methtuple/issues/60)

- Test whether read contains complicated cigar and skip read if it does [\#57](https://github.com/PeteHaitch/methtuple/issues/57)

- Add Travis CI integration [\#56](https://github.com/PeteHaitch/methtuple/issues/56)

- Allow multiple methylation types [\#54](https://github.com/PeteHaitch/methtuple/issues/54)

- Write `QNAME` of reads that failed QC to a separate file [\#50](https://github.com/PeteHaitch/methtuple/issues/50)

- Helpful error messages [\#49](https://github.com/PeteHaitch/methtuple/issues/49)

- Include helper scripts [\#46](https://github.com/PeteHaitch/methtuple/issues/46)

- Extend methodology from 2-tuples to m-tuples \(m = 1, 2, ...\) [\#39](https://github.com/PeteHaitch/methtuple/issues/39)

- Write histogram of CpGs per read to a separate file [\#37](https://github.com/PeteHaitch/methtuple/issues/37)

- Add option to filter reads by `MAPQ` in `comethylation.py` [\#36](https://github.com/PeteHaitch/methtuple/issues/36)

- Rewrite `correct\_Bismark\_PE\_SAM.py` [\#34](https://github.com/PeteHaitch/methtuple/issues/34)

- Check that `--oldBismark` flags works correctly [\#24](https://github.com/PeteHaitch/methtuple/issues/24)

- Include `--mTuple` and `--methylationType` parameters in output filenames [\#23](https://github.com/PeteHaitch/methtuple/issues/23)

- Choose license [\#22](https://github.com/PeteHaitch/methtuple/issues/22)

- Add `--strandCollapse` or `--noStrandCollapse` option [\#18](https://github.com/PeteHaitch/methtuple/issues/18)

- Parse SAM/BAM header to extract useful information [\#14](https://github.com/PeteHaitch/methtuple/issues/14)

- Add a --bismark option to directly process SAM/BAM files created by Bismark, i.e. allowing for the non-compliant FLAG and read name suffixes imposed by Bismark. [\#11](https://github.com/PeteHaitch/methtuple/issues/11)

- Refine process for dealing with overlapping read-pairs [\#8](https://github.com/PeteHaitch/methtuple/issues/8)

- Add count of number of reads skipped due to non-identical overlapping sequence for paired-end data [\#7](https://github.com/PeteHaitch/methtuple/issues/7)

- List of key features for test data [\#4](https://github.com/PeteHaitch/methtuple/issues/4)

- Documentation [\#2](https://github.com/PeteHaitch/methtuple/issues/2)

- Overview of plan for initial public release [\#1](https://github.com/PeteHaitch/methtuple/issues/1)

**Fixed bugs:**

- python setup.py test has exit status 1 [\#63](https://github.com/PeteHaitch/methtuple/issues/63)

- Changing name from Comethylation to comethylation [\#58](https://github.com/PeteHaitch/methtuple/issues/58)

- Test against non-directional data [\#53](https://github.com/PeteHaitch/methtuple/issues/53)

- Simplify command line arguments [\#52](https://github.com/PeteHaitch/methtuple/issues/52)

- Ensure scripts fail gracefully if packages aren't available [\#51](https://github.com/PeteHaitch/methtuple/issues/51)

- Decide on names of concepts [\#48](https://github.com/PeteHaitch/methtuple/issues/48)

- Change all references to n-tuple to m-tuple [\#47](https://github.com/PeteHaitch/methtuple/issues/47)

- `chromosome\_index` is hard-coded [\#43](https://github.com/PeteHaitch/methtuple/issues/43)

- Develop consistent order of column names in output of `comethylation.py` [\#40](https://github.com/PeteHaitch/methtuple/issues/40)

- Check all scripts for compatibility with latest version of pysam [\#38](https://github.com/PeteHaitch/methtuple/issues/38)

- Gracefully fail for reads without XM-, XR- or XG-tags [\#17](https://github.com/PeteHaitch/methtuple/issues/17)

- Allow ignoring of a different number of bases for read\_1 and read\_2 in PE data [\#13](https://github.com/PeteHaitch/methtuple/issues/13)

- Read-pairs with overlapping sequence are not being properly handled [\#10](https://github.com/PeteHaitch/methtuple/issues/10)

- Refine process for dealing with overlapping read-pairs [\#8](https://github.com/PeteHaitch/methtuple/issues/8)

- Remove hardcoding from `create\_chromosome\_index\(\)` [\#6](https://github.com/PeteHaitch/methtuple/issues/6)

- Make all function variables explicitly scoped [\#5](https://github.com/PeteHaitch/methtuple/issues/5)

- Warning about constructing bookended CpG-n-tuples from paired-end reads [\#3](https://github.com/PeteHaitch/methtuple/issues/3)

**Closed issues:**

- Acknowledgements [\#66](https://github.com/PeteHaitch/methtuple/issues/66)

- Remove unused files for initial public release [\#59](https://github.com/PeteHaitch/methtuple/issues/59)

- Deprecate `--pairChoice` option [\#41](https://github.com/PeteHaitch/methtuple/issues/41)

- Make all scripts executable [\#28](https://github.com/PeteHaitch/methtuple/issues/28)

- Add proper version numbers to all scripts [\#26](https://github.com/PeteHaitch/methtuple/issues/26)

**Merged pull requests:**

- Reinstated tests folder on add\_tests branch [\#55](https://github.com/PeteHaitch/methtuple/pull/55) ([PeteHaitch](https://github.com/PeteHaitch))

- Begin adding support for reads containing INDELs to both bismarkify.py a... [\#16](https://github.com/PeteHaitch/methtuple/pull/16) ([PeteHaitch](https://github.com/PeteHaitch))

- Master [\#21](https://github.com/PeteHaitch/methtuple/pull/21) ([PeteHaitch](https://github.com/PeteHaitch))



\* *This Change Log was automatically generated by [github_changelog_generator](https://github.com/skywinder/Github-Changelog-Generator)*