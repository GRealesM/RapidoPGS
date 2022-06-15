# RapidoPGS 2.2.0

RapidoPGS 2.2.0.9000
2022-06-15
- We updated hg38 LD blocks from liftovered Berisa & Pickrell to recomputed MacDonald et al. 2022 (https://github.com/jmacdon/LDblocks_GRCh38, https://www.biorxiv.org/content/10.1101/2022.03.04.483057v1), which offer a more accurate representation of LD structure in hg38.

RapidoPGS 2.1.0.9009
2021-09-28
- Fixed a bug derived from automatically supplying "nref" to SuSIE in `rapidopgs_multi()`, which is no longer required.

RapidoPGS 2.1.0.9008
2021-08-13
- Fixed a change in the hard-coded url to download 1000G panel from the official server to meet a change in versioning (v5a -> v5b) at the source.

RapidoPGS 2.1.0.9007
2021-07-13
- Fixed a bug caused by changes in behavior of `runsusie()` in `rapidopgs_multi()` that used to supply an extra zero element which is not supplied anymore.

RapidoPGS 2.1.0.9006 
2021-07-12
- Fixed an error popping up when `rapidopgs_multi()` is not supplied input of data.table class.
- Removed a deprecated argument in `runsusie()` internal function that was preventing `rapidopgs_multi()` to run properly.

RapidoPGS 2.1.0.9005 
2021-07-09
- Fixed a change in the hard-coded url to download 1000G panel from the official server to meet a change in versioning for sex chromosomes at the source.

RapidoPGS 2.1.0.9004
2021-06-17
- Tiny updates in README and vignettes concerning dependencies and versions.

RapidoPGS 2.1.0.9003
2021-06-09
- Removed need for N_LD in `rapidopgs_multi()`.


RapidoPGS 2.1.0.9002
2021-06-08
- Updated dependency versions.
- Change some stuff for CRAN submission.

RapidoPGS 2.1.0.9001
2021-04-14
- New vignette for `rapidopgs_multi()` and minor changes in `rapidopgs_single()` vignette.
- New dataset `michailidou19`, analogous to `michailidou` but in hg19 build, required for `rapidopgs_multi()` examples.


2021-03-08
- New implementation of `rapidopgs_multi()`, which now allows to use pre-computed LD matrices.
- `rapidopgs_multi()` now uses a different susie implementation, as implemented in `coloc` package.
- `rapidopgs_multi()` added parameters to select trait type, as well as path to LD matrices, and sample size for the dataset used for LD matrix computation.
- `coloc` added as a dependency.


# RapidoPGS 2.0.0 

RapidoPGS 2.0.0.9001 
2021-01-27
- Fixed bug in `rapidopgs_single()` that prevented it for running for quantitative traits.
- Fixed missing space in message for users in `rapidopgs_single()`.

2020-12-17
- Changed function name `computePGS()` to `rapidopgs_single()`.
- `rapidopgs_single()` dropped the requirement for MAF for case-control traits.
- `rapidopgs-single()` added a parameter to select type of trait (case-control or quantitative) of the dataset used in the call.
- `rapidopgs_multi()` function added, which applies SuSiE method to PGS computing.
- `susieR` added as a dependency.
- `create_1000G()` added, which automatically download and sets up a reference panel for `rapidopgs_multi()` from 1000 Genomes Phase III.
- `sd.prior.est()` function added to use heritability estimates to create an informed SD prior.

# RapidoPGS 1.0.2
2020-08-03
- Corrected description paragraph to not mention the package itself. Now it starts as "Quickly computes...".
- Included package reference (authors (year), <doi:...>) in DESCRIPTION.
- Olly Burren added as a contributor
- Included reference to preprint in README.
- Included executable examples in the exported functions.
- Updated the website to include links to GitHub on the side.

# RapidoPGS 1.0.1
2020-07-24
- First package version and attempt to submission to CRAN.
  
