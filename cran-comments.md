## Test environments
* local R installation, R 4.1.0
* ubuntu 16.04 (on travis-ci), R 4.1.0
* win-builder (devel)

## R CMD check results

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

## Resubmission

In my previous submission the problem was that one dependency (coloc 5.1.0) was not yet available in CRAN. It is now, so hopefully it will be alright.

When checking for R-hub, Windows Server 2008 R2 SP1, R-devel, 32/64 bit, and 	Fedora Linux, R-devel, got PREPERROR: "Bioconductor does not yet build and check packages for R version 4.2". This is a known issue (https://github.com/r-hub/rhub/issues/471).

winbuilder finished with one note, related to (possibly) mis-spelled words in the DESCRIPTION file. I checked those and can confirm they're correct.