# RapidoPGS 2.0.0 
2020-12-17
- Changed function name `computePGS()` to `rapidopgs_single()`.
- `rapidopgs_single()` dropped the requirement for MAF for case-control traits.
- `rapidopgs-single()` added a parameter to select type of trait (case-control or quantitative) of the dataset used in the call.
- `rapidopgs_multi()` function added, which applies SuSiE method to PGS computing.
- `susieR` as a dependency.
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
  