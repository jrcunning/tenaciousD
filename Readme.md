This repository contains all of the raw data and analysis scripts to accompany the manuscript (in prep):

# Tenacious D: bleaching-resistant Symbiodinium in clade D are retained in reef corals at both high and low temperature extremes despite photophysiological stress

##by R Silverstein, R Cunning, and AC Baker

### Repository contents:

**tenaciousD_analysis.Rmd:** R Markdown file containing code to perform all analyses and generate all figures used in the study, along with detailed description and rationale for statistical approach

[**tenaciousD_analysis.html:**](tenaciousD_analysis.html) output of R Markdown file

**setup.R:** R script that imports data and performs initial data carpentry to prepare for analyses.

**supp_analysis.R:** Supplementary R script that tests the relative impact of thermal history vs. dominant symbiont clade on coral responses.

**tenaciousD.Rproj:** R project file used by RStudio to manage repository.

**tenaciousD_data.csv:** CSV file containing all raw data. Column headings are as follows:

- sample: coral core unique identifier
- mother: mother colony from which core was taken
- history: thermal history of core in previous study [(Silverstein et al. 2015)](http://dx.doi.org/10.1111/gcb.12706)

| history | bleach1 | recovtemp1 | bleach2 | recovtemp2 |
|---------|---------|------------|---------|------------|
| A'      | control | 24         | control | 24         |
| A*      | control | 24         | heat    | 24         |
| B'      | control | 29         | control | 29         |
| B*      | control | 29         | heat    | 29         |
| C'      | heat    | 24         | control | 24         |
| C*      | heat    | 24         | heat    | 24         |
| D'      | heat    | 29         | control | 29         |
| D*      | heat    | 29         | heat    | 29         |
| E'      | DCMU    | 24         | control | 24         |
| E*      | DCMU    | 24         | heat    | 24         |
| F'      | DCMU    | 29         | control | 29         |
| F*      | DCMU    | 29         | heat    | 29         |

- ramp: whether heating or cooling treatment was applied
- time: day on which datum or sample was collected
- fvfm: maximum quantum yield of photosystem II
- C.SH: symbiont to host ratio of clade C *Symbiodinium*
- D.SH: symbiont to host ratio of clade D *Symbiodinium*
- tot.SH: total symbiont to host ratio (sum of C and D ratios)


