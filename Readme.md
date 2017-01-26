[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.260218.svg)](https://doi.org/10.5281/zenodo.260218)

This repository contains all of the raw data and analysis scripts to accompany the manuscript:

# Tenacious D: *Symbiodinium* in clade D remain in reef corals at both high and low temperature extremes despite impairment
### Authors: Rachel N. Silverstein, Ross Cunning, and Andrew C. Baker
### Journal: _Journal of Experimental Biology_
### Link: [doi:10.1242/jeb.148239](http://dx.doi.org/10.1242/jeb.148239)

-----

### Description:
With this work, we challenge the paradigm that temperature-induced photodamage in coral symbionts leads directly to coral ‘bleaching’ (loss of symbionts), and that clade D symbionts are able to resist bleaching by resisting temperature stress. In fact, we demonstrate that clade D experiences equivalent or greater photodamage relative to clade C in response to experimental heating and cooling, yet remains at higher abundance within the coral host. Therefore, the bleaching-resistance of clade D may actually be due to its tenacity in intracellular occupation rather than its sensitivity to temperature stress.

### Repository contents:
#### Data:
**tenaciousD_data.csv:** CSV file containing all raw data. Column headings are as follows:

- sample: coral core unique identifier
- mother: mother colony from which core was taken
- history: thermal history of core in previous study [(Silverstein et al. 2015)](http://dx.doi.org/10.1111/gcb.12706)
    - see table below for description of codes used to indicate thermal history
- ramp: whether heating or cooling treatment was applied
- time: day on which datum or sample was collected
- fvfm: maximum quantum yield of photosystem II
- C.SH: symbiont to host ratio of clade C *Symbiodinium*
- D.SH: symbiont to host ratio of clade D *Symbiodinium*
- tot.SH: total symbiont to host ratio (sum of C and D ratios)

#### Code:
**setup.R:** R script that imports data and performs initial data carpentry to prepare for analyses.

**tenaciousD_analysis.Rmd:** R Markdown file containing code to perform all analyses and generate all figures used in the study, along with detailed description and rationale for statistical approach

#### Output:
[**tenaciousD_analysis.html:**](https://cdn.rawgit.com/jrcunning/tenaciousD/b1fc6c17/tenaciousD_analysis.html) output of R Markdown file

-----

**Thermal history:** 'bleach1' and 'bleach2' indicate whether a coral was bleached by heat or DCMU or neither (control) during the first and second experimental bleaching events in Silverstein et al. 2015. 'recovtemp1' and 'recovtemp2' indicate the temperature (°C) that the corals were held at after application of the first and second experimental bleaching events.

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

