# geoBAMr

<img src="https://raw.githubusercontent.com/markwh/mcfli-swotr/master/logos/bamr/logo.png" width=200 alt="bamr Logo"/>

## Overview
An update to the [**bamr R package**](https://github.com/markwh/bamr) (built by [**Mark Hagemann**](https://scholar.google.com/citations?user=_-XH9u4AAAAJ&hl=en&oi=ao) while at UMass) that uses more geomorphically-informed prior knowledge for discharge inversion.

From bamr: "The bamr package facilitates Bayesian AMHG + Manning discharge estimation using stream slope, width, and partial cross-section area. It includes functions to preprocess and visualize data, perform Bayesian inference using Hamiltonian Monte Carlo (via models pre-written in the Stan language), and analyze the results."

geoBAMr expands upon this project by defining prior river knowledge using river classification frameworks.  Internal to geoBAMr, geomorphic river types are assigned to rivers using stream widths, which in turn determine which priors are fed into the BAM algorithm.  geoBAMr uses the identical Bayesian model as used in bamr.


## Installation

First, you need to have installed rstan from source on your local machine. To do that, follow the directions at [**this link**](https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-source-on-Windows) verbatim. Otherwise, an error will be thrown during package installation. This only needs to be done the first time you wish to install geoBAMr.

Following that, you can install geoBAMr:

```
 #First get devtools package
if (!require("devtools")) {
  install.packages("devtools")
  library("devtools")
}

#Then install from github
devtools:: install_github("craigbrinkerhoff/geoBAMr", force=TRUE)
```
Note that if you need to reinstall this package, uninstall the current version first.  Also, check if R has created the folder '~R/win-library/3.6/00LOCK-geoBAMr'. If so, delete it and reinstall.

## Usage
The best way to get started is to follow the examples in the included vignettes, now located at the [**bamr** website](https://markwh.github.io/bamr/index.html)

geoBAMr is packaged with two approaches to defining prior river knowledge: 1) an expert geomorphology framework or 2) an unsupervised one.  The assigned priors will be different, depending on which is chosen. 'Expert' is default and accounts for both really big rivers and 'highly width variable' ones.  The unsupervised approach is purely statistical and implicitly accounts for edge cases like these, but was not explictly designed to do so.

Regardless of classification chosen, resulting river types are accessible using the bam_priors() object$river_type.  For a read on what these river types qualitvately represent, follow this link: ().

## Notes

1) The Sacramento test case in the bamr package is not included with geoBAMr.

2) geoBAMr is far more memory intensive than bamr.  If you have large time-series of stream widths and/or slopes and geoBAMr crashes, this is likely a memory issue.  One solution is to run geoBAMr on a computing cluster, the other is to use bamr!

3) If both bamr and geoBAMr are installed, make sure to explictly call functions by package as they have the same names. Otherwise, chaotic confusion will ensue!

For example:

```
geoBAMr:: bam_estimate()
bamr:: bam_estimate()
```

## Contact
For any questions regarding this package, I am reachable at cbrinkerhoff@umass.edu
