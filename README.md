## betaclust
The goal of betaclust is to appositely model the beta-valued cytosine-guanine dinucleotide (CpG) sites, to objectively identify methylation state thresholds and to identify the differentially methylated CpG (DMC) sites using a model-based clustering approach. The family of BMMs employs different parameter constraints applicable to different study settings. The EM algorithm is used for parameter estimation, with a novel approximation during the M-step providing tractability and ensuring computational feasibility.

## Installation

You can install the development version of betaclust like so:

``` r
library(devtools)
install_github('koyelucd/betaclust',force = TRUE)
library(betaclust)
```
