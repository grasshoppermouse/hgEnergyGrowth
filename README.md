
# hgEnergyGrowth

<!-- badges: start -->
<!-- badges: end -->

This package provides most of the functionality used for the paper:

*Menopause averted a midlife energetic crisis with help from older children and parents: A simulation study*

The paper repo is available here: <https://github.com/grasshoppermouse/menopause>

## Installation

You can install hgEnergyGrowth from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("grasshoppermouse/hgEnergyGrowth")
```

## Example

This is a basic example which shows you how to create a data frame where each row is one year 
in the life of a hunter-gatherer family, and the columns include variables such as `wife_age`,
`husband_age`, `child_age`, `wife_survival`, `wife_consumption`, `wife_production`, and so forth.

``` r
library(hgEnergyGrowth)
out <- hg_lifecourse
```

