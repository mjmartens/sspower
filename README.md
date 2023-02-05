
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sspower

<!-- badges: start -->
<!-- badges: end -->

The goal of sspower is to provide functions to compute sample size and
power for studies that will be analyzed using generalized linear models
and from time-to-event regression models, including the Cox proportional
hazards and Fine-Gray models.

## Installation

You can install the released version of sspower from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("sspower")
```

## Example

Suppose we want to estimate the sample size needed for a joint test of
the effects of two main variable in a logistic regression model while
adjusting for other covariates. The test will have a 5% significance
level and 80% power. Log odds ratios targeted for the variables are 0.5
and 0.3. The main variables are assumed to be independent with variances
0.5 and 0.25. The coefficient of determination matrix (R2) is assumed to
be diag(c(0.,0.5)). The expected residual variance E\[Var(Y\|X)\] is
specified as 0.2.

Sample size calculation proceeds as follows from the ssGLM function:

``` r
library(sspower)
level = 0.05
pow = 0.80
del = c(0.5,0.3)
vz = diag(c(0.5,0.25))
rsq = diag(c(0,0.5))
f0 = 0.20

ssGLM(level,pow,del,vz,rsq,f0)
#> [1] 353.5656
```

Approximately 354 patients in total are required for this study.
