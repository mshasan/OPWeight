[![Build Status](https://travis-ci.org/mshasan/OPWeight.svg?branch=master)](https://travis-ci.org/mshasan/OPWeight)
[![codecov.io](https://codecov.io/github/mshasan/OPWeight/coverage.svg?branch=master)](https://codecov.io/github/mshasan/OPWeight?branch=master)

# OPWeight
An R package to perform optimal p-value weighting with indpendent information

The large scale multiple testing inherent to high throughput biological data necessitates very high statistical stringency and thus true effects in data are difficult to detect unless they have high effect sizes. One solution to this problem is to use independent information to prioritize the most promising features of the data and thus increase the power to detect them. Weighted p-values provide a general framework for doing this in a statistically rigorous fashion.  However, calculating weights that incorporate the independent information and optimize statistical power remains a challenging problem despite recent advances in this area. Existing methods tend to perform poorly in the common situation that true positive features are rare and/or of low effect size. We introduce a method for calculating approximate optimal weights conditioned on the ranking of tests by an external covariate. This approach uses the probabilistic relationship between covariate ranking and test effect size to calculate more informative weights that are not diluted by null effects as is common with group-based methods. This relationship can be calculated theoretically for normally distributed covariates or estimated empirically in other cases. We show via simulations and applications to data that this method outperforms existing methods by a large margin in the rare/low effect size scenario and has at least comparable performance in all scenarios. Use of this approach can greatly enhance discovery probabilities for high throughput data.


You can install the package as follows:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("OPWeight")
```


