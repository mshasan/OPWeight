# OPWeight
An R package to perform Optimal Pvalue Weighting for High-Throughput Data

1) A novel pvalue weighting method for high throughput data. 
2) The uniqueness of our method is that it does not require estimating true effect sizes.
3) This method is able to produce a distinct weight for individual tests, and hence the true effect can stand out from the false effect. 
4) The method is based on a probabilistic relationship between the true effect size and the ranking of tests by an 
independent covariate. 
5) Perform better when the number of true alternate tests are only a small fraction of all tests and have low effect sizes. 
6) The proposed method shows significant improvement in power compared to other methods while maintaining 
the family wise error rate as well as false discovery rate. 
7) The Proposed method is not sensitive to correlations between the test 
statistics and does not significantly affect by the variation of the effect sizes. 


You can install the package as follows:

```{r}
library("devtools")

# Bioconductor prerequisites
source("http://bioconductor.org/biocLite.R")
biocLite(c("genefilter","DESeq2","qvalue","Biobase",
            "BiocParallel","airway","pasilla", "BiocStyle"))
    
# install OPWeight
install_github("vladchimescu/lpsymphony", subdir="lpsymphony")
install_github("mshasan/OPWeight")
```
