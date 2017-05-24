## ----setup, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, autodep = TRUE)

## ----dataAnalysis--------------------------------------------------------
# initial stage--------
pvals = de_res_air$pvalue
tests = qnorm(pvals/2, lower.tail = FALSE)
filters = de_res_air$baseMean + .0001    # to ensure filters are postive for the box-cox

# formulate a data set-------------
Data = tibble(pvals, filters)
OD <- Data[order(Data$filters, decreasing=T), ]
Ordered.pvalue <- OD$pvals


# estimate the true null and alternative test sizes------
m = length(Data$pvals); m
nullProp = qvalue(Data$pvals, pi0.method="bootstrap")$pi0; nullProp
m0 = ceiling(nullProp*m); m0
m1 = m - m0; m1

# fit box-cox regression
#--------------------------------

bc <- boxcox(filters ~ tests)
lambda <- bc$x[which.max(bc$y)]; lambda
model <- lm(filters^lambda ~ tests)

# If lambda = 0. use log-transformation
# model <- lm(log(Data$filters) ~ Data$tests)

# etimated test efects of the true altrnatives------------
test_effect = if(m1 == 0) {0
} else {sort(tests, decreasing = T)[1:m1]}		# two-tailed test

# for the continuous effects etimated mean effects
mean_testEffect = mean(test_effect, na.rm = T)
mean_testEffect
mean_filterEffect = model$coef[[1]] + model$coef[[2]]*mean_testEffect
mean_filterEffect

# # for the binary effects estiamted median effects 
# mean_testEffect = median(test_effect, na.rm = T)
# mean_testEffect
# mean_filterEffect = model$coef[[1]] + model$coef[[2]]*mean_testEffect
# mean_filterEffect

