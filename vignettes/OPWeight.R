## ----setup, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, autodep = TRUE)

## ----ranksProb-----------------------------------------------------------
set.seed(123)
prob_cont <- sapply(1:m, prob_rank_givenEffect, et = mean_filterEffect,
                   ey = mean_filterEffect, nrep = 10000, m0 = m0, m1 = m1)

## ----weights-------------------------------------------------------------
wgt <- weight_continuous(alpha = .1, et = mean_testEffect, m = m, ranksProb = prob_cont)

## ----results-------------------------------------------------------------
alpha = .1
padj <- p.adjust(Ordered.pvalue/wgt, method = "BH")
rejections_list = OD[which((padj <= alpha) == TRUE), ]
n_rejections = dim(rejections_list)[1]
n_rejections

## ------------------------------------------------------------------------
opw_results2 <- opw(pvalue = pvals, filter = filters, weight = wgt, 
                     effectType = "continuous", alpha = .1, method = "BH")
opw_results2$rejections

