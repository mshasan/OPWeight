## ----setup, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, autodep = TRUE)

## ----loadOPWeight, message=FALSE, warning=FALSE--------------------------
library("OPWeight")
opw_results <- opw(pvalue = de_res_air$pvalue, filter = de_res_air$baseMean, 
                  ranks = TRUE, test = de_res_air$stat, alpha = .05, tail = 2, 
                  effectType = "continuous", method = "BH")

## ----outputsName---------------------------------------------------------
names(opw_results)

## ----nullProp------------------------------------------------------------
opw_results$propNulls

## ----no.OfRejections-----------------------------------------------------
opw_results$rejections

## ----probVsWeight--------------------------------------------------------
m = opw_results$totalTests
ranks = 1:m
probs = opw_results$probGivenEffect
weights = opw_results$weight
dat = data.frame(ranks, probs, weights)
p_prob = ggplot(dat, aes(x = ranks, y = probs)) + geom_point(size = .5, col="firebrick4") + 
    labs(y = "p(rank | effect)")
p_weight = ggplot(dat, aes(x = ranks, y = weights)) + geom_point(size = .5, col="firebrick4")
grid.arrange(p_prob, p_weight, ncol = 2)

