## ----setup, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, autodep = TRUE)

## ----loadlibs, message=FALSE, warning=FALSE------------------------------
library("ggplot2")
library("airway")
library("DESeq2")
library(gridExtra)
library(cowplot)
library(tibble)
library(qvalue)
library(MASS)

## ----dataProcessing, message=FALSE, warning=FALSE------------------------
data("airway", package = "airway")
dds <- DESeqDataSet(se = airway, design = ~ cell + dex)
dds <- DESeq(dds)
de_res_air <- as.data.frame(results(dds))

## ----colnames_de_res_air-------------------------------------------------
colnames(de_res_air)
dim(de_res_air)

## ----loadOPWeight, message=FALSE, warning=FALSE--------------------------
library("OPWeight")
filters = de_res_air$baseMean + .0001  # add a small constant to make all values positive
opw_results <- opw(pvalue = de_res_air$pvalue, filter = filters, 
                  ranks = TRUE, test = de_res_air$stat, alpha = .1, tail = 2, 
                  effectType = "continuous", method = "BH")

## ----outputsName---------------------------------------------------------
names(opw_results)

## ----nullProp------------------------------------------------------------
opw_results$propNulls

## ----no.OfRejections-----------------------------------------------------
opw_results$rejections

## ----probVsWeight--------------------------------------------------------
m = opw_results$totalTests
testRanks = 1:m
probs = opw_results$probGivenEffect
weights = opw_results$weight
dat = data.frame(testRanks, probs, weights)
p_prob = ggplot(dat, aes(x = testRanks, y = probs)) + geom_point(size = .5, col="firebrick4") + 
    labs(y = "p(rank | effect)")
p_weight = ggplot(dat, aes(x = testRanks, y = weights)) + geom_point(size = .5, col="firebrick4")
plot_grid(p_prob, p_weight, labels = c("A", "B"))

## ----sumPlots------------------------------------------------------------
Data <- tibble(pval=de_res_air$pvalue, filter=de_res_air$baseMean)

barlines <- "#1F3552"

hist_pval <- ggplot(Data, aes(x = Data$pval)) +
        geom_histogram(colour = barlines, fill = "#4281AE")+
		labs(x = "pvalues")

hist_filter <- ggplot(Data, aes(x = Data$filter)) +
        geom_histogram( colour = barlines, fill = "#4274AE") +
		labs(x = "filter statistics")
        
pval_filter <- ggplot(Data, aes(x = rank(-Data$filter), y = -log10(pval))) +
		geom_point()+
		labs(x = "ranks of filters", y = "-log(pvalue)")
		
p_ecdf <- ggplot(Data, aes(x = pval)) +
			stat_ecdf(geom = "step")+
			labs(x = "pvalues", title="empirical cumulative distribution")+
			theme(plot.title = element_text(size = rel(.7)))

p <- plot_grid(hist_pval, hist_filter, pval_filter, p_ecdf, 
               labels = c("A", "B", "c", "D"), ncol = 2)
# now add the title
title <- ggdraw() + draw_label("Airway: data summary")
plot_grid(title, p, ncol = 1, rel_heights=c(0.1, 1)) 


## ----dataAnalysis--------------------------------------------------------
# initial stage--------
pvals = de_res_air$pvalue
tests = qnorm(pvals/2, lower.tail = FALSE)
filters = de_res_air$baseMean + .0001
Data <- tibble(pvals, tests, filters)


# estimate the true null and alternative test sizes------
m = length(Data$pvals); m
null = qvalue(Data$pvals, pi0.method="bootstrap")$pi0; null
m0 = ceiling(null*m); m0
m1 = m - m0; m1

# fit box-cox regression
#--------------------------------

bc <- boxcox(Data$filters ~ Data$tests)
lambda <- bc$x[which.max(bc$y)]; lambda
model <- lm(Data$filters^lambda ~ Data$tests)

# If lambda = 0. use log-transformation
# model <- lm(log(Data$filters) ~ Data$tests)

# etimated test efects of the true altrnatives------------
test_effect = if(m1 == 0) {0
} else {sort(Data$tests, decreasing = T)[1:m1]}		# two-tailed test

# for the continuous effects etimated mean effects
mean_testEffect = mean(test_effect, na.rm = T)
mean_testEffect
mean_filterEffect = model$coef[[1]] + model$coef[[2]]*mean_testEffect
mean_filterEffect

# for the binary effects estiamted median effects 
mean_testEffect = median(test_effect, na.rm = T)
mean_testEffect
mean_filterEffect = model$coef[[1]] + model$coef[[2]]*mean_testEffect
mean_filterEffect

