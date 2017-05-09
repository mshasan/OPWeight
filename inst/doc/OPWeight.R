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

data("airway", package = "airway")
dds <- DESeqDataSet(se = airway, design = ~ cell + dex)
dds <- DESeq(dds)
de_res_air <- as.data.frame(results(dds))

## ----colnames_de_res_air-------------------------------------------------
colnames(de_res_air)
dim(de_res_air)

## ----loadihw, message=FALSE, warning=FALSE-------------------------------
library("OPWeight")
opw_results <- opw(pvalue = de_res_air$pvalue, filter = de_res_air $baseMean, 
                  test = de_res_air$stat, alpha = .05, tail = 2, 
                  effectType = "continuous", method = "BH")

## ------------------------------------------------------------------------
names(opw_results)

## ------------------------------------------------------------------------
opw_results$propNulls

## ------------------------------------------------------------------------
opw_results$rejections

## ------------------------------------------------------------------------
m = opw_results$totalTests
ranks = 1:m
probs = opw_results$ranksProb
weights = opw_results$weight
dat = data.frame(ranks, probs, weights)
p_prob = ggplot(dat, aes(x = ranks, y = probs)) + geom_point(size = .5, col="firebrick4") + 
    labs(y = "p(rank | effect)")
p_weight = ggplot(dat, aes(x = ranks, y = weights)) + geom_point(size = .5, col="firebrick4")
grid.arrange(p_prob, p_weight, ncol = 2)

## ------------------------------------------------------------------------
Data <- tibble(test=de_res_air$stat, pval=de_res_air$pvalue, filter=de_res_air$baseMean)

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

grid.arrange(hist_pval, hist_filter, pval_filter, p_ecdf, 
		ncol = 2, heights = c(7, 7), top = "Airway: data summary")


## ------------------------------------------------------------------------
# estimate the true null and alternative test sizes------
m = length(Data$pval); m
null = qvalue(Data$pval, pi0.method="bootstrap")$pi0; null
m0 = ceiling(null*m); m0
m1 = m - m0; m1

# fit box-cox regression
#--------------------------------
bc <- boxcox(Data$filter ~ Data$test, plotit = FALSE)
lambda <- bc$x[which.max(bc$y)]; lambda
model <- lm(Data$filter^lambda ~ Data$test)

# lambda = 0. Therefore, we use log-transformation

model <- lm(log(Data$filter) ~ Data$test)

# etimated test efects of the true altrnatives------------
test_effect = if(m1 == 0) {0
} else {sort(abs(Data$test), decreasing = T)[1:m1]}		# two-tailed test

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

