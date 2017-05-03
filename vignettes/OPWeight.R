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
#data("airway", package = "airway")
#dds <- DESeqDataSet(se = airway, design = ~ cell + dex)
#dds <- DESeq(dds)
#de_res_air <- as.data.frame(results(dds))
de_res_air <- data.frame(pvalue=runif(100), stat = rnorm(100), baseMean = sample(1:100))

## ----colnames_de_res_air-------------------------------------------------
colnames(de_res_air)
dim(de_res_air)

## ----loadihw, message=FALSE, warning=FALSE-------------------------------
library("OPWeight")
opw_results <- opw(pvalue = de_res_air$pvalue, filter = de_res_air $baseMean, 
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
p_prob = ggplot(dat, aes(x = ranks, y = probs)) + geom_line(size=1.5, col="firebrick4") + 
    labs(y = "p(rank | effect)")
p_weight = ggplot(dat, aes(x = ranks, y = weights)) + geom_line(size=1.5, col="firebrick4")
grid.arrange(p_prob, p_weight, ncol = 2)

## ------------------------------------------------------------------------
Data <- tibble(test=de_res_air$stat, pval=de_res_air$pvalue, filter=de_res_air$baseMean)

barlines <- "#1F3552"

hist_test <- ggplot(Data, aes(x = Data$test)) +
        geom_histogram(aes(y = ..density..), binwidth = 1,
	  colour = barlines, fill = "#4271AE") +
		labs(x = "test statistics")

hist_pval <- ggplot(Data, aes(x = Data$pval)) +
        geom_histogram(aes(y = ..density..),
	  colour = barlines, fill = "#4281AE")+
		labs(x = "pvalues")

hist_filter <- ggplot(Data, aes(x = Data$filter)) +
        geom_histogram(aes(y = ..density..),
	  colour = barlines, fill = "#4274AE") +
		labs(x = "filter statistics")

test_filter <- ggplot(Data, aes(x = Data$test, y = Data$filter)) +
		geom_point() + labs(x = "test statstics", y = "filter statistics")
		#scale_x_continuous(limits = c(0, 25000), breaks=seq(0, 25000, 10000))

pval_filter <- ggplot(Data, aes(x = rank(-Data$filter), y = -log10(pval))) +
		geom_point()+
		labs(x = "ranks of filters", y = "-log(pvalue)")
		#scale_x_continuous(limits = c(0, 25000), breaks=seq(0, 25000, 10000))

p_ecdf <- ggplot(Data, aes(x = pval)) +
			stat_ecdf(geom = "step")+
			labs(x = "pvalues", title="empirical cumulative distribution")+
			theme(plot.title = element_text(size = rel(.7)))


grid.arrange(hist_test,  hist_pval, hist_filter, test_filter , pval_filter, p_ecdf, 
		ncol = 3, heights = c(7, 7), top = "Airway: data summary")


## ------------------------------------------------------------------------
# fite box-cox regression
#--------------------------------
m = length(Data$pval)
null = qvalue(Data$pval, pi0.method="bootstrap")$pi0
m0 = ceiling(null*m)
m1 = m - m0

bc <- boxcox(Data$filter ~ Data$test, plotit = FALSE)
trans <- bc$x[which.max(bc$y)]
model <- lm(Data$filter^trans ~ Data$test)

test_effect = if(m1 == 0) {0
} else {sort(abs(Data$test), decreasing = T)[1:m1]}		# two-tailed test

# for the continuous effects 
mean_testEffect = mean(test_effect, na.rm = T)
mean_filterEffect = model$coef[[1]] + model$coef[[2]]*mean_testEffect

# for the binary effects 
mean_testEffect = median(test_effect, na.rm = T)
mean_filterEffect = model$coef[[1]] + model$coef[[2]]*mean_testEffect

