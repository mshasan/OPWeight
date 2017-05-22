#' @title Perform Optimal Pvalue Weighting
#'
#' @description A function to perform weighted pvalue multiple hypothesis test.
#' This function compute either the probabilities of the tests' ranks or the
#' probabilities of the filters given the effect sizes, and consequently
#' the weights, then provides the number of rejected null hypothesis and the list of
#' the rejected pvalues as well as the corresponing filter statistics.
#'
#' @param  pvalue vector of pvalues of the test statistics
#' @param filter vector of filter statistics
#' @param prob_givenEffect the probabilities of the filters or filters' ranks given
#' the mean of the filter effects
#' @param ranks determine what type of probabilities \code{prob_givenEffect} is used
#' @param mean_filterEffect mean filter effect of the true alternatives
#' @param mean_testEffect mean test effect of the true alterantives
#' @param effectType type of effect sizes; c("continuous", "binary")
#' @param alpha significance level of the hypothesis test
#' @param nrep number of replications for importance sampling, default value is 10,000,
#' can be increased to obtain smoother probability curves
#' @param tail right-tailed or two-tailed hypothesis test. default is right-tailed test.
#' @param delInterval interval between the \code{delta} values of a sequence. Note that,
#' \code{delta} is a LaGrange multiplier, necessary to normalize the weight
#' @param method type of methods is used to obtain the results; c("BH", "BON"),
#' Benjemini-Hochberg or Bonferroni
#' @param monitor to observe the progress of the computation
#' @param ... Arguments passed to internal functions
#'
#' @details If one wants to test \deqn{H_0: epsilon_i = 0 vs. H_a: epsilon_i > 0,}
#' then the \code{mean_testEffect}  and \code{mean_filterEffect} should be mean of the test
#' and filter effect sizes, respectively. This is called hypothesis testing for
#' the continuous effect sizes.\cr
#'
#' If one wants to test \deqn{H_0: epsilon_i = 0 vs. H_a: epsilon_i = epsilon,}
#' then \code{mean_testEffect} and \code{mean_filterEffect} should be median or
#' any discrete value of the test and filter effect sizes. This is called hypothesis
#' testing for the Binary effect sizes, where \code{epsilon} refers to a fixed value.\cr
#'
#' Internally, \code{opw} function compute the \code{prob_givenEffect} and consequently
#' the weights, then uses the pvalues to make conclusions about hypotheses.
#' Therefore, if \code{prob_givenEffect} is given then \code{mean_filterEffect}
#' and are redundant, and should not be provided to the funciton.
#' Although \code{prob_givenEffect} is not required to the function,
#' One can compute \code{prob_givenEffect} by using either the function
#' \code{\link{prob_rank_givenEffect}} if \code{ranks == TRUE} or
#' \code{\link{dnorm}} if \code{ranks == FALSE}.\cr
#'
#' The function internally compute \code{mean_filterEffect} and \code{mean_testEffect}
#' from a simple linear regression with box-cox transformation between the test
#' and filter statistics, where the filters are regressed on the test statistics.
#' Thus, filters need to be positive to apply \code{boxcox} from the \code{R}
#' library \code{MASS}. Then the estimated \code{mean_filterEffect} and
#' \code{mean_testEffect} are used to obtian the \code{prob_givenEffect} and the weights.
#' Thus, in order to apply the function properly, it is crucial to understand the
#' uses \code{mean_filterEffect} and \code{mean_testEffect}. If \code{mean_filterEffect} and
#' \code{mean_testEffect} are not provided then the test statistics computed from
#' the pvalues will be used to compute the relationship between the filter
#' statistics and the test statistics.\cr
#'
#' If one of the mean effects \code{mean_filterEffect} and \code{mean_testEffect}
#' are not provided then the missing mean effect will be computed internally.
#'
#' If \code{monitor = TRUE} then a window will open to see the progress of the
#' computation. It is useful for a large number of tests
#'
#' @author Mohamad S. Hasan and Paul Schliekelman
#'
#' @export
#'
#' @import qvalue qvalue
#' @import tibble tibble
#' @import MASS boxcox
#'
#' @seealso \code{\link{prob_rank_givenEffect}} \code{\link{weight_binary}}
#' \code{\link{weight_continuous}} \code{\link{qvalue}} \code{\link{dnorm}}
#'
#'
#' @return \code{totalTests} total number of hypothesis tests evaluated
#' @return \code{propNulls} estimated propotion of the true null hypothesis
#' @return \code{probGivenEffect} probability of the ranks given the mean filter effect,
#' p(rank | ey = mean_filterEffect)
#' @return \code{weight} normalized weight
#' @return \code{rejections} total number of rejections
#' @return \code{rejections_list} list of rejected pvalues and the corresponding
#' filter statistics
#'
#'
#' @examples
#' # generate pvalues and filter statistics
#' m = 1000
#' set.seed(3)
#' filters = runif(m, min = 0, max = 2.5)          # filter statistics
#' H = rbinom(m, size = 1, prob = 0.1)             # hypothesis true or false
#' tests = rnorm(m, mean = H * filters)            # Z-score
#' pvals = 1 - pnorm(tests)                        # pvalue
#'
#' # ranks = FLASE to use probability from the normal density
#' results <- opw(pvalue = pvals, filter = filters, ranks = FALSE,
#'                      effectType = "continuous", method = "BH")
#'
#' # ranks = TRUE to use ranks probability
#' results2 <- opw(pvalue = pvals, filter = filters, ranks = TRUE,
#'                effectType = "continuous", tail = 2, method = "BH", monitor = TRUE)
#'
#' # supply the mean effects for both the filters and the tests externally as follows:
#' mod <- lm(log(filters) ~ tests)
#' et = mean(tests)
#' ey = mod$coef[[1]] + mod$coef[[2]]*et
#' results3 <- opw(pvalue = pvals, filter = filters, ranks = FALSE,
#'                mean_filterEffect = ey, mean_testEffect = et, tail = 2,
#'                effectType = "continuous", method = "BH")
#'
#' # supply the probabilities externally as follows:
#' library(qvalue)
#' ranks <- 1:m
#' nullProp = qvalue(p = pvals, pi0.method = "bootstrap")$pi0
#' m0 = ceiling(nullProp*m)
#' m1 = m - m0
#' probs <- sapply(ranks, prob_rank_givenEffect, et = ey, ey = ey,
#'                                         nrep = 10000, m0 = m0, m1 = m1)
#' results4 <- opw(pvalue = pvals, filter = filters, prob_givenEffect = probs,
#'                      effectType = "continuous", tail = 2, method = "BH")
#'
#' probs2 <- dnorm(filters, mean = 0, sd = 1)
#' results5 <- opw(pvalue = pvals, filter = filters, prob_givenEffect = probs2,
#'                      effectType = "continuous", tail = 2, method = "BH")
#'
#============================================================================
# function to apply opw methods on data
#---------------------------------------------------
# Input:
#----------------------------
# pvalue = vector of pvalues
# filter = vector of filter statistics
# prob_givenEffect = the probabilities of the filters or filters' ranks given
# the mean of the filter effects
# ranks = determine what type of probabilities \code{prob_givenEffect} is used
# mean_filterEffect = filter effect size
# mean_testEffect = test effect size
# effectType = type of effect size c("binary","continuous")
# alpha = significance level of the hypotheis test
# nrep = number of replications for importance sampling
# tail = right-tailed or two-tailed hypothesis test. default is two-tailed test
# delInterval = interval between the delta values of a sequence
# method = Benjamini_HOchberg (BH) or Bonferroni (BON)

# internal parameters:-----
# m = number of hypothesis test
# nullProp = proportion of true null hypothesis
# m0 =  number of the true null tests
# m1 = number of the true alternative tests
# test =  compute test statistics from the pvalues if not given
# test_effect_vec = estiamted number of the true alternaitve test statistics
# mean_testEffect = mean test effect sizes of the true alternaive hypotheis
# mean_filterEffect = mean filter effect sizes of the true alternaive hypotheis
# prob = probailities of the ranks given the mean effect size
# wgt = weights
# Data = create a data set
# OD = odered by covariate
# odered.pvalues = odered pvalues for all tests
# padj = adjusted pvalues for FDR uses

# Output:
#-------------------------
# totalTests = total number of hypothesis tests
# propNulls = estimated propotion of the true null hypothesis
# probGivenEffect = probability of the ranks given the mean filter effect,
#                                           p(rank | ey = mean_filterEffet)
# weight = normalized weight
# rejections = total number of rejections
# rejections_list = list of rejected pvalues and corresponding filter statistics
#-------------------------------------------------------------------------------

opw <- function(pvalue, filter, prob_givenEffect = NULL, ranks = FALSE,
                mean_filterEffect = NULL, mean_testEffect = NULL,
                effectType = c("continuous", "binary"), alpha = .05, nrep = 10000,
                tail = 1L, delInterval = .0001, method = c("BH", "BON"), monitor = FALSE, ... )
    {
        # compute the number of tests------------
        m = length(pvalue)
        nullProp = qvalue(p = pvalue, pi0.method = "bootstrap")$pi0
        m0 = ceiling(nullProp*m)
        m1 = m - m0

        # determine the side of the tests-------------
        if(any(filter <= 0)){
             stop("filter statistics need to be positive")
        }

        # compute test statistics from the pvalues---------
        if(tail == 1){
            test <- qnorm(pvalue, lower.tail = FALSE)
        } else {
            test <- qnorm(pvalue/2, lower.tail = FALSE)
        }

        test[which(!is.finite(test))] <- NA

        # estimate the true alterantive test effect sizes----------------
        if(m1 == 0){
            test_effect_vec <- 0
        } else {
            test_effect_vec <-  sort(test, decreasing = TRUE)[1:m1]
        }

        # estimate the mean test effect size-------------
        if(!is.null(mean_testEffect)){
            mean_testEffect <- mean_testEffect
        } else {
            if(effectType == "continuous"){
                mean_testEffect <- mean(test_effect_vec, na.rm = TRUE)
            } else {
                mean_testEffect <- median(test_effect_vec, na.rm = TRUE)
            }
        }

        # estimate lambda from the box-cox transformation----------------
        bc <- boxcox(filter ~ test)
        lambda <- bc$x[which.max(bc$y)]

        # estimate the mean filter effect size------------------
        if(!is.null(mean_filterEffect)){
            mean_filterEffect <- mean_filterEffect
        } else {
            if(lambda == 0){
                model <- lm(log(filter + .0001) ~ test)
            } else {
                model <- lm(filter**lambda ~ test)
            }
            mean_filterEffect <- model$coef[[1]] + model$coef[[2]]*mean_testEffect
        }

        # compute the probability of the filter given the mean filter effect
        message("computing probabilities")
        if(!is.null(prob_givenEffect)){
            prob <- prob_givenEffect
        } else {
            if(ranks == FALSE){
                if(lambda == 0){
                    prob <- dnorm(log(filter + .0001), mean = 0, sd = 1)
                } else {
                    prob <- dnorm(filter**lambda, mean = 0, sd = 1)
                }
            } else {
                prob <- sapply(1:m, prob_rank_givenEffect, et = mean_filterEffect,
                            ey = mean_filterEffect, nrep = nrep, m0 = m0, m1 = m1,
                            monitor = monitor)
            }
        }
        message("finished computing the probabilities")

        # compute the weights (always right-tailed)------------
        message("computing weights")
        if(effectType == "continuous"){
            wgt = weight_continuous(alpha = alpha, et = mean_testEffect, m = m,
                            tail = 1, delInterval = delInterval, prob = prob)
        } else {
            wgt = weight_binary(alpha = alpha, et = mean_testEffect, m = m, m1 = m1,
                                tail = 1, delInterval = delInterval, prob = prob)
        }
        message("finished computing the weights")

        # formulate a data set-------------
        Data = tibble(pvalue, filter)
        OD <- Data[order(Data$filter, decreasing=T), ]
        Ordered.pvalue <- OD$pvalue

        if(method == "BH"){
            padj <- p.adjust(Ordered.pvalue/wgt, method = "BH")
            rejections_list = OD[which((padj <= alpha) == TRUE), ]
        } else {
            rejections_list = OD[which((Ordered.pvalue <= alpha*wgt/m) == TRUE), ]
        }

        # outputs--------------
        n_rejections = dim(rejections_list)[1]

        return(list(totalTests = m, propNulls = nullProp,
                    probGivenEffect = prob, weight = wgt,
                    rejections = n_rejections, rejections_list = rejections_list))
    }





