#' @title Perform Optimal Pvalue Weighting
#'
#' @description A function to perform weighted pvalue multiple hypothesis test.
#' This function compute the probabilities of the ranks of the filter statistics
#' given the effect sizes, and consequently the weights if neighter the weights
#' nor the probabilities are given. Then provides the number of rejected null
#' hypothesis and the list of the rejected pvalues as well as the corresponing
#' filter statistics.
#'
#' @param pvalue Numeric vector of pvalues of the test statistics
#' @param filter Numeric vector of filter statistics
#' @param weight An optional numeric weight vector not required
#' @param ranksProb An optional numeric vector of the ranks probability of the
#' filters given the mean effect
#' @param mean_filterEffect Numeric, value of the mean filter effect
#' of the true alternatives
#' @param mean_testEffect Numeric, value of the mean test effect
#' of the true alterantives
#' @param effectType Character ("continuous" or "binary"), type of effect sizes
#' @param alpha Numeric, significance level of the hypothesis test
#' @param nrep Integer, number of replications for importance sampling, default
#' value is 10,000, can be increased to obtain smoother probability curves
#' @param tail Integer (1 or 2), right-tailed or two-tailed hypothesis test.
#' default is right-tailed test.
#' @param delInterval Numeric, interval between the \code{delta} values of a
#' sequence. Note that, \code{delta} is a LaGrange multiplier, necessary to
#' normalize the weight
#' @param method Character ("BH" or "BON"), type of methods is used to obtain
#' the results; Benjemini-Hochberg or Bonferroni
#' @param ... Arguments passed to internal functions
#'
#' @details If one wants to test \deqn{H_0: epsilon_i = 0 vs. H_a: epsilon_i > 0,}
#' then the \code{mean_testEffect}  and \code{mean_filterEffect} should be mean
#' of the test and filter effect sizes, respectively. This is called hypothesis
#' testing for the continuous effect sizes.\cr
#'
#' If one wants to test \deqn{H_0: epsilon_i = 0 vs. H_a: epsilon_i = epsilon,}
#' then \code{mean_testEffect} and \code{mean_filterEffect} should be median or
#' any discrete value of the test and filter effect sizes. This is called hypothesis
#' testing for the Binary effect sizes, where \code{epsilon} refers to a fixed value.\cr
#'
#' The main goal of the function is to compute the probabilities of the ranks from
#' the pvalues and the filter statistics, consequently the weights. Although \code{weights}
#' \code{ranksProb} are optional, \code{opw} has the options so that one can compute
#' the probabilities and the weights externally if necessary (see examples).\cr
#'
#' Internally, \code{opw} function compute the \code{ranksProb} and consequently
#' the weights, then uses the pvalues to make conclusions about hypotheses.
#' Therefore, if \code{ranksProb} is given then \code{mean_filterEffect}
#' and are redundant, and should not be provided to the funciton.
#' Although \code{ranksProb} is not required to the function,
#' One can compute \code{ranksProb} by using the function
#' \code{\link{prob_rank_givenEffect}}.\cr
#'
#' The function internally compute \code{mean_filterEffect} and \code{mean_testEffect}
#' from a simple linear regression with box-cox transformation between the test
#' and filter statistics, where the filters are regressed on the test statistics.
#' Thus, filters need to be positive to apply \code{boxcox} from the \code{R}
#' library \code{MASS}. Then the estimated \code{mean_filterEffect} and
#' \code{mean_testEffect} are used to obtian the \code{ranksProb} and the weights.
#' Thus, in order to apply the function properly, it is crucial to understand the
#' uses \code{mean_filterEffect} and \code{mean_testEffect}. If \code{mean_filterEffect} and
#' \code{mean_testEffect} are not provided then the test statistics computed from
#' the pvalues will be used to compute the relationship between the filter
#' statistics and the test statistics.\cr
#'
#' If one of the mean effects \code{mean_filterEffect} and \code{mean_testEffect}
#' are not provided then the missing mean effect will be computed internally.
#'
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
#' @return \code{totalTests} Integer, total number of hypothesis tests evaluated
#' @return \code{nullProp} Numeric, estimated propotion of the true null
#' hypothesis
#' @return \code{ranksProb} Numeric vector of ranks probability given the
#' mean filter effect, p(rank | ey = mean_filterEffect)
#' @return \code{weight} Numeric vector of normalized weight
#' @return \code{rejections} Integer, total number of rejections
#' @return \code{rejections_list} data frame, list of rejected pvalues and the
#' corresponding filter statistics
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
#' # general use
#' results <- opw(pvalue = pvals, filter = filters, effectType = "continuous",
#'                                               method = "BH")
#'
#' # supply the mean effects for both the filters and the tests externally
#' mod <- lm(log(filters) ~ tests)
#' et = mean(tests)
#' ey = mod$coef[[1]] + mod$coef[[2]]*et
#' results2 <- opw(pvalue = pvals, filter = filters,
#'                mean_filterEffect = ey, mean_testEffect = et, tail = 2,
#'                effectType = "continuous", method = "BH")
#'
#' # supply the rank probabilities externally
#' library(qvalue)
#' ranks <- 1:m
#' nullProp = qvalue(p = pvals, pi0.method = "bootstrap")$pi0
#' m0 = ceiling(nullProp*m)
#' m1 = m - m0
#' probs <- sapply(ranks, prob_rank_givenEffect, et = ey, ey = ey,
#'                                         nrep = 10000, m0 = m0, m1 = m1)
#' results3 <- opw(pvalue = pvals, filter = filters, ranksProb = probs,
#'                  effectType = "continuous", tail = 2, method = "BH")
#'
#' # supply weight externally
#' wgt <- weight_continuous(alpha = .05, et = et, m = m, ranksProb = probs)
#' results4 <- opw(pvalue = pvals, filter = filters, weight = wgt,
#'                         effectType = "continuous", alpha = .05, method = "BH")
#'
#===============================================================================
# # function to apply opw methods on data
#---------------------------------------------------
# Input:
#----------------------------
# pvalue = vector of pvalues
# filter = vector of filter statistics
# ranksProb = the probabilities of the filters or filters' ranks given
# the mean of the filter effects
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
# ranksProb = probailities of the ranks given the mean effect size
# wgt = weights
# Data = create a data set
# OD = odered by covariate
# odered.pvalues = odered pvalues for all tests
# padj = adjusted pvalues for FDR uses

# Output:
#-------------------------
# totalTests = total number of hypothesis tests
# nullProp = estimated propotion of the true null hypothesis
# ranksProb = probability of the ranks given the mean filter effect,
#                                           p(rank | ey = mean_filterEffet)
# weight = normalized weight
# rejections = total number of rejections
# rejections_list = list of rejected pvalues and corresponding filter statistics
#-------------------------------------------------------------------------------

opw <- function(pvalue, filter, weight = NULL, ranksProb = NULL, mean_filterEffect = NULL,
                mean_testEffect = NULL, effectType = c("continuous", "binary"),
                alpha = .05, nrep = 10000, tail = 1L, delInterval = .0001,
                method = c("BH", "BON"), ... )
    {
        # compute the number of tests------------
        m = length(pvalue)
        nullProp = qvalue(p = pvalue, pi0.method = "bootstrap")$pi0
        m0 = ceiling(nullProp*m)
        m1 = m - m0


        # formulate a data set-------------
        Data = tibble(pvalue, filter)
        OD <- Data[order(Data$filter, decreasing=TRUE), ]
        Ordered.pvalue <- OD$pvalue


        #check whether weight is provided------------
        if(!is.null(weight)){
            wgt <- weight
        } else {

            # compute test statistics from the pvalues---------
            test <- qnorm(pvalue/tail, lower.tail = FALSE)
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


            #check whether filter ranks probability is provided------------
            if(!is.null(ranksProb)){
                ranksProb <- ranksProb
            } else {
                # estimate the mean filter effect size------------------
                if(!is.null(mean_filterEffect)){
                    mean_filterEffect <- mean_filterEffect
                } else {
                    # check whether the filters are positive for the box-cox--------
                    if(any(filter <= 0)){
                        stop("filter statistics need to be positive")
                    }

                    bc <- boxcox(filter ~ test)
                    lambda <- bc$x[which.max(bc$y)]

                    if(lambda == 0){
                        model <- lm(log(filter + .0001) ~ test)
                    } else {
                        model <- lm(filter**lambda ~ test)
                    }
                    mean_filterEffect <- model$coef[[1]] + model$coef[[2]]*mean_testEffect
                }

                message("computing rank probabilities")
                # compute the probability of the rank of the filter given the mean effect
                ranksProb <- sapply(1:m, prob_rank_givenEffect, et = mean_filterEffect,
                               ey = mean_filterEffect, nrep = nrep, m0 = m0, m1 = m1)

                message("finished computing the rank probabilities")
            }

            # compute the weights (always right-tailed)------------
            message("computing weights")
            if(effectType == "continuous"){
                wgt = weight_continuous(alpha = alpha, et = mean_testEffect, m = m,
                                        tail = 1, delInterval = delInterval,
                                        ranksProb = ranksProb)
            } else {
                wgt = weight_binary(alpha = alpha, et = mean_testEffect, m = m,
                                    m1 = m1, tail = 1, delInterval = delInterval,
                                    ranksProb = ranksProb)
            }
            message("finished computing the weights")
        }

        message("comparing pvalues with thresholds")
        if(method == "BH"){
            padj <- p.adjust(Ordered.pvalue/wgt, method = "BH")
            rejections_list = OD[which((padj <= alpha) == TRUE), ]
        } else {
            rejections_list = OD[which((Ordered.pvalue <= alpha*wgt/m) == TRUE), ]
        }


        # outputs--------------
        n_rejections = dim(rejections_list)[1]

        return(list(totalTests = m, nullProp = nullProp,
                    ranksProb = ranksProb, weight = wgt,
                    rejections = n_rejections, rejections_list = rejections_list))
    }





