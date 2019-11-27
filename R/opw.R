#' @title Perform Optimal Pvalue Weighting
#'
#' @description A function to perform weighted pvalue multiple hypothesis test.
#' This function compute the probabilities of the ranks of the covariate statistics
#' given the effect sizes, and consequently the weights if neighter the weights
#' nor the probabilities are given. Then provides the number of rejected null
#' hypothesis and the list of the rejected pvalues as well as the corresponing
#' covariate statistics.
#'
#' @param pvalue Numeric vector of pvalues of the test statistics
#' @param covariate Numeric vector of covariate statistics
#' @param weight An optional numeric weight vector not required
#' @param ranksProb An optional numeric vector of the ranks probability of the
#' covariates given the mean effect
#' @param mean_covariateEffect Numeric, value of the mean covariate effect
#' of the true alternatives
#' @param mean_testEffect Numeric, value of the mean test effect
#' of the true alterantives
#' @param effectType Character ("continuous" or "binary"), type of effect sizes
#' @param x0 Numeric, a initial value for the Newton-Raphson mehod.
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
#' then the \code{mean_testEffect}  and \code{mean_covariateEffect} should be mean
#' of the test and covariate effect sizes, respectively. This is called hypothesis
#' testing for the continuous effect sizes.\cr
#'
#' If one wants to test \deqn{H_0: epsilon_i = 0 vs. H_a: epsilon_i = epsilon,}
#' then \code{mean_testEffect} and \code{mean_covariateEffect} should be median or
#' any discrete value of the test and covariate effect sizes. This is called hypothesis
#' testing for the Binary effect sizes, where \code{epsilon} refers to a fixed value.\cr
#'
#' The main goal of the function is to compute the probabilities of the ranks from
#' the pvalues and the covariate statistics, consequently the weights. Although \code{weights}
#' \code{ranksProb} are optional, \code{opw} has the options so that one can compute
#' the probabilities and the weights externally if necessary (see examples).\cr
#'
#' Internally, \code{opw} function compute the \code{ranksProb} and consequently
#' the weights, then uses the pvalues to make conclusions about hypotheses.
#' Therefore, if \code{ranksProb} is given then \code{mean_covariateEffect}
#' and are redundant, and should not be provided to the funciton.
#' Although \code{ranksProb} is not required to the function,
#' One can compute \code{ranksProb} by using the function
#' \code{\link{prob_rank_givenEffect}}.\cr
#'
#' The function internally compute \code{mean_covariateEffect} and \code{mean_testEffect}
#' from a simple linear regression with box-cox transformation between the test
#' and covariate statistics, where the covariates are regressed on the test statistics.
#' Thus, covariates need to be positive to apply \code{boxcox} from the \code{R}
#' library \code{MASS}. Then the estimated \code{mean_covariateEffect} and
#' \code{mean_testEffect} are used to obtian the \code{ranksProb} and the weights.
#' Thus, in order to apply the function properly, it is crucial to understand the
#' uses \code{mean_covariateEffect} and \code{mean_testEffect}. If \code{mean_covariateEffect} and
#' \code{mean_testEffect} are not provided then the test statistics computed from
#' the pvalues will be used to compute the relationship between the covariate
#' statistics and the test statistics.\cr
#'
#' If one of the mean effects \code{mean_covariateEffect} and \code{mean_testEffect}
#' are not provided then the missing mean effect will be computed internally.
#'
#'
#' @author Mohamad S. Hasan, shakilmohamad7@gmail.com
#'
#' @export
#'
#' @import qvalue qvalue
#' @import tibble tibble
#' @import MASS
#'
#' @seealso \code{\link{prob_rank_givenEffect}} \code{\link{weight_binary}}
#' \code{\link{weight_continuous}} \code{\link{qvalue}} \code{\link{dnorm}}
#'
#'
#' @return \code{totalTests} Integer, total number of hypothesis tests evaluated
#' @return \code{nullProp} Numeric, estimated propotion of the true null
#' hypothesis
#' @return \code{rejections} Integer, total number of rejections
#' @return \code{mean_testEffect} Numeric, mean of the alternative test effects
#' @return \code{mean_covariateEffect} Numeric, mean of the covariate effects
#' corresponding to the mean of the alternative test effects
#' @return \code{rejections_list} Data frame, list of rejected p-values and the
#' corresponding covariate statistics and the adjusted p-values if method = "BH" used.
#' @return \code{dataOut} Data frame, ordered covariate and the corresponding pvalues,
#' ranks probabilities, weights, weight-adjusted pvalues, and decision of null rejection.
#'
#' @examples
#' # generate pvalues and covariate statistics
#' m = 1000
#' set.seed(3)
#' covariates = runif(m, min = 1, max = 3)          # covariate statistics
#' H = rbinom(m, size = 1, prob = 0.1)             # hypothesis true or false
#' tests = rnorm(m, mean = H * covariates)            # Z-score
#' pvals = 1 - pnorm(tests)                        # pvalue
#'
#' # general use
#' results <- opw(pvalue = pvals, covariate = covariates, effectType = "continuous",
#'                                               method = "BH")
#'
#' # supply the mean effects for both the covariates and the tests externally
#' library(qvalue)
#' ranks <- 1:m
#' nullProp = qvalue(p = pvals, pi0.method = "bootstrap")$pi0
#' m0 = ceiling(nullProp*m)
#' m1 = m - m0
#' mod <- lm(log(covariates) ~ tests)
#' et = mean(sort(tests, decreasing = TRUE)[1:m1], na.rm = TRUE)
#' ey = mod$coef[[1]] + mod$coef[[2]]*et
#' results2 <- opw(pvalue = pvals, covariate = covariates,
#'                mean_covariateEffect = ey, mean_testEffect = et, tail = 2,
#'                effectType = "continuous", method = "BH")
#'
#' # supply the rank probabilities externally
#' probs <- sapply(ranks, prob_rank_givenEffect, et = ey, ey = ey,
#'                                         nrep = 10000, m0 = m0, m1 = m1)
#' results3 <- opw(pvalue = pvals, covariate = covariates, ranksProb = probs,
#'                  effectType = "continuous", tail = 2, method = "BH")
#'
#' # supply weight externally
#' wgt <- weight_continuous(alpha = .05, et = et, m = m, ranksProb = probs)
#' results4 <- opw(pvalue = pvals, covariate = covariates, weight = wgt,
#'                         effectType = "continuous", alpha = .05, method = "BH")
#'
#===============================================================================
# # function to apply opw methods on data
#---------------------------------------------------
# internal parameters:-----
# m = number of hypothesis test
# nullProp = proportion of true null hypothesis
# m0 =  number of the true null tests
# m1 = number of the true alternative tests
# test =  compute test statistics from the pvalues if not given
# test_effect_vec = estiamted number of the true alternaitve test statistics
# mean_testEffect = mean test effect sizes of the true alternaive hypotheis
# mean_covariateEffect = mean covariate effect sizes of the true alternaive hypotheis
# ranksProb = probailities of the ranks given the mean effect size
# wgt = weights
# Data = create a data set
# OD = odered by covariate
# odered.pvalues = odered pvalues for all tests
# padj = adjusted pvalues for FDR uses
#===============================================================================

opw <- function(pvalue, covariate, weight = NULL, ranksProb = NULL,
            mean_covariateEffect = NULL, mean_testEffect = NULL,
            effectType = c("continuous", "binary"), x0 = NULL,
            alpha = .05, nrep = 10000, tail = 2L, delInterval = .001,
            method = c("BH", "BON"), ... )
    {
        # compute the number of tests------------
        m = length(pvalue)
        nullProp = qvalue(p = pvalue, pi0.method = "bootstrap")$pi0
        m0 = ceiling(nullProp*m)
        m1 = m - m0


        # formulate a data set and order by covariate-------------
        Data = tibble(pvalue, covariate)
        OD <- Data[order(Data$covariate, decreasing=TRUE), ]
        Ordered.pvalue <- OD$pvalue


        #check whether weight is provided------------
        if(!is.null(weight)){
            wgt <- weight
        } else {

            # compute test statistics from the pvalues and
            # keep top m1 tests for regression ---------
            test <- qnorm(pvalue/tail, lower.tail = FALSE)
            test[which(!is.finite(test))] <- NA
            Data2 = add_column(Data, test)
            OD2 <- Data2[order(Data2$test, decreasing=TRUE), ][1:m1, ]

            # estimate the true alterantive test effect sizes----------------
            if(m1 == 0){
                test_effect_vec <- 0
            } else {
                test_effect_vec <-  OD2$test
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


            #check whether covariate ranks probability is provided------------
            if(!is.null(ranksProb)){
                ranksProb <- ranksProb
            } else {
                # estimate the mean covariate effect size------------------
                if(!is.null(mean_covariateEffect)){
                    mean_covariateEffect <- mean_covariateEffect
                } else {
                    # check whether the covariates are positive for the box-cox--------
                    if(any(covariate <= 0)){
                        stop("covariate statistics need to be positive")
                    }

                    bc <- MASS::boxcox(covariate ~ test, data = OD2)
                    lambda <- bc$x[which.max(bc$y)]

                    if(lambda == 0){
                        model <- lm(log(covariate + .0001) ~ test, data = OD2)
                    } else {
                        model <- lm(covariate**lambda ~ test, data = OD2)
                    }
                    mean_covariateEffect <- model$coef[[1]] + model$coef[[2]]*mean_testEffect
                }

                message("computing rank probabilities")
                # compute the probability of the rank of the covariate given the mean effect
                ranksProb <- sapply(1:m, prob_rank_givenEffect, et = mean_covariateEffect,
                               ey = mean_covariateEffect, nrep = nrep, m0 = m0, m1 = m1)

                message("finished computing the rank probabilities")
            }

            # compute the weights (always right-tailed)------------
            message("computing weights")
            if(effectType == "continuous"){
                if(mean_testEffect > 1){
                    wgt = weight_continuous_nwt(alpha = alpha, et = mean_testEffect,
                                                ranksProb = ranksProb)$w
                } else {
                    wgt = weight_continuous(alpha = alpha, et = mean_testEffect, m = m,
                                            tail = tail, delInterval = delInterval,
                                            ranksProb = ranksProb)
                }

            } else {
                if(mean_testEffect > 1){
                    wgt = weight_binary_nwt(alpha = alpha, m1 = m1,
                                et = mean_testEffect, ranksProb = ranksProb)$w
                } else {
                    wgt = weight_binary(alpha = alpha, et = mean_testEffect, m = m,
                                        m1 = m1, tail = tail, delInterval = delInterval,
                                        ranksProb = ranksProb)
                }
            }
            message("finished computing the weights")
        }

        message("comparing pvalues with thresholds")
        if(method == "BH"){
            padj <- p.adjust(Ordered.pvalue/wgt, method = "BH")
            #OD <- add_column(OD, padj)
            reject <- (padj <= alpha)
            #rejections_list = OD[which(reject == TRUE), ]
        } else {
            reject <- Ordered.pvalue <= alpha*wgt/m
            #rejections_list = OD[which(reject == TRUE), ]
        }

        # outputs--------------
        if(!is.null(weight)){
          dataOut <- data.frame(OD, weight=wgt, adjPvalue=padj, nullReject=reject)
        } else {
          dataOut <- data.frame(OD, ranksProb, weight=wgt, adjPvalue=padj, nullReject=reject)
        }

        rejections_list = dataOut[which(dataOut$nullReject == TRUE), ]
        n_rejections = length(which(dataOut$nullReject == TRUE))

        return(list(totalTests = m,
                    nullProp = nullProp,
                    mean_testEffect = mean_testEffect,
                    mean_covariateEffect = mean_covariateEffect,
                    rejections = n_rejections,
                    rejections_list = rejections_list,
                    dataOut = dataOut))
    }





