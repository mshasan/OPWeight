
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


    # formulate a data set-------------
    Data = tibble(pvalue, covariate)
    OD <- Data[order(Data$covariate, decreasing=TRUE), ]
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

                bc <- MASS::boxcox(covariate ~ test)
                lambda <- bc$x[which.max(bc$y)]

                if(lambda == 0){
                    model <- lm(log(covariate + .0001) ~ test)
                } else {
                    model <- lm(covariate**lambda ~ test)
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
    n_rejections = sum(dataOut$nullReject == TRUE)

    return(list(totalTests = m,
                nullProp = nullProp,
                rejections = n_rejections,
                rejections_list = rejections_list,
                dataOut = dataOut))
}
