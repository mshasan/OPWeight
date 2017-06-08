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
        OD <- Data[order(Data$filter, decreasing=T), ]
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
                                        tail = 1, delInterval = delInterval, ranksProb = ranksProb)
            } else {
                wgt = weight_binary(alpha = alpha, et = mean_testEffect, m = m, m1 = m1,
                                    tail = 1, delInterval = delInterval, ranksProb = ranksProb)
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





