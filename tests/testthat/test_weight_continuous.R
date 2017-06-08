
# function to compute weight continuous case
#--------------------------------------
weight_continuous <- function(alpha, et, m, tail = 1L, delInterval = .0001, ranksProb)
{
    prob <- ranksProb/sum(ranksProb, na.rm = T)
    delta <- seq(0, 1, delInterval)

    weightSumVec <- sapply(delta, weight_by_delta, alpha = alpha, et = et, m = m,
                           m1 = NULL, tail = tail, ranksProb = prob,
                           effectType = "continuous")

    deltaOut <- delta[min(abs(weightSumVec - m)) == abs(weightSumVec - m)]
    deltaOut <- ifelse(length(deltaOut) > 1, .0001, deltaOut)
    weight.out <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(deltaOut/(alpha*prob)),
                                       lower.tail=FALSE)
    sumWeight <- sum(weight.out, na.rm = TRUE)
    normWeight <- if(sumWeight == 0) {rep(1, m)} else {weight.out/sumWeight*m}
    return(normWeight)
}
