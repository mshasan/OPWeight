#' @title Weight for the continuous effect sizes
#'
#' @description This is a Newton Raphson based algorithm to compute weight from
#' the ranks probability for the continuous effect sizes.
#' @param alpha Numeric, significance level of the hypothesis test
#' @param et Numeric, mean effect size of the test statistics
#' @param ranksProb Numeric vector of ranks probability of the tests given
#' the effect size
#' @param x0 Numeric, a initial value for the Newton-Raphson mehod.
#'
#' @details
#' If one wants to test \deqn{H_0: epsilon_i = 0 vs. H_a: \epsilon_i > 0,}
#' then \code{et} should be mean of the test effect sizes.
#' This is called hypothesis testing for the continuous effect sizes.\cr
#'
#' For the Newton-Raphson mehthod the initial value \code{x0} is very important.
#' One may need to supply externally if the function shows error. Newton-Raphson
#' method may not work in some cases. In that situations, use \code{weight_continuous}
#' function, which is a general but slow approach to compute weights beased on the
#' grid search algorithm.
#'
#' @author Mohamad S. Hasan, shakilmohamad7@gmail.com
#'
#' @export
#'
#' @import OPWeight prob_rank_givenEffect
#'
#' @seealso \code{\link{prob_rank_givenEffect}} \code{\link{weight_continuous}}
#'
#' @return \code{w} Numeric vector of weights of the tests for the continuous case
#' @return \code{lambda} Numeric value of the LaGrange multiplier
#' @return \code{n} Integer, number of iteration needed to obtain the initial value
#' @return \code{k} Integer, number of iteration needed to obtain the optimal
#' value of \code{lambda}
#'
#' @examples
#'
#' # compute the probabilities of the ranks of a test being rank 1 to 100 if the
#' # targeted test effect is 2 and the overall mean covariate effect is 1.
#' ranks <- 1:100
#' prob <- sapply(ranks, prob_rank_givenEffect, et = 2, ey = 1, nrep = 10000,
#'                               m0 = 50, m1 = 50)
#'
#' # compute weight for the continuous case
#' results = weight_continuous_nwt(alpha = .05, et = 2,
#'                                              ranksProb = prob)
#'
#===============================================================================
# function to compute weight from p(rank=k|covariateEffect=ey)
#-----------------------------------------------------
# internal parameters:-----
# m = total numebr of tests
# f = function to be optimized
# df = first of derivative of the function f
#===============================================================================

weight_continuous_nwt <- function(alpha, et, ranksProb, x0 = NULL)
{
    m = length(ranksProb)
    prob <- ranksProb/sum(ranksProb, na.rm = TRUE)

    # define the functions used by Newton's method---------
    # the goal is to find the zero of the function f below
    f <- function(c)
    {
        sum(pnorm(et/2 + c/et - log(alpha*prob)/et, lower.tail = FALSE)) - alpha
    }

    # first derivative of f---------
    df <- function(c)
    {
        sum(-dnorm(et/2 + c/et - log(alpha*prob)/et)/et)
    }

    # applying Newton's method--------------
    # get appropriate initial value automatically
    # need to make sure that f(x0) > 0
    nmax <- 1000
    n = 1
    if(is.null(x0)) {x0 = 0} else {x0 = x0}
    while(f(x0) < 0 && (n <= nmax))
    {
        if(f(0) > f(.5)){x0 = x0 - .5} else {x0 = x0 + .5}
        n = n + 1
    }

    if(f(x0) < 0 && (n >= nmax)) {
        cat("f(x0) = ", f(x0))
        stop("f(x0) is negative, modify initial guess x0 so that f(x0) >= 0")
    }

    #the vector x will hold the trajectory
    x <- rep(0, nmax)
    #ex will hold the proportion of error
    ex <- rep(0, nmax)

    # Newton's method iterations----------
    # Compute first step, and error
    x[1] <- x0 - f(x0) / df(x0)
    ex[1] <- abs((x[1] - x0)/x[1])
    k <- 2

    while ((ex[k - 1] >= .0001) && (k <= nmax))
    {
        x[k] <- x[k - 1] - (f(x[k - 1]) / df(x[k - 1]))
        ex[k] <- abs((x[k] - x[k - 1])/x[k])
        k <- k + 1
    }

    #store the last point in the trajectory
    c <- x[k - 1] # could return this

    #compute weights and set output
    w <- (m/alpha)*pnorm(et/2 + c/et - log(alpha*prob)/et, lower.tail = FALSE)
    lambda <- exp(c)

    # error checking---------
    epsi <- 1e-2
    if (abs(sum(w, na.rm = TRUE) - m) > epsi) {
        cat("Warning: Weights do not average to 1. Mean weight = ",
                                        mean(w, na.rm = TRUE), "\n")
    }

    if (any(w < 0)) {
        cat("Some weights are negative\n")
    }

    results <- list(w = w, lambda = lambda, n = n, k = k)
    return(results)
}

