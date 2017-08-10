m = 100
ranks = 1:m
prob <- sapply(ranks, prob_rank_givenEffect, et = 2, ey = 1, nrep = 10000,
               m0 = 50, m1 = 50)
prob <- prob/sum(prob, na.rm = TRUE)

w = weight_continuous_nwt(alpha = .05, et = 2, ranksProb = prob)$w

test_that("sum of ranks prob. is 1", {
    expect_equal(sum(prob, na.rm = TRUE), 1.0)
})

test_that("sum of weights is m", {
    expect_equal(sum(w, na.rm = TRUE), m)
})

if (any(w < 0)) {
    cat("Some weights are negative\n")
}

