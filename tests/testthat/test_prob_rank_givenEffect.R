prob_rank_givenEffect <- function(k, et, ey, nrep = 10000, m0, m1)
	{
        m = m0 + m1
		t <- rnorm(nrep, et, 1)
		p0 <- pnorm(-t)
		p1 <- pnorm(ey - t)

		mean0 <- (m0 - 1)*p0 + m1*p1 + 1
		mean1 <- m0*p0 + (m1 - 1)*p1 + 1

		var0 <- (m0 - 1)*p0*(1 - p0) + m1*p1*(1 - p1)
		var1 <- m0*p0*(1 - p0) + (m1 - 1)*p1*(1 - p1)

		prob <- ifelse(et == 0, mean(dnorm(k, mean0, sqrt(var0))),
					   mean(dnorm(k, mean1, sqrt(var1))))

		return(prob)
	}






