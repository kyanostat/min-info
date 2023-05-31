source("min-info.R")

demo_paper_ex1_2 = function(n=50, d=4, rho=0.5, LOOP=200, fname="demo_paper_ex1_2_time.R"){
	set.seed(1)
	S0 = toeplitz(rho^(0:(d-1)))  # AR model
	Sinv0 = solve(S0)
	th0 = c(-Sinv0[lower.tri(Sinv0)])  # true value
	C = chol(S0)

	time_ple = numeric(LOOP)
	time_cle = numeric(LOOP)

	K = d*(d-1)/2
	th_s = matrix(0, LOOP, K)
	th_mle_s = matrix(0, LOOP, K)
	th_ple_s = matrix(0, LOOP, K)
	for(Li in 1:LOOP){
		U = matrix(rnorm(n*d), n, d)
		X = U %*% C

		S = cov(X)
		Sinv = solve(S)
		th_mle = c(-Sinv[lower.tri(Sinv)])  # MLE based on Gaussian
		th_mle_s[Li,] = th_mle

		time_ple[Li] = system.time(result <- min_info_Besag(X, h_2dim, Mar_2dim, thin=n))[1]
		th_ple = result$th
		th_ple_s[Li,] = th_ple
		
		time_cle[Li] = system.time(result <- min_info_mle(X, h_2dim, Mar_2dim, th=th_ple))[1]
		th = result$th
		th_s[Li,] = th

		cat(Li, "-th turn\n")
		cat("th", th, "\n")
		cat("th_mle", th_mle, "\n")
		cat("th_ple", th_ple, "\n")
	}
#	cat("CLE bias", colMeans(th_s)-th0, "standard deviation", apply(th_s,2,sd), "\n")
#	cat("MLE bias", colMeans(th_mle_s)-th0, "standard deviation", apply(th_mle_s,2,sd), "\n")
#	cat("PLE bias", colMeans(th_ple_s)-th0, "standard deviation", apply(th_ple_s,2,sd), "\n")
	th0_s = matrix(th0, LOOP, K, byrow=TRUE)
	cat("CLE RMSE", sqrt(colMeans((th_s-th0_s)^2)), "\n")
	cat("MLE RMSE", sqrt(colMeans((th_mle_s-th0_s)^2)), "\n")
	cat("PLE RMSE", sqrt(colMeans((th_ple_s-th0_s)^2)), "\n")
	result = list(th_s = th_s, th_mle_s = th_mle_s, th_ple_s = th_ple_s, th0=th0, time_cle=time_cle, time_ple=time_ple)
	dump("result", fname)
	invisible(result)
}
