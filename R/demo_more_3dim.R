source("min-info.R")

# 3-dim
demo_3dim_repeat = function(th_s = c(0, 1, 2), LOOP = 200, n=50, N=1000){
  th_n = length(th_s)
  result = list()
  for(i in 1:th_n){
    show(start <- date())
    # th = c(0, 0, 0, th_s[i])
    th = th_s[i] * c(1, 0, 0, -1)
    result[[i]] = demo_3dim_main(n=n, N=N, d=3, th=th, fname=paste0("demo_3dim_result-", i, ".R"), LOOP=LOOP)
    result[[i]]$start = start
    result[[i]]$goal = date()
  }
  invisible(result)
}

### main body

demo_3dim_main = function(n=50, N=1000, d=3, th0=c(0,0,0,1), LOOP=200, fname="result.R"){
	set.seed(1)
	
	K = d*(d-1)/2 + d*(d-1)*(d-2)/6
	th_s = matrix(0, LOOP, K)
	th_mle_s = matrix(0, LOOP, K)
	th_ple_s = matrix(0, LOOP, K)
	for(Li in 1:LOOP){
		# marginal sampling
		X_pop = matrix(rnorm(N*d), N, d)

		# exchange algorithm
		X_pop = min_info_sampling(X_pop, th0, h_3dim, Mar_3dim)  # "population"
		X = X_pop[sample(1:N, n), ]  # sample

		# PLE
		result = min_info_Besag(X, h_3dim, Mar_3dim, thin=n)
		th_ple = result$th
		th_ple_s[Li,] = th_ple
		
		# conditional MLE
		result = min_info_mle(X, h_3dim, Mar_3dim, thin=n, th=th_ple)
		th = result$th
		th_s[Li,] = th

		cat(Li, "-th turn\n")
		cat("th", th, "\n")
		cat("th_ple", th_ple, "\n")
	}
	th0_s = matrix(th0, LOOP, K, byrow=TRUE)
	cat("CLE RMSE", sqrt(colMeans((th_s-th0_s)^2)), "\n")
	cat("PLE RMSE", sqrt(colMeans((th_ple_s-th0_s)^2)), "\n")
	result = list(th_s = th_s, th_ple_s = th_ple_s, th0=th0)
	dump("result", fname)
	invisible(result)
}
