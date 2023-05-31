source("min-info.R")

# "hybrid data"
# beta-Poisson, h(x,y) = x / (1 + y)

demo_hybrid_repeat = function(th_s = c(0, 10, 100), LOOP = 200, n = 50, N=1000){
  th_n = length(th_s)
  result = list()
  for(i in 1:th_n){
    show(start <- date())
    th = th_s[i]
    result[[i]] = demo_hybrid_main(n=n, N=N, th=th, fname=paste0("demo_hybrid_result-", i, ".R"), LOOP=LOOP)
    result[[i]]$start = start
    result[[i]]$goal = date()
  }
  invisible(result)
}

h_hybrid = function(Z, X_aux){
	if(is.vector(Z)) Z = matrix(Z, 1, length(Z))
	H = matrix(0, nrow(Z), 1)
	X = X_aux$X
	H[,1] = X[Z[,1], 1] / (1 + X[Z[,2], 2])
	return(H)
}

Mar_hybrid = function(X){
	X_aux = list(X=X)
	return(X_aux)
}

### main body

demo_hybrid_main = function(n=50, N=1000, th0=1, LOOP=200, fname="result.R"){
	set.seed(2)

	K = 1	
	th_s = matrix(0, LOOP, K)
	th_mle_s = matrix(0, LOOP, K)
	th_ple_s = matrix(0, LOOP, K)
	for(Li in 1:LOOP){
		# marginal sampling
		X_pop = matrix(0, N, 2)
		X_pop[, 1] = rbeta(N, 10, 10)
		X_pop[, 2] = rpois(N, 3)

		# exchange algorithm
		X_pop = min_info_sampling(X_pop, th0, h_hybrid, Mar_hybrid) # "population"
		X = X_pop[sample(1:N, n), ] # sample

		# PLE
		result = min_info_Besag(X, h_hybrid, Mar_hybrid, thin=n, eps=1e-4)
		th_ple = result$th
		th_ple_s[Li,] = th_ple

		# conditional MLE
		result = min_info_mle(X, h_hybrid, Mar_hybrid, th=th_ple, thin=n, eps=1e-4)
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
