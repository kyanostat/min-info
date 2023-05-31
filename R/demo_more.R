source("min-info.R")

# Gaussian AR(1) with parameter rho
demo_AR = function(rho_s = c(0, 1/4, 1/2, 3/4)){
  rho_n = length(rho_s)
  result = list()
  for(i in 1:rho_n){
    result[[i]] = demo_paper_ex1_1(rho=rho_s[i])
  }
  dump("result", "demo_AR_result.R")
  invisible(result)
}

demo_AR_repeat = function(rho_s = c(0, 1/4, 1/2, 3/4), LOOP = 200){
  rho_n = length(rho_s)
  result = list()
  for(i in 1:rho_n){
    show(start <- date())
    result[[i]] = demo_paper_ex1_2(rho=rho_s[i], fname=paste0("demo_AR_result2-", i, ".R"), LOOP=LOOP)
    result[[i]]$start = start
    result[[i]]$goal = date()
  }
  invisible(result)
}

# Gaussian Exchangeable with common correlation rho
demo_EX = function(rho_s = c(0, 1/4, 1/2, 3/4)){
	rho_n = length(rho_s)
	result = list()
	for(i in 1:rho_n){
		result[[i]] = demo_paper_EX_1(rho=rho_s[i])
	}
	dump("result", "demo_EX_result.R")
	invisible(result)
}

demo_EX_repeat = function(rho_s = c(0, 1/4, 1/2, 3/4), LOOP = 200){
  rho_n = length(rho_s)
  result = list()
  for(i in 1:rho_n){
    show(start <- date())
    result[[i]] = demo_paper_EX_2(rho=rho_s[i], fname=paste0("demo_EX_result2-", i, ".R"), LOOP=LOOP)
    result[[i]]$start = start
    result[[i]]$goal = date()
  }
  invisible(result)
}

### main body

demo_paper_ex1_1 = function(n=50, d=4, rho=0.5){
	set.seed(1)
	
	S0 = toeplitz(rho^(0:(d-1)))  # AR model
	C = chol(S0)
	U = matrix(rnorm(n*d), n, d)
	X = U %*% C

	# CMLE
	result = min_info_mle(X, h_2dim, Mar_2dim)
	th = result$th
	th_se = result$th_se

	# MLE
	S = cov(X)
	Sinv = solve(S)
	th_mle = c(-Sinv[lower.tri(Sinv)])  # MLE based on Gaussian
	th_mle_se = numeric(d*(d-1)/2)
	n = nrow(X)
	cnt = 0
	for(i in 1:(d-1)){
		for(j in (i+1):d){
			cnt = cnt + 1
			th_mle_se[cnt] = sqrt(Sinv[i,i]*Sinv[j,j]+Sinv[i,j]^2)/sqrt(n)
			# the inverse of Fisher information is g^{ij,kl} = (S^{ii}*S^{jj} - S^{ij}*S^{ij})
		}
	}
	
	# true
	S0inv = solve(S0)
	th0 = c(-S0inv[lower.tri(S0inv)])

	K = length(th)
	par(cex=1.3, lwd=2)
	plot(1:K, th, xlab=expression(k), ylab=expression(theta[k]), ylim=range(c(th-2*th_se, th+2*th_se, th_mle-2*th_mle_se, th_mle+2*th_mle_se)), xlim=c(0.5, K+0.5))
	segments(-1, 0, K+1, 0)
	for(i in 1:K){
		arrows(i, th[i] - 1.96*th_se[i], i, th[i] + 1.96*th_se[i], angle=90, code=3, length=0.1)
	}
	points(1:K+0.2, th_mle, pch=4)
	for(i in 1:K){
		arrows(i+0.2, th_mle[i] - 1.96*th_mle_se[i], i+0.2, th_mle[i] + 1.96*th_mle_se[i], angle=90, code=3, length=0.1)
	}
#	dev.copy2pdf(file="2-dim.pdf")
	result = list(th=th, th_se=th_se, th_mle=th_mle, th_mle_se=th_mle_se, th0=th0)
#	dump("result", "demo_paper_ex1_1.R")
	invisible(result)
}

demo_paper_ex1_2 = function(n=50, d=4, rho=0.5, LOOP=200, fname="demo_paper_ex1_2.R"){
	set.seed(1)
	S0 = toeplitz(rho^(0:(d-1)))  # AR model
	Sinv0 = solve(S0)
	th0 = c(-Sinv0[lower.tri(Sinv0)])  # true value
	C = chol(S0)

	K = d*(d-1)/2
	th_s = matrix(0, LOOP, K)
	th_mle_s = matrix(0, LOOP, K)
	th_ple_s = matrix(0, LOOP, K)
	for(Li in 1:LOOP){
		U = matrix(rnorm(n*d), n, d)
		X = U %*% C

		# MLE
		S = cov(X)
		Sinv = solve(S)
		th_mle = c(-Sinv[lower.tri(Sinv)])  # MLE based on Gaussian
		th_mle_s[Li,] = th_mle
		
		# PLE
		result = min_info_Besag(X, h_2dim, Mar_2dim, thin=n)
		th_ple = result$th
		th_ple_s[Li,] = th_ple
		
		# CLE
		result = min_info_mle(X, h_2dim, Mar_2dim, thin=n, th=th_ple)
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
	result = list(th_s = th_s, th_mle_s = th_mle_s, th_ple_s = th_ple_s, th0=th0)
	dump("result", fname)
	invisible(result)
}

demo_paper_EX_1 = function(n=50, d=4, rho=0.5){
	set.seed(1)
	
	S0 = matrix(rho, d, d)  # exchange model
	diag(S0) = 1
#	show(S0)
	C = chol(S0)
	U = matrix(rnorm(n*d), n, d)
	X = U %*% C

	# CMLE
	result = min_info_mle(X, h_2dim, Mar_2dim)
	th = result$th
	th_se = result$th_se

	# MLE
	S = cov(X)
	Sinv = solve(S)
	th_mle = c(-Sinv[lower.tri(Sinv)])  # MLE based on Gaussian
	th_mle_se = numeric(d*(d-1)/2)
	n = nrow(X)
	cnt = 0
	for(i in 1:(d-1)){
		for(j in (i+1):d){
			cnt = cnt + 1
			th_mle_se[cnt] = sqrt(Sinv[i,i]*Sinv[j,j]+Sinv[i,j]^2)/sqrt(n)
			# the inverse of Fisher information is g^{ij,kl} = (S^{ii}*S^{jj} - S^{ij}*S^{ij})
		}
	}
	
	# true
	S0inv = solve(S0)
	th0 = c(-S0inv[lower.tri(S0inv)])

	K = length(th)
	par(cex=1.3, lwd=2)
	plot(1:K, th, xlab=expression(k), ylab=expression(theta[k]), ylim=range(c(th-2*th_se, th+2*th_se, th_mle-2*th_mle_se, th_mle+2*th_mle_se)), xlim=c(0.5, K+0.5))
	segments(-1, 0, K+1, 0)
	for(i in 1:K){
		arrows(i, th[i] - 1.96*th_se[i], i, th[i] + 1.96*th_se[i], angle=90, code=3, length=0.1)
	}
	points(1:K+0.2, th_mle, pch=4)
	for(i in 1:K){
		arrows(i+0.2, th_mle[i] - 1.96*th_mle_se[i], i+0.2, th_mle[i] + 1.96*th_mle_se[i], angle=90, code=3, length=0.1)
	}
#	dev.copy2pdf(file="2-dim.pdf")
	result = list(th=th, th_se=th_se, th_mle=th_mle, th_mle_se=th_mle_se, th0=th0)
#	dump("result", "demo_paper_ex1_1.R")
	invisible(result)
}

demo_paper_EX_2 = function(n=50, d=4, rho=0.5, LOOP=200, fname="demo_paper_EX_2.R"){
	set.seed(1)
	S0 = matrix(rho, d, d)  # exchange model
	diag(S0) = 1
	Sinv0 = solve(S0)
	th0 = c(-Sinv0[lower.tri(Sinv0)])  # true value
	C = chol(S0)

	K = d*(d-1)/2
	th_s = matrix(0, LOOP, K)
	th_mle_s = matrix(0, LOOP, K)
	th_ple_s = matrix(0, LOOP, K)
	for(Li in 1:LOOP){
		U = matrix(rnorm(n*d), n, d)
		X = U %*% C

		# MLE
		S = cov(X)
		Sinv = solve(S)
		th_mle = c(-Sinv[lower.tri(Sinv)])  # MLE based on Gaussian
		th_mle_s[Li,] = th_mle
		
		# PLE
		result = min_info_Besag(X, h_2dim, Mar_2dim, thin=n)
		th_ple = result$th
		th_ple_s[Li,] = th_ple

		# CLE
		result = min_info_mle(X, h_2dim, Mar_2dim, thin=n, th=th_ple)
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
	result = list(th_s = th_s, th_mle_s = th_mle_s, th_ple_s = th_ple_s, th0=th0)
	dump("result", fname)
	invisible(result)
}
