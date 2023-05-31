# Minimum Information Dependence modeling
# last updated on 2021-01-20 (specification of "Mar" is changed)
#
# Usage:
#   min_info_Besag(X, h, Mar)
#   min_info_one_step(X, h, Mar)
#   min_info_mle(X, h, Mar)
#   min_info_mcmc(Z, th, h, L, ...)
#
# Arguments:
#   X: a matrix, a data frame or a list
#   h: a function of permutation Z and marginal information Mar (see h_2dim)
#   Mar: marginal information (function or matrix) (see Mar_2dim)
#   Z: permutation matrix
#   th: parameter
#   L: the number of MCMC steps
#
# Details:
#   Fit a minimum information dependence model via MCMC
#    p(x;th) = exp(sum_k(th_k * h_k(x)) - sum_i a_i(x_i;th) - psi(th))
#
# Value:
#   th: estimated parameter
#   th_se: standard error of estimated parameter

# sampling
min_info_mcmc = function(Z, th, h, L,...){
	# Z is permutation matrix (n*d-dim)
	# th is parameter (K-dim)
	# h is canonical statistic (function)
	# L is number of MCMC steps
	# ... are arguments for function h

	d = ncol(Z)
	M = nrow(Z)
	if(L <= 0) return(Z)
	M_2 = M %/% 2
	for(Li in 1:L){
		for(i in 1:d){
			st_pairs = sample(M) # batch sampling
			s = st_pairs[1:M_2]
			t = st_pairs[(M_2+1):(2*M_2)]
			z = zp = Z[s,]
			w = wp = Z[t,]
			zp[,i] = w[,i]  # swap
			wp[,i] = z[,i]  # swap
			del = c((h(z, ...) + h(w, ...) - h(zp, ...) - h(wp, ...)) %*% th)
			r = 1 / (1 + exp(-del))
			u = runif(M_2)
			Z[s,i] = ifelse(u > r, zp[,i], z[,i])
			Z[t,i] = ifelse(u > r, wp[,i], w[,i])
		}
	}
	invisible(Z)
}

min_info_sampling = function(X, th, h, Mar, L_b=10, L_n=300, monitor=FALSE, ...){
	# X is assumed to be a numeric matrix for simplicity
	# L_b: burn-in
	# if monitor = TRUE, return the whole samples and draw the autocorrelation function

	n = nrow(X)	
	d = ncol(X)
	X_aux = Mar(X)  # marginal information
	
	Z = Z0 = matrix(1:n, n, d, byrow=FALSE)  # identity permutation
	Z = min_info_mcmc(Z, th, h, L_b, X_aux, ...)  # burn-in
	if(!monitor){
		Z = min_info_mcmc(Z, th, h, L_n, X_aux, ...)
		for(i in 1:d){
			X[,i] = X[Z[,i], i]
		}
		return(X)
	}else{
		Z_s = array(0, c(L_n, n, d))
		X_s = array(0, c(L_n, n, d))
		for(L_i in 1:L_n){
			Z = min_info_mcmc(Z, th, h, 1, X_aux, ...)
			Z_s[L_i,,] = Z
			for(i in 1:d){
				X_s[L_i, , i] = X[Z[,i], i]
			}
		}
		t = sample(1:n, 1)
		acf(X_s[,t,])
		return(invisible(X_s))
	}
}

# Besag pseudo-MLE
min_info_Besag = function(X, h, Mar, th=NULL, eps=1e-3, max_iter=100, thin=5, ...){
	# X is data frame (n*d-dim) or list (length d, each element is a vector or a data frame)
	# h is canonical statistic (function)
	# Mar is marginal information (function or matrix)
	# thin is a thin-out parameter (expected number of edges per vertex)

	if(is.matrix(X)) X = as.data.frame(X)
	d = length(X)  # X is interpreted as a list
	if(is.matrix(X[[1]])){
		n = nrow(X[[1]])
	}else{
		n = length(X[[1]])
	}

	if(is.function(Mar)){  # marginal information
		X_aux = Mar(X)
	}else{
		X_aux = Mar
	}
	
	Z = matrix(1:n, n, d, byrow=FALSE)  # identity permutation
	h_bar = colMeans(h(Z, X_aux, ...))  # canonical statistic
	K = length(h_bar)  # dimension of parameter
	if(is.null(th)) th = rep(0, K)  # initial value

	err = Inf
	iter = 1
	t_s = list()
	for(s in 1:(n-1)){
		smp = sample(1:n, thin)
		smp = smp[smp > s]
		t_s[[s]] = smp  # sampled pairs
	}
	while(err > eps && iter <= max_iter){
		obj = 0
		gr = rep(0, K)
		hess = matrix(0, K, K)
		for(i in 1:d){
			for(s in 1:(n-1)){
				for(t in t_s[[s]]){
					z = zp = Z[s,]
					w = wp = Z[t,]
					zp[i] = w[i] # swap
					wp[i] = z[i] # swap
					H = c(h(z, X_aux, ...) + h(w, X_aux, ...) - h(zp, X_aux, ...) - h(wp, X_aux, ...))
					del = sum(th * H)
					obj = obj - log(1/(1+exp(-del)))  # negative sign for minimization
					gr = gr - exp(-del)/(1+exp(-del)) * H
					hess = hess + exp(-del)/(1+exp(-del))^2 * outer(H, H)
				}
			}
		}
		err = max(abs(gr))
		th = th - solve(hess, gr)  # Newton method
		cat("err (Besag) =", err, "  obj =", obj, "  max(abs(th)) =", max(abs(th)), "\n")
#		show(th)
		iter = iter + 1
	}
	if(iter > max_iter) warning("not converged")
	list(th=th, obj=obj, gr=gr, hess=hess)
}

# one-step estimator
min_info_one_step = function(X, h, Mar, th=NULL, L_b=10, L_n=300, return_detail=FALSE, only_res=FALSE, pseudo=FALSE, Z0=NULL, ...){
	# L_b: burn-in
	# L_n: number of MCMC samples

	if(is.matrix(X)) X = as.data.frame(X)
	d = length(X)  # X is interpreted as a list
	if(is.matrix(X[[1]])){
		n = nrow(X[[1]])
	}else{
		n = length(X[[1]])
	}

	if(is.function(Mar)){  # marginal information
		X_aux = Mar(X)
	}else{
		X_aux = Mar
	}

	Z = matrix(1:n, n, d, byrow=FALSE)  # identity permutation
	h_bar = colMeans(h(Z, X_aux, ...))  # canonical statistic
	K = length(h_bar)  # dimension of parameter
	if(is.null(th)) th = rep(0, K)  # initial value

	if(is.null(Z0)){
		Z0 = Z  # initialize by identity permutation (this is "cheat", but work)
	}
	H_s = array(0, c(L_n, n, K))
	Z_s = array(0, c(L_n, n, d))
	
	Z = min_info_mcmc(Z0, th, h, L_b, X_aux, ...)  # burn-in (batch sampling)
	for(L_i in 1:L_n){
		Z = min_info_mcmc(Z, th, h, 1, X_aux, ...)  # (batch sampling)
		H = h(Z, X_aux, ...)
		Z_s[L_i,,] = Z
		H_s[L_i,,] = H
	}
	eta = apply(H_s, 3, mean)  # mean parameter (length K)
	cat("err (MLE) =", max(abs(h_bar-eta)), "  L_n =", L_n, "  max(abs(th)) =", max(abs(th)), "\n")
	if(only_res) return(h_bar-eta)   # if only residual is evaluated
	FI = min_info_fisher(H_s)   # Fisher information
	if(pseudo){
		SI = min_info_pseudo_fisher(H_s)   # pseudo Fisher information (time consuming)
		FI = SI
	}
#	FIback = min_info_fisher_back(Z_s, H_s)   # Fisher information (by backfitting)
#	cat("check backfitting:", max(abs(FI-FIback)), "\n")

	th_v = solve(FI, h_bar - eta)  # Newton direction (Fisher)
	for(k in 1:5){  # line search
		LLR = min_info_LLR(H_s, th_v, h_bar)  # approximate log-likelihood ratio
		if(LLR > 0) break
		th_v = th_v / 2
		cat("line search:", k, "\n")
	}
	th = th + th_v  # new estimate

	V = solve(FI)  # asymptotic covariance estimate (valid only if pseudo==FALSE)
	th_se = sqrt(diag(V)) / sqrt(n)  # standard error

	if(return_detail) return(list(th=th, th_se=th_se, X_aux=X_aux, Z_s=Z_s, H_s=H_s, h_bar=h_bar, FI=FI, res=h_bar-eta))
	list(th=th, th_se=th_se, res=h_bar-eta)
}

# MLE
min_info_mle = function(X, h, Mar, th=NULL, L_b=10, L_n=300, eps=1e-2, max_iter=50, max_L=1e6, draw=TRUE, thin=5, pseudo_fix=FALSE, return_detail=FALSE, ...){
	# L_b: burn-in
	# L_n: number of MCMC samples
	
	if(is.matrix(X)) X = as.data.frame(X)
	d = length(X)  # X is interpreted as a list
	if(is.matrix(X[[1]])){
		n = nrow(X[[1]])
	}else{
		n = length(X[[1]])
	}

	if(is.function(Mar)){  # marginal information
		X_aux = Mar(X)
	}else{
		X_aux = Mar
	}
	
	# error estimate
	res = min_info_one_step(X, h, X_aux, th=th, L_b=L_b, L_n=L_n, only_res=TRUE, ...)
	K = length(res)
	if(draw){
		plot(c(0,max_iter), range(c(-res,res)), type="n", xlab="iteration", ylab="residual", col=1:K)
		segments(-1, 0, max_iter+1, 0)  # zero line
	}

	# Besag
	result = min_info_Besag(X, h, X_aux, th=th, thin=thin, ...)
	th = th0 = result$th

	# error estimate
	res_next = min_info_one_step(X, h, X_aux, th=th, L_b=L_b, L_n=L_n, only_res=TRUE, ...)
	if(draw){
		for(k in 1:K) segments(0, res[k], 1, res_next[k], col=k, lty=k)
	}
	res = res_next
	
	# MLE
	err = Inf
	iter = 2  # Besag is counted as 1 iteration
	while(err > eps && iter <= max_iter && L_n <= max_L){
		pseudo = ifelse(pseudo_fix || err > 1, TRUE, FALSE)  # adaptive
		result = min_info_one_step(X, h, X_aux, th=th, L_b=L_b, L_n=L_n, return_detail=return_detail, pseudo=pseudo, ...)
		th = result$th
		res_next = result$res
		if(draw){
			for(k in 1:K) segments(iter-1, res[k], iter, res_next[k], col=k, lty=k)
		}
		res = res_next
		err_next = max(abs(res))
		if(err_next > 1.1 * err) L_n = floor(1.5 * L_n)  # adaptive
		err = err_next
		if(max(abs(th) / (1 + abs(th0))) > 1e2){
		  th = th0
		  message("parameter reset.")
		}
		iter = iter + 1
	}
	if(iter > max_iter) warning("not converged")

	result
}

# MLE (calibration)
min_info_calibration = function(result, L_n=1e4, ...){
	# result is the result of min_info_mle() with return_detail=TRUE
	# result = list(th=th, th_se=th_se, X_aux=X_aux, Z_s=Z_s, H_s=H_s, h_bar=h_bar, FI=FI, res=h_bar-eta)
	th = result$th
	X_aux = X_aux
	L_n = max(L_n, dim(Z_s)[1])
	n = dim(Z_s)[2]
	d = dim(Z_s)[3]
	K = dim(H_s)[3]
	
	H_s = array(0, c(L_n, n, K))
	Z_s = array(0, c(L_n, n, d))
	
	Z = min_info_mcmc(Z0, th, h, L_b, X_aux, ...)  # burn-in (batch sampling)
	for(L_i in 1:L_n){
		Z = min_info_mcmc(Z, th, h, 1, X_aux, ...)  # (batch sampling)
		H = h(Z, X_aux, ...)
		Z_s[L_i,,] = Z
		H_s[L_i,,] = H
	}

	# under construction
}

# Fisher information
min_info_fisher = function(H_s){
	L = dim(H_s)[1]  # number of MCMC samples
	n = dim(H_s)[2]  # sample size
	K = dim(H_s)[3]  # dimension of parameter

	FI = matrix(0, K, K)   # Fisher information
	HH = matrix(0, L, K)
	for(k in 1:K){
		HH[,k] = rowSums(H_s[,,k])
	}
	for(k in 1:K){
		for(l in 1:K){
			FI[k,l] = cov(HH[,k], HH[,l]) * (L-1) / L / n
		}
	}
	FI
}

# Fisher information by backfitting algorithm
min_info_fisher_back = function(Z_s, H_s, eps=1e-3, max_iter=100){
	L = dim(H_s)[1]  # number of MCMC samples
	n = dim(H_s)[2]  # sample size
	K = dim(H_s)[3]  # dimension of parameter
	d = dim(Z_s)[3]  # dimension of data

	R_s = H_s   # residual
	for(k in 1:K){
		R_s[,,k] = R_s[,,k] - mean(R_s[,,k])
	}
	err = Inf
	iter = 0
	while(err > eps && iter < max_iter){
		R_s_pre = R_s
		for(i in 1:d){
			a = matrix(0, n, K)  # conditional mean given x_i
			for(L_i in 1:L){
				a[rank(Z_s[L_i,,i]),] = a[rank(Z_s[L_i,,i]),] + R_s[L_i,,]
			}
			a = a / L
			for(L_i in 1:L){
				R_s[L_i,,] = R_s[L_i,,] - a[rank(Z_s[L_i,,i]),]  # projection to i-th plane
			}
		}
		err = max(abs(R_s - R_s_pre))
		# show(fisher_residual_check(Z_s, R_s))  # To do
		iter = iter + 1
	}

	FI = matrix(0, K, K)   # Fisher information
	for(k in 1:K){
		for(l in 1:K){
			FI[k,l] = cov(c(R_s[,,k]), c(R_s[,,l])) * (L*n-1) / (L*n)
		}
	}
	FI
}

# pseudo Fisher information (greater than Fisher information)
min_info_pseudo_fisher = function(H_s){
	L = dim(H_s)[1]
	n = dim(H_s)[2]
	K = dim(H_s)[3]
	SI = matrix(0, K, K)
	H_vec = matrix(0, L*n, K)
	for(k in 1:K){
		H_vec[,k] = c(H_s[,,k])
	}
	for(k in 1:K){
		for(l in 1:K){
			SI[k,l] = cov(H_vec[,k], H_vec[,l]) * (L*n-1) / (L*n)  # covariance
		}
	}
	SI
}

# Log-likelihood ratio
min_info_LLR = function(H_s, th_v, h_bar){
	# th_v: parameter difference
	
	L_n = dim(H_s)[1]  # number of MCMC samples
	n = dim(H_s)[2]   # sample size
	K = dim(H_s)[3]   # dimension of parameter

	log_ratio = rep(0, L_n)
	for(k in 1:K){
		log_ratio = log_ratio + rowSums(th_v[k]*H_s[,,k])
	}
	baseline = max(log_ratio)
	LLR = n * sum(th_v * h_bar) - baseline - log(mean(exp(log_ratio - baseline)))
#	show(LLR)
	LLR
}

# likelihood ratio test
min_info_LRT = function(X, h, Mar, nzs, L_n=1e5, ...){
	# nzs specifies non-zero parameters in null hypothesis
	my_h = function(Z, X_aux, nzs){
		H = h(Z, X_aux)
		H[, nzs, drop=FALSE]  # nzs specifies non-zero parameters
	}
	result1 = min_info_mle(X, h, Mar, return_detail=TRUE, ...)  # alternative hypothesis
	th1 = result1$th
	if(any(nzs)){
		result0 = min_info_mle(X, my_h, Mar, return_detail=TRUE, nzs=nzs, ...)  # null hypothesis
		th0 = numeric(length(result1$th))
		th0[nzs] = result0$th
	}else{
		th0 = 0 * result1$th  # zero vector
	}
	# sampling at th = 0
	n = dim(result1$Z_s)[2]
	d = dim(result1$Z_s)[3]
	K = dim(result1$H_s)[3]
	H_s = array(0, c(L_n, n, K))
	Z_s = array(0, c(L_n, n, d))
	X_aux = Mar(X)
	Z = matrix(1:n, n, d)
	Z = min_info_mcmc(Z, th0, h, L_n/10, X_aux)  # burn-in
	for(L_i in 1:L_n){
#		for(i in 1:d){
#			Z[,i] = sample(1:n, n)
#		}
		Z = min_info_mcmc(Z, th0, h, 1, X_aux)
		H = h(Z, X_aux)
		Z_s[L_i,,] = Z
		H_s[L_i,,] = H
	}
	h_bar = result1$h_bar
	result = list()
	result$LLR = min_info_LLR(H_s, th1-th0, h_bar)
	result$th0 = th0
	result$th0_se = result0$th_se
	result$th1 = th1
	result$th1_se = result1$th_se
	result
}

# Wald test
#min_info_Wald = function(X, h, Mar, nzs, ...){
min_info_Wald = function(result, nzs, ...){
	my_h = function(Z, X_aux, nzs){
		H = h(Z, X_aux)
		H[, nzs, drop=FALSE]  # nzs specifies non-zero parameters
	}
#	result1 = min_info_mle(X, h, Mar, return_detail=TRUE, ...)  # alternative hypothesis
	n = dim(result$Z_s)[2]
	th1 = result$th
	FI = result$FI
	Wald1 = n * sum(th1 * (FI %*% th1))
	if(any(nzs)){
		V1 = solve(FI)
		th0 = th1[nzs]
		V0 = V1[nzs,nzs,drop=FALSE]
		FI0 = solve(V0)
		Wald0 = n * sum(th0 * (FI0 %*% th0))
	}else{
		Wald0 = 0
	}
	Wald1 - Wald0
}

# model selection
min_info_selection = function(X, h, Mar, max_iter=50, draw=TRUE, return_detail=FALSE, method="AIC", strct=NULL, ...){
	my_h = function(Z, X_aux, nzs, ...){
		H = h(Z, X_aux, ...)
		H[, nzs, drop=FALSE]  # nzs specifies non-zero parameters
	}
	if(is.matrix(X)) X = as.data.frame(X)
	result = min_info_mle(X, h, Mar, max_iter=max_iter, draw=draw, return_detail=TRUE, ...)
	Kmax = length(result$th)
	if(!is.null(strct)){
		X_aux = Mar(X)
		X_str = strct(X_aux)  # structure of canonical statistic (group and/or hierachy)
		G_s = X_str$G_s  # a list specifying groups of parameter components
		higher = X_str$higher  # a boolean matrix specifying higher interactions
	}else{
		G_s = as.list(1:Kmax)  # default
		higher = matrix(FALSE, Kmax, Kmax)
	}
	Gmax = length(G_s)  # number of groups
	nzs = rep(TRUE, Kmax)  # non-zero parameters
	nzG = rep(TRUE, Gmax)  # non-zero groups
	n = dim(result$Z_s)[2]
	while(any(nzG)){
		th = result$th
		th_se = result$th_se
		FI = result$FI  # Fisher information matrix
		V = solve(FI)
		G = sum(nzG)  # number of non-zero groups
		K = sum(nzs)  # number of non-zero parameters
		if(!is.null(strct)) cat("G=", G, "\n")  # number of non-zero groups
		cat("K=", K, "\n")  # number of non-zero parameters
		AIC_s = numeric(G)
		cnt = 0
		for(g in 1:Gmax){
			if(!nzG[g]) next
			cnt = cnt + 1
			group = G_s[[g]]
			nzs_v = rep(FALSE, Kmax)
			nzs_v[group] = TRUE   # new zero parameters
			th_large = rep(0, Kmax)
			th_large[nzs] = th
			th_v = th_large[nzs_v]   # difference for Wald test
			V_large = matrix(0, Kmax, Kmax)
			V_large[nzs, nzs] = V
			V_v = V_large[nzs_v, nzs_v, drop=FALSE]
			FI_v = solve(V_v)   # partial Fisher information for Wald test
			if(method=="AIC") d_n = 2*length(group)
			if(method=="BIC") d_n = log(nrow(X))*length(group)
			AIC_s[cnt] = n*sum(th_v * c(FI_v%*%th_v)) - d_n  # based on Wald's test statistic (fixed 2022/04/24)
			if(any(nzG[higher[g,]])) AIC_s[cnt] = Inf  # do not remove if higher interaction exists
		}
		cat("AIC diff:", AIC_s, "\n")
		if(min(AIC_s) > 0) break
		w = which.min(AIC_s)[1]
		cnt = 0
		for(g in 1:Gmax){
			if(!nzG[g]) next
			cnt = cnt + 1
			if(cnt == w){
				nzG[g] = FALSE
				nzs[G_s[[g]]] = FALSE
			}
		}
		show(nzs)
		result = min_info_mle(X, my_h, Mar, max_iter=max_iter, draw=draw, return_detail=TRUE, nzs=nzs,...)
	}
	
	if(return_detail){
		result$nzs = nzs
		return(result)
	}
  	list(th=th, th_se=th_se, nzs=nzs)
}

# example: 2-dim interaction for quantitative/qualitative data
h_2dim = function(Z, X_aux){
	if(is.vector(Z)) Z = matrix(Z, 1, length(Z))
	
	n = nrow(Z)
	d = ncol(Z)
	Y = X_aux$Y
	d_s = X_aux$d_s

	H = matrix(0, n, 0)
	d_cum = c(0, cumsum(d_s))
	for(i in 1:(d-1)){
		for(j in (i+1):d){
			for(k in d_cum[i]+(1:d_s[i])){
				for(l in d_cum[j]+(1:d_s[j])){
					H = cbind(H, Y[Z[,i],k]*Y[Z[,j],l]) # 2-dim interaction
				}
			}
		}
	}
	H
}

Mar_2dim = function(X){
	n = nrow(X)
	d = ncol(X)
	Y = matrix(0, n, 0)
	d_s = numeric(d)
	for(i in 1:d){
		if(class(X[,i]) == "angle"){
			Y = cbind(Y, cos(2*pi*X[,i]/360), sin(2*pi*X[,i]/360))
			d_s[i] = 2
			next
		}
		if(is.numeric(X[,i])){
			Y = cbind(Y, X[,i])
			d_s[i] = 1
		}
		if(is.factor(X[,i])){
			m = length(levels(X[,i]))
			X_num = as.numeric(X[,i])
			for(j in 2:m){
				Y = cbind(Y, as.numeric(X_num == j))  # dummy variable
			}
			d_s[i] = m-1
		}
	}
	X_aux = list(Y=Y, d_s=d_s)
	return(X_aux)
}

strct_2dim = function(X_aux){
	d_s = X_aux$d_s
	d = length(d_s)
	G_s = list()
	d_cum = c(0, cumsum(d_s))
	g_cnt = 0
	th_cnt = 0
	for(i in 1:(d-1)){
		for(j in (i+1):d){
			g_cnt = g_cnt + 1
			group = c()
			for(k in d_cum[i]+(1:d_s[i])){
				for(l in d_cum[j]+(1:d_s[j])){
					th_cnt = th_cnt + 1
					group = c(group, th_cnt)
				}
			}
			G_s[[g_cnt]] = group
		}
	}
	if(g_cnt != d*(d-1)/2) stop("bug.")
	higher = matrix(FALSE, g_cnt, g_cnt)
	list(G_s = G_s, higher = higher)
}

# example: 3-dim interaction
h_3dim = function(Z, X_aux){
	if(is.vector(Z)) Z = matrix(Z, 1, length(Z))
	
	n = nrow(Z)
	d = ncol(Z)
	Y = X_aux$Y
	d_s = X_aux$d_s

	H = matrix(0, n, 0)
	d_cum = c(0, cumsum(d_s))
	for(i in 1:(d-1)){
		for(j in (i+1):d){
			for(k in d_cum[i]+(1:d_s[i])){
				for(l in d_cum[j]+(1:d_s[j])){
					H = cbind(H, Y[Z[,i],k]*Y[Z[,j],l]) # 2-dim interaction
				}
			}
		}
	}
	for(i in 1:(d-2)){
		for(j in (i+1):(d-1)){
			for(k in (j+1):d){
				for(a in d_cum[i]+(1:d_s[i])){
					for(b in d_cum[j]+(1:d_s[j])){
						for(c in d_cum[k]+(1:d_s[k])){
							H = cbind(H, Y[Z[,i],a]*Y[Z[,j],b]*Y[Z[,k],c]) # 3-dim interaction
						}
					}
				}
			}
		}
	}
	H
}

Mar_3dim = function(X){
	Mar_2dim(X)
}

strct_3dim = function(X_aux){
	d_s = X_aux$d_s
	d = length(d_s)
	G_s = list()
	d_cum = c(0, cumsum(d_s))
	g_cnt = 0
	th_cnt = 0
	g_num = matrix(0, d, d)  # group number
	for(i in 1:(d-1)){
		for(j in (i+1):d){
			g_cnt = g_cnt + 1
			group = c()
			for(k in d_cum[i]+(1:d_s[i])){
				for(l in d_cum[j]+(1:d_s[j])){
					th_cnt = th_cnt + 1
					group = c(group, th_cnt)
				}
			}
			G_s[[g_cnt]] = group
			g_num[i,j] = g_num[j,i] = g_cnt
		}
	}
	for(i in 1:(d-2)){
		for(j in (i+1):(d-1)){
			for(k in (j+1):d){
				g_cnt = g_cnt + 1
				group = c()
				for(a in d_cum[i]+(1:d_s[i])){
					for(b in d_cum[j]+(1:d_s[j])){
						for(c in d_cum[k]+(1:d_s[k])){
							th_cnt = th_cnt + 1
							group = c(group, th_cnt)
						}
					}
				}
				G_s[[g_cnt]] = group
			}
		}
	}
	if(g_cnt != d*(d-1)/2 + d*(d-1)*(d-2)/6) stop("bug.")
	higher = matrix(FALSE, g_cnt, g_cnt)
	g_cnt = d*(d-1)/2
	for(i in 1:(d-2)){
		for(j in (i+1):(d-1)){
			for(k in (j+1):d){
				g_cnt = g_cnt + 1
				higher[g_num[i,j], g_cnt] = TRUE
				higher[g_num[i,k], g_cnt] = TRUE
				higher[g_num[j,k], g_cnt] = TRUE
			}
		}
	}
	list(G_s = G_s, higher = higher)
}

# EOF
