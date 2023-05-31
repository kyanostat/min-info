# possible correlation coefficients for Poisson marginals

lam1 = 1
lam2_s = 2^((-2):2)
lam2_n = length(lam2_s)

range_s = matrix(0, lam2_n, 2)

SK_algo = function(A, r1, r2, eps=1e-4, max_iter=1000){  # Sinkhorn-Knopp algorithm
	n1 = nrow(A)
	d1 = rep(1, n1)
	err = Inf
	iter = 0
	while(err > eps && iter < max_iter){
		A = diag(d1) %*% A
		p2 = colSums(A)
		d2 = r2 / p2
		A = A %*% diag(d2)
		p1 = rowSums(A)
		d1 = r1 / p1
		err = max(abs(d1 - 1))
		iter = iter + 1
	}
	list(d1=d1, d2=d2, A=A, conv=ifelse(iter < max_iter, 0, 1))
}

for(lam2_i in 1:lam2_n){
	lam2 = lam2_s[lam2_i]
	cat(lam2, "\n")

	x1_s = 0:20
	x2_s = 0:20
	x1_n = length(x1_s)
	x2_n = length(x2_s)
	r1 = dpois(x1_s, lam1)
	r2 = dpois(x2_s, lam2)
	r1 = r1 / sum(r1)
	r2 = r2 / sum(r2)

#	th_s = seq(-10, 10, len = 201)
	th_s = seq(-9, 9, len = 201)
	th_n = length(th_s)
	result = numeric(th_n)

	for(th_i in 1:th_n){
		th = th_s[th_i]
		cat(th_i, "")
		if(th < 0){
			A = A0 = exp(th*(outer(x1_s, x2_s)))
		}else{
			A = A0 = exp(th*(outer(x1_s, x2_s) - outer(x1_s^2/2, rep(1,x2_n)) - outer(rep(1,x1_n), x2_s^2/2)))
		}

		# Sinkhorn algorithm
		SK_result = SK_algo(A, r1, r2)
		A = SK_result$A
		if(SK_result$conv == 1) warning(paste("not converge", "lam2=", lam2, "th=", th))
		v12 = sum(outer(x1_s, x2_s) * A) - sum(x1_s * r1) * sum(x2_s * r2)
		v1 = sum(x1_s^2 * r1) - sum(x1_s * r1)^2
		v2 = sum(x2_s^2 * r2) - sum(x2_s * r2)^2
		rho = v12 / sqrt(v1*v2)
		result[th_i] = rho
	}
	plot(th_s, result)
	show(range(result))
	range_s[lam2_i,] = range(result)
}
