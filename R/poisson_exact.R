# possible correlation coefficients for Poisson marginals
# using the optimal transport

lam1 = 1
lam2_s = 2^((-2):2)
lam2_n = length(lam2_s)

range_s = matrix(0, lam2_n, 2)

exact_OT = function(p_s, q_s){  # exact optimal transport
	m = length(p_s)
	n = length(q_s)
	i = j = 1
	r_mat = matrix(0, n, n)
	p_now = q_now = 0
	while(i <= m && j <= n){
		p = p_s[i]
		q = q_s[j]
		if(p-p_now > q-q_now){
			r_mat[i,j] = q - q_now
			q_now = 0
			p_now = p_now + r_mat[i,j]
			j = j + 1
		}else{
			r_mat[i,j] = p - p_now
			p_now = 0
			q_now = q_now + r_mat[i,j]
			i = i + 1
		}
	}
	r_mat
}

get_cor = function(r_mat){ # find correlation (only for count data)
	d = dim(r_mat)
	p = rowSums(r_mat)
	q = colSums(r_mat)
	x = 0:(d[1]-1)
	y = 0:(d[2]-1)
	mu1 = sum(p * x)
	mu2 = sum(q * y)
	var1 = sum(p * x^2) - mu1^2
	var2 = sum(q * y^2) - mu2^2
	my_cov = sum(r_mat * outer(x, y)) - mu1 * mu2
	my_cov / sqrt(var1 * var2)
}

for(lam2_i in 1:lam2_n){
	lam2 = lam2_s[lam2_i]
	
	# max
	p_s = dpois(0:20, lam1)
	q_s = dpois(0:20, lam2)
	r_mat = exact_OT(p_s, q_s)
	cor_max = get_cor(r_mat)
	
	# min
	p_s = dpois(0:20, lam1)
	q_s = dpois(20:0, lam2)
	r_mat = exact_OT(p_s, q_s)
	r_mat = r_mat[,21:1]
	cor_min = get_cor(r_mat)

	range_s[lam2_i,] = c(cor_min, cor_max)

	cat(lam2, round(range_s[lam2_i,],2), "\n")
}
