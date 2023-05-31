source("min-info.R")
X = X_ori = read.csv("Mechsol_format.csv")

set.seed(1)

# preprocess 1
lat1 = as.numeric(substr(X[,2], 1, 2))
lat2 = as.numeric(substr(X[,2], 4, 7))
lat = lat1 * 60 + lat2
lon1 = as.numeric(substr(X[,3], 1, 3))
lon2 = as.numeric(substr(X[,3], 5, 8))
lon = lon1 * 60 + lon2
dep = as.numeric(sub("km", "", X[,4]))
X[,2] = lat
X[,3] = lon
X[,4] = dep

# preprocess 2
# Reference:
#   https://www.data.jma.go.jp/svd/eqev/data/mech/kaisetu/mechkaisetu2.html
angle2vec = function(strike, dip, slip){
	phi = strike * (2*pi)/360
	del = dip * (2*pi)/360
	lam = slip * (2*pi)/360
	v1 = cos(lam)*cos(phi) + sin(lam)*cos(del)*sin(phi)
	v2 = cos(lam)*sin(phi) - sin(lam)*cos(del)*cos(phi)
	v3 = -sin(lam)*sin(del)
	v = cbind(v1, v2, v3)
	v
}
vec1 = angle2vec(X$strike1, X$dip1, X$slip1)
vec2 = angle2vec(X$strike2, X$dip2, X$slip2)
vec2 = vec2 - rowSums(vec2*vec1) / rowSums(vec1^2) * vec1  # orghotonalize
vec2 = vec2 / sqrt(rowSums(vec2^2))  # normalize

# vector product
vec3 = cbind(vec1[,2]*vec2[,3] - vec1[,3]*vec2[,2],
			vec1[,3]*vec2[,1] - vec1[,1]*vec2[,3],
			vec1[,1]*vec2[,2] - vec1[,2]*vec2[,1])  # orthogonalize
vec3 = vec3 / sqrt(rowSums(vec3^2))  # normalize


### h function and Mar function

# Bingham type (not used)

h_Bingham = function(Z, X_aux){
	if(is.vector(Z)) Z = matrix(Z, 1, length(Z))
	n = nrow(Z)
	d = ncol(Z)
	Y = X_aux
	H = matrix(0, n, 0)
	Y_B = 5  # dimension of Bingham part
	for(a in 1:(d-1)){
		for(Y_i in 1:Y_B){
			H = cbind(H, Y[Z[,1],Y_i]*Y[Z[,1+a],Y_B+a])
		}
	}
	H
}

Mar_Bingham = function(X){
	# X is assumed to be a list consisting of
	#   X[[1]]: n * 3 matrix (axial data)
	#   X[[i]]: numeric for i >= 2
	d = length(X)
	n = nrow(X[[1]])

	Y_B = 5  # dimension of Bingham part	
	Y = matrix(0, n, Y_B+(d-1))
	idx = 0
	for(j in 1:3){
		for(k in 1:3){
			if(j > k) next
			if(j == k && k == 3) next
			idx = idx + 1
			Y[,idx] = X[[1]][,j] * X[[1]][,k]
		}
	}
	for(a in 1:(d-1)){
		Y[,Y_B+a] = X[[1+a]]
	}
	Y
}

# frame-Bingham type

h_frame = function(Z, X_aux){
	if(is.vector(Z)) Z = matrix(Z, 1, length(Z))
	n = nrow(Z)
	d = ncol(Z)
	Y = X_aux
	H = matrix(0, n, 0)
	Y_f = 10  # dimension of frame part (see Arnold and Jupp (2013), Table 1)
	for(a in 1:(d-1)){
		for(Y_i in 1:Y_f){
			H = cbind(H, Y[Z[,1],Y_i]*Y[Z[,1+a],Y_f+a])
		}
	}
	H
}

Mar_frame = function(X){
	# X is assumed to be a list consisting of
	#   X[[1]]: frame n * (3*3) matrices
	#   X[[i]]: numeric for i >= 2
	d = length(X)
	n = nrow(X[[1]])
	
	Y_f = 10  # dimension of frame part (see Arnold and Jupp (2013), Table 1)
	Y = matrix(0, n, Y_f+(d-1))
	idx = 0
	for(i in 1:2){
		for(j in 1:3){
			for(k in 1:3){
				if(j > k) next
				if(j == k && k == 3) next
				idx = idx + 1
				Y[,idx] = X[[1]][,(i-1)*3+j] * X[[1]][,(i-1)*3+k]
			}
		}
	}
	for(a in 1:(d-1)){
		Y[,Y_f+a] = X[[a+1]]
	}
	Y
}

strct_frame = function(X_aux){
	Y = X_aux
	Y_f = 10  # dimension of frame part
	m = ncol(Y) - Y_f  # number of numeric variables
	G_s = list()
	for(a in 1:m) G_s[[a]] = (a-1)*Y_f + (1:Y_f)
	higher = matrix(FALSE, m, m)
	list(G_s = G_s, higher = higher)
}


# main

X = list()
# compressional, tensional and null axes (see Arnold and Jupp (2013), Sec 10)
X[[1]] = cbind((vec1-vec2)/sqrt(2), (vec1+vec2)/sqrt(2), vec3)
X[[2]] = scale(dep)


#w = c(dep <= 100)
#X[[1]] = X[[1]][w,]
#X[[2]] = X[[2]][w]

show(lapply(X, head))

#result = min_info_mle(X, h_frame, Mar_frame, thin=nrow(X[[1]]), eps=0.0001, return_detail=TRUE, max_iter=100)
result = min_info_mle(X, h_frame, Mar_frame, thin=nrow(X[[1]]), eps=0.001, return_detail=TRUE)

th = result$th
th_se = result$th_se
K = length(th)
#plot(1:K, th, ylab="theta", ylim=range(c(th-2*th_se, th+2*th_se)))
#segments(-1, 0, K+1, 0)
#for(i in 1:K){
#	arrows(i, th[i] - 2*th_se[i], i, th[i] + 2*th_se[i], angle=90, code=3, length=0.1)
#}
plot(1:K, th/th_se)


# interpretation
A = matrix(0, 3, 3)  # compressional part
A[1,1] = th[1]
A[1,2] = th[2]
A[1,3] = th[3]
A[2,2] = th[4]
A[2,3] = th[5]
A = (A + t(A))/2
A_lam = eigen(A)$val
A_adjust = A - A_lam[2]*diag(1,3,3)

B = matrix(0, 3, 3)  # tensional part
B[1,1] = th[6]
B[1,2] = th[7]
B[1,3] = th[8]
B[2,2] = th[9]
B[2,3] = th[10]
B = (B + t(B))/2
B_lam = eigen(B)$val
B_adjust = B - B_lam[2]*diag(1,3,3)

show(eigen(A_adjust))
show(eigen(B_adjust))

##### eps=0.01
#$values
#[1]  4.132908e+00 -4.440892e-16 -1.090449e+00
#
#$vectors
#            [,1]       [,2]       [,3]
#[1,]  0.06006996 0.88259858  0.4662739
#[2,] -0.26790209 0.46423445 -0.8442244
#[3,]  0.96157167 0.07420323 -0.2643365
#
#eigen() decomposition
#$values
#[1]  1.045246e+00 -4.440892e-16 -1.225560e+00
#
#$vectors
#           [,1]      [,2]       [,3]
#[1,]  0.1244186 0.9424927  0.3102055
#[2,] -0.6150571 0.3185828 -0.7212557
#[3,]  0.7786044 0.1010565 -0.6193245
#
#
##### eps=0.001
#> eigen(A_adjust)
#eigen() decomposition
#$values
#[1]  3.915036  0.000000 -1.018220
#
#$vectors
#           [,1]       [,2]       [,3]
#[1,]  0.0708846 0.90282040  0.4241352
#[2,] -0.2734219 0.42650490 -0.8621682
#[3,]  0.9592788 0.05485339 -0.2770835
#
#> eigen(B_adjust)
#eigen() decomposition
#$values
#[1]  1.002095  0.000000 -1.178124
#
#$vectors
#           [,1]       [,2]       [,3]
#[1,]  0.1624190 0.92787719  0.3356546
#[2,] -0.5906066 0.36392419 -0.7202381
#[3,]  0.7904454 0.08125945 -0.6071186


# Wald test
WT = numeric(3)
#WT[1] = min_info_Wald(X, h_frame, Mar_frame, rep(FALSE, 10), thin=nrow(X[[1]]))
#WT[2] = min_info_Wald(X, h_frame, Mar_frame, c(rep(FALSE, 5), rep(TRUE,5)), thin=nrow(X[[1]]))
#WT[3] = min_info_Wald(X, h_frame, Mar_frame, c(rep(TRUE, 5), rep(FALSE,5)), thin=nrow(X[[1]]))
WT[1] = min_info_Wald(result, rep(FALSE, 10), thin=nrow(X[[1]]))
WT[2] = min_info_Wald(result, c(rep(FALSE, 5), rep(TRUE,5)), thin=nrow(X[[1]]))
WT[3] = min_info_Wald(result, c(rep(TRUE, 5), rep(FALSE,5)), thin=nrow(X[[1]]))
df = c(10,5,5)
pval_WT = 1 - pchisq(WT, df)
cat("Wald test:", WT, "\n")
cat("p-value:", pval_WT, "\n")

# likelihood ratio test
LRT = numeric(3)
LRT[1] = min_info_LRT(X, h_frame, Mar_frame, rep(FALSE, 10), thin=nrow(X[[1]]), eps=0.001)
LRT[2] = min_info_LRT(X, h_frame, Mar_frame, c(rep(FALSE, 5), rep(TRUE,5)), thin=nrow(X[[1]]), eps=0.001)
LRT[3] = min_info_LRT(X, h_frame, Mar_frame, c(rep(TRUE, 5), rep(FALSE,5)), thin=nrow(X[[1]]), eps=0.001)
df = c(10,5,5)
pval_LRT = 1 - pchisq(LRT, df)
cat("likelihood ratio test:", LRT, "\n")
cat("p-value:", pval_LRT, "\n")

dump(c("WT", "pval_WT", "LRT", "pval_LRT"), "test_stat.R")

dump(c("A_adjust", "B_adjust"), "AB-estimate.R")
