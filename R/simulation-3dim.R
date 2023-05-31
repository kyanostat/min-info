source("min-info.R")

set.seed(1)

th0 = c(1,0,0,-1)

N = 1e3

### simulation (1 time)

X = matrix(rnorm(N*3), N, 3)
X_pop = min_info_sampling(X, th0, h_3dim, Mar_3dim)

par(mfrow=c(1,2), pty="s")
plot(X_pop[X_pop[,3]>0,1:2], xlim=c(-3,3), ylim=c(-3,3), xlab="x", ylab="y", main="z>0")
plot(X_pop[X_pop[,3]<0,1:2], xlim=c(-3,3), ylim=c(-3,3), xlab="x", ylab="y", main="z<0")
dev.copy2pdf(file="simulation-3dim-1.pdf")

par(mfrow=c(1,1))
n = 1e2
X = X_pop[sample(1:N, n),]  # sample

result = min_info_mle(X, h_3dim, Mar_3dim, eps=0.01)
th = result$th
th_se = result$th_se
show(th)
show(th_se)

dev.copy2pdf(file="simulation-3dim-2.pdf")


### simulation (many times)

date0 = date()

LOOP = 1e3
K = length(th0)
th_s = matrix(0, LOOP, K)
th_se_s = matrix(0, LOOP, K)
th_ple_s = matrix(0, LOOP, K)
CI_in = matrix(0, LOOP, K)

for(Li in 1:LOOP){
	X = matrix(rnorm(N*3), N, 3)
	X_pop = min_info_sampling(X, th0, h_3dim, Mar_3dim)
	X = X_pop[sample(1:N, n), ]
	result = min_info_mle(X, h_3dim, Mar_3dim, eps=0.01)
	th = result$th
	th_se = result$th_se

	th_s[Li,] = th
	th_se_s[Li,] = th_se
	CI_in[Li,] = (abs(th - th0) < qnorm(0.975) * th_se)  # confidence intervals

	result = min_info_Besag(X, h_3dim, Mar_3dim, thin=n)
	th_ple = result$th
	th_ple_s[Li,] = th_ple

	cat(Li, "-th turn\n")
	cat("th", th, "\n")
	cat("th_ple", th_ple, "\n")
}

date1 = date()

th0_s = matrix(th0, LOOP, K, byrow=TRUE)
cat("CLE RMSE", sqrt(colMeans((th_s-th0_s)^2)), "\n")
cat("CLE confidence level", colMeans(CI_in), "\n")
cat("PLE RMSE", sqrt(colMeans((th_ple_s-th0_s)^2)), "\n")
result = list(th_s = th_s, th_se_s = th_se_s, th_ple_s = th_ple_s, th0=th0, N=N, date0=date0, date1=date1)
dump("result", "simulation-3dim-dump.R")
