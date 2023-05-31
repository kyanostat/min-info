source("min-info.R")
library(palmerpenguins)
data(penguins)

set.seed(1)
#set.seed(20220424)

X_ori = as.data.frame(penguins)
w_na = which(is.na(X_ori$sex))  # all missing cases are removed
X = X_ori[-w_na,]
# X = X[,-2]   # 2列目を除くだけで推定は可能（1列目と2列目に structural zero があるため、全変数だとエラーが出る）
#	X = X[,3:7]
#	X = X[,c(1,3:7)]
# X = X[,c(1,3,5,7)]  # 4 variables
#X = X[X[,1]=="Adelie", c(3:7)]  # 5 variables
X = X[X[,1]=="Adelie", c(3:7)]  # 5 variables
w_num = unlist(lapply(X, is.numeric))
X[,w_num] = scale(X[,w_num])  # 量的データは標準化

result = min_info_mle(X, h_2dim, Mar_2dim, eps=0.01, thin=nrow(X))
th = result$th
th_se = result$th_se

dev.copy2pdf(file="penguin-1.pdf")

par(bg="white", cex=1.3)

# 95% confidence interval
show(th)
K = length(th)
plot(1:K, th, ylab="theta", ylim=range(c(th-2*th_se, th+2*th_se)))
segments(-1, 0, K+1, 0)
for(i in 1:K){
	arrows(i, th[i] - 2*th_se[i], i, th[i] + 2*th_se[i], angle=90, code=3, length=0.1)
}

dev.copy2pdf(file="penguin-2.pdf")

# z-value
show(th / th_se)

result_2dim_full = result


d = ncol(X)
X_aux = Mar_2dim(X)
Y = X_aux$Y
d_s = X_aux$d_s

cat("th:", th, "\n")
cat("th_se:", th_se, "\n")
cnt = 0
for(i in 1:(d-1)){
	for(j in (i+1):d){
		for(k in 1:d_s[i]){
			for(l in 1:d_s[j]){
			  cnt = cnt + 1
				cat("edge",cnt,": ",i,k,";",j,l,"\n")
			}
		}
	}
}

# relation to regression analysis
show(th[c(4,7,9,10)])  # related to sex
result_glm = glm(sex ~ ., family="binomial", data=X)
show(result_glm$coef[-1])
#show(summary(result_lm1)$coef)

show(th[c(1,2,3,4)])  # related to bill_length
result_lm1 = lm(bill_length_mm ~ ., data=X)
show(result_lm1$coef[-1] / summary(result_lm1)$sigma^2)

show(th[c(1,5,6,7)])  # related to bill_depth
result_lm2 = lm(bill_depth_mm ~ ., data=X)
show(result_lm2$coef[-1] / summary(result_lm2)$sigma^2)

show(th[c(2,5,8,9)])  # related to flipper_length
result_lm3 = lm(flipper_length_mm ~ ., data=X)
show(result_lm3$coef[-1] / summary(result_lm3)$sigma^2)

show(th[c(3,6,8,10)])  # related to body_mass
result_lm4 = lm(body_mass_g ~ ., data=X)
show(result_lm4$coef[-1] / summary(result_lm4)$sigma^2)



result_2dim = min_info_selection(X, h_2dim, Mar_2dim, strct=strct_2dim, eps=0.01, thin=30)

# AIC diff: 3.236658 13.0839 6.340738 9.427295 23.85622 16.9435 
# $th
# [1] 0.3500276 1.3402972 0.4575508 1.3815040 0.5978457 2.9286544
# $th_se
# [1] 0.1529590 0.3450995 0.1584299 0.4086772 0.1175727 0.6728806
# $nzs
# [1] FALSE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE FALSE  TRUE
### 1-4, 1-5, 2-4, 2-5, 3-4, 4-5

show(result_2dim$th[c(1,2)])  # related to bill_length
result_lm1 = lm(bill_length_mm ~ body_mass_g + sex, data=X)
show(result_lm1$coef[-1] / summary(result_lm1)$sigma^2)

show(result_2dim$th[c(3,4)])  # related to bill_depth
result_lm2 = lm(bill_depth_mm ~ body_mass_g + sex, data=X)
show(result_lm2$coef[-1] / summary(result_lm2)$sigma^2)

show(result_2dim$th[c(5)])  # related to flipper_length
result_lm3 = lm(flipper_length_mm ~ body_mass_g, data=X)
show(result_lm3$coef[-1] / summary(result_lm3)$sigma^2)

show(result_2dim$th[c(1,3,5,6)])  # related to body_mass
result_lm4 = lm(body_mass_g ~ ., data=X)
show(result_lm4$coef[-1] / summary(result_lm4)$sigma^2)

show(result_2dim$th[c(2,4,6)])  # related to sex
result_glm = glm(sex ~ bill_length_mm + bill_depth_mm + body_mass_g, family="binomial", data=X)
show(result_glm$coef[-1])


result_3dim = min_info_selection(X, h_3dim, Mar_3dim, strct=strct_3dim, eps=0.01, thin=30)

#my_strct_3dim = function(X_aux){
#  a = strct_3dim(X_aux)
#  a$higher = a$higher & FALSE
#  a
#}
#result_3dim_not_hier = min_info_selection(X, h_3dim, Mar_3dim, strct=my_strct_3dim, eps=0.01, thin=30)

show(result_2dim$th[c(3,4)])  # related to bill_depth
result_lm2a = lm(bill_depth_mm ~ body_mass_g * sex, data=X)
show(result_lm2a$coef[-1] / summary(result_lm2a)$sigma^2)


LRT_2dim = min_info_LRT(X, h_3dim, Mar_3dim, nzs=result_2dim$nzs, thin=30)
LRT_3dim = min_info_LRT(X, h_3dim, Mar_3dim, nzs=result_3dim$nzs, thin=30)

#> LRT_2dim
#[1] 7.405898
#> LRT_3dim
#[1] 0.6988411

# figure
library(igraph)
# set.seed(3)  # for random layout
g1 = graph(edges=c(1,4,1,5,2,4,2,5,3,4,4,5), directed=F)
#plot(g1, edge.color=ifelse(th>0, "black", "red"), edge.width=abs(th/th_se))
par(cex=1.5, mar=c(1,1,1,1))
#g1_layout = matrix(c(2,0,2,2,-1,1,1,1,3,1), 5, 2, byrow=TRUE)
g1_layout = matrix(c(1,2,1,0,4,1,2,1,0,1), 5, 2, byrow=TRUE)
#
v1_coef = format(round(result_lm1$coef[-1] / summary(result_lm1)$sigma^2, 2), width=4)
v2_coef = format(round(result_lm2$coef[-1] / summary(result_lm2)$sigma^2, 2), width=4)
v3_coef = format(round(result_lm3$coef[-1] / summary(result_lm3)$sigma^2, 2), width=4)
v4_coef = format(round(result_lm4$coef[-1] / summary(result_lm4)$sigma^2, 2), width=4)
v5_coef = format(round(result_glm$coef[-1], 2), width=4)
#
plot(g1, vertex.size=20, vertex.color="white", edge.width=1.5, edge.color="black", layout=g1_layout)
text(c(-0.1,-0.9,-0.1,-0.9,0.5,-0.5), c(0.55,0.55,-0.45,-0.45,-0.07,-0.07), format(round(result_2dim$th, 2), width=4))
text(c(-0.1,-0.9,-0.1,-0.9,0.5,-0.5), c(0.45,0.45,-0.55,-0.55,-0.17,-0.17), paste("(",format(round(result_2dim$th_se, 2), width=4),")",sep=""))
text_col = "skyblue3" # "grey40"
text(c(-0.32,-0.68), c(0.85,0.85), v1_coef, cex=0.7, col=text_col)
text(c(-0.32,-0.68), c(-0.85,-0.85), v2_coef, cex=0.7, col=text_col)
text(c(0.8), c(-0.05), v3_coef, cex=0.7, col=text_col)
text(c(0,0,0.2,-0.2), c(0.15,-0.15,-0.05,-0.05), v4_coef, cex=0.7, col=text_col)
text(c(-1,-1,-0.8), c(0.15,-0.15,-0.05), v5_coef, cex=0.7, col=text_col)
dev.copy2pdf(file="penguin-3.pdf")

nzs1 = c(result_2dim$nzs, rep(FALSE, 10))
result_model_1 = min_info_LRT(X, h_3dim, Mar_3dim, nzs1, thin=30)
nzs2 = nzs3 = nzs4 = nzs1
nzs2[16] = TRUE
# 123, 124, 125, 134, 135, 145, 234, 235, 245, 345
result_model_2 = min_info_LRT(X, h_3dim, Mar_3dim, nzs2, thin=30)
nzs3[19] = TRUE
result_model_3 = min_info_LRT(X, h_3dim, Mar_3dim, nzs3, thin=30)
nzs4[16] = nzs4[19] = TRUE
result_model_4 = min_info_LRT(X, h_3dim, Mar_3dim, nzs4, thin=30)

dump(c("result_2dim_full", "result_2dim", "result_3dim", "LRT_2dim", "LRT_3dim", "result_model_1", "result_model_2", "result_model_3", "result_model_4"), file="penguin-dump.R")
