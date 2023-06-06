source("min-info.R")

#Example 1: continuous data in the unit interval and count data
#Marginal Beta distribution for one varible $X$ in $[0,1]$
#Marginal Poisson distribution for the other variable $Y$ in $\mathbb{N}$
#Joint structure $h$ is given by $h(X,Y)=X/(Y+1)$ with canonical parameter $\theta=100$
########################################################################################
########################################################################################


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

n   = 100
N   = 10000
th0 = 100

# marginal sampling
X_pop = matrix(0, N, 2)
X_pop[, 1] = rbeta(N, 10, 10)
X_pop[, 2] = rpois(N, 3)

# exchange algorithm
X_pop = min_info_sampling(X_pop, th0, h_hybrid, Mar_hybrid) # "population"
X = X_pop[sample(1:N, n), ] # sample

print("+++++++Mixed variables (continuous and discrete variables)+++++++")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
result_ple = min_info_Besag(X, h_hybrid, Mar_hybrid, thin=n, eps=1e-4)
th_ple = result_ple$th
result_cle = min_info_mle(X, h_hybrid, Mar_hybrid, th=th_ple, thin=n, eps=1e-4)
th_cle = result_cle$th
th_se_cle = result_cle$th_se
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("Estimation result")
print(paste0("Target theta : ",th0))
print(paste0("PLE          : ",th_ple))
print(paste0("CLE (sd)     : ",th_cle," (",th_se_cle,")"))
plot(X,xlab=expression(italic(X)),ylab=expression(italic(Y)))
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


########################################################################################
########################################################################################


#Example 2: Penguins data (Mixed graphical model with three-dimensional interactions)
#
#We use the palmerpenguins dataset from Horst, Hill and Gorman (2020).
#
#- The dataset contains several qualitative and quantitative variables
#- We use the variables of (1: continuous) Bill length, (2: continuous) Bill depth, (3: continuous) flipper length, (4: continuous) body mass, (5: discrete) sex
#- We test $h(x)=\sum \theta_{ij}x_i x_j$ (two-dimensional effect)
#and $h(x)=\sum \theta_{ij}x_i x_j + \sum \theta_{ijk}x_i x_j x_k$ (three-dimensional effect) after one hot encoding for sex
#
#
###References
#- K. Gorman, T. Williams, W. Fraser (2014). Ecological sexual dimorphism and environmental variability within a community of Antarctic penguins (genus Pygoscelis). PLoS ONE 9(3):e90081. https://doi.org/10.1371/journal.pone.0090081
#
#- A. Horst, A. Hill, K. Gorman (2020). palmerpenguins: Palmer Archipelago (Antarctica) penguin data. R package version 0.1.0. https://allisonhorst.github.io/palmerpenguins/. doi: 10.5281/zenodo.3960218.


library(palmerpenguins)
data(penguins)

X_penguins_pre = as.data.frame(penguins)
w_na = which(is.na(X_penguins_pre$sex))  # all missing cases are removed
X_penguins = X_penguins_pre[-w_na,]

X_penguins = X_penguins[X_penguins[,1]=="Adelie", c(3:7)]  # 5 variables
w_num = unlist(lapply(X_penguins, is.numeric))
X_penguins[,w_num] = scale(X_penguins[,w_num])

print("++++++Penguins data (Mixed graphical model)++++++++++++++++++++++")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

#h_2dim and Mar_2dim are defined in min-info.R
result_penguins_ple = min_info_Besag(X_penguins, h_2dim, Mar_2dim, thin=n, eps=1e-4)
th_penguins_ple = result_penguins_ple$th
result_penguins_cle = min_info_mle(X_penguins, h_2dim, Mar_2dim, th=th_penguins_ple, eps=0.01, thin=nrow(X_penguins))
th_penguins_cle = result_penguins_cle$th
th_penguins__cle_se = result_penguins_cle$th_se
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("Estimation result")
print(paste0("Interaction, PLE, CLE"))
print("1-2,    1-3,   1-4,   1-5,   2-3,   2-4,   2-5,   3-4,   3-5,   4-5")
print(round(th_penguins_ple,3))
print(round(th_penguins_cle,3))
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


########################################################################################
########################################################################################


#Example 3: Earthquake catalog data (orthogonal frame data and real valued data)
#
#- We analyze the dependence between mechanism solution and depth in the earthquake catalog
#- We use the 158 earthquake data that occurred in Japan during the period from January 1st, 2021 to December 8th, 2021 
#- The data we use are provided by the Japan Meteorological Agency and contain
#the information about (0) time (1) latitude (2) longitude (3) depth (4) Magnitude (5) JMA Magnitude (6) strike 1 (7) dip 1 (8) slip 1 (9) strike 2 (10) dip 2 (11) slip 2 of the events
#- We use the columns of (3) depth ,(6)-(11) mechanism solution
#- The mechanism solution of an earthquake event is given in the form of two orthogonal axes in the 3D space (Arnold and Jupp, 2013): $\{v_P,v_T\in \mathbb{R}^{3}: (v_P)^{\top}v_T=0, \|v_P\|=1, \|v_T\|=1\}$
#  - To analyze the dependence between mechanism solution and depth, we consider the minimum information dependence model with the form
#$\exp(\mathrm{tr}(A(z v_{P}v_{P}^{\top})) + \mathrm{tr}(B(z v_{T}v_{T}^{\top})) +\text{terms related to marginals})$, where $z$ is the depth.
#
###References
#- JAPAN METEOROLOGICAL AGENCY (2022). The seismological bulletin of Japan. https://www.data.jma.go.jp/svd/eqev/data/bulletin/index_e.html.
#- R. Arnold and P. Jupp (2013). Statistics of orthogonal axial frames. Biometrika 100 571–586.


print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("++++++Earthquake data (manifold and continuous data)+++++++++++++")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


X_eq = X_eq_ori = read.csv("Mechsol_format.csv")

# preprocess 1
lat1 = as.numeric(substr(X_eq[,2], 1, 2))
lat2 = as.numeric(substr(X_eq[,2], 4, 7))
lat = lat1 * 60 + lat2
lon1 = as.numeric(substr(X_eq[,3], 1, 3))
lon2 = as.numeric(substr(X_eq[,3], 5, 8))
lon = lon1 * 60 + lon2
dep = as.numeric(sub("km", "", X_eq[,4]))
X_eq[,2] = lat
X_eq[,3] = lon
X_eq[,4] = dep

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

vec1 = angle2vec(X_eq$strike1, X_eq$dip1, X_eq$slip1)
vec2 = angle2vec(X_eq$strike2, X_eq$dip2, X_eq$slip2)
vec2 = vec2 - rowSums(vec2*vec1) / rowSums(vec1^2) * vec1  # orghotonalize
vec2 = vec2 / sqrt(rowSums(vec2^2))  # normalize

# vector product
vec3 = cbind(vec1[,2]*vec2[,3] - vec1[,3]*vec2[,2],
             vec1[,3]*vec2[,1] - vec1[,1]*vec2[,3],
             vec1[,1]*vec2[,2] - vec1[,2]*vec2[,1])  # orthogonalize
vec3 = vec3 / sqrt(rowSums(vec3^2))  # normalize


### h function and Mar function

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


# main

X_eq_mechdep = list()
# compressional, tensional and null axes (see Arnold and Jupp (2013), Sec 10)
X_eq_mechdep[[1]] = cbind((vec1-vec2)/sqrt(2), (vec1+vec2)/sqrt(2), vec3)
X_eq_mechdep[[2]] = scale(dep)


result_mechdep = min_info_mle(X_eq_mechdep, h_frame, Mar_frame, thin=nrow(X_eq_mechdep[[1]]), eps=0.001, return_detail=TRUE)

print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("Estimation result")


th_mechdep = result_mechdep$th
th_mechdep_se = result_mechdep$th_se
K = length(th_mechdep)

# interpretation
print("Interpretable matrix representation")
A = matrix(0, 3, 3)  # compressional part
A[1,1] = th_mechdep[1]
A[1,2] = th_mechdep[2]
A[1,3] = th_mechdep[3]
A[2,2] = th_mechdep[4]
A[2,3] = th_mechdep[5]
A = (A + t(A))/2
A_lam = eigen(A)$val
A_adjust = A - A_lam[2]*diag(1,3,3)


show(eigen(A_adjust))


print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
