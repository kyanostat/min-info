# -*- coding: utf-8 -*-
#Import libraries

import numpy as np
import gc
from sklearn.linear_model import LogisticRegression
from scipy import linalg
import random

## MCMC sampler for the minimum information dependence model with canonical statistics h and canonical parameter theta

def min_info_mcmc(X,h,loop,theta,printop=False):

###############################################################################
  #INPUT
  ## X: initial sample list [[multivariate sample 1],[multivariate sample 2],...]
  ## h   : canonical statistics function (from X (d-dim list) to R^K K-dim array)
  ## loop: number of mcmc steps
  ## theta: canonical parameter (np.array)
  ## printop: print option (if True then print the acceptance rate every 1000 steps)

  #OUTPUT
  ## list of Markov chains (list of Xs)
  ## output[0]  yields the input
  ## output[-1] yields the final output
###############################################################################

  if len(X) ==0:
    print("NULL list")
  else:
    n = len(X)
    d = len(X[0])
    output_list = [X]
    acceptancenum = 0
    for Li in range(loop):
      i  = random.sample(range(d),1)[0]
      s,t = random.sample(range(n),2)[0:2]
      X_iswap_s = output_list[Li][s].copy()
      X_iswap_t = output_list[Li][t].copy()

      X_iswap_s[i] = output_list[Li][t][i]
      X_iswap_t[i] = output_list[Li][s][i]

      h_diff = h(output_list[Li][s])+h(output_list[Li][t])-h(X_iswap_s)-h(X_iswap_t)
      log_rho_tilde  = np.min([np.log(1),-theta.dot(h_diff)])
      u  = np.random.uniform(low=0,high=1,size=1)
      if np.log(u) < log_rho_tilde:
        X_new = output_list[Li].copy()
        X_new[s] = X_iswap_s
        X_new[t] = X_iswap_t
        output_list.append(X_new)
        del X_new

        acceptancenum = acceptancenum + 1
      else:
        output_list.append(output_list[Li])
      del X_iswap_s
      del X_iswap_t
      if (Li % 10000 == 0):
        gc.collect()
      if (Li % 1000 == 0) and (printop == True):
        print("Step/Total : "+str(Li)+"/"+str(loop)+" Acceptance rate : "+str(np.round(acceptancenum/(Li+1),2)))
    return output_list

## Besag pseudo-likelihood estimation 

def min_info_Besag(X,h,max_iter,tol,sparseop=False,sparsepen=1.):

###############################################################################
  #INPUT
  ## X        : initial sample list [[multivariate sample 1],[multivariate sample 2],...]
  ## h        : canonical statistics (from X (d-dim list) to R^K)
  ## max_iter : number of iteration
  ## tol      : tolerance value for the optimizer
  ## sparseop : if True, add l1 penalty
  ## sparsepen: penalty strength valid only if sparseop is True. 

  #OUTPUT
  ## estimate (K-dim array)
###############################################################################

  if len(X) ==0:
    print("NULL list")
  else:
    n = len(X)
    d = len(X[0])
    Us = []
    Ys = []
    for i in range(d):
          for t in range(n):
            for s in range(t):

              X_iswap_s = X[s][:]
              X_iswap_t = X[t][:]
              X_iswap_s[i] = X[t][i]
              X_iswap_t[i] = X[s][i]
              Us.append(h(X[s])+h(X[t])-h(X_iswap_s)-h(X_iswap_t))
              if (t==(n-1))&(s == (n-2)):
                Ys.append(0)
              else:
                Ys.append(1)
    if sparseop==False:
      model = LogisticRegression(solver='liblinear',random_state=0,C=1e10,fit_intercept=False,tol=tol,max_iter=max_iter).fit(np.array(Us), Ys)
    else:
      model = LogisticRegression(solver='saga',random_state=0,C=sparsepen,penalty="l1",fit_intercept=False,tol=tol,max_iter=max_iter).fit(np.array(Us), Ys)
    return model.coef_[0]

## Bootstrap for Besag pseudo-likelihood estimation

def min_info_Besag_bootstrap(X,h,tol,max_iter,B,sparseop=False,sparsepen=1.):

###############################################################################
  #INPUT
  ## X       : initial sample list [[multivariate sample 1],[multivariate sample 2],...]
  ## h       : canonical statistics (from X (d-dim list) to R^K)
  ## max_iter: size of iteration
  ## tol     : tolerance value for the optimizer
  ## B       : size of bootstrap samples 
  ## sparseop : if True, add l1 penalty
  ## sparsepen: penalty strength valid only if sparseop is True. 

  #OUTPUT
  ## bootstrap samples of estimates
  ## mean of bootstraped estimates
  ## covariance of bootstraped estimates
###############################################################################

  n = len(X)
  bootstrap_Besag=[]
  for b in range(B):
    boot_sample=[X[i] for i in random.choices(range(n),k=n)]
    bootstrap_Besag.append(min_info_Besag(X=boot_sample,h=h,tol=tol,max_iter=max_iter,sparseop=sparseop,sparsepen=sparsepen))
  boot_mean = sum(bootstrap_Besag)/B
  boot_cov_1   = sum([np.tensordot(bootstrap_Besag[l],bootstrap_Besag[l],axes=0) for l in range(B)])/B
  boot_cov_2   = np.tensordot(boot_mean,boot_mean,axes=0)
  boot_cov     = boot_cov_1-boot_cov_2
  return bootstrap_Besag, boot_mean, boot_cov

## Calculation of score vector & Hessian of Besag pseudo likelihood
## Calculation of asymptotic covariance of Besag pseudo likelihood estimates

def min_info_Besag_score(x1,x2,h,theta):
###############################################################################
  #INPUT
  ## x1      : multivariate sample
  ## x2      : multivariate sample
  ## h       : canonical statistics (from X (d-dim list) to R^K)
  ## theta   : dependence parameter 

  #OUTPUT
  ## score of Pseudo likelihood

  d = len(x1)
  K = theta.shape[0]
  score = 0*h(x1)

  for i in range(d):
    x_iswap_1 = x1[:]
    x_iswap_2 = x2[:]
    x_iswap_1[i] = x2[i]
    x_iswap_2[i] = x1[i]
    
    ui = h(x1)+h(x2)-h(x_iswap_1)-h(x_iswap_2)
    score = score + (-1)*ui/(1+np.exp(theta.dot(ui)))

  return score

def min_info_Besag_Hessian(x1,x2,h,theta):
###############################################################################
  #INPUT
  ## x1      : multivariate sample
  ## x2      : multivariate sample
  ## h       : canonical statistics (from X (d-dim list) to R^K)
  ## theta   : dependence parameter 

  #OUTPUT
  ## Hessian of Pseudo likelihood

  d = len(x1)
  K = theta.shape[0]
  H = np.zeros(K*K).reshape((K,K))

  for i in range(d):
    x_iswap_1 = x1[:]
    x_iswap_2 = x2[:]
    x_iswap_1[i] = x2[i]
    x_iswap_2[i] = x1[i]
    
    ui = h(x1)+h(x2)-h(x_iswap_1)-h(x_iswap_2)
    Mi = np.tensordot(ui,ui,axes=0)
    H = H + Mi/(1+np.exp(theta.dot(ui)))**2

  return H

def min_info_Besag_asymptoticCovariance(X,h,tol,max_iter):

###############################################################################
  #INPUT
  ## X       : initial sample list [[multivariate sample 1],[multivariate sample 2],...]
  ## h       : canonical statistics (from X (d-dim list) to R^K)
  ## max_iter: size of iteration
  ## tol     : tolerance value for the optimizer

  #OUTPUT
  ## Asymptotic Covariance of PLE
###############################################################################

  n = len(X)
  theta_PLE=min_info_Besag(X=X,h=h,tol=tol,max_iter=max_iter,sparseop=False)
  d = theta_PLE.shape[0]
  I = np.zeros(d*d).reshape((d,d))
  J = np.zeros(d*d).reshape((d,d))

  for s in range(n):
    conditional_score = np.zeros(d)
    for t in range(n):
      conditional_score = conditional_score + (1/n)* min_info_Besag_score(x1=X[s],x2=X[t],h=h,theta=theta_PLE)
      if (t<s):
        B=min_info_Besag_Hessian(x1=X[s],x2=X[t],h=h,theta=theta_PLE)
        J=J+B/(n*(n-1)/2)
    I = I + np.tensordot(conditional_score ,conditional_score ,axes=0) / n
  
  asymptoticCov = (linalg.inv(J) @ I @ linalg.inv(J)) * (4/n)

  return asymptoticCov

## conditional maximum likelihood estimation 

def min_info_CLE(X,h,L,burnin,thin,max_L,max_iter,tol,detailop=False):

###############################################################################
  #INPUT	
  ## X       : initial sample list [[multivariate sample 1],[multivariate sample 2],...]
  ## h       : canonical statistics (from X (d-dim list) to R^K array)
	## L       : MCMC iteration (length of chain)
	## burnin  : burnin of MCMC
	## thin    : thinning
  ## max_iter: number of iteration
  ## tol     : tolerance value for the optimizer
  ##detailop : if True, print (iteration,L,current estimate, current residual for optimization) at each step

	#OUTPUT
	## list of estimates (by iteration)
	## list of residuals for optimization (by iteration)
	## list of L (by iteration)
	## list of asymptotic covariances (by iteration)
###############################################################################

	if len(X) ==0:
		print("NULL list")
	else:
		n = len(X)
		d = len(X[0])
		init_theta = min_info_Besag(X=X,h=h,tol=10**(-2),max_iter=1000)
		init_res = 10**3
		Li         = 0
		theta_list = [init_theta]
		res_list = [init_res]
		L_list   = [L]
		AsymptoticCov_list = [0]
		current_theta = init_theta #array
		current_res = init_res #real number

		while ((current_res>tol) and (Li<max_iter) and (L<max_L)):
			Li = Li + 1
			current_theta = theta_list[Li-1]
			current_res = res_list[Li-1]
			current_X_MCMC= min_info_mcmc(X=X,h=h,loop=L,theta=current_theta,printop=False)[burnin:L:thin]
			total_MCMC_num=len(current_X_MCMC)
			current_hsum             = sum(list(map(h,X)))
			current_hsum_MCMCsamples = [sum(list(map(h, current_X_MCMC[l]))) for l in range(total_MCMC_num)]
			current_hsum_mean        = sum(current_hsum_MCMCsamples)/total_MCMC_num
			current_hsum_cov_1         = sum([np.tensordot(current_hsum_MCMCsamples[l],current_hsum_MCMCsamples[l],axes=0) for l in range(total_MCMC_num)])/total_MCMC_num
			current_hsum_cov_2         = np.tensordot(current_hsum_mean,current_hsum_mean,axes=0)
			current_hsum_cov           = current_hsum_cov_1-current_hsum_cov_2
			score = current_hsum-current_hsum_mean
			modified_score            = linalg.solve(current_hsum_cov, score)
			theta_list.append(current_theta+modified_score)
			res_list.append(max(abs(score)))
			L_list.append(L)
			AsymptoticCov_list.append( linalg.pinv(current_hsum_cov) )
			##for debug
			if (detailop==True):
				print(Li,L,theta_list[Li],res_list[Li])
			##adaptive scaling of L on the basis of the residual update
			if(current_res>1.1*res_list[Li]):
				L = int(np.floor(1.5*L))
			else:
				L = L
			##restart if estimates are very far from PLE
			if(max(abs(theta_list[Li]) / (1 + abs(theta_list[0]))) > 100):
				theta_list[Li] = theta_list[0]
				print("parameter reseted.")
	      
		theta_list.pop(-1)
		res_list.pop(-1)
		AsymptoticCov_list.pop(-1)
	return theta_list,res_list,L_list,AsymptoticCov_list
