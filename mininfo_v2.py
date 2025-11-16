# -*- coding: utf-8 -*-

import numpy as np
import gc
import copy
import matplotlib.pyplot as plt
import matplotlib
from sklearn.linear_model import LogisticRegression
from scipy import linalg
import pandas as pd
import random
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Tuple, Optional
import time

#Definition of Functions We Use
#* min_info_mcmc : MCMC sampler for the model.
#* min_info_Besag_fullbatch : Pseudo-likelihood estimation (PLE) using the full dataset.
#* min_info_Besag_SGD : Pseudo-likelihood estimation via Stochastic Gradient Descent (SGD).
#* min_info_Besag_bootstrap_fullbatch : Empirical bootstrap for pseudo-likelihood estimation using the full dataset.
#* min_info_Besag_bootstrap_SGD_Multicore : Empirical bootstrap for pseudo-likelihood estimation via SGD with parallel computation.
#* min_info_Besag_score / Hessian / asymptoticCovariance_fullbatch : Option to calculate the asymptotic covariance of the pseudo-likelihood estimator (full-batch version).
#* min_info_Besag_score / Hessian / asymptoticCovariance_minibatch : Option to calculate the asymptotic covariance of the pseudo-likelihood estimator via SGD (mini-batch version).
#* min_info_CLE : Conditional maximum likelihood estimation (CLE).

def min_info_mcmc(X, h, loop, theta, printop=False):
###############################################################################
## MCMC sampler for the minimum information dependence model with canonical statistics h and canonical parameter theta
###############################################################################
  #INPUT
  ## X        : initial sample list [[multivariate sample 1],[multivariate sample 2],...]
  ## h        : canonical statistics (from X (d-dim list) to R^K)
  ## max_iter : number of iteration
  ## tol      : tolerance value for the optimizer
  ## sparseop : if True, add l1 penalty
  ## sparsepen: penalty strength valid only if sparseop is True.
###############################################################################
  #OUTPUT
  ## estimate (K-dim array)
###############################################################################
    if len(X) == 0:
        print("NULL list")
        return None

    n = len(X)
    d = len(X[0])
    output_list = [copy.deepcopy(X)]
    acceptancenum = 0

    for Li in range(loop):
        i = random.randint(0, d - 1)
        s, t = random.sample(range(n), 2)

        X_iswap_s = list(output_list[Li][s])
        X_iswap_t = list(output_list[Li][t])

        X_iswap_s[i] = output_list[Li][t][i]
        X_iswap_t[i] = output_list[Li][s][i]

        h_diff = h(output_list[Li][s]) + h(output_list[Li][t]) - h(X_iswap_s) - h(X_iswap_t)
        log_rho_tilde = min(0, -theta.dot(h_diff))
        u = np.random.uniform()

        if np.log(u) < log_rho_tilde:
            X_new = [list(row) for row in output_list[Li]]
            X_new[s] = X_iswap_s
            X_new[t] = X_iswap_t
            output_list.append(X_new)
            acceptancenum += 1
        else:
            output_list.append(output_list[Li])

        if (Li % 10000 == 0):
            gc.collect()

        if (Li % 1000 == 0) and printop:
            rate = acceptancenum / (Li + 1)
            print(f"Step/Total : {Li}/{loop} Acceptance rate : {rate:.2f}")

    return output_list

def min_info_Besag_fullbatch(X,h,max_iter,tol,sparseop=False,sparsepen=1.):
###############################################################################
## Besag pseudo-likelihood estimation (full batch implementation)
###############################################################################
  #INPUT
  ## X        : initial sample list [[multivariate sample 1],[multivariate sample 2],...]
  ## h        : canonical statistics (from X (d-dim list) to R^K)
  ## max_iter : number of iteration
  ## tol      : tolerance value for the optimizer
  ## sparseop : if True, add l1 penalty
  ## sparsepen: penalty strength valid only if sparseop is True.
###############################################################################
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


def min_info_Besag_bootstrap_fullbatch(X,h,tol,max_iter,B,sparseop=False,sparsepen=1.):
###############################################################################
## Bootstrap for Besag pseudo-likelihood estimation (full batch)
###############################################################################
  #INPUT
  ## X       : initial sample list [[multivariate sample 1],[multivariate sample 2],...]
  ## h       : canonical statistics (from X (d-dim list) to R^K)
  ## max_iter: size of iteration
  ## tol     : tolerance value for the optimizer
  ## B       : size of bootstrap samples
  ## sparseop : if True, add l1 penalty
  ## sparsepen: penalty strength valid only if sparseop is True.
###############################################################################
  #OUTPUT
  ## bootstrap samples of estimates
  ## mean of bootstraped estimates
  ## covariance of bootstraped estimates
###############################################################################
  n = len(X)
  bootstrap_Besag=[]
  for b in range(B):
    boot_sample=[X[i] for i in random.choices(range(n),k=n)]
    bootstrap_Besag.append(min_info_Besag_fullbatch(X=boot_sample,h=h,tol=tol,max_iter=max_iter,sparseop=sparseop,sparsepen=sparsepen))
  boot_mean = sum(bootstrap_Besag)/B
  boot_cov_1   = sum([np.tensordot(bootstrap_Besag[l],bootstrap_Besag[l],axes=0) for l in range(B)])/B
  boot_cov_2   = np.tensordot(boot_mean,boot_mean,axes=0)
  boot_cov     = boot_cov_1-boot_cov_2
  return bootstrap_Besag, boot_mean, boot_cov

def min_info_Besag_SGD(
    X,
    h,
    max_iter=10000,
    tol=10**(-5),
    learningrate=0.1,
    stepsize_decay = 0.00,
    init = None,
    batch_size=1,
    *,
    # Optional arguments
    ema_beta=0.9,
    patience=5,
    clip_norm=None,
    verbose_every=0,
    seed=None,
    # --- Ruppert--Polyak--Juditsky option ---
    use_rpj=True,
    avg_start=1000,
    return_avg=True,
    # --- Adam option ---
    use_adam=True,
    beta1=0.9,
    beta2=0.999,
    eps=1e-8
):
###############################################################################
## Besag pseudo-likelihood estimation via SGD
###############################################################################
  #INPUT
  ## X             : initial sample list [[multivariate sample 1],[multivariate sample 2],...]
  ## h             : canonical statistics (from X (d-dim list) to R^K)
  ## max_iter      : size of iteration
  ## tol           : tolerance value for the optimizer
  ## learningrate  : base constant step size.
  ## stepsize_decay: decay rate of step size
  ## init          : initial value for the optimization
  ## batch_size    : size of minibatch 
  ## ema_beta      : float value in [0,1) for exponential moving average (EMA) of gradient norm and parameter change.
  ## patience      : Number (int) of consecutive iterations where both EMA(gradient norm) <= tol and EMA(parameter change) <= tol * step size must hold to stop early.
  ## clip_norm     : float or None. If set, clip gradient g by global L2 norm to improve stability.
  ## verbose_every : Logging frequency. Set 0/False to silence.
  ## seed          : seed for reproducibility.
  ## use_rpj       : True (Ruppert--Polyak--Juditsky averaging)
  ## avg_start     : Burn-in for the Ruppert--Polyak--Juditsky averaging
  ## return_avg    : True (use the averaged value for the final parameter vector )
  ## use_adam      : If True, use Adam update instead of plain SGD (default False).
  ## beta1, beta2  : Adam momentum and RMS coefficients.
  ## eps           : Adam numerical stabilizer.
###############################################################################
  #OUTPUT
  ## coef : np.ndarray, shape (p,). The final parameter vector.
  ## RPJ_count: size of the Ruppert--Polyak--Juditsky averaging.
###############################################################################
# --- Random number generation setup ---
  if seed is not None:
    random.seed(seed)
    np.random.seed(seed)

  # --- Shapes & casting ---
  n = len(X)
  d = len(X[0])
  p = len(h(X[0]))
  bs = batch_size

  # --- initialization ---
  if init is None:
    coef = np.zeros(p, dtype=float)
  else:
    coef = np.asarray(init, dtype=float)

  # --- Adam state ---
  if use_adam:
    m = np.zeros_like(coef)
    v = np.zeros_like(coef)

  # --- Numerically stable logistic sigmoid via tanh ---
  def safe_sigmoid(x: float) -> float:
    return 0.5 * (1.0 + np.tanh(0.5 * x))

  # --- Exponential Moving Averages and early-stop bookkeeping ---
  ema_grad = None
  ema_dtheta = None
  consec_small = 0

  # --- Ruppert--Polyak--Juditsky averaging state ---
  if use_rpj:
    avg_coef = np.zeros_like(coef)
    avg_count = 0

  for step in range(1, int(max_iter) + 1):

    i = random.randrange(0, d)
    s = np.empty(bs, dtype=int)
    t = np.empty(bs, dtype=int)
    
    for k in range(bs):
        s[k] = random.randrange(0, n - 1)
        t[k] = random.randrange(s[k] + 1, n)

    
    Xs = np.asarray([X[idx] for idx in s], dtype=float)  # (bs, d)
    Xt = np.asarray([X[idx] for idx in t], dtype=float)  # (bs, d)
    
    Xs_sw = Xs.copy()
    Xt_sw = Xt.copy()
    Xs_sw[:, i] = Xt[:, i]
    Xt_sw[:, i] = Xs[:, i]
    
    hs    = np.vstack([np.asarray(h(x),    dtype=float) for x in Xs])     # (bs, p)
    ht    = np.vstack([np.asarray(h(x),    dtype=float) for x in Xt])     # (bs, p)
    hs_sw = np.vstack([np.asarray(h(x),    dtype=float) for x in Xs_sw])  # (bs, p)
    ht_sw = np.vstack([np.asarray(h(x),    dtype=float) for x in Xt_sw])  # (bs, p)
    
    Usti = (hs_sw + ht_sw) - (hs + ht)          # (bs, p)
    inner = Usti @ np.asarray(coef, dtype=float)  # (bs,)
    
    x = inner.astype(float, copy=False)
    sig = np.empty_like(x, dtype=float)
    pos = x >= 0
    neg = ~pos
    sig[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    expx = np.exp(x[neg])
    sig[neg] = expx / (1.0 + expx)
    
    g_batch = -(sig[:, None]) * Usti    # (bs, p)
    g = g_batch.mean(axis=0)            # (p,) 

    # --- Optional gradient clipping ---
    if clip_norm is not None and clip_norm > 0:
      gnorm = float(np.linalg.norm(g))
      if gnorm > clip_norm:
          g = g * (clip_norm / (gnorm + 1e-12))

    # --- Parameter update (SGD or Adam) ---
    prev = coef.copy()
    if use_adam:
      # Adam moments
      m = beta1 * m + (1.0 - beta1) * g
      v = beta2 * v + (1.0 - beta2) * (g * g)
      # Bias correction
      m_hat = m / (1.0 - (beta1 ** step))
      v_hat = v / (1.0 - (beta2 ** step))
      # Update (note: original SGD uses "+ learningrate * g")
      coef = coef + learningrate * (m_hat / (np.sqrt(v_hat) + eps))
    else:
      coef = coef + (learningrate/(step**stepsize_decay)) * g

    # === Update Ruppert--Polyak--Juditsky running average ===
    if use_rpj and step >= int(avg_start):
      avg_count += 1
      avg_coef += (coef - avg_coef) / avg_count
    # --- Monitoring ---
    grad_norm = float(np.linalg.norm(g))
    dtheta = float(np.linalg.norm(coef - prev))

    # --- Update EMAs ---
    if ema_grad is None:
      ema_grad = grad_norm
      ema_dtheta = dtheta
    else:
      ema_grad = ema_beta * ema_grad + (1 - ema_beta) * grad_norm
      ema_dtheta = ema_beta * ema_dtheta + (1 - ema_beta) * dtheta

    # --- Logging ---
    if verbose_every and (step % verbose_every == 0 or step == 1):
      print(
          f"[{step:6d}/{max_iter}] "
          f"||grad||={grad_norm:.3e} EMA(||grad||)={ema_grad:.3e} "
          f"||update||={dtheta:.3e} EMA(||update||)={ema_dtheta:.3e} "
          f"coef[:3]={coef[:3]}"
      )

    # --- SGD stopping rule ---
    # Require BOTH EMAs to be small, for 'patience' consecutive iterations:
    #   EMA(||g||) ≤ tol
    #   EMA(||update||) ≤ tol * learningrate
    small_grad = ema_grad <= tol
    small_step = ema_dtheta <= (tol * learningrate)

    if small_grad and small_step:
      consec_small += 1
      if consec_small >= patience:
        if verbose_every:
          print(
                          f"Early stop at step {step}: "
                          f"EMA(||g||)={ema_grad:.3e} ≤ {tol:.3e} and "
                          f"EMA(||update||)={ema_dtheta:.3e} ≤ {tol*learningrate:.3e} "
                          f"for {patience} consecutive steps."
          )
          break
    else:
      consec_small = 0

  # --- Return the Ruppert--Polyak--Juditsky averagingif requested ---
  if use_rpj and return_avg and (avg_count > 0):
    return avg_coef,avg_count
  return coef,np.inf # in this case, the asysmptotic covariance evaluation does not work well

def min_info_Besag_bootstrap_SGD(X,h,B,
    max_iter=10000,
    tol=10**(-5),
    learningrate=0.1,
    stepsize_decay = 0.0,
    init = None,
    batch_size=1,
    *,
    # Optional arguments
    ema_beta=0.9,
    patience=5,
    clip_norm=None,
    verbose_every=0,
    seed=None,
    # --- Ruppert--Polyak--Juditsky option ---
    use_rpj=True,
    avg_start=1000,
    return_avg=True,
    # --- Adam option ---
    use_adam=True,
    beta1=0.9,
    beta2=0.999,
    eps=1e-8
):
###############################################################################
## Bootstrap for Besag pseudo-likelihood estimation with SGD
###############################################################################
  #INPUT
  ## X             : initial sample list [[multivariate sample 1],[multivariate sample 2],...]
  ## h             : canonical statistics (from X (d-dim list) to R^K)
  ## B             : size of bootstrap samples
  ## max_iter      : size of iteration
  ## tol           : tolerance value for the optimizer
  ## learningrate  : base constant step size.
  ## stepsize_decay: decay rate of step size
  ## init          : initial value for the optimization
  ## batch_size    : size of minibatch
  ## ema_beta      : float value in [0,1) for exponential moving average (EMA) of gradient norm and parameter change.
  ## patience      : Number (int) of consecutive iterations where both EMA(gradient norm) <= tol and EMA(parameter change) <= tol * step size must hold to stop early.
  ## clip_norm     : float or None. If set, clip gradient g by global L2 norm to improve stability.
  ## verbose_every : Logging frequency. Set 0/False to silence.
  ## seed          : seed for reproducibility.
  ## use_rpj       : True (Ruppert--Polyak--Juditsky averaging)
  ## avg_start     : Burn-in for the Ruppert--Polyak--Juditsky averaging
  ## return_avg    : True (use the averaged value for the final parameter vector )
  ## use_adam      : If True, use Adam update instead of plain SGD (default False).
  ## beta1, beta2  : Adam momentum and RMS coefficients.
  ## eps           : Adam numerical stabilizer.
###############################################################################
  #OUTPUT
  ## coef : np.ndarray, shape (p,). The final parameter vector.
###############################################################################
  n = len(X)
  bootstrap_Besag=[]
  for b in range(B):
    boot_sample=[X[i] for i in random.choices(range(n),k=n)]
    bootstrap_Besag.append(min_info_Besag_SGD(X=boot_sample,h=h,tol=tol,max_iter=max_iter,
                                               learningrate=learningrate,stepsize_decay = stepsize_decay,
                                              init=init,
                                              batch_size = batch_size,
                                              ema_beta=ema_beta,patience=patience,clip_norm=clip_norm,verbose_every=verbose_every,
                                              seed=seed,use_rpj=use_rpj,avg_start=avg_start,return_avg=return_avg,
                                              use_adam=use_adam,beta1=beta1,beta2=beta2,eps=eps)[0])
  boot_mean = sum(bootstrap_Besag)/B
  boot_cov_1   = sum([np.tensordot(bootstrap_Besag[l],bootstrap_Besag[l],axes=0) for l in range(B)])/B
  boot_cov_2   = np.tensordot(boot_mean,boot_mean,axes=0)
  boot_cov     = boot_cov_1-boot_cov_2
  return bootstrap_Besag, boot_mean, boot_cov

def _one_bootstrap_run(
    indices: np.ndarray,
    X: list,
    h,
    tol: float,
    max_iter: int,
    learningrate: float,
    stepsize_decay: float,
    init,
    batch_size: int,
    ema_beta: float,
    patience: int,
    clip_norm: Optional[float],
    verbose_every: int,
    seed: Optional[int],
    use_rpj: bool,
    avg_start: int,
    return_avg: bool,
    use_adam: bool,
    beta1: float,
    beta2: float,
    eps: float
):
###############################################################################
#####The function required for multicore bootstrap
#####Perform a single bootstrap replicate (executed in a worker process).
###############################################################################
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed % (2**32 - 1))

    boot_sample = [X[i] for i in indices]
    theta = min_info_Besag_SGD(
        X=boot_sample, h=h, tol=tol, max_iter=max_iter,
        learningrate=learningrate,
    stepsize_decay = stepsize_decay,
        init = init,
        batch_size=batch_size,
        ema_beta=ema_beta,
        patience=patience, clip_norm=clip_norm, verbose_every=verbose_every,
        seed=seed, use_rpj=use_rpj, avg_start=avg_start, return_avg=return_avg,
        use_adam=use_adam, beta1=beta1, beta2=beta2, eps=eps
    )[0]
    return theta


def min_info_Besag_bootstrap_SGD_Multicore(
    X,
    h,
    B,
    max_iter=10000,
    tol=10**(-5),
    learningrate=0.1,
    stepsize_decay = 0.00,
    init = None,
    batch_size=1,
    *,
    # Optional arguments
    ema_beta=0.9,
    patience=5,
    clip_norm=None,
    verbose_every=0,
    seed=None,
    # --- Ruppert--Polyak--Juditsky option ---
    use_rpj=True,
    avg_start=1000,
    return_avg=True,
    # --- Adam option ---
    use_adam=True,
    beta1=0.9,
    beta2=0.999,
    eps=1e-8,
    # --- Parallel option ---
    n_jobs: int = -1,
    chunksize: int = 1,
    progress: bool = False
):
###############################################################################
## Multicore bootstrap for Besag pseudo-likelihood estimation with SGD
###############################################################################
  #INPUT
  ## X             : initial sample list [[multivariate sample 1],[multivariate sample 2],...]
  ## h             : canonical statistics (from X (d-dim list) to R^K)
  ## B             : size of bootstrap samples
  ## max_iter      : size of iteration
  ## tol           : tolerance value for the optimizer
  ## learningrate  : base step size
  ## stepsize_decay: decay rate of step size
  ## init          : initial value for the optimization
  ## batch_size    : size of minibatch
  ## ema_beta      : float value in [0,1) for exponential moving average (EMA) of gradient norm and parameter change.
  ## patience      : Number (int) of consecutive iterations where both EMA(gradient norm) <= tol and EMA(parameter change) <= tol * step size must hold to stop early.
  ## clip_norm     : float or None. If set, clip gradient g by global L2 norm to improve stability.
  ## verbose_every : Logging frequency. Set 0/False to silence.
  ## seed          : seed for reproducibility.
  ## use_rpj       : True (Ruppert--Polyak--Juditsky averaging)
  ## avg_start     : Burn-in for the Ruppert--Polyak--Juditsky averaging
  ## return_avg    : True (use the averaged value for the final parameter vector )
  ## use_adam      : If True, use Adam update instead of plain SGD (default False).
  ## beta1, beta2  : Adam momentum and RMS coefficients.
  ## eps           : Adam numerical stabilizer.
  ## n_jobs        : Number of jobs to run in parallel
  ## chunksize     : Number of samples per chunk
  ## progress      : If True, print progress
###############################################################################
  #OUTPUT
  ## coef : np.ndarray, shape (p,). The final parameter vector.
###############################################################################
    n = len(X)
    if B <= 0:
        raise ValueError("B must be >= 1")

    # Auto-detect CPU count if n_jobs == -1
    if n_jobs == -1:
        import os
        n_jobs = max(1, os.cpu_count() or 1)

    # Precompute all bootstrap indices and seeds for reproducibility
    rng = np.random.default_rng(seed)
    all_indices = rng.integers(0, n, size=(B, n))
    worker_seeds = rng.integers(0, 2**32 - 1, size=B) if seed is not None else [None] * B

    # Sequential fallback for debugging
    if n_jobs == 1:
        bootstrap_Besag = []
        for b in range(B):
            if progress and (b % max(1, B // 10) == 0):
                print(f"[{b}/{B}] running…")
            theta = _one_bootstrap_run(
                indices=all_indices[b], X=X, h=h, tol=tol, max_iter=max_iter,
                learningrate=learningrate,
                stepsize_decay = stepsize_decay,
                init = init,
                batch_size=batch_size,
                ema_beta=ema_beta,
                patience=patience, clip_norm=clip_norm, verbose_every=verbose_every,
                seed=int(worker_seeds[b]) if worker_seeds[b] is not None else None,
                use_rpj=use_rpj, avg_start=avg_start, return_avg=return_avg,
                use_adam=use_adam, beta1=beta1, beta2=beta2, eps=eps
            )
            bootstrap_Besag.append(np.asarray(theta))
    else:
        bootstrap_Besag = [None] * B
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            futures = []
            for b in range(B):
                fut = ex.submit(
                    _one_bootstrap_run,
                    all_indices[b], X, h, tol, max_iter, learningrate, stepsize_decay,
                    init,
                    batch_size,
                    ema_beta, patience, clip_norm, verbose_every,
                    int(worker_seeds[b]) if worker_seeds[b] is not None else None,
                    use_rpj, avg_start, return_avg, use_adam, beta1, beta2, eps
                )
                futures.append((b, fut))

            done_cnt = 0
            for b, fut in futures:
                theta = fut.result()
                bootstrap_Besag[b] = np.asarray(theta)
                done_cnt += 1
                if progress and (done_cnt % max(1, B // 10) == 0):
                    print(f"[{done_cnt}/{B}] finished")

    # Compute statistics
    theta_arr = np.stack(bootstrap_Besag, axis=0)
    boot_mean = theta_arr.mean(axis=0)
    boot_cov = np.cov(theta_arr, rowvar=False, bias=True)

    return bootstrap_Besag, boot_mean, boot_cov

def min_info_Besag_score_fullbatch(x1,x2,h,theta):
###############################################################################
## Calculation of score vector (full batch)
###############################################################################
  #INPUT
  ## x1      : multivariate sample
  ## x2      : multivariate sample
  ## h       : canonical statistics (from X (d-dim list) to R^K)
  ## theta   : dependence parameter
###############################################################################
  #OUTPUT
  ## score of Pseudo likelihood (full batch)
###############################################################################

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

def min_info_Besag_Hessian_fullbatch(x1,x2,h,theta):
###############################################################################
## Calculation of Hessian of Besag pseudo likelihood (full batch)
###############################################################################
  #INPUT
  ## x1      : multivariate sample
  ## x2      : multivariate sample
  ## h       : canonical statistics (from X (d-dim list) to R^K)
  ## theta   : dependence parameter
###############################################################################
  #OUTPUT
  ## Hessian of Pseudo likelihood (full batch)
###############################################################################
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

def min_info_Besag_asymptoticCovariance_fullbatch(X,h,theta_PLE):
###############################################################################
## Calculation of full batched asymptotic covariance at theta_PLE
###############################################################################
  #INPUT
  ## X         : initial sample list [[multivariate sample 1],[multivariate sample 2],...]
  ## h         : canonical statistics (from X (d-dim list) to R^K)
  ## theta_PLE : the value of PLE
###############################################################################
  #OUTPUT
  ## Asymptotic Covariance of PLE
###############################################################################

  n = len(X)
  d = theta_PLE.shape[0]
  I = np.zeros(d*d).reshape((d,d))
  J = np.zeros(d*d).reshape((d,d))

  for s in range(n):
    conditional_score = np.zeros(d)
    for t in range(n):
      conditional_score = conditional_score + (1/n)* min_info_Besag_score_fullbatch(x1=X[s],x2=X[t],h=h,theta=theta_PLE)
      if (t<s):
        B=min_info_Besag_Hessian_fullbatch(x1=X[s],x2=X[t],h=h,theta=theta_PLE)
        J=J+B/(n*(n-1)/2)
    I = I + np.tensordot(conditional_score ,conditional_score ,axes=0) / n

  asymptoticCov = (linalg.inv(J) @ I @ linalg.inv(J)) * (4/n)

  return asymptoticCov


def min_info_Besag_asymptoticCovariance_minibatch(
    X, h, theta_PLE, RPJ_count,*,
    num_s=None,
    num_t_per_s=10,
    num_pairs_J=100,
    seed=None,
    replace_s=True
):
###############################################################################
    ## Monte Carlo / minibatch approximation of the asymptotic covariance at theta_PLE
###############################################################################
    #INPUT
    ##   X           : list of samples [[x^(1)], [x^(2)], ...], each x is d-dimensional
    ##   h           : canonical statistics; maps a d-dim vector to R^K
    ##   theta_PLE   : PLE parameter vector (shape (d,))
    ##   RPJ_count    : sample sizes in the Ruppert--Polyak--Juditsky averaging
    ##   num_s       : number of s indices drawn to approximate the average over s
    ##   num_t_per_s : number of t indices drawn for each s to average the conditional score
    ##   num_pairs_J : number of unordered distinct pairs (s, t) to approximate J
    ##   seed        : RNG seed
    ##   replace_s   : whether to sample s with replacement
###############################################################################
    # OUTPUT
    ##   Asymptotic covariance (Monte Carlo estimate)
###############################################################################

    n = len(X)
    d = theta_PLE.shape[0]
    I = np.zeros((d, d), dtype=float)
    J = np.zeros((d, d), dtype=float)

    # RNG setup
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    rng = np.random.default_rng(seed)

    # Default sample sizes
    if num_s is None:
        num_s = n

    # --- Minibatch approximation of I ---
    # I is the average over s of the outer product of the conditional score (averaged over t).
    # We sample s (with or without replacement) and, for each s, sample t's with replacement.
    if replace_s:
        s_indices = rng.integers(0, n, size=num_s)           # with replacement
    else:
        s_indices = rng.choice(n, size=min(num_s, n), replace=False)  # without replacement

    Gs = []
    for s in s_indices:
      t_idx = rng.integers(0, n, size=num_t_per_s)
      gbar = np.zeros(d)
      x_s = X[s]
      for t in t_idx:
        gbar += min_info_Besag_score_fullbatch(x1=x_s, x2=X[t], h=h, theta=theta_PLE)
        Gs.append(gbar / num_t_per_s)
    Gbar = np.mean(Gs, axis=0)
    I = sum(np.outer(g - Gbar, g - Gbar) for g in Gs) / len(Gs)

    # --- Minibatch approximation of J ---
    # Sample ordered pairs with t < s to mirror the full-batch loop.
    for _ in range(int(num_pairs_J)):
      s, t = rng.choice(n, size=2, replace=False)
      J += min_info_Besag_Hessian_fullbatch(x1=X[s], x2=X[t], h=h, theta=theta_PLE)
    J /= float(num_pairs_J)

    # === Asymptotic covariance (same scaling as in the full-batch code) ===
    samplesize_effective = 1/(1/RPJ_count+1/n)
    asymptoticCov = (linalg.inv(J) @ I @ linalg.inv(J)) * (4.0 / samplesize_effective)
    return asymptoticCov

def min_info_CLE_fullbatch(X,h,L,burnin,thin,max_L,max_iter,tol,detailop=False):
###############################################################################
	## conditional maximum likelihood estimation
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
###############################################################################
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
		init_theta = min_info_Besag_fullbatch(X=X,h=h,tol=10**(-2),max_iter=1000)
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