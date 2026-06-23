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
#* mindemo_mcmc : MCMC sampler for the model.
#* mindemo_Besag_fullbatch : Pseudo-likelihood estimation (PLE) using the full dataset.
#* mindemo_Besag_SGD : Pseudo-likelihood estimation via Stochastic Gradient Descent (SGD).
#* mindemo_Besag_bootstrap_fullbatch : Empirical bootstrap for pseudo-likelihood estimation using the full dataset.
#* mindemo_Besag_bootstrap_SGD : Empirical bootstrap for pseudo-likelihood estimation via SGD.
#* mindemo_Besag_bootstrap_SGD_Multicore : Empirical bootstrap for pseudo-likelihood estimation via SGD with parallel computation.
#* mindemo_bh_qvalues : Benjamini-Hochberg adjustment of p-values, also called q-values.
#* mindemo_bootstrapinference_summary : Summary table for bootstrap inference, including confidence intervals, p-values, q-values, and FDR selection.
#* mindemo_plot_bootstrap_forest : Forest plot visualization for bootstrap inference summaries.
#* mindemo_Besag_score / Hessian / asymptoticCovariance_fullbatch : Option to calculate the asymptotic covariance of the pseudo-likelihood estimator (full-batch version).
#* mindemo_Besag_score / Hessian / asymptoticCovariance_minibatch : Option to calculate the asymptotic covariance of the pseudo-likelihood estimator via SGD (mini-batch version).
#* mindemo_CLE_fullbatch : Conditional maximum likelihood estimation (CLE).
#* mindemo_CLE_damped_fullbatch : Stabilized conditional maximum likelihood estimation with damping, ridge regularization, and step-size control.

    
def mindemo_mcmc(X, h, loop, theta, burnin=0, thin=1, printop=False):
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

    rng = np.random.default_rng()
    theta = np.asarray(theta, dtype=float)

    current = [list(row) for row in X]
    n = len(current)
    d = len(current[0])

    samples = []
    acceptancenum = 0

    for Li in range(loop):
        i = rng.integers(0, d)
        s, t = rng.choice(n, size=2, replace=False)

        x_s = current[s]
        x_t = current[t]

        x_s_sw = x_s.copy()
        x_t_sw = x_t.copy()
        x_s_sw[i] = x_t[i]
        x_t_sw[i] = x_s[i]

        h_diff = h(x_s) + h(x_t) - h(x_s_sw) - h(x_t_sw)
        log_rho = min(0.0, -float(theta.dot(h_diff)))

        if np.log(rng.random()) < log_rho:
            current[s] = x_s_sw
            current[t] = x_t_sw
            acceptancenum += 1

        if Li >= burnin and ((Li - burnin) % thin == 0):
            samples.append([row.copy() for row in current])

        if printop and Li % 1000 == 0:
            rate = acceptancenum / (Li + 1)
            print(f"Step/Total : {Li}/{loop} Acceptance rate : {rate:.2f}")

    return samples

def mindemo_Besag_fullbatch(X,h,max_iter,tol,sparseop=False,sparsepen=1.):
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


def mindemo_Besag_bootstrap_fullbatch(X,h,tol,max_iter,B,sparseop=False,sparsepen=1.):
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
    bootstrap_Besag.append(mindemo_Besag_fullbatch(X=boot_sample,h=h,tol=tol,max_iter=max_iter,sparseop=sparseop,sparsepen=sparsepen))
  boot_mean = sum(bootstrap_Besag)/B
  boot_cov_1   = sum([np.tensordot(bootstrap_Besag[l],bootstrap_Besag[l],axes=0) for l in range(B)])/B
  boot_cov_2   = np.tensordot(boot_mean,boot_mean,axes=0)
  boot_cov     = boot_cov_1-boot_cov_2
  return bootstrap_Besag, boot_mean, boot_cov

def mindemo_Besag_SGD(
    X,
    h,
    max_iter=10000,
    tol=10**(-5),
    learningrate=0.1,
    stepsize_decay = 0.00,
    init = None,
    batch_size=1,
    *,
    # --- Stopping criteria option ---
    ema_beta=0.9,
    patience=5,
    # --- norm stabilization option ---
    clip_norm=None,
    # --- verbose option ---
    verbose_every=0,
    # --- seed option ---
    seed=None,
    # --- Ruppert--Polyak--Juditsky option ---
    use_rpj=True,
    avg_start=1000,
    return_avg=True,
    # --- Optimizer option ---
    # "sgd", "adam", "adamw", "rmsprop", or "adagrad"
    optimizer=None,
    # Backward-compatible optimizer flags. Used only when optimizer is None.
    use_adam=True,
    use_adamw=False,
    use_rmsprop=False,
    use_adagrad=False,
    # --- Regularization option ---
    # We maximize the penalized objective: pseudo-log-likelihood - (l2_penalty / 2) * ||coef||^2.
    l2_penalty=0.0,
    # --- Adam / AdamW option ---
    beta1=0.9,
    beta2=0.999,
    eps=1e-8,
    adamw_weight_decay=0.0,
    # --- RMSProp option ---
    rmsprop_beta=0.9,
    rmsprop_eps=1e-8,
    rmsprop_bias_correction=False,
    # --- AdaGrad option ---
    adagrad_eps=1e-8,
    adagrad_initial_accumulator=0.0
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
  ## optimizer     : Optional string; "sgd", "adam", "adamw", "rmsprop", or "adagrad".
  ##                 If None, use_adagrad/use_rmsprop/use_adamw/use_adam flags are used.
  ## use_adam      : If True, use Adam update when optimizer is None. Kept for backward compatibility.
  ## use_adamw     : If True, use AdamW update when optimizer is None.
  ## use_rmsprop   : If True, use RMSProp update when optimizer is None.
  ## use_adagrad   : If True, use AdaGrad update when optimizer is None.
  ## l2_penalty    : L2 penalty strength. In this code's sign convention, g is replaced by g - l2_penalty * coef.
  ## beta1, beta2  : Adam/AdamW momentum and RMS coefficients.
  ## eps           : Adam/AdamW numerical stabilizer.
  ## adamw_weight_decay : Decoupled weight decay strength for AdamW.
  ## rmsprop_beta  : RMSProp exponential moving average coefficient for squared gradients.
  ## rmsprop_eps   : RMSProp numerical stabilizer.
  ## rmsprop_bias_correction : If True, correct the initial bias of RMSProp's second moment estimate.
  ## adagrad_eps   : AdaGrad numerical stabilizer.
  ## adagrad_initial_accumulator : Initial value for the AdaGrad squared-gradient accumulator.
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

  # --- Optimizer selection ---
  if optimizer is None:
    # Backward-compatible precedence: explicit adaptive flags override use_adam.
    if use_adagrad:
      optimizer = "adagrad"
    elif use_rmsprop:
      optimizer = "rmsprop"
    elif use_adamw:
      optimizer = "adamw"
    elif use_adam:
      optimizer = "adam"
    else:
      optimizer = "sgd"

  optimizer = str(optimizer).lower()
  optimizer_aliases = {"adam_w": "adamw"}
  optimizer = optimizer_aliases.get(optimizer, optimizer)

  if optimizer not in {"sgd", "adam", "adamw", "rmsprop", "adagrad"}:
    raise ValueError(
        'optimizer must be one of "sgd", "adam", "adamw", '
        '"rmsprop", or "adagrad"'
    )

  use_adam = (optimizer == "adam")
  use_adamw = (optimizer == "adamw")
  use_rmsprop = (optimizer == "rmsprop")
  use_adagrad = (optimizer == "adagrad")

  if l2_penalty < 0:
    raise ValueError("l2_penalty must be non-negative")
  if adamw_weight_decay < 0:
    raise ValueError("adamw_weight_decay must be non-negative")
  if not (0.0 <= beta1 < 1.0):
    raise ValueError("beta1 must satisfy 0 <= beta1 < 1")
  if not (0.0 <= beta2 < 1.0):
    raise ValueError("beta2 must satisfy 0 <= beta2 < 1")

  # --- Optimizer state ---
  if use_adam or use_adamw:
    m = np.zeros_like(coef)
    v = np.zeros_like(coef)
  if use_rmsprop:
    if not (0.0 <= rmsprop_beta < 1.0):
      raise ValueError("rmsprop_beta must satisfy 0 <= rmsprop_beta < 1")
    rmsprop_v = np.zeros_like(coef)
  if use_adagrad:
    if adagrad_initial_accumulator < 0:
      raise ValueError("adagrad_initial_accumulator must be non-negative")
    adagrad_accumulator = np.full_like(coef, float(adagrad_initial_accumulator))

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

    # --- Optional L2 penalty ---
    # The code updates coef by adding g. Therefore, for the penalized objective
    # pseudo-log-likelihood - (l2_penalty / 2) * ||coef||^2, L2 contributes
    # -l2_penalty * coef to the update direction.
    if l2_penalty != 0.0:
      g = g - float(l2_penalty) * coef

    # --- Optional gradient clipping ---
    if clip_norm is not None and clip_norm > 0:
      gnorm = float(np.linalg.norm(g))
      if gnorm > clip_norm:
          g = g * (clip_norm / (gnorm + 1e-12))

    # --- Parameter update (SGD, Adam, AdamW, RMSProp, or AdaGrad) ---
    prev = coef.copy()
    eta_t = learningrate / (step ** stepsize_decay)
    if use_adam or use_adamw:
      # AdamW uses decoupled weight decay: do not add the decay term to g,
      # and do not accumulate it in m or v.
      if use_adamw and adamw_weight_decay != 0.0:
        coef = coef * (1.0 - eta_t * float(adamw_weight_decay))

      # Adam/AdamW moments
      m = beta1 * m + (1.0 - beta1) * g
      v = beta2 * v + (1.0 - beta2) * (g * g)
      # Bias correction
      m_hat = m / (1.0 - (beta1 ** step))
      v_hat = v / (1.0 - (beta2 ** step))
      # Update (note: this code uses "+ eta_t * g" sign convention)
      coef = coef + eta_t * (m_hat / (np.sqrt(v_hat) + eps))
    elif use_rmsprop:
      # RMSProp: coordinate-wise scaling by EMA of squared gradients
      rmsprop_v = rmsprop_beta * rmsprop_v + (1.0 - rmsprop_beta) * (g * g)
      if rmsprop_bias_correction:
        rmsprop_v_eff = rmsprop_v / (1.0 - (rmsprop_beta ** step))
      else:
        rmsprop_v_eff = rmsprop_v
      coef = coef + eta_t * (g / (np.sqrt(rmsprop_v_eff) + rmsprop_eps))
    elif use_adagrad:
      # AdaGrad: coordinate-wise scaling by cumulative squared gradients
      adagrad_accumulator += g * g
      coef = coef + eta_t * (g / (np.sqrt(adagrad_accumulator) + adagrad_eps))
    else:
      coef = coef + eta_t * g

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
          f"optimizer={optimizer} "
          f"||grad||={grad_norm:.3e} "
          f"||update||={dtheta:.3e} "
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
                          f"EMAed (||grad||)={ema_grad:.3e} ≤ {tol:.3e} and "
                          f"EMAed (||update||)={ema_dtheta:.3e} ≤ {tol*learningrate:.3e} "
                          f"for {patience} consecutive steps."
          )
          break
    else:
      consec_small = 0

  # --- Return the Ruppert--Polyak--Juditsky averagingif requested ---
  if use_rpj and return_avg and (avg_count > 0):
    return avg_coef,avg_count
  return coef,np.inf # in this case, the asysmptotic covariance evaluation does not work well

def mindemo_Besag_bootstrap_SGD(X,h,B,
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
    # --- Optimizer option ---
    optimizer=None,
    # --- Adam / AdamW option ---
    use_adam=True,
    use_adamw=False,
    beta1=0.9,
    beta2=0.999,
    eps=1e-8,
    # --- Regularization / decoupled weight decay option ---
    l2_penalty=0.0,
    adamw_weight_decay=0.0,
    # --- RMSProp option ---
    use_rmsprop=False,
    rmsprop_beta=0.9,
    rmsprop_eps=1e-8,
    rmsprop_bias_correction=False,
    # --- AdaGrad option ---
    use_adagrad=False,
    adagrad_eps=1e-8,
    adagrad_initial_accumulator=0.0
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
  ## optimizer     : Optional string; "sgd", "adam", "adamw", "rmsprop", or "adagrad".
  ##                 If None, use_adagrad/use_rmsprop/use_adamw/use_adam flags are used.
  ## use_adam      : If True, use Adam update when optimizer is None. Kept for backward compatibility.
  ## use_adamw     : If True, use AdamW update when optimizer is None.
  ## beta1, beta2  : Adam/AdamW momentum and RMS coefficients.
  ## eps           : Adam/AdamW numerical stabilizer.
  ## l2_penalty    : L2 penalty strength. In this code's sign convention, g is replaced by g - l2_penalty * coef.
  ## adamw_weight_decay : Decoupled weight decay strength for AdamW.
  ## use_rmsprop   : If True, use RMSProp update. Takes precedence over use_adamw/use_adam when optimizer is None.
  ## rmsprop_beta  : RMSProp exponential moving average coefficient for squared gradients.
  ## rmsprop_eps   : RMSProp numerical stabilizer.
  ## rmsprop_bias_correction : If True, correct the initial bias of RMSProp's second moment estimate.
  ## use_adagrad   : If True, use AdaGrad update. Takes precedence over use_rmsprop/use_adamw/use_adam when optimizer is None.
  ## adagrad_eps   : AdaGrad numerical stabilizer.
  ## adagrad_initial_accumulator : Initial value for the AdaGrad squared-gradient accumulator.
###############################################################################
  #OUTPUT
  ## coef : np.ndarray, shape (p,). The final parameter vector.
###############################################################################
  n = len(X)
  bootstrap_Besag=[]
  for b in range(B):
    boot_sample=[X[i] for i in random.choices(range(n),k=n)]
    bootstrap_Besag.append(mindemo_Besag_SGD(X=boot_sample,h=h,tol=tol,max_iter=max_iter,
                                               learningrate=learningrate,stepsize_decay = stepsize_decay,
                                              init=init,
                                              batch_size = batch_size,
                                              ema_beta=ema_beta,patience=patience,clip_norm=clip_norm,verbose_every=verbose_every,
                                              seed=seed,use_rpj=use_rpj,avg_start=avg_start,return_avg=return_avg,
                                              optimizer=optimizer,
                                              use_adam=use_adam,
                                              use_adamw=use_adamw,
                                              beta1=beta1,beta2=beta2,eps=eps,
                                              l2_penalty=l2_penalty,
                                              adamw_weight_decay=adamw_weight_decay,
                                              use_rmsprop=use_rmsprop,
                                              rmsprop_beta=rmsprop_beta,
                                              rmsprop_eps=rmsprop_eps,
                                              rmsprop_bias_correction=rmsprop_bias_correction,
                                              use_adagrad=use_adagrad,
                                              adagrad_eps=adagrad_eps,
                                              adagrad_initial_accumulator=adagrad_initial_accumulator)[0])
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
    optimizer: Optional[str],
    use_adam: bool,
    use_adamw: bool,
    beta1: float,
    beta2: float,
    eps: float,
    l2_penalty: float,
    adamw_weight_decay: float,
    use_rmsprop: bool,
    rmsprop_beta: float,
    rmsprop_eps: float,
    rmsprop_bias_correction: bool,
    use_adagrad: bool,
    adagrad_eps: float,
    adagrad_initial_accumulator: float
):
###############################################################################
#####The function required for multicore bootstrap
#####Perform a single bootstrap replicate (executed in a worker process).
###############################################################################
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed % (2**32 - 1))

    boot_sample = [X[i] for i in indices]
    theta = mindemo_Besag_SGD(
        X=boot_sample, h=h, tol=tol, max_iter=max_iter,
        learningrate=learningrate,
    stepsize_decay = stepsize_decay,
        init = init,
        batch_size=batch_size,
        ema_beta=ema_beta,
        patience=patience, clip_norm=clip_norm, verbose_every=verbose_every,
        seed=seed, use_rpj=use_rpj, avg_start=avg_start, return_avg=return_avg,
        optimizer=optimizer,
        use_adam=use_adam,
        use_adamw=use_adamw,
        beta1=beta1, beta2=beta2, eps=eps,
        l2_penalty=l2_penalty,
        adamw_weight_decay=adamw_weight_decay,
        use_rmsprop=use_rmsprop,
        rmsprop_beta=rmsprop_beta,
        rmsprop_eps=rmsprop_eps,
        rmsprop_bias_correction=rmsprop_bias_correction,
        use_adagrad=use_adagrad,
        adagrad_eps=adagrad_eps,
        adagrad_initial_accumulator=adagrad_initial_accumulator
    )[0]
    return theta


def mindemo_Besag_bootstrap_SGD_Multicore(
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
    # --- Optimizer option ---
    optimizer=None,
    # --- Adam / AdamW option ---
    use_adam=True,
    use_adamw=False,
    beta1=0.9,
    beta2=0.999,
    eps=1e-8,
    # --- Regularization / decoupled weight decay option ---
    l2_penalty=0.0,
    adamw_weight_decay=0.0,
    # --- RMSProp option ---
    use_rmsprop=False,
    rmsprop_beta=0.9,
    rmsprop_eps=1e-8,
    rmsprop_bias_correction=False,
    # --- AdaGrad option ---
    use_adagrad=False,
    adagrad_eps=1e-8,
    adagrad_initial_accumulator=0.0,
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
  ## optimizer     : Optional string; "sgd", "adam", "adamw", "rmsprop", or "adagrad".
  ##                 If None, use_adagrad/use_rmsprop/use_adamw/use_adam flags are used.
  ## use_adam      : If True, use Adam update when optimizer is None. Kept for backward compatibility.
  ## use_adamw     : If True, use AdamW update when optimizer is None.
  ## beta1, beta2  : Adam/AdamW momentum and RMS coefficients.
  ## eps           : Adam/AdamW numerical stabilizer.
  ## l2_penalty    : L2 penalty strength. In this code's sign convention, g is replaced by g - l2_penalty * coef.
  ## adamw_weight_decay : Decoupled weight decay strength for AdamW.
  ## use_rmsprop   : If True, use RMSProp update. Takes precedence over use_adamw/use_adam when optimizer is None.
  ## rmsprop_beta  : RMSProp exponential moving average coefficient for squared gradients.
  ## rmsprop_eps   : RMSProp numerical stabilizer.
  ## rmsprop_bias_correction : If True, correct the initial bias of RMSProp's second moment estimate.
  ## use_adagrad   : If True, use AdaGrad update. Takes precedence over use_rmsprop/use_adamw/use_adam when optimizer is None.
  ## adagrad_eps   : AdaGrad numerical stabilizer.
  ## adagrad_initial_accumulator : Initial value for the AdaGrad squared-gradient accumulator.
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
                optimizer=optimizer,
                use_adam=use_adam,
                use_adamw=use_adamw,
                beta1=beta1, beta2=beta2, eps=eps,
                l2_penalty=l2_penalty,
                adamw_weight_decay=adamw_weight_decay,
                use_rmsprop=use_rmsprop,
                rmsprop_beta=rmsprop_beta,
                rmsprop_eps=rmsprop_eps,
                rmsprop_bias_correction=rmsprop_bias_correction,
                use_adagrad=use_adagrad,
                adagrad_eps=adagrad_eps,
                adagrad_initial_accumulator=adagrad_initial_accumulator
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
                    use_rpj, avg_start, return_avg, optimizer,
                    use_adam, use_adamw, beta1, beta2, eps,
                    l2_penalty, adamw_weight_decay,
                    use_rmsprop, rmsprop_beta, rmsprop_eps, rmsprop_bias_correction,
                    use_adagrad, adagrad_eps, adagrad_initial_accumulator
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

def mindemo_bh_qvalues(pvals):
###############################################################################
## Benjamini-Hochberg adjusted p-values, also called q-values.
###############################################################################
  # INPUT
  ## pvals : array-like, shape (m,). Raw p-values.
###############################################################################
  # OUTPUT
  ## qvals : np.ndarray, shape (m,). Benjamini-Hochberg adjusted p-values.
###############################################################################

    pvals = np.asarray(pvals, dtype=float)

    if pvals.ndim != 1:
        raise ValueError("pvals must be a 1D array.")

    if np.any((pvals < 0.0) | (pvals > 1.0)):
        raise ValueError("pvals must be between 0 and 1.")

    m = pvals.size
    if m == 0:
        return np.asarray([], dtype=float)

    order = np.argsort(pvals)
    sorted_p = pvals[order]

    q_sorted = sorted_p * m / np.arange(1, m + 1)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.minimum(q_sorted, 1.0)

    qvals = np.empty(m, dtype=float)
    qvals[order] = q_sorted

    return qvals

def mindemo_bootstrapinference_summary(
    coef,
    boot,
    labels=None,
    fdr_alpha=0.05,
    ci_percentiles=(2.5, 97.5),
    p_floor=None
):

###############################################################################
## Bootstrap inference summary for coefficient estimates.
###############################################################################
  # INPUT
  ## coef          : array-like, shape (p,). Point estimates of coefficients.
  ## boot          : array-like, shape (B, p). Bootstrap estimates.
  ## labels        : array-like or None. Feature / coefficient labels. If None, labels are generated automatically.
  ## fdr_alpha     : FDR threshold for Benjamini-Hochberg selection.
  ## ci_percentiles: tuple(float, float). Percentiles used for bootstrap confidence intervals.
  ## p_floor       : float or None. Lower bound used only for -log10(q) visualization.
###############################################################################
  # OUTPUT
  ## result_df : pandas.DataFrame. 
  ## Summary table with coefficients, bootstrap confidence intervals,
  ## p-values, Benjamini-Hochberg q-values, -log10(q), FDR selection,
  ## and confidence-interval exclusion indicator.
###############################################################################

    coef = np.asarray(coef, dtype=float)
    boot = np.asarray(boot, dtype=float)

    if boot.ndim != 2:
        raise ValueError("boot must be a 2D array with shape (B, p).")

    B, p = boot.shape

    if coef.shape != (p,):
        raise ValueError(f"coef must have shape ({p},), but got {coef.shape}.")

    if labels is None:
        feature_labels = np.asarray([f"term_{j}" for j in range(p)])
    else:
        feature_labels = np.asarray(labels)
        if feature_labels.shape[0] != p:
            raise ValueError(
                f"labels must have length {p}, but got {feature_labels.shape[0]}."
            )

    if not (0.0 < fdr_alpha < 1.0):
        raise ValueError("fdr_alpha must be between 0 and 1.")

    if p_floor is None:
        p_floor = 1.0 / (B + 1.0)

    # Bootstrap two-sided p-values for H0: coefficient = 0
    pvals = 2.0 * np.minimum(
        np.mean(boot <= 0.0, axis=0),
        np.mean(boot >= 0.0, axis=0),
    )
    pvals = np.minimum(pvals, 1.0)

    # BH-FDR correction
    qvals = mindemo_bh_qvalues(pvals)

    # Stable score for visualization / ranking
    qvals_for_score = np.maximum(qvals, p_floor)
    importance = -np.log10(qvals_for_score)

    # Bootstrap confidence intervals
    ci_low = np.percentile(boot, ci_percentiles[0], axis=0)
    ci_high = np.percentile(boot, ci_percentiles[1], axis=0)

    selected_fdr = qvals < fdr_alpha
    ci_excludes_zero = (ci_low > 0.0) | (ci_high < 0.0)

    result_df = pd.DataFrame({
        "term": feature_labels,
        "coef": coef,
        f"ci_low_{ci_percentiles[0]}%": ci_low,
        f"ci_high_{ci_percentiles[1]}%": ci_high,
        "p_value": pvals,
        "q_value_BH": qvals,
        "-log10(q)": importance,
        f"FDR_q<{fdr_alpha}": selected_fdr,
        "CI_excludes_0": ci_excludes_zero,
    })

    return result_df
    
def mindemo_plot_bootstrap_forest(
    result_df,
    *,
    top_n=None,
    sort_by="-log10(q)",
    ascending=False,
    fdr_alpha=0.05,
    title="Bootstrap inference summary",
):

###############################################################################
## Forest plot for bootstrap inference summary.
###############################################################################
  # INPUT
  ## result_df : pandas.DataFrame. Output of mindemo_bootstrapinference_summary.
  ## top_n     : int or None. Number of top rows to plot after sorting. If None, plot all rows.
  ## sort_by   : str or None. Column name used for sorting. If None, keep the original order.
  ## ascending : bool. Sorting order.
  ## fdr_alpha : FDR threshold used to identify selected terms.
  ## title     : title of the plot.
###############################################################################
  # OUTPUT
  ## fig : matplotlib.figure.Figure. Figure object.
  ## ax  : matplotlib.axes.Axes. Axes object.
###############################################################################

    df = result_df.copy()

    if sort_by is not None:
        df = df.sort_values(sort_by, ascending=ascending)

    if top_n is not None:
        df = df.head(top_n)

    df = df.iloc[::-1].reset_index(drop=True)

    ci_low_col = [c for c in df.columns if c.startswith("ci_low_")][0]
    ci_high_col = [c for c in df.columns if c.startswith("ci_high_")][0]

    terms = df["term"].astype(str).values
    coef = df["coef"].values
    ci_low = df[ci_low_col].values
    ci_high = df[ci_high_col].values

    xerr = np.vstack([coef - ci_low, ci_high - coef])

    sig_col = f"FDR_q<{fdr_alpha}"
    if sig_col in df.columns:
        sig = df[sig_col].values
    else:
        sig = df["q_value_BH"].values < fdr_alpha

    fig_height = max(4, 0.4 * len(df))
    fig, ax = plt.subplots(figsize=(9, fig_height))

    y = np.arange(len(df))

    ax.errorbar(
        coef[~sig],
        y[~sig],
        xerr=xerr[:, ~sig],
        fmt="o",
        capsize=3,
        label="Not FDR-selected",
    )

    ax.errorbar(
        coef[sig],
        y[sig],
        xerr=xerr[:, sig],
        fmt="o",
        capsize=3,
        label="FDR-selected",
    )

    ax.axvline(0.0, linestyle="--", linewidth=1)

    ax.set_yticks(y)
    ax.set_yticklabels(terms)
    ax.set_xlabel("PLE coefficient")
    ax.set_title(title)
    ax.legend()

    fig.tight_layout()
    return fig, ax

def mindemo_Besag_score_fullbatch(x1,x2,h,theta):
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

def mindemo_Besag_Hessian_fullbatch(x1,x2,h,theta):
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

def mindemo_Besag_asymptoticCovariance_fullbatch(X,h,theta_PLE):
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
      conditional_score = conditional_score + (1/n)* mindemo_Besag_score_fullbatch(x1=X[s],x2=X[t],h=h,theta=theta_PLE)
      if (t<s):
        B=mindemo_Besag_Hessian_fullbatch(x1=X[s],x2=X[t],h=h,theta=theta_PLE)
        J=J+B/(n*(n-1)/2)
    I = I + np.tensordot(conditional_score ,conditional_score ,axes=0) / n

  asymptoticCov = (linalg.inv(J) @ I @ linalg.inv(J)) * (4/n)

  return asymptoticCov


def mindemo_Besag_asymptoticCovariance_minibatch(
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
        gbar += mindemo_Besag_score_fullbatch(x1=x_s, x2=X[t], h=h, theta=theta_PLE)
        Gs.append(gbar / num_t_per_s)
    Gbar = np.mean(Gs, axis=0)
    I = sum(np.outer(g - Gbar, g - Gbar) for g in Gs) / len(Gs)

    # --- Minibatch approximation of J ---
    # Sample ordered pairs with t < s to mirror the full-batch loop.
    for _ in range(int(num_pairs_J)):
      s, t = rng.choice(n, size=2, replace=False)
      J += mindemo_Besag_Hessian_fullbatch(x1=X[s], x2=X[t], h=h, theta=theta_PLE)
    J /= float(num_pairs_J)

    # === Asymptotic covariance (same scaling as in the full-batch code) ===
    samplesize_effective = 1/(1/RPJ_count+1/n)
    asymptoticCov = (linalg.inv(J) @ I @ linalg.inv(J)) * (4.0 / samplesize_effective)
    return asymptoticCov

def mindemo_CLE_fullbatch(
    X,
    h,
    L=1000,
    burnin=100,
    thin=10,
    max_L=100000,
    max_iter=20,
    tol=1e-4,
    detailop=False,
    init=None
):
###############################################################################
    ## conditional maximum likelihood estimation
###############################################################################
  #INPUT
  ## X       : initial sample list [[multivariate sample 1],[multivariate sample 2],...]
  ## h       : canonical statistics (from X (d-dim list) to R^K array)
  ## L       : MCMC iteration (length of chain)
  ## burnin  : burnin of MCMC
  ## thin    : thinning
  ## max_L   : upper bound for adaptive MCMC length
  ## max_iter: number of iterations
  ## tol     : tolerance value for the optimizer
  ## detailop: if True, print (iteration, L, current estimate, current residual)
  ## init    : initial value of theta for warm start. If None, use Besag estimator.
###############################################################################
    #OUTPUT
    ## list of estimates (by iteration)
    ## list of residuals for optimization (by iteration)
    ## list of L (by iteration)
    ## list of asymptotic covariances (by iteration)
###############################################################################

    if len(X) == 0:
        print("NULL list")
    else:
        n = len(X)
        d = len(X[0])

        # --- warm start option ---
        if init is None:
            init_theta = mindemo_Besag_fullbatch(
                X=X, h=h, tol=10**(-2), max_iter=1000
            )
        else:
            init_theta = np.asarray(init, dtype=float)
            p = len(h(X[0]))
            if init_theta.shape != (p,):
                raise ValueError(f"init must have shape ({p},), but got {init_theta.shape}")

        init_res = 10**3
        Li = 0
        theta_list = [init_theta]
        res_list = [init_res]
        L_list = [L]
        AsymptoticCov_list = [0]
        current_theta = init_theta
        current_res = init_res

        while ((current_res > tol) and (Li < max_iter) and (L < max_L)):
            Li += 1
            current_theta = theta_list[Li - 1]
            current_res = res_list[Li - 1]

            current_X_MCMC = mindemo_mcmc(
                X=X, h=h, loop=L, theta=current_theta, printop=False
            )[burnin:L:thin]

            total_MCMC_num = len(current_X_MCMC)
            current_hsum = sum(list(map(h, X)))
            current_hsum_MCMCsamples = [
                sum(list(map(h, current_X_MCMC[l])))
                for l in range(total_MCMC_num)
            ]
            current_hsum_mean = sum(current_hsum_MCMCsamples) / total_MCMC_num
            current_hsum_cov_1 = sum([
                np.tensordot(
                    current_hsum_MCMCsamples[l],
                    current_hsum_MCMCsamples[l],
                    axes=0
                )
                for l in range(total_MCMC_num)
            ]) / total_MCMC_num
            current_hsum_cov_2 = np.tensordot(current_hsum_mean, current_hsum_mean, axes=0)
            current_hsum_cov = current_hsum_cov_1 - current_hsum_cov_2

            score = current_hsum - current_hsum_mean
            modified_score = linalg.solve(current_hsum_cov, score)

            theta_list.append(current_theta + modified_score)
            res_list.append(max(abs(score)))
            L_list.append(L)
            AsymptoticCov_list.append(linalg.pinv(current_hsum_cov))

            if detailop:
                print(Li, L, theta_list[Li], res_list[Li])

            if current_res > 1.1 * res_list[Li]:
                L = int(np.floor(1.5 * L))

            if max(abs(theta_list[Li]) / (1 + abs(theta_list[0]))) > 100:
                theta_list[Li] = theta_list[0]
                print("parameter reseted.")

        theta_list.pop(-1)
        res_list.pop(-1)
        AsymptoticCov_list.pop(-1)

    return theta_list, res_list, L_list, AsymptoticCov_list
    
def mindemo_CLE_damped_fullbatch(
    X,
    h,
    L=1000,
    burnin=100,
    thin=10,
    max_L=100000,
    max_iter=20,
    tol=1e-4,
    detailop=False,
    init=None,
    step_size=0.2,
    ridge=1e-6,
    max_step_norm=1.0,
    use_pinv_fallback=True
):
###############################################################################
    ## conditional maximum likelihood estimation (damped / stabilized version)
###############################################################################
  #INPUT
  ## X                 : initial sample list [[multivariate sample 1],[...],...]
  ## h                 : canonical statistics (from X (d-dim list) to R^K array)
  ## L                 : MCMC iteration (length of chain)
  ## burnin            : burnin of MCMC
  ## thin              : thinning
  ## max_L             : upper bound for adaptive MCMC length
  ## max_iter          : number of outer iterations
  ## tol               : tolerance value for the optimizer
  ## detailop          : if True, print iteration info
  ## init              : initial value of theta for warm start. If None, use Besag estimator.
  ## step_size         : damping factor for CLE update
  ## ridge             : ridge term added to covariance matrix for numerical stability
  ## max_step_norm     : clip update vector by its L2 norm if needed
  ## use_pinv_fallback : if True, use pseudo-inverse when linear solve fails
###############################################################################
    #OUTPUT
    ## list of estimates (by iteration)
    ## list of residuals for optimization (by iteration)
    ## list of L (by iteration)
    ## list of asymptotic covariances (by iteration)
###############################################################################

    if len(X) == 0:
        print("NULL list")
        return None

    n = len(X)
    d = len(X[0])
    p = len(h(X[0]))

    # --- warm start option ---
    if init is None:
        init_theta = mindemo_Besag_fullbatch(
            X=X,
            h=h,
            tol=10**(-2),
            max_iter=1000
        )
    else:
        init_theta = np.asarray(init, dtype=float)
        if init_theta.shape != (p,):
            raise ValueError(f"init must have shape ({p},), but got {init_theta.shape}")

    init_res = 10**3
    Li = 0
    theta_list = [init_theta.copy()]
    res_list = [init_res]
    L_list = [L]
    AsymptoticCov_list = [np.zeros((p, p))]
    current_theta = init_theta.copy()
    current_res = init_res

    while ((current_res > tol) and (Li < max_iter) and (L < max_L)):
        Li += 1
        current_theta = theta_list[Li - 1]
        current_res = res_list[Li - 1]

        current_X_MCMC = mindemo_mcmc(
            X=X,
            h=h,
            loop=L,
            theta=current_theta,
            printop=False
        )[burnin:L:thin]

        total_MCMC_num = len(current_X_MCMC)
        if total_MCMC_num == 0:
            raise ValueError("No MCMC samples left after burnin/thin. Increase L or reduce burnin/thin.")

        current_hsum = sum(list(map(h, X)))
        current_hsum_MCMCsamples = [
            sum(list(map(h, current_X_MCMC[l])))
            for l in range(total_MCMC_num)
        ]

        current_hsum_mean = sum(current_hsum_MCMCsamples) / total_MCMC_num
        current_hsum_cov_1 = sum([
            np.tensordot(
                current_hsum_MCMCsamples[l],
                current_hsum_MCMCsamples[l],
                axes=0
            )
            for l in range(total_MCMC_num)
        ]) / total_MCMC_num
        current_hsum_cov_2 = np.tensordot(current_hsum_mean, current_hsum_mean, axes=0)
        current_hsum_cov = current_hsum_cov_1 - current_hsum_cov_2

        # ridge stabilization
        stabilized_cov = current_hsum_cov + ridge * np.eye(p)

        score = current_hsum - current_hsum_mean

        # stable linear solve
        try:
            raw_step = linalg.solve(stabilized_cov, score)
        except Exception:
            if use_pinv_fallback:
                raw_step = linalg.pinv(stabilized_cov) @ score
            else:
                raise

        # step clipping
        raw_step_norm = np.linalg.norm(raw_step)
        if (max_step_norm is not None) and (raw_step_norm > max_step_norm):
            raw_step = raw_step * (max_step_norm / raw_step_norm)

        # damped update
        theta_candidate = current_theta + step_size * raw_step

        theta_list.append(theta_candidate)
        res_list.append(np.max(np.abs(score)))
        L_list.append(L)
        AsymptoticCov_list.append(linalg.pinv(stabilized_cov))

        if detailop:
            print(Li, L, theta_list[Li], res_list[Li])

        # If residual improves sufficiently, modestly increase MCMC length
        if current_res > 1.1 * res_list[Li]:
            L = int(np.floor(1.5 * L))

        # reset if estimate moves too far from initial point
        if np.max(np.abs(theta_list[Li]) / (1 + np.abs(theta_list[0]))) > 100:
            theta_list[Li] = theta_list[0].copy()
            print("parameter reseted.")

    theta_list.pop(-1)
    res_list.pop(-1)
    AsymptoticCov_list.pop(-1)

    return theta_list, res_list, L_list, AsymptoticCov_list