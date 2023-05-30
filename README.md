# Minimum information dependence modeling

(Written by Tomonari Sei and Keisuke Yano)

This pape provides R and python codes for minimum information dependence models developed in the following paper:

==========================================================================

Paper inforation: arXiv:2206.06792

Title: Minimum information dependence modeling

Authors: Tomonari Sei (The university of Tokyo), Keisuke Yano (The institute of statistical mathematis)

Abstract: We propose a method of constructing a joint statistical model for mixed-domain data to analyze their dependence. Multivariate Gaussian and log-linear models are particular examples of the proposed model. It is shown that the functional equation defining the model has a unique solution under fairly weak conditions. The model is characterized by two orthogonal sets of parameters: the dependence parameter and the marginal parameter. To estimate the dependence parameter, a conditional inference together with a sampling procedure is established and is shown to provide a consistent estimator of the dependence parameter. Illustrative examples of data analyses involving penguins and earthquakes are presented.

==========================================================================


# Codes 

- min-info.R : definition of the functions used
- min-info.py : definition of the functions used
- min-info.ipynb: the instruction of the usage (created by using Google Colaboratory; filepath should be adequately changed)
- Mechsol_format.csv: Mechanism solution catalog :The original catalog is in the web page of Japan Meteological Agency and processed by the authors.

# Brief summary of the minimum information dependence model

## 1. Minimum information dependence model

The minimum information dependence model is a joint model for a mixed-domain data proposed by Sei and Yano, Minimum information dependence modeling:

$p(x; \theta, \nu)=\exp(\theta^{\top}h(x)-\sum_{j=1,\ldots,d}a_{j}(x_j;\theta,\nu)-\psi(\theta,\nu))\prod_{j=1,\ldots,d}r_{j}(x_{j};\nu)$

such that $a_{j}(x_j;\theta,\nu)$ and $\psi(\theta,\nu)$ are determined by

--mariginal condition $\int p(x;\theta,\nu)dx_{-j}=r_{j}(x_{j};\nu)$

--identifiability condition $\int \sum_{j=1,\ldots,d} a_{j}(x_{j};\theta,\nu) p(x;\theta,\nu)dx=0$.

The model admits various types of the domains of variables (say, continuou/categorical/manifold/etc...) and various types of dependence (higher-order interaction/negative interaction/etc...)

An example of the minimum information dependence model of Poisson and Beta marginals with a negative interaction is given in the following figure:

![Sampling from minimum information dependence model with Poisson /Beta marginals](img/Figure_PoissonBeta.png "Poisson Beta marginals")


## 2. Inference 





References:
-  JAPAN METEOROLOGICAL AGENCY (2022). The seismological bulletin of Japan. https://www.data.jma.go.jp/svd/eqev/data/bulletin/index_e.html.
-  A. Horst, A. Hill, K. Gorman (2020). palmerpenguins: Palmer Archipelago (Antarctica) penguin data. R package version 0.1.0. https://allisonhorst.github.io/palmerpenguins/. doi: 10.5281/zenodo.3960218.
