# Minimum information dependence modeling

(Written by Tomonari Sei and Keisuke Yano)

This pape provides R and python codes for minimum information dependence models developed in the following paper.

==========================================================================

Paper inforation: arXiv:2206.06792

Title: Minimum information dependence modeling

Authors: Tomonari Sei (The university of Tokyo), Keisuke Yano (The institute of statistical mathematis)

Abstract: We propose a method of constructing a joint statistical model for mixed-domain data to analyze their dependence. Multivariate Gaussian and log-linear models are particular examples of the proposed model. It is shown that the functional equation defining the model has a unique solution under fairly weak conditions. The model is characterized by two orthogonal sets of parameters: the dependence parameter and the marginal parameter. To estimate the dependence parameter, a conditional inference together with a sampling procedure is established and is shown to provide a consistent estimator of the dependence parameter. Illustrative examples of data analyses involving penguins and earthquakes are presented.

==========================================================================

The codes 

- min-info.R : definition of the functions used
- min-info.py : definition of the functions used
- min-info.ipynb: the instruction of the usage
- Mechsol_format.csv: Mechanism solution catalog :The original catalog is in the web page of Japan Meteological Agency and processed by the authors.

References:
-  JAPAN METEOROLOGICAL AGENCY (2022). The seismological bulletin of Japan. https://www.data.jma.go.jp/svd/eqev/data/bulletin/index_e.html.
-  A. Horst, A. Hill, K. Gorman (2020). palmerpenguins: Palmer Archipelago (Antarctica) penguin data. R package version 0.1.0. https://allisonhorst.github.io/palmerpenguins/. doi: 10.5281/zenodo.3960218.
