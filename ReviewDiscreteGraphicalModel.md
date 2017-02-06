## Review on Discrete Graphical Model Structure Estimation

- *S.-I. Lee, V. Ganapathi, and D. Koller. Efficient structure learning of markov networks using l1-regularization. NIPS. 2007.*

  Log-linear model parameterization with L1 penalty on canonical association parameters. Dense graph with large number of nodes can't be optimized directly. Node introduction procedure and approximate inference (loopy belief propagation) is used. 

  ​

- *P. Ravikumar, M. J. Wainwright, and J. Lafferty. High-dimensional ising model selection using ℓ1-regularized logistic regression. Annals of Statistics. 2010*

  Auto-logistic regrerssion with l1-penalty for each node, proof of sparsistency, i.e. the correct recovery of signed neighborhood set with probability 1. Asymptotic learning only. No discussion on finite sample conflict resolve.

  ​

- *Jalali, Ali, et al. On Learning Discrete Graphical Models using Group-Sparse Regularization. AISTATS. 2011.*

  ​

- *Loh, Po-Ling, and Martin J. Wainwright. Structure estimation for discrete graphical models: Generalized covariance matrices and their inverses. NIPS. 2012.*

  The inverse of the covariance matrix of an augmented discrete graph (by triangulation) is block-structured. Tree structure is already triangulated, its inverse covariance matrix is thus block-structured. Log-determinant method (Graphical Lasso) with thresholding  can estimate the covariance matrix of a tree structured graph consistently.