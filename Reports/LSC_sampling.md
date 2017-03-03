# LSC sampling algorithm

## Latent variable sampling strategy

### Sample exactly from multinomial distribution

- Time complexity

$$
O(K^22^K + n)
$$



### Sample sequentially from auto-logistic conditional distribution

- Time complexity
  $$
  O(K^2n)
  $$








## Regression coefficients sampling strategy

### Metropolis-Hastings algorithm with pseudolikelihood

$$
pL = \prod_{k=1}^K p(L_k | L_{-k}; \beta_k, \theta^{(2)})
$$

---



#Latent Gaussian with INLA

Measurement $M$, latent Gaussian variable $Z$, latent state $L = h(Z)$, parameters  $\theta = (\delta, \gamma, \mu_z, \rho_z)$

Bronze standard only  $[M, Z, \theta] =\{\prod_{k=1}^K[M_k|Z_k, \theta]\} [Z|\theta] [\theta]$
$$
[M, Z, \theta] \propto \exp \Big\{ \sum_{k=1}^K \log(1-\delta_k) + \sum_{k=1}^K M_k\log \frac{\delta_k}{1-\delta_k} + \sum_{k=1}^K h(Z_k)\log \frac{1-\gamma_k}{1-\delta_k} + \\ \sum_{k=1}^K h(Z_k)M_k \log \frac{\gamma_k(1-\delta_k)}{\delta_k(1-\gamma_k)} - \frac{1}{2}(Z - \mu_z)^T \Sigma^{-1} (Z-\mu_z)\Big\} [\theta]
$$


Gaussian approximation of the posterior distribution of $Z$, by using a quadratic approximation of $\sum_{k=1}^Kg(M_k, Z_k, \theta) = \log \pi(M|Z, \theta)$ w.r.t $Z$.
$$
\tilde{\pi}_G(Z|\theta, M) \propto [Z|\theta] [M|Z,\theta] \\
\propto \exp \Big\{ - \frac{1}{2}(Z - \mu_z)^T \Sigma^{-1} (Z-\mu_z) + \sum_{k=1}^K g(M_k, Z_k, \theta) \Big\}\\
\propto  \exp \Big\{ - \frac{1}{2} [Z - \mu_{z|m}(\theta)]^T \Lambda(\theta) [Z-\mu_{z|m}(\theta)] \Big\}
$$


The approximated posterior distribution of $\theta$, where $Z^*$ is the mode of $\tilde{\pi}_G(Z|\theta,M)$.
$$
\tilde{\pi}(\theta|M) = \frac{[M, Z, \theta]}{\tilde{\pi}_G(Z|\theta, M)} \Big|_{Z = Z^*} \propto \frac{[Z,\theta| M]}{[Z|\theta, M]}
$$
This approximation is not a proper density function, but it helps with the exploration of the posterior space and the numeric integration for getting the marginal density.
$$
\hat{\pi}(\theta_j|M) = \int \tilde{\pi}(\theta|M) d\theta_{-j} \\
\hat{\pi}(Z_k|M) = \int \tilde{\pi}(Z_k|\theta, M) \tilde{\pi}(\theta|M) d\theta
$$
$\tilde{\pi}(Z_k|\theta, M) $ can be approximated in three ways, where simplified Laplacian approximation is recommended.

---

# Latent Discrete Graphical Model

Measurement $M$, latent state L, parameters  $\Theta = (\delta, \gamma, \theta^{(1)}, \theta^{(2)})$

Bronze standard only  $[M, L, \Theta] =\{\prod_{k=1}^K[M_k|L_k, \Theta]\} [L|\Theta][\Theta]$
$$
[M, L, \Theta] \propto \exp \Big\{ \sum_{k=1}^K \log(1-\delta_k) + \sum_{k=1}^K M_k\log \frac{\delta_k}{1-\delta_k} + \sum_{k=1}^K L_k \big(\theta_k^{(1)} + \log \frac{1-\gamma_k}{1-\delta_k} \big) + \\ \sum_{k=1}^K L_k M_k \log \frac{\gamma_k(1-\delta_k)}{\delta_k(1-\gamma_k)} - \sum_{k \neq k'}L_kL_{k'} \theta^{(2)}_{kk'} \Big\} [\Theta]
$$



### May use mean field approximation of $[L|M,\Theta]$ to do:

- i. Variational EM to find MAP of $\Theta$.


- ii. Approximate posterior distribution of $\Theta$. 
- ​






