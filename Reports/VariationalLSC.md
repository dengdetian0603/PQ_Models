# Latent Discrete Graphical Model

Measurement $M$, latent state L, parameters  $\Theta = (\delta, \gamma, \theta, \rho, D)$, case/control indicator $Y$.

Bronze standard only  $[M, L, \Theta] =\{\prod_{k=1}^K[M_k|L_k, \Theta]\} [L|\Theta][\Theta]$
$$
[M, L, \Theta] \propto \exp \Big\{ \sum_{k=1}^K \Big [\log(1-\delta_k) + M_k\log \frac{\delta_k}{1-\delta_k} + YL_k \log \frac{1-\gamma_k}{1-\delta_k} + YL_k M_k \log \frac{\gamma_k(1-\delta_k)}{\delta_k(1-\gamma_k)} \Big] + \\
Y\Big[\sum_{k=1}^K (L_k\theta_k +\rho\sum_{k \neq k'}L_kL_{k'}D_{kk'}) - A(\theta, \rho, D)\Big]\Big\} [\Theta]
$$

Silver standard data likelihood 
$$
[M^{SS}|L, \Theta] \propto \exp \Big \{ \sum_{k=1}^K L_kM_k^{SS} \log \frac{\eta_k}{1 - \eta_k} + L_k \log(1 - \eta_k) + M_k^{SS}\log L_k \Big \}\\
$$

### Variational approximate posterior distribution of $\Theta$.

Assume the approximate posterior $Q(\Theta, L)$ to be fully factorized, i.e.
$$
Q(\Theta, L) = q(\delta, \gamma)q(L)q(\theta)q(\rho)q(D) = \prod_{\eta \in \{\Theta, L\}}q(\eta)
$$
We maximize the lower bound of the joint density (minimize the KL divergence) to find $Q$.
$$
Q(\Theta, L) = \text{argmax} \{ \mathbb{E}_Q\Big[\log\frac{[M, L, \Theta]}{Q(\Theta, L)}\Big]\} \\
q(\eta) \propto \exp\{\mathbb{E}_{Q(-\eta)}(\log [M, L, \Theta])\} \\
$$
__For q$(\gamma, \delta)$__:
$$
q(\delta, \gamma) \propto \exp\Big\{ \mathbb{E}_L\big [\sum_{k=1}^K \log [M_k|L_k, \delta, \gamma]+\log[\delta, \gamma]\big]\Big \} \\
\propto \exp \Big\{ \sum_{k=1}^K \Big[\log(1-\delta_k) + M_k\log \frac{\delta_k}{1-\delta_k} +
Yq_{l_k} \big( \log \frac{1-\gamma_k}{1-\delta_k} \big) +
Yq_{l_k} M_k \log \frac{\gamma_k(1-\delta_k)}{\delta_k(1-\gamma_k)} \Big] \Big \}[\gamma, \delta]\\
= \prod_{k=1}^K \Big[ (1-\gamma_k)^{Yq_{l_k}(1-M_k)} \gamma_k^{Yq_{l_k}M_k}\Big]\Big[ (1-\delta_k)^{1+Yq_{l_k}M_k-Yq_{l_k}-M_k}\delta_k^{M_k(1-Yq_{l_k})}\Big] [\gamma, \delta]
$$
Where $q_{l_k} = \mathbb{E}_Q(L_k|Y=1)$. Let $a_{k}^*,b_{k}^*$ be the hyper-parametes for $\gamma_k$, and let $c^*_k,d^*_k$ be the hyper-parameters for $\delta_k$. Then the posterior distribution is: $\gamma_k \sim \text{Beta}(A^*_k, B^*_k)$ and $\delta_k \sim \text{Beta}(A'_k, B'_k)$ where
$$
A^*_k = Yq_{l_k}M_k + a^*_k \\
B^*_k = Yq_{l_k}(1-M_k) + b^*_k\\
A'_K = M_k(1-Yq_{l_k}) + c^*_K\\
B'_K = 1+Yq_{l_k}M_k-Yq_{l_k}-M_k + d^*_K
$$
**With $n​$ independent samples**:
$$
A^*_{kn} = \sum_{i=1}^n Y_iq_{l_{ik}}M_{ik} + a^*_k \\
B^*_{kn} = \sum_{i=1}^n Y_iq_{l_{ik}}(1-M_{ik}) + b^*_k\\
A'_{kn} = \sum_{i=1}^n M_{ik}(1-Y_iq_{l_{ik}}) + c^*_K\\
B'_{kn} = n+\sum_{i=1}^n (Y_iq_{l_{ik}}M_{ik}-Y_iq_{l_{ik}}-M_{ik}) + d^*_K\\
\tilde{A}_{kn} = \sum_{i=1}^{n_{case}} q_{l_{ik}}M^{SS}_{ik} + \tilde{a}_k\\
\tilde{B}_{kn} = \sum_{i=1}^{n_{case}} q_{l_{ik}}(1-M^{SS}_{ik}) + \tilde{b}_k\\
$$
__For $q(L|Y=1)$__:
$$
q(L|Y=1) = \prod_k q(L_k)\propto \exp\Big\{ \mathbb{E}_{\Theta}\big [\sum_{k=1}^K \log [M_k|L_k, \delta, \gamma] +
\log[L|\theta, \rho,D]\big]\Big \} \\
= \exp \Big\{ \mathbb{E} \sum_{k=1}^K \Big[ L_k (\theta_k + \log \frac{1-\gamma_k}{1-\delta_k} +
M_k\log \frac{\gamma_k(1-\delta_k)}{\delta_k(1-\gamma_k)}) + \rho \sum_{k'\neq k}L_kL_{k'}D_{kk'}\Big]\Big\}
$$
Let $h_k = \theta_k + \log \frac{1-\gamma_k}{1-\delta_k} +M_k\log \frac{\gamma_k(1-\delta_k)}{\delta_k(1-\gamma_k)}​$, and let $d_{kk'} = \mathbb{E}_Q(D_{kk'})​$, then
$$
q(L_k) \propto \exp \Big\{ Y\Big[L_k\mathbb{E}(h_k) + \rho L_k\sum_{k'\neq k}q_{l_{k'}}d_{kk'}\Big]\Big\}\\
q_{l_k} = \frac{\exp(H_k + \mu_\rho\sum_{k'\neq k}q_{l_{k'}}d_{kk'})}{1 + \exp(H_k + \mu_\rho\sum_{k'\neq k}q_{l_{k'}}d_{kk'})}
$$
Where with $\psi(x) = \frac{\partial \log \Gamma(x)}{\partial x} \text{ and } \psi(x+y) = \frac{\partial \log \Gamma(x+y)}{\partial x}​$
$$
H^{BS}_k = \mathbb{E}_Q(h_k) =  M_k[\psi(A^*_k) - \psi(A'_k)] + (1 - M_k)[\psi(B^*_k) - \psi(B'_k)] \\
+ \mu_{\theta_k} - \psi(A^*_k + B^*_k) + \psi(A'_k + B'_k)\\
$$
**with $n$ independent samples, indexed by $i$**:
$$
q_{l_{ik}} = \frac{\exp(H_{ik} + \mu_\rho\sum_{k'\neq k}q_{l_{ik'}}d_{kk'})}{1 +
\exp(H_{ik} + \mu_\rho\sum_{k'\neq k}q_{l_{ik'}}d_{kk'})}\\
H^{BS}_{ik} =  M_{ik}[\psi(A^*_{kn}) - \psi(A'_{kn})] + (1 - M_{ik})[\psi(B^*_{kn}) - \psi(B'_{kn})] \\
+ \mu_{\theta_k} - \psi(A^*_{kn} + B^*_{kn}) + \psi(A'_{kn} + B'_{kn})\\
H^{SS}_{ik}|_{M_{ik}^{SS}=0} = \psi(\tilde{B}_{kn}) - \psi(\tilde{A}_{kn} + \tilde{B}_{kn})\\
H^{BS+SS}_{ik}|_{M_{ik}^{SS}=0} = H^{BS}_{ik} + \psi(\tilde{B}_{kn}) - \psi(\tilde{A}_{kn} + \tilde{B}_{kn})\\
q_{l_{ik}} = 1 \text{ if } M_{ik}^{SS} = 1
$$
__For $q(\theta)​$__:
$$
q(\theta) = \prod_kq(\theta_k) \propto \exp \Big\{ Y\mathbb{E}_{(L,\rho,D)}\big[\sum_{k=1}^K L_k \theta_k - A(\theta, \rho, D)\big] \Big\}[\theta]\\
\approx \exp \Big\{ Y\sum_{k=1}^K \big[q_{l_k} \theta_k -
\mathbb{E}\log \big(1 + \exp(\theta_k +\rho\sum_{k' \neq k}L_{k'}D_{kk'})\big)\big] \Big\}[\theta] 
\text{ by pseudo-likelihood}
$$
By taylor expansion, 
$$
\mathbb{E}_y\log(1+e^{x+y}) = \log(1+e^{x + \mu_y}) +
\frac{e^{x + \mu_y}}{2(1 + e^{x + \mu_y})^2}\text{Var}(y) + \mathbb{E}_y[o(y-\mu_y)^3]
$$
Let $R_k = \sum_{k' \neq k}L_{k'}D_{kk'}$, then 
$$
\mathbb{E}(R_k) = \sum_{k' \neq k}q_{l_{k'}}d_{kk'}\\
\mathbb{Var}(R_k) = \sum_{k' \neq k}q_{l_{k'}}d_{kk'}(1 - q_{l_{k'}}d_{kk'})\\
\mathbb{Var}(\rho R_k) = \mathbb{E}(\rho^2)\mathbb{E}(R_k^2) - \mathbb{E}(\rho)^2\mathbb{E}(R_k)^2\\
= (\sigma^2_{\rho} +  \mu_{\rho}^2)[\mathbb{Var}(R_k)  +  \mathbb{E}(R_k)^2] - \mu_{\rho}^2\mathbb{E}(R_k) ^2\\
= (\sigma^2_{\rho} +  \mu_{\rho}^2)\sum_{k' \neq k}q_{l_{k'}}d_{kk'}(1 - q_{l_{k'}}d_{kk'}) +
\sigma^2_{\rho} [\sum_{k' \neq k}q_{l_{k'}}d_{kk'}]^2
$$
Thus, $q(\theta_k)$ is further approximated by
$$
q(\theta_k) \propto \exp \Big \{ Y \big [q_{l_k}\theta_k - \log(1 + e^{\theta_k + \mu_{\rho}\sum_{k' \neq k}q_{l_{k'}}d_{kk'}}) - 
\frac{ e^{\theta_k + \mu_{\rho}\sum_{k' \neq k}q_{l_{k'}}d_{kk'}}}{2(1 + e^{\theta_k + \mu_{\rho}\sum_{k' \neq k}q_{l_{k'}}d_{kk'}})^2}\mathbb{Var}(\rho R_k)\big ] -\frac{1}{2}(\theta_k - \mu^*_{\theta})^2\tau^*_{\theta}\Big\}
$$


**with $n$ independent samples**:
$$
q(\theta_k) \propto \exp \Big \{ \sum_{i=1}^n Y_i \Big [q_{l_{ik}} \theta_k - 
\log(1 + e^{\theta_k + \mu_{\rho}\sum_{k' \neq k}q_{l_{ik'}}d_{kk'}}) 
- \frac{ e^{\theta_k + \mu_{\rho}\sum_{k' \neq k}q_{l_{ik'}}d_{kk'}}}{2(1 + e^{\theta_k + 
\mu_{\rho}\sum_{k' \neq k}q_{l_{ik'}}d_{kk'}})^2} \mathbb{Var}(\rho R_{ik}) \Big ]
- \frac{1}{2}(\theta_k - \mu^*_{\theta})^2\tau^*_{\theta}\Big\}\\
\mathbb{Var}(\rho R_{ik}) = (\sigma^2_{\rho} +  \mu_{\rho}^2)\sum_{k' \neq k}q_{l_{ik'}}d_{kk'}(1 - q_{l_{ik'}}d_{kk'}) +
\sigma^2_{\rho} [\sum_{k' \neq k}q_{l_{ik'}}d_{kk'}]^2
$$
Then $\mu_{\theta_k} \approx \text{argmax}\ q(\theta_k)$, and $\sigma^2_{\theta_k} \approx - \Big[ \frac{\partial^2 \log q(\theta_k)}{\partial \theta_k^2} \Big|_{\theta_k = \mu_{\theta_k}} \Big]^{-1}$.



__For $q(\rho)$__:
$$
q(\rho) \propto \exp \Big\{ Y\mathbb{E}_{(\theta,L,D)}\big[\rho\sum_{k=1}^K \sum_{k'\neq k} L_kL_{k'}D_{kk'}
- A(\theta, \rho, D)\big] \Big\}[\rho]\\
\approx  \exp \Big\{ Y \sum_{k=1}^K \big[ \rho \sum_{k'\neq k} q_{l_k}q_{l_{k'}}d_{kk'} -
\mathbb{E}\log \big(1 + \exp(\theta_k +\rho\sum_{k' \neq k}L_{k'}D_{kk'})\big)\big] \Big\} [\rho]
$$
By Le Cam's theorem, $R_k$ can be approximated by Poisson$(\sum_{k' \neq k} q_{l_{k'}}d_{kk'})$, thus
$$
q(\rho) \propto \exp \Big\{ Y \sum_{k=1}^K \big[ \rho \sum_{k'\neq k} q_{l_k}q_{l_{k'}}d_{kk'} -
\sum_{j=0}^{K-1} \text{Poi}(j;\lambda = \sum_{k' \neq k} q_{l_{k'}}d_{kk'})\Big ( \log (1 + e^{\mu_{\theta_k} + j \rho}) +
\frac{e^{\mu_{\theta_k} + j\rho}}{2(1 + e^{\mu_{\theta_k} + j\rho})^2}\sigma^2_{\theta_k}\Big) \big]\Big\} [\rho]
$$


**With $n$ independent samples**:
$$
q(\rho) \propto \exp \Big\{ \sum_{i=1}^n Y_i \sum_{k=1}^K \big[ \rho\sum_{k'\neq k} q_{l_{ik}}q_{l_{ik'}}d_{kk'} -
\sum_{j=0}^{K-1} \text{Poi}(j;\lambda = \sum_{k' \neq k} q_{l_{ik'}}d_{kk'})\Big ( \log (1 + e^{\mu_{\theta_k} + j \rho}) +
\frac{e^{\mu_{\theta_k} + j\rho}}{2(1 + e^{\mu_{\theta_k} + j\rho})^2}\sigma^2_{\theta_k}\Big) \big]\Big\} [\rho]
$$
Then $\mu_{\rho} \approx \text{argmax}\ q(\rho)$, and $\sigma^2_{\rho} \approx - \Big[ \frac{\partial^2 \log q(\rho)}{\partial \rho^2} \Big|_{\rho = \mu_{\rho}} \Big]^{-1}$.



__For $q(D)$__:
$$
q(D) \propto \exp \Big\{ Y\mathbb{E}_{(\theta,L,\rho)}\big[\rho\sum_{k=1}^K \sum_{k'\neq k} L_kL_{k'}D_{kk'}
- A(\theta, \rho, D)\big] \Big\}[D] \\
\approx  \exp \Big\{ Y \sum_{k=1}^K \big[ \mu_\rho\sum_{k'\neq k} q_{l_k}q_{l_{k'}}d_{kk'} -
\mathbb{E}\log \big(1 + \exp(\theta_k +\rho\sum_{k' \neq k}L_{k'}D_{kk'})\big)\big] \Big\} [D]\\
q(D_{k_1k_2}) \propto \exp \Big\{ 2Y\mu_{\rho} q_{l_{k_1}}q_{l_{k_2}}D_{k_1k_2} -
Y\mathbb{E} \log (1 + e^{\theta_{k_1} + \rho\sum_{k \neq k_1, k\neq k_2}L_kD_{kk_1} + \rho L_{k_2}D_{k_1k_2}}) \\
- Y\mathbb{E} \log (1 + e^{\theta_{k_2} + \rho\sum_{k \neq k_1, k\neq k_2}L_kD_{kk_2} + \rho L_{k_1}D_{k_1k_2}}) +
D_{k_1k_2}\mathbb{E}_{\pi_d}(\log\frac{\pi_d}{1 - \pi_d})\Big\}
$$
Let $S_{k_1k_2} = \theta_{k_1} + \rho\sum_{k \neq k_1, k\neq k_2}L_kD_{kk_1}$, $T_{k_1k_2} = \mathbb{E} \log ( 1 + e^{S_{k_1k_2} + \rho L_{k_2}}), T'_{k_1k_2} = \mathbb{E} \log ( 1 + e^{S_{k_1k_2}})$, and let $A_{\pi_d}, B_{\pi_d}$ be the Beta hyper-parameters for $\pi_d$, then for all $k_1 < k_2$:
$$
\mathbb{E}_{\pi_d} (\log \frac{\pi_d}{1-\pi_d}) = \psi(A_{\pi_d}) - \psi(B_{\pi_d})\\
\mathbb{Var}(S_{k_1k_2}) = \sigma^2_{\theta_{k_1}} + (\sigma^2_{\rho} + \mu^2_{\rho})
\sum_{k\neq k_1,k \neq k_2} q_{l_k}d_{kk_1}(1 - q_{l_k}d_{kk_1}) +
\sigma^2_{\rho}[\sum_{k\neq k_1,k \neq k_2} q_{l_k}d_{kk_1}]^2\\
\mathbb{Var}(\rho L_{k}) =  (\sigma^2_{\rho} + \mu^2_{\rho}) q_{l_k} - \mu^2_{\rho}q^2_{l_k}\\
$$
Thus
$$
T_{k_1k_2} \approx \log ( 1 + e^{ \mu_{\theta_{k_1}} + \mu_{\rho}\sum_{k \neq k_1, k\neq k_2}q_{l_k}d_{kk_1} +
\mu_{\rho} q_{l_{k_2}}}) + \frac{e^{ \mu_{\theta_{k_1}} + \mu_{\rho}\sum_{k \neq k_1, k\neq k_2}q_{l_k}d_{kk_1} +
\mu_{\rho} q_{l_{k_2}}}}{2(1 + e^{ \mu_{\theta_{k_1}} + \mu_{\rho}\sum_{k \neq k_1, k\neq k_2}q_{l_k}d_{kk_1} +
\mu_{\rho} q_{l_{k_2}}})^2} [\mathbb{Var}(S_{k_1k_2}) + \mathbb{Var}(\rho L_{k_2})]\\
T'_{k_1k_2} \approx \log ( 1 + e^{ \mu_{\theta_{k_1}} + \mu_{\rho}\sum_{k \neq k_1, k\neq k_2}q_{l_k}d_{kk_1}}) +
\frac{e^{ \mu_{\theta_{k_1}} + \mu_{\rho}\sum_{k \neq k_1, k\neq k_2}q_{l_k}d_{kk_1}}}{2(1 + e^{ \mu_{\theta_{k_1}} +
\mu_{\rho}\sum_{k \neq k_1, k\neq k_2}q_{l_k}d_{kk_1} })^2} \mathbb{Var}(S_{k_1k_2}) \\
d_{k_1k_2} = \frac{\exp \big \{ 2Y\mu_{\rho}q_{l_{k_1}}q_{l_{k_2}} - YT_{k_1k_2} -
YT_{k_2k_1} + \psi(A_{\pi_d}) - \psi(B_{\pi_d}) \big \}}
{\exp \big \{ 2Y\mu_{\rho}q_{l_{k_1}}q_{l_{k_2}} - YT_{k_1k_2} -
YT_{k_2k_1} + \psi(A_{\pi_d}) - \psi(B_{\pi_d}) \big \} +
\exp \big \{  - YT'_{k_1k_2} - YT'_{k_2k_1}\big \}}
$$


**With $n$ independent samples**:
$$
d_{k_1k_2} = \frac{\exp \big \{ \sum_{i=1}^n Y_i [ 2\mu_{\rho}q_{l_{ik_1}}q_{l_{ik_2}} - T^{(i)}_{k_1k_2} -
T^{(i)}_{k_2k_1}] + \psi(A_{\pi_d}) - \psi(B_{\pi_d}) \big \}}
{\exp \big \{ \sum_{i=1}^n Y_i [ 2\mu_{\rho}q_{l_{ik_1}}q_{l_{ik_2}} - T^{(i)}_{k_1k_2} -
T^{(i)}_{k_2k_1}] + \psi(A_{\pi_d}) - \psi(B_{\pi_d}) \big \} +
\exp \big \{  - \sum_{i=1}^n Y_i [T'^{(i)}_{k_1k_2} + T'^{(i)}_{k_2k_1}]\big \}}
$$





### Covariance Matrix Calibration

$$
\mathbb{E}_Q\log[M, L, \Theta] = \sum_{k=1}^K \Big \{ (\sum_{i=1}^n Y_iq_{l_{ik}}M_{ik} + a^*_k  ) \mathbb{E}_Q (\log \gamma_k) +
\big [\sum_{i=1}^n Y_iq_{l_{ik}}(1-M_{ik}) + b^*_k \big ] \mathbb{E}_Q [\log(1-\gamma_k)] +\\
\big [\sum_{i=1}^n M_{ik}(1-Y_iq_{l_{ik}}) + c^*_K \big] \mathbb{E}_Q (\log \delta_k) + 
\big [ n+\sum_{i=1}^n (Y_iq_{l_{ik}}M_{ik}-Y_iq_{l_{ik}}-M_{ik}) + d^*_K\big] \mathbb{E}_Q [\log (1- \delta_k)] +\\
(\sum_{i=1}^{n_{case}} q_{l_{ik}}M^{SS}_{ik} + \tilde{a}_k) \mathbb{E}_Q (\log \eta_k) +
\big[ \sum_{i=1}^{n_{case}} q_{l_{ik}}(1-M^{SS}_{ik}) + \tilde{b}_k \big] \mathbb{E}_Q[\log(1 - \eta_k)] + 
\sum_{i=1}^{n_{case}} M^{SS}_{ik} (1 - q_{ik}) (-\infty) \Big \} + \\
\sum_{i=1}^{n_{case}} \sum_{k=1}^K \Big \{ q_{l_{ik}} X_i^T \hat{\beta}_{\cdot k} +
\hat{\rho} \sum_{k'\neq k} q_{l_{ik}}q_{l_{ik'}}d_{kk'}  -
\mathbb{E}_Q\log \big [1 + \exp(X_i^T\beta_{\cdot k}+\rho\sum_{k' \neq k}L_{ik'}D_{kk'})\big]\Big \}  +
\mathbb{E}_Q \log [\beta, \rho, D]\\
$$

where we denote $\hat{\theta}_{ik} = X_i^T \hat{\beta}_{\cdot k}$ and 
$$
\mathbb{E}_Q\log \big [1 + \exp(X_i^T\beta_{\cdot k}+\rho\sum_{k' \neq k}L_{ik'}D_{kk'})\big] \\
\approx \sum_{j=0}^{K-1} \text{Poi}(j;\lambda = \sum_{k' \neq k} q_{l_{ik'}}d_{kk'})
\Big [ \log (1 + e^{\hat{\theta}_{ik} + j \hat{\rho}}) +
\frac{e^{\hat{\theta}_{ik} + j\hat{\rho}}}{2(1 + e^{\hat{\theta}_{ik} + j\hat{\rho}})^2}
(\sigma^2_{\theta_{ik}} + j^2\sigma^2_{\rho})\Big]
$$

Denote $\mathbb{E}_Q \log[M, L, \Theta]$ as $\mathcal{L}$, and let $\lambda_{ik} = \sum_{k' \neq k} q_{l_{ik'}}d_{kk'}$ then

qL:
$$
\frac{\partial \mathcal{L}}{\partial q_{l_{ik}}} = 
M_{ik}^{SS} \mathbb{E}_Q(\log \eta_k) + (1 - M_{ik}^{SS}) \mathbb{E}_Q[\log(1 - \eta_k)] + \\
M_{ik}^{BS} \mathbb{E}_Q(\log \gamma_k) + (1 - M_{ik}^{BS}) \mathbb{E}_Q[\log(1 - \gamma_k)] - \\
M_{ik}^{BS} \mathbb{E}_Q(\log \delta_k) - (1 - M_{ik}^{BS}) \mathbb{E}_Q[\log(1 - \delta_k)] +\\
X_i^T\beta_{\cdot k} + 2 \rho \lambda_{ik} - \\
\sum_{k' \neq k} \sum_{j=0}^{K-1} \frac{\lambda_{ik'}^{j - 1}e^{-\lambda_{ik'}}}{j\text{!}} (j - \lambda_{ik'}) d_{kk'}
\Big [ \log (1 + e^{\hat{\theta}_{ik'} + j \hat{\rho}}) +
\frac{e^{\hat{\theta}_{ik'} + j\hat{\rho}}}{2(1 + e^{\hat{\theta}_{ik'} + j\hat{\rho}})^2}
(\sigma^2_{\theta_{ik'}} + j^2\sigma^2_{\rho})\Big]\\
$$
D:
$$
\text{for } k_1 < k_2, \;\frac{\partial \mathcal{L}}{\partial d_{k_1k_2}} = 
2  \sum_{i=1}^{n_{case}} \Big \{ \hat{\rho}\sum_{k_1 <k_2} q_{l_{ik_1}} q_{l_{ik_2}} - \\
\sum_{j=0}^{K-1} \frac{\lambda_{ik_1}^{j - 1}e^{-\lambda_{ik_1}}}{j\text{!}} (j - \lambda_{ik_1}) q_{l_{ik_2}}
\Big [ \log (1 + e^{\hat{\theta}_{ik_1} + j \hat{\rho}}) +
\frac{e^{\hat{\theta}_{ik_1} + j\hat{\rho}}}{2(1 + e^{\hat{\theta}_{ik_1} + j\hat{\rho}})^2}
(\sigma^2_{\theta_{ik_1}} + j^2\sigma^2_{\rho})\Big]  -\\
\sum_{j=0}^{K-1} \frac{\lambda_{ik_2}^{j - 1}e^{-\lambda_{ik_2}}}{j\text{!}} (j - \lambda_{ik_2}) q_{l_{ik_1}}
\Big [ \log (1 + e^{\hat{\theta}_{ik_2} + j \hat{\rho}}) +
\frac{e^{\hat{\theta}_{ik_2} + j\hat{\rho}}}{2(1 + e^{\hat{\theta}_{ik_2} + j\hat{\rho}})^2}
(\sigma^2_{\theta_{ik_2}} + j^2\sigma^2_{\rho})\Big] \Big \}\\
$$
$\mathbb{\beta}_{\cdot k}$:
$$
\frac{\partial \mathcal{L}}{\partial \hat{\beta}_{\cdot k}} = 
\sum_{i=1}^{n_{case}} q_{l_{ik}}X_i - 
\sum_{i=1}^{n_{case}} \sum_{j=0}^{K-1}\text{Poi}(j; \lambda_{ik})
\Big[ \frac{e^{\hat{\theta}_{ik} + j \hat{\rho}}}{1 + e^{\hat{\theta}_{ik} + j \hat{\rho}}}+
\frac{e^{\hat{\theta}_{ik} + j\hat{\rho}} - e^{3(\hat{\theta}_{ik} + j\hat{\rho})}}
{2(1 + e^{\hat{\theta}_{ik} + j\hat{\rho}})^4} 
(\sigma^2_{\theta_{ik}} + j^2\sigma^2_{\rho})
\Big] X_i \\
$$
$\rho$:
$$
\frac{\partial \mathcal{L}}{\partial \hat{\rho}} = 
2 \sum_{i=1}^{n_{case}}\sum_{k_1 <k_2} q_{l_{ik_1}} q_{l_{ik_2}} d_{k_1k_2} - 
\sum_{i=1}^{n_{case}} \sum_{k=1}^K\sum_{j=0}^{K-1}\text{Poi}(j; \lambda_{ik})
\Big[ \frac{e^{\hat{\theta}_{ik} + j \hat{\rho}}}{1 + e^{\hat{\theta}_{ik} + j \hat{\rho}}}+
\frac{e^{\hat{\theta}_{ik} + j\hat{\rho}} - e^{3(\hat{\theta}_{ik} + j\hat{\rho})}}
{2(1 + e^{\hat{\theta}_{ik} + j\hat{\rho}})^4} 
(\sigma^2_{\theta_{ik}} + j^2\sigma^2_{\rho})
\Big] j \\
$$


Second order partial derivatives for qL:

vs. qL:
$$
\frac{\partial^2 \mathcal{L}}{\partial q_{l_{ik}} \partial q_{l_{i'k'}}} =  0\\
\frac{\partial^2 \mathcal{L}}{\partial q_{l_{ik_1}} \partial q_{l_{ik_2}}} =  2 \hat{\rho} d_{k_1k_2} - \\
\sum_{k' \neq k_1, k' \neq k_2}\sum_{j = 0}^{K-1} \frac{1}{j\text{!}} \lambda_{ik'}^{j - 2}e^{-\lambda_{ik'}}[j(j-1) -2 j\lambda_{ik'} +
\lambda_{ik'}^2]d_{k_1k'}d_{k_2k'} 
\Big [ \log (1 + e^{\hat{\theta}_{ik'} + j \hat{\rho}}) +
\frac{e^{\hat{\theta}_{ik'} + j\hat{\rho}}}{2(1 + e^{\hat{\theta}_{ik'} + j\hat{\rho}})^2}
(\sigma^2_{\theta_{ik'}} + j^2\sigma^2_{\rho})\Big]\\
$$


vs. D:
$$
\frac{\partial^2 \mathcal{L}}{\partial q_{l_{ik_1}} \partial d_{k_1k_2}} =  2 \hat{\rho} q_{l_{ik_2}} - \\
\sum_{j=0}^{K-1} \Big \{ \frac{1}{j\text{!}} \Big [ \lambda_{ik_2}^{j - 2}e^{-\lambda_{ik_2}}[j(j-1) -2 j\lambda_{ik_2} +
\lambda_{ik_2}^2]q_{l_{ik_1}}^2 + \lambda_{ik_2}^{j - 1}e^{-\lambda_{ik_2}} (j - \lambda_{ik_2}) \Big ]\\
\Big [ \log (1 + e^{\hat{\theta}_{ik_2} + j \hat{\rho}}) +
\frac{e^{\hat{\theta}_{ik_2} + j\hat{\rho}}}{2(1 + e^{\hat{\theta}_{ik_2} + j\hat{\rho}})^2}
(\sigma^2_{\theta_{ik_2}} + j^2\sigma^2_{\rho})\Big] \Big \}\\
\\ %----------------------------------------------------------------------------------------------------------------------------
\frac{\partial^2 \mathcal{L}}{\partial q_{l_{ik_2}} \partial d_{k_1k_2}} =  2 \hat{\rho} q_{l_{ik_1}} - \\
\sum_{j=0}^{K-1} \Big \{ \frac{1}{j\text{!}} \Big [ \lambda_{ik_1}^{j - 2}e^{-\lambda_{ik_1}}[j(j-1) -2 j\lambda_{ik_1} +
\lambda_{ik_1}^2]q_{l_{ik_2}}^2 + \lambda_{ik_1}^{j - 1}e^{-\lambda_{ik_1}} (j - \lambda_{ik_1}) \Big ]\\
\Big [ \log (1 + e^{\hat{\theta}_{ik_1} + j \hat{\rho}}) +
\frac{e^{\hat{\theta}_{ik_1} + j\hat{\rho}}}{2(1 + e^{\hat{\theta}_{ik_1} + j\hat{\rho}})^2}
(\sigma^2_{\theta_{ik_1}} + j^2\sigma^2_{\rho})\Big] \Big \}\\
$$
​	for $k_3 \neq k_2$ and $k_3 \neq k1$,

​	
$$
\frac{\partial^2 \mathcal{L}}{\partial q_{l_{ik_3}} \partial d_{k_1k_2}} = -
\sum_{j=0}^{K-1} \Big \{ \frac{1}{j\text{!}} \lambda_{ik_2}^{j - 2}e^{-\lambda_{ik_2}}[j(j-1) -2 j\lambda_{ik_2} +
\lambda_{ik_2}^2]q_{l_{ik_1}}q_{l_{ik_3}}\\
\Big [ \log (1 + e^{\hat{\theta}_{ik_2} + j \hat{\rho}}) +
\frac{e^{\hat{\theta}_{ik_2} + j\hat{\rho}}}{2(1 + e^{\hat{\theta}_{ik_2} + j\hat{\rho}})^2}
(\sigma^2_{\theta_{ik_2}} + j^2\sigma^2_{\rho})\Big] \Big \} -\\
\sum_{j=0}^{K-1} \Big \{ \frac{1}{j\text{!}} \lambda_{ik_1}^{j - 2}e^{-\lambda_{ik_1}}[j(j-1) -2 j\lambda_{ik_1} +
\lambda_{ik_1}^2]q_{l_{ik_2}}q_{l_{ik_3}}\\
\Big [ \log (1 + e^{\hat{\theta}_{ik_1} + j \hat{\rho}}) +
\frac{e^{\hat{\theta}_{ik_1} + j\hat{\rho}}}{2(1 + e^{\hat{\theta}_{ik_1} + j\hat{\rho}})^2}
(\sigma^2_{\theta_{ik_1}} + j^2\sigma^2_{\rho})\Big] \Big \}\\
$$
vs. $\beta$:
$$
\frac{\partial^2 \mathcal{L}}{\partial q_{l_{ik}} \partial \hat{\beta}_{\cdot k}} = X_i \\
\frac{\partial^2 \mathcal{L}}{\partial q_{l_{ik}} \partial \hat{\beta}_{\cdot k'}} =
\sum_{j=0}^{K-1} \frac{\lambda_{ik'}^{j - 1}e^{-\lambda_{ik'}}}{j\text{!}} (j - \lambda_{ik'}) d_{kk'}
\Big[ \frac{e^{\hat{\theta}_{ik'} + j \hat{\rho}}}{1 + e^{\hat{\theta}_{ik'} + j \hat{\rho}}}+
\frac{e^{\hat{\theta}_{ik'} + j\hat{\rho}} - e^{3(\hat{\theta}_{ik'} + j\hat{\rho})}}
{2(1 + e^{\hat{\theta}_{ik'} + j\hat{\rho}})^4} 
(\sigma^2_{\theta_{ik'}} + j^2\sigma^2_{\rho})
\Big] X_i \\
$$
vs. $\rho$:
$$
\frac{\partial^2 \mathcal{L}}{\partial q_{l_{ik}} \partial \hat{\rho}} = 
2 \lambda_{ik} - \\
\sum_{k' \neq k} \sum_{j=0}^{K-1} \frac{\lambda_{ik'}^{j - 1}e^{-\lambda_{ik'}}}{j\text{!}} (j - \lambda_{ik'}) d_{kk'}
\Big[ \frac{e^{\hat{\theta}_{ik'} + j \hat{\rho}}}{1 + e^{\hat{\theta}_{ik'} + j \hat{\rho}}}+
\frac{e^{\hat{\theta}_{ik'} + j\hat{\rho}} - e^{3(\hat{\theta}_{ik'} + j\hat{\rho})}}
{2(1 + e^{\hat{\theta}_{ik'} + j\hat{\rho}})^4} 
(\sigma^2_{\theta_{ik'}} + j^2\sigma^2_{\rho})
\Big]j
$$
