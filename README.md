# Generalized Finite Difference Method (GFDM) Formulation

## 1. Problem Statement
Consider a 2D Poisson equation:
$$
\nabla^2 u = f(x,y) \quad \text{in} \quad \Omega
$$
with boundary conditions:
$$
u = g(x,y) \quad \text{on} \quad \partial\Omega
$$

## 2. GFDM Approximation
For each point $i$ and its neighboring points $j$, we use a Taylor series expansion:

### Taylor Series Expansion
$$
u(x_j, y_j) = u(x_i, y_i) + \frac{\partial u}{\partial x}\Delta x_j + \frac{\partial u}{\partial y}\Delta y_j + \frac{1}{2}\frac{\partial^2 u}{\partial x^2}\Delta x_j^2 + \frac{\partial^2 u}{\partial x\partial y}\Delta x_j\Delta y_j + \frac{1}{2}\frac{\partial^2 u}{\partial y^2}\Delta y_j^2 + O(h^3)
$$
where:
- $\Delta x_j = x_j - x_i$
- $\Delta y_j = y_j - y_i$

### Matrix Form
For each point $i$ with $n$ neighbors, we form:
$$
\mathbf{M}\boldsymbol{\alpha} = \mathbf{b}
$$
where:
$$
\mathbf{M} = \begin{bmatrix}
\Delta x_1 & \Delta y_1 & \frac{1}{2}\Delta x_1^2 & \Delta x_1\Delta y_1 & \frac{1}{2}\Delta y_1^2 \\
\vdots & \vdots & \vdots & \vdots & \vdots \\
\Delta x_n & \Delta y_n & \frac{1}{2}\Delta x_n^2 & \Delta x_n\Delta y_n & \frac{1}{2}\Delta y_n^2
\end{bmatrix}
$$

$$
\boldsymbol{\alpha} = \begin{bmatrix}
\frac{\partial u}{\partial x} \\
\frac{\partial u}{\partial y} \\
\frac{\partial^2 u}{\partial x^2} \\
\frac{\partial^2 u}{\partial x\partial y} \\
\frac{\partial^2 u}{\partial y^2}
\end{bmatrix}
$$

### GFDM Weights
The weights $\gamma_{ij}$ are computed by:
$$
\boldsymbol{\gamma} = \mathbf{M}^+ \mathbf{L}
$$
where:
- $\mathbf{M}^+$ is the pseudo-inverse of $\mathbf{M}$
- $\mathbf{L}$ is the vector representing the Laplacian operator: $\mathbf{L} = [0,0,1,0,1]^T$

### Discrete Laplacian
The Laplacian at point $i$ is approximated as:
$$
\nabla^2 u_i \approx \sum_{j=1}^n \gamma_{ij}(u_j - u_i)
$$

## 3. Solution Process
The discrete system becomes:
$$
\sum_{j=1}^n \gamma_{ij}(u_j - u_i) = f_i
$$
for each interior point $i$, with:
$$
u_i = g_i \quad \text{for boundary points}
$$

## 4. Iterative Solution
Using an iterative scheme:
$$
u_i^{k+1} = u_i^k + \omega \left(\frac{f_i - \sum_{j=1}^n \gamma_{ij}(u_j^k - u_i^k)}{\sum_{j=1}^n \gamma_{ij}}\right)
$$
where:
- $k$ is the iteration number
- $\omega$ is a relaxation parameter