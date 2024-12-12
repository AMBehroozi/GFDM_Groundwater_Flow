# Generalized Finite Difference Method (GFDM) for Groundwater Flow

This repository contains the implementation of the Generalized Finite Difference Method (GFDM) for solving partial differential equations, specifically focusing on groundwater flow applications.

## 1. Problem Statement
The GFDM framework is applied to solve the 2D Poisson equation:

$\nabla^2 u = f(x, y) \quad \text{in} \quad \Omega$

with boundary conditions:

$u = g(x, y) \quad \text{on} \quad \partial\Omega$

## 2. GFDM Approximation
For each point $i$ and its neighboring points $j$, the Taylor series expansion provides a local approximation.

### Taylor Series Expansion
$u(x_j, y_j) = u(x_i, y_i) + \frac{\partial u}{\partial x}\Delta x_j + \frac{\partial u}{\partial y}\Delta y_j + \frac{1}{2}\frac{\partial^2 u}{\partial x^2}\Delta x_j^2 + \frac{\partial^2 u}{\partial x\partial y}\Delta x_j\Delta y_j + \frac{1}{2}\frac{\partial^2 u}{\partial y^2}\Delta y_j^2 + O(h^3)$

where:
- $\Delta x_j = x_j - x_i$
- $\Delta y_j = y_j - y_i$

### Matrix Form
For a point $i$ with $n$ neighbors:

$\mathbf{M}\boldsymbol{\alpha} = \mathbf{b}$

where:

M = | Δx₁    Δy₁    ½Δx₁²   Δx₁Δy₁   ½Δy₁² |
    | ...    ...    ...     ...      ...   |
    | Δxₙ    Δyₙ    ½Δxₙ²   ΔxₙΔyₙ   ½Δyₙ² |

α = | ∂u/∂x           |
    | ∂u/∂y           |
    | ∂²u/∂x²         |
    | ∂²u/∂x∂y        |
    | ∂²u/∂y²         |

    
## 3. GFDM Weights
The weights $\gamma_{ij}$ are computed using:

$\boldsymbol{\gamma} = \mathbf{M}^+ \mathbf{L}$

where:
- $\mathbf{M}^+$ is the pseudo-inverse of $\mathbf{M}$
- $\mathbf{L}$ is the Laplacian operator vector: $\mathbf{L} = [0, 0, 1, 0, 1]^T$

### Discrete Laplacian
The Laplacian at point $i$ is approximated as:

$\nabla^2 u_i \approx \sum_{j=1}^n \gamma_{ij}(u_j - u_i)$

## 4. Solution Process
The discrete system of equations for interior points becomes:

$\sum_{j=1}^n \gamma_{ij}(u_j - u_i) = f_i$

with boundary conditions:

$u_i = g_i \quad \text{on boundary points}$

## 5. Iterative Solution
Using an iterative scheme:

$u_i^{k+1} = u_i^k + \omega \left(\frac{f_i - \sum_{j=1}^n \gamma_{ij}(u_j^k - u_i^k)}{\sum_{j=1}^n \gamma_{ij}}\right)$

where:
- $k$ is the iteration number
- $\omega$ is a relaxation parameter for convergence control

## 6. Applications
This GFDM implementation is designed for groundwater flow problems but can be extended to other domains requiring meshless numerical solutions for PDEs.
