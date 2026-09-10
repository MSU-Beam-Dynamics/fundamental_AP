---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Linear Algebra

## Matrices and linear relations

A linear relation can be expressed by an $m\times n$ matrix connecting an input
of dimension $m$ and a response of dimension $n$:

$$
A=\begin{pmatrix}
a_{11} & a_{12} & \cdots & a_{1n}\\
a_{21} & a_{22} & \cdots & a_{2n}\\
\vdots & \vdots & \ddots & \vdots\\
a_{m1} & a_{m2} & \cdots & a_{mn}
\end{pmatrix}
$$

### Linear properties

If $A$ and $B$ are two matrices of the same size $m\times n$, their sum
$C=A+B=B+A$ has the same size, with $C_{ij}=A_{ij}+B_{ij}$. A matrix can be
scaled by a constant number $k$ (scalar multiplication):

$$
C=kA,\qquad C_{ij}=k A_{ij}
$$

The transpose of a matrix $A$ of size $m\times n$ is $A^T$ with

$$
(A^{T})_{ij}=A_{ji}
$$

### Matrix multiplication

Two matrices, $A$ with size $m\times n$ and $B$ with size $n\times p$, can be
multiplied to give $C=AB$ with size $m\times p$:

$$
C_{ij}=\sum_{k=1}^{n}A_{ik}B_{kj}
$$

```{figure} ../images/matrix_multiplication.svg
:name: fig:matmul
Matrix multiplication as a combination of rows and columns.
```

In general $AB\ne BA$, even when both products are defined.

### Square matrix

When both dimensions have the same size $N$, the matrix is called square.

### Diagonal matrix

If the elements $a_{ij}$ of a square matrix are zero for all $i\ne j$, $A$ is a
diagonal matrix. The product of two diagonal matrices commutes,
$AB=BA$. A special diagonal matrix is the unit (identity) matrix $I$, the
matrix version of the number 1:

$$
AI = IA = A
$$

### Trace of a square matrix

The trace of a square matrix is the sum of its diagonal entries,

$$
\mathrm{Tr}\,A=\sum_{i=1}^{n}A_{ii},
$$

which is linear, $\mathrm{Tr}(\alpha A+\beta B)=\alpha\,\mathrm{Tr}A+\beta\,\mathrm{Tr}B$.
Most importantly for accelerator physics, the trace is invariant under cyclic
permutations:

$$
\mathrm{Tr}(A_1A_2\cdots A_n)=\mathrm{Tr}(A_nA_1A_2\cdots A_{n-1})=\cdots=\mathrm{Tr}(A_2A_3\cdots A_nA_1)
$$

### Determinant of a square matrix

A square matrix has a determinant. If $\det A\ne 0$ the matrix is invertible
with inverse $A^{-1}$ satisfying $AA^{-1}=A^{-1}A=I$. For a $2\times2$ matrix,

$$
\det\begin{pmatrix}a & b\\ c & d\end{pmatrix}=ad-bc
$$

Useful properties on an $n\times n$ matrix:

$$
\det A=\det A^T,\qquad
\det(kA)=k^n\det A,\qquad
\det(AB)=\det A \det B
$$

The condition $\det M = 1$ for transfer matrices will return as the
**symplectic condition** in the Transverse Dynamics chapter.

### Eigensystems of a square matrix

A number $\lambda$ and a non-zero vector $X$ satisfying

$$
AX=\lambda X
\quad\Leftrightarrow\quad
(A-\lambda I)X=0
$$

require $\det(A-\lambda I)=0$, which is an $n^{\text{th}}$-order polynomial in
$\lambda$. Its roots $\lambda_i$ are the eigenvalues. If there are $N$
independent eigenvectors $X_i$, the eigenvalue decomposition reads

$$
Q^{-1}AQ=\Lambda,\qquad
Q=(X_1\ X_2\ \cdots\ X_N),\quad
\Lambda=\operatorname{diag}(\lambda_1,\dots,\lambda_N)
$$

## Numerical examples in Julia

Eigenvalues of a one-turn map determine stability. A matrix is easiest to
understand as a picture: it sends the unit circle to an ellipse, and the only
directions it does not also rotate are its eigenvectors, which it simply
stretches by the corresponding eigenvalue. Drag the off-diagonal entry of a
symmetric $2\times2$ matrix and watch the ellipse and eigenvectors respond
together:

```{code-cell} julia
:tags: [hide-input]

using LinearAlgebra, TrackPadWidgets

θc = range(0, 2π, length=121)
cxs, cys = cos.(θc), sin.(θc)

explorer(
    title   = "A 2×2 matrix as a linear map",
    sliders = [Knob("off-diagonal b", range(-2.5, 2.5, length=41);
                    fmt = b -> string(round(b; digits=2)), init = 24)],
    panels  = [Panel(xlabel="x", ylabel="y", title="unit circle → A·(circle)",
                     xlim=(-4.0,4.0), ylim=(-4.0,4.0), equal=true, height=280, legend=:bottomleft)],
    statics = [line(cxs, cys; color="#9aa4b2", dash=true, label="unit circle")],
    note = "The eigenvectors are the only directions the map does not rotate — it just "*
           "stretches them by their eigenvalue. Every other radius of the circle both "*
           "rotates and stretches, tracing the ellipse.",
) do b
    A = [2.0 b; b 3.0]
    pts = [A*[cxs[i], cys[i]] for i in eachindex(cxs)]
    ex = [p[1] for p in pts]; ey = [p[2] for p in pts]
    λ, X = eigen(A)
    arrows = Any[]
    for k in 1:2
        v = λ[k]*X[:,k]
        push!(arrows, line([0.0, v[1]], [0.0, v[2]]; color=PALETTE[k+2], width=2.5,
                           label="λ$k · eigenvector $k"))
    end
    (series = vcat([line(ex, ey; color=PALETTE[1], label="A·(circle)")], arrows),
     readouts = ["A" => "[2 $(round(b;digits=2)); $(round(b;digits=2)) 3]",
                 "λ₁" => round(λ[1]; digits=3),
                 "λ₂" => round(λ[2]; digits=3),
                 "Tr A" => round(tr(A); digits=3),
                 "det A" => round(det(A); digits=3)])
end
```

The same picture decides stability once a matrix is applied over and over. If
$\det M=1$ its two eigenvalues satisfy $\lambda_1\lambda_2=1$, so they are
either a complex-conjugate pair on the unit circle or a reciprocal pair on the
real axis. The characteristic equation $\lambda^2-\mathrm{Tr}(M)\lambda+1=0$
puts them on the unit circle exactly when $|\mathrm{Tr}\,M|<2$, and then
$M^{k}$ stays bounded for all $k$ — this is the long-term stability condition we
will derive for periodic lattices in the Transverse Dynamics chapter.
