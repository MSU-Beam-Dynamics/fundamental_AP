---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Linear Transverse Motion and Transfer Matrices

## Motion in a quadrupole

Consider a normal quadrupole with field $B_x=Gy$, $B_y=Gx$, where
$G=B_0b_1$ is the gradient. The quadrupole is positioned so its zero-field axis
lies on the ideal trajectory — the ideal particle feels no field.

From the Lorentz force, $\gamma m \frac{d\mathbf{v}}{dt} = q \mathbf{v}\times \mathbf{B}$.
Expanding the cross product,

$$
\gamma m\frac{d^{2}x}{dt^{2}}=-qGxv_{z},\qquad
\gamma m\frac{d^{2}y}{dt^{2}}=qGyv_{z}
$$

Using the **paraxial approximation** ($v_x,v_y \ll |\mathbf v|$, hence
$v_z\simeq v$) we replace $\frac{d}{dt}=v_s\frac{d}{ds}\sim v\frac{d}{ds}$ and
obtain

$$
\frac{d^{2}x}{ds^{2}}=-\frac{G}{B\rho}x,
\qquad
\frac{d^{2}y}{ds^{2}}=\frac{G}{B\rho}y .
$$

With the normalized quadrupole strength $k=G/(B\rho)$:

$$
x''+kx=0,\qquad y''-ky=0
$$

(prime denotes $d/ds$ from now on). Depending on the sign of $k$:

$$
x(s)=
\begin{cases}
a \cos (\sqrt{k}s)+b \sin (\sqrt{k}s)&k>0\\
a s + b &k=0\\
a \cosh (\sqrt{-k}s)+b\sinh (\sqrt{-k}s)&k<0
\end{cases}
$$

The three cases correspond to a focusing quadrupole, a drift space and a
defocusing quadrupole: a quadrupole always focuses one plane while defocusing
the other.

## Motion in a dipole

```{figure} ../images/dipole_focusing.png
:width: 300px
:name: fig:dipolefocusing
Sketch of the horizontal focusing effect of a dipole.
```

A dipole only bends the beam; however, the choice of the Frenet–Serret system
creates geometric focusing. A particle injected a distance $x>0$ off the ideal
trajectory travels longer in the dipole and bends more. For small bending
angle $\theta$:

$$
\Delta x'= -\frac{x}{\rho} \theta = -\frac{x}{\rho^2} \Delta s,
\qquad \Delta y' = 0
$$

so inside a dipole

$$
x''+\frac{x}{\rho^2}=0,\qquad y''=0 .
$$

## Hill's equation

Combining dipole and quadrupole effects gives **Hill's equation**:

$$
x''+\left(\frac{1}{\rho^2(s)}  + k(s)\right)x =0,
\qquad
y''-k(s)\, y=0
$$

In general $\rho$ and $k$ are functions of the longitudinal coordinate $s$.
If they are piecewise constant we can solve section by section.  This is exactly how
tracking codes integrate piecewise-hard-edge lattices.

## Transfer matrices

The solution of Hill's equation in each element can be written in matrix
form — the **betatron transfer matrix** $M$ acting on $(x, x')$:

$$
\begin{pmatrix}x\\x'\end{pmatrix}_{\text{exit}}
=M\,\begin{pmatrix}x\\x'\end{pmatrix}_{\text{entrance}},
\qquad
M=\begin{pmatrix}m_{11} & m_{12}\\ m_{21} & m_{22}\end{pmatrix}
$$

### Drift space ($k=0$, $1/\rho=0$)

$$
M_{\text{drift}}(l)=\begin{pmatrix}1 & l\\ 0 & 1\end{pmatrix}
$$

### Dipole (horizontal plane)

$$
M_{\text{dipole}}(l)=\begin{pmatrix}
\cos(l/\rho) & \rho\sin(l/\rho)\\
-\sin(l/\rho)/\rho & \cos(l/\rho)
\end{pmatrix}
\;\xrightarrow{\;l/\rho\ll1\;}\;
\begin{pmatrix}1 & l\\ 0 & 1\end{pmatrix}
$$

### Quadrupoles (thick and thin)

$$
M_{\text{quad}}=
\begin{cases}
\begin{pmatrix}
\cos(\sqrt{k}l) & \sin(\sqrt{k}l)/\sqrt{k}\\
-\sqrt{k}\sin(\sqrt{k}l) & \cos(\sqrt{k}l)
\end{pmatrix}&k>0\\[2ex]
\begin{pmatrix}
\cosh(\sqrt{-k}l) & \sinh(\sqrt{-k}l)/\sqrt{-k}\\
\sqrt{-k}\sinh(\sqrt{-k}l) & \cosh(\sqrt{-k}l)
\end{pmatrix}&k<0
\end{cases}
$$

For a short, strong quad the **thin-lens limit** uses the focal length
$f=1/(|k|l)$ as $l\to0$:

$$
M_{\text{thin quad}}=
\begin{cases}
\begin{pmatrix}1 & 0\\ -1/f & 1\end{pmatrix}&\text{focusing}\\[1ex]
\begin{pmatrix}1 & 0\\ 1/f & 1\end{pmatrix}&\text{defocusing}
\end{cases}
$$

### Chain of linear elements

```{figure} ../images/chain_of_elements.png
:width: 400px
:name: fig:chain
Chain of linear elements — a lattice.
```

For a sequence of magnets (a ***lattice***) the total transfer matrix is the
product of individual matrices, in *reverse* order of placement:

$$
M(s_n,s_0)=M(s_n,s_{n-1})\cdots M(s_2,s_1)\,M(s_1,s_0)
$$

### Properties of transfer matrices

**Symplectic condition.** The matrix represents energy-conserving dynamics.
In the $2\times2$ case this reduces to $\det(M)=1$.

**Long-term stability condition.** A periodic chain (cell) with matrix $M$ is
traversed $k$ times, so the total effect is $M^k$. Requiring bounded motion for
$k\to\infty$ means eigenvalues $\lambda$ of $M$ must satisfy $|\lambda|\le1$.
From $\lambda^2-\mathrm{Tr}(M)\lambda+1=0$ (using $\det M=1$), stability demands

$$
\left|\mathrm{Tr}(M)\right|=\left|m_{11}+m_{22}\right|\le 2
$$

## Numerical exploration

### Mapping the unit circle

First, verify the matrix properties directly. All three named elements are
just 2×2 matrices acting on the unit circle in $(x,x')$: pick one, drag its
length or strength, and watch the image ellipse and its determinant:

```{code-cell} julia
:tags: [hide-input]

using LinearAlgebra, TrackPadWidgets

drift(L) = [1.0 L; 0.0 1.0]
dipole(L, ρ) = [cos(L/ρ) ρ*sin(L/ρ); -sin(L/ρ)/ρ cos(L/ρ)]
function quadrupole(L, k)
    if k > 0
        q = sqrt(k); return [cos(q*L) sin(q*L)/q; -q*sin(q*L) cos(q*L)]
    elseif k < 0
        q = sqrt(-k); return [cosh(q*L) sinh(q*L)/q; q*sinh(q*L) cosh(q*L)]
    else
        return drift(L)
    end
end

θc = range(0, 2π, length=91)
cxs, cys = cos.(θc), sin.(θc)
KINDS = [("drift", 1), ("dipole", 2), ("quadrupole", 3)]

explorer(
    title   = "Drift, dipole and quadrupole are all just 2×2 matrices",
    sliders = [Knob("element", 1:3; fmt = i -> KINDS[i][1]),
               Knob("length / strength", range(0.1, 3.0, length=25);
                    fmt = u -> string(round(u; digits=2)), init = 13)],
    panels  = [Panel(xlabel="x", ylabel="x′", title="unit circle → M·(circle)",
                     xlim=(-3.2,3.2), ylim=(-3.2,3.2), equal=true, height=300)],
    statics = [line(cxs, cys; color="#9aa4b2", dash=true, label="unit circle")],
    note = "det M = 1 for all three — every one of these elements is symplectic.",
) do kind, u
    M = kind == 1 ? drift(u) : kind == 2 ? dipole(u, 1.5) : quadrupole(0.3, u*4 - 2)
    pts = [M*[cxs[i], cys[i]] for i in eachindex(cxs)]
    ex = [p[1] for p in pts]; ey = [p[2] for p in pts]
    label = kind == 1 ? "drift(L=$(round(u;digits=2)))" :
            kind == 2 ? "dipole(L=$(round(u;digits=2)), ρ=1.5)" :
                        "quadrupole(L=0.3, k=$(round(u*4-2;digits=2)))"
    (series = [line(ex, ey; color=PALETTE[kind], label=label)],
     readouts = ["m₁₁" => round(M[1,1]; digits=3), "m₁₂" => round(M[1,2]; digits=3),
                 "m₂₁" => round(M[2,1]; digits=3), "m₂₂" => round(M[2,2]; digits=3),
                 "det M" => round(det(M); digits=6)])
end
```

### Tracking through a drift
Now reproduce the same physics with TrackPad's tracking engine: launch an ensemble of particles through a drift, watch the phase-space distribution shear:

```{code-cell} julia
:tags: [hide-input]

using StaticArrays, TrackPad, CairoMakie, Random, Statistics, TrackPadWidgets

# Reference beam for every tracking cell and widget on this page. Defined at top
# level, in the first cell that needs it, so that later cells share it and the
# page keeps working whatever order the cells are executed in.
beam = Beam(3.0e9)

NPSHEAR   = 1000                              # particles in the shear demo
entrance0 = optics4DUC(1.0, 0.0, 1.0, 0.0)
coords0   = matched_gaussian(
    MersenneTwister(1234), NPSHEAR, entrance0;
    emitx=20e-9,
    emity=20e-9,
    emitz=1e-6,
    betaz=0.2,
)

explorer(
    title   = "A drift shears the phase-space distribution",
    sliders = [Knob("drift length [m]", 1.0:1.0:5.0;
                    fmt = L -> string(round(Int, L)), init = 2)],
    panels  = [Panel(xlabel="x [mm]", ylabel="pₓ [mrad]", title="Phase space",
                     height=380, legend=:bottomright)],
    # The entrance cloud does not depend on the knob, so it is drawn once.
    statics = [points(coords0[:, 1] .* 1e3, coords0[:, 2] .* 1e3;
                      color="#808000", size=2.2, alpha=0.3,
                      label="initial distribution")],
    note = "x grows by x′·L while pₓ is untouched: the cloud shears, but its " *
           "area — the emittance — never changes.",
) do L
    c = copy(coords0)
    linepass!(c, Lattice(AbstractElement[Drift(L)]), beam, zeros(Int, NPSHEAR))
    x, px = c[:, 1], c[:, 2]
    ε = sqrt(var(x) * var(px) - cov(x, px)^2)          # rms emittance, invariant
    (series = [points(x .* 1e3, px .* 1e3; color="#1e90ff", size=2.2, alpha=0.3,
                      label="after $(round(Int, L)) m")],
     readouts = ["σₓ"            => string(round(std(x) * 1e3;  sigdigits=3), " mm"),
                 "σₚₓ"           => string(round(std(px) * 1e3; sigdigits=3), " mrad"),
                 "⟨x pₓ⟩"        => string(round(cov(x, px) * 1e6; sigdigits=3), " mm·mrad"),
                 "rms emittance" => string(round(ε * 1e9; sigdigits=4), " nm·rad")])
end
```
Now we can see how particle envelope evolves in a drift space.  You will note how different initial condition matters in their evolution:

```{code-cell} julia
:tags: [hide-input]

using Random, Statistics, StaticArrays, TrackPad, TrackPadWidgets

const Ldrift = 2.0                            # drift length [m]
const NP     = 300                            # particles per frame
const σ0     = sqrt(100.0 * 20e-9)            # 1.414 mm — entrance rms size, held fixed
BETAS  = 10.0 .^ range(2, 0, length=5)        # 100 → 1 m, geometric
ALPHAS = collect(range(100.0, 0.0, length=5)) # 100 → 0
sgrid  = collect(range(0.0, Ldrift, length=41))

explorer(
    title   = "A drift acting on a matched bunch of fixed entrance size",
    sliders = [Knob("β₀ [m]", BETAS; fmt = b -> string(round(b; sigdigits=3)), init = 1),
               Knob("α₀", ALPHAS; fmt = a -> string(round(a; digits=1)), init = 1)],
    panels  = [Panel(xlabel="s [m]", ylabel="x [mm]", title="trajectories through the drift",
                     autoscale=:frame, height=270, legend=:bottomleft),
               Panel(xlabel="x [mm]", ylabel="pₓ [mrad]", title="phase space",
                     autoscale=:frame, height=270, legend=:bottomleft)],
    note = "The meaning of beta and alpha function will be explained later.  The beam size at begining is fixed.  as beta function changes, the 'beam quality' changes accordingly.",
) do β0, α0
    ε = σ0^2 / β0                             # holds σₓ(0) = √(β₀ε) at σ₀ for every β₀
    entrance = optics4DUC(β0, α0, β0, α0)     # (βₓ, αₓ, βᵧ, αᵧ) — the two planes match
    c0 = matched_gaussian(MersenneTwister(1234), NP, entrance;
                          emitx=ε, emity=ε, emitz=1e-9, betaz=0.2)
    cL = copy(c0)
    linepass!(cL, Lattice(AbstractElement[Drift(Ldrift)]), beam, zeros(Int, NP))

    # In a drift x(s) = x₀ + s·pₓ₀, so one segment per particle draws the whole fan.
    tx = Union{Float64,Nothing}[]; ty = Union{Float64,Nothing}[]
    for i in 1:NP
        push!(tx, 0.0);     push!(ty, c0[i,1]*1e3)
        push!(tx, Ldrift);  push!(ty, cL[i,1]*1e3)
        push!(tx, nothing); push!(ty, nothing)      # lift the pen between particles
    end

    # rms envelope measured from the ensemble itself — no optics functions needed
    env = [std(c0[:,1] .+ s .* c0[:,2]) * 1e3 for s in sgrid]
    iw  = argmin(env)

    (series = [line(tx, ty; panel=1, color=PALETTE[1], width=0.8, alpha=0.18,
                    label="particles"),
               line(sgrid,  env; panel=1, color=PALETTE[4], width=2.0, label="±σₓ(s)"),
               line(sgrid, -env; panel=1, color=PALETTE[4], width=2.0),
               points(c0[:,1].*1e3, c0[:,2].*1e3; panel=2, color=PALETTE[1],
                      size=2.2, alpha=0.45, label="entrance"),
               points(cL[:,1].*1e3, cL[:,2].*1e3; panel=2, color=PALETTE[2],
                      size=2.2, alpha=0.45, label="after $(Ldrift) m")],
     readouts = ["ε = σ₀²/β₀"    => string(round(ε; sigdigits=3), " m·rad"),
                 "σₓ at entrance" => string(round(env[1];   sigdigits=3), " mm"),
                 "σₓ at exit"     => string(round(env[end]; sigdigits=3), " mm"),
                 "waist σₓ"       => string(round(env[iw];  sigdigits=3), " mm"),
                 "σ₀/√(1+α₀²)"    => string(round(σ0*1e3/sqrt(1+α0^2); sigdigits=3), " mm"),
                 "waist at s"     => string(round(sgrid[iw]; digits=2), " m"),
                 "α₀β₀/(1+α₀²)"   => string(round(α0*β0/(1+α0^2); digits=2), " m"),
                 "σₓ′ = √(γ₀ε)"   => string(round(σ0*1e3*sqrt(1+α0^2)/β0; sigdigits=3), " mrad")])
end
```

A drift does not change any particle's angle, so every trajectory is a straight line.  But the envelope is not a straight line!

### Tracking through a quadrupole.

A quadrupole tends to rotate the particle distribution — the essence of
alternating-gradient transport. 

Before assembling cells, look at a single quadrupole doing its job. The rays
below enter parallel to the axis, so wherever they cross is by definition the
focal point, and the same magnet that focuses in $x$ defocuses in $y$:

```{code-cell} julia
:tags: [hide-input]

using TrackPadWidgets

"Track a bunch element by element, recording x and y at every boundary."
function track_s(pieces, beam, coords0)
    c = copy(coords0); flags = zeros(Int, size(c, 1))
    S = Float64[0.0]; X = [copy(c[:, 1])]; Y = [copy(c[:, 3])]
    for e in pieces
        linepass!(c, Lattice(AbstractElement[e]), beam, flags)
        push!(S, S[end] + get_length(e))
        push!(X, copy(c[:, 1])); push!(Y, copy(c[:, 3]))
    end
    S, reduce(hcat, X), reduce(hcat, Y)
end

Lq = 0.3
xr = collect(range(-4e-3, 4e-3, length=9))

explorer(
    title   = "A quadrupole focuses one plane and defocuses the other",
    sliders = [Knob("k₁ [m⁻²]", range(-4.0, 4.0, length=41);
                      fmt = k -> string(round(k; digits=2)), init = 31)],
    panels  = [Panel(xlabel="s [m]", ylabel="x [mm]", title="horizontal",
                     ylim=(-9.0, 9.0), height=240, legend=:bottomleft),
               Panel(xlabel="s [m]", ylabel="y [mm]", title="vertical",
                     ylim=(-9.0, 9.0), height=240, legend=:bottomleft)],
    note = "The thin-lens estimate f = 1/(k₁ℓ) and the thick-lens value −1/M₂₁ agree only "*
           "while f is much longer than the magnet itself.",
) do k
    pieces = AbstractElement[Drift(0.1) for _ in 1:5]
    append!(pieces, [Quadrupole(Lq/6, k) for _ in 1:6])
    append!(pieces, [Drift(0.15)         for _ in 1:20])
    c0 = zeros(length(xr), 6); c0[:, 1] .= xr; c0[:, 3] .= xr
    S, X, Y = track_s(pieces, beam, c0)

    M  = transfer_map(Lattice(AbstractElement[Quadrupole(Lq, k)]), beam)
    marker = ([0.5, 0.5, nothing, 0.8, 0.8], [-9.0, 9.0, nothing, -9.0, 9.0])

    (series = vcat(
        [line(S, X[j, :] .* 1e3; panel=1, color=PALETTE[1], alpha=0.85, width=1.4)
         for j in eachindex(xr)],
        [line(S, Y[j, :] .* 1e3; panel=2, color=PALETTE[2], alpha=0.85, width=1.4)
         for j in eachindex(xr)],
        [line(marker...; panel=1, color="#9aa4b2", dash=true, width=1, label="quad"),
         line(marker...; panel=2, color="#9aa4b2", dash=true, width=1, label="quad")]),
     readouts = ["k₁ℓ [m⁻¹]" => round(k*Lq; sigdigits=4),
                 "thin lens f = 1/(k₁ℓ)" => k == 0 ? "∞" : string(round(1/(k*Lq); sigdigits=4), " m"),
                 "thick lens f = −1/M₂₁" => M[2,1] == 0 ? "∞" : string(round(-1/M[2,1]; sigdigits=4), " m"),
                 "det M (x block)" => round(M[1,1]*M[2,2] - M[1,2]*M[2,1]; digits=9)])
end
```

Now we can insert a short dipole in the middle of a 2 meter drfit space. Launch a 2 mm bunch a metre upstream and collect them at downstream: the trajectory fan and
the rms envelope give the absolute picture, and the three insets show the
phase-space cloud at entrance, at the magnet and at exit.

```{code-cell} julia
:tags: [hide-input]

using Random, Statistics, StaticArrays, TrackPad, TrackPadWidgets

"""
Track a bunch through `pieces`, returning the s grid and the full coordinate
array recorded at every element boundary.
"""
function track_states(pieces, beam, coords0)
    c = copy(coords0); flags = zeros(Int, size(c, 1))
    S = Float64[0.0]; C = [copy(c)]
    for e in pieces
        linepass!(c, Lattice(AbstractElement[e]), beam, flags)
        push!(S, S[end] + get_length(e))
        push!(C, copy(c))
    end
    S, C
end

const LQ   = 0.1        # quadrupole length [m], centred on s = 0
const LARM = 0.95       # drift on each side, so s runs −1 m → +1 m
const NDR  = 10         # slices per drift arm (for a smooth envelope)
const NSQ  = 4          # slices inside the quadrupole
const NPQ  = 500        # particles in the cloud
const NRAY = 50         # trajectories drawn individually
const SIGQ = 2.0e-3     # entrance rms beam size, held at 2 mm for every setting

# s = −1 m (entrance), 0 (quad centre) and +1 m (exit) in the recorded grid
const ISNAP = (1, 1 + NDR + NSQ ÷ 2, 1 + 2NDR + NSQ)
# A ray is straight in each drift, so five points draw it exactly; only the
# envelope, a hyperbola, needs the full grid.
const IRAY  = (1, 1 + NDR, ISNAP[2], 1 + NDR + NSQ, 1 + 2NDR + NSQ)

# Plot coordinates are in mm/mrad; four significant digits is far finer than a
# pixel and keeps the data embedded in the page small.
r4(v) = round.(v; sigdigits=4)

explorer(
    title   = "A short quadrupole (ℓ = 10 cm) seen by a 2 mm beam",
    sliders = [Knob("β₀ [m]", [10.0, 5.0, 2.0, 1.0];
                    fmt = b -> string(round(Int, b)), init = 2),
               Knob("α₀", [-4.0, 0.0, 4.0];
                    fmt = a -> string(round(Int, a)), init = 2),
               Knob("k₁ [m⁻²]", [-20.0, -10.0, 0.0, 10.0, 20.0];
                    fmt = k -> string(round(Int, k)), init = 4)],
    panels  = [Panel(xlabel="s [m]", ylabel="x [mm]",
                     title="trajectories through the quadrupole",
                     ylim=(-12.0, 12.0), height=300, legend=:bottomleft,
                     basis="100%"),
               Panel(title="(x, pₓ) at s = −1 m", autoscale=:frame, share="xpx",
                     height=185, basis="30%", minwidth=165, ticklabels=false),
               Panel(title="(x, pₓ) at s = 0", autoscale=:frame, share="xpx",
                     height=185, basis="30%", minwidth=165, ticklabels=false),
               Panel(title="(x, pₓ) at s = +1 m", autoscale=:frame, share="xpx",
                     height=185, basis="30%", minwidth=165, ticklabels=false)],
    note = "The entrance size is pinned at 2 mm, so β₀ trades spot size against "*
           "divergence and α₀ says whether the bunch arrives converging or diverging. "*
           "The three insets share one x/pₓ box, sized to hold every particle at all "*
           "three stations, so the clouds can be compared directly; the box itself is "*
           "redrawn for each setting of the knobs.",
) do β0, α0, k
    ε  = SIGQ^2 / β0                       # holds σₓ(−1 m) = 2 mm for every β₀
    c0 = matched_gaussian(MersenneTwister(1234), NPQ, optics4DUC(β0, α0, β0, α0);
                          emitx=ε, emity=ε, emitz=1e-9, betaz=0.2)

    pieces = AbstractElement[]
    append!(pieces, [Drift(LARM / NDR)  for _ in 1:NDR])
    append!(pieces, [Quadrupole(LQ / NSQ, k) for _ in 1:NSQ])
    append!(pieces, [Drift(LARM / NDR)  for _ in 1:NDR])
    S, C = track_states(pieces, beam, c0)
    s   = S .- 1.0                         # put the quadrupole centre at s = 0
    env = [std(c[:, 1]) for c in C] .* 1e3

    sray = r4([s[j] for j in IRAY])
    rays = [line(sray, r4([C[j][i, 1] * 1e3 for j in IRAY]); panel=1,
                 color=PALETTE[1], alpha=0.4, width=1.0) for i in 1:NRAY]
    insets = [points(r4([c[i, 1] * 1e3 for i in 1:NPQ]), r4([c[i, 2] * 1e3 for i in 1:NPQ]);
                     panel=p, color=PALETTE[p - 1], size=2.6, alpha=0.5)
              for (p, c) in zip(2:4, (C[ISNAP[1]], C[ISNAP[2]], C[ISNAP[3]]))]

    (series = vcat(rays,
        [line(r4(s), r4(env);  panel=1, color=PALETTE[4], width=2.0, label="±σₓ(s)"),
         line(r4(s), r4(-env); panel=1, color=PALETTE[4], width=2.0),
         line([-LQ/2, -LQ/2, nothing, LQ/2, LQ/2], [-12.0, 12.0, nothing, -12.0, 12.0];
              panel=1, color="#9aa4b2", dash=true, width=1, label="quadrupole")],
        insets),
     readouts = ["ε = σ₀²/β₀"  => string(round(ε * 1e9;   sigdigits=3), " nm·rad"),
                 "σₓ′(−1 m)"   => string(round(std(C[ISNAP[1]][:, 2]) * 1e3; sigdigits=3), " mrad"),
                 "k₁ℓ [m⁻¹]"   => string(round(k * LQ; sigdigits=3)),
                 "f = 1/(k₁ℓ)" => k == 0 ? "∞" : string(round(1 / (k * LQ); sigdigits=3), " m"),
                 "σₓ(0)"       => string(round(env[ISNAP[2]]; sigdigits=3), " mm"),
                 "σₓ(+1 m)"    => string(round(env[ISNAP[3]]; sigdigits=3), " mm"),
                 "waist σₓ"    => string(round(minimum(env);  sigdigits=3), " mm"),
                 "waist at s"  => string(round(s[argmin(env)]; digits=2), " m")])
end
```

### Tracking through a cell of quadrupoles

One quadrupole cannot confine both planes at once, but two of opposite sign
can. Put a focusing and a defocusing quadrupole of equal strength at
$s=\pm1$ m and launch a matched bunch at $s=-2$ m. Watching the horizontal
plane alone is enough to see the mechanism: the bunch arrives at each magnet
with a different width, so the same $|k_1|$ produces a different kick at each,
and the pair nets out to focusing.

```{code-cell} julia
:tags: [hide-input]

const FLQ  = 0.2                    # quadrupole length [m]
const FN1, FNQ, FN2 = 6, 4, 12      # slices: outer drift, quadrupole, middle drift
const NPF  = 200                    # particles in the cloud
const NRF  = 10                     # trajectories drawn individually
const SIGF = 1.5e-3                 # entrance rms beam size, held at 1.5 mm

# Boundary indices of the two quadrupoles in the recorded grid, and of the exit.
const IQ1 = (1 + FN1,             1 + FN1 + FNQ)
const IQ2 = (1 + FN1 + FNQ + FN2, 1 + FN1 + 2FNQ + FN2)
const NFP = 1 + 2FN1 + 2FNQ + FN2
# A ray is straight in every drift, so it only needs the quadrupole boundaries
# and the two ends; the envelope, a hyperbola between the magnets, needs them all.
const IFR = (1, IQ1[1]:IQ1[2]..., IQ2[1]:IQ2[2]..., NFP)
# The four phase-space stations: s = −2, −1 (QF centre), +1 (QD centre), +2.
const IST = (1, 1 + FN1 + FNQ ÷ 2, 1 + FN1 + FNQ + FN2 + FNQ ÷ 2, NFP)

"The −2 m → +2 m cell: drift, QF at −1 m, drift, QD at +1 m, drift."
function fodo_pieces(k)
    p = AbstractElement[]
    append!(p, [Drift(0.9 / FN1)  for _ in 1:FN1])
    append!(p, [Quadrupole(FLQ / FNQ, +k) for _ in 1:FNQ])
    append!(p, [Drift(1.8 / FN2)  for _ in 1:FN2])
    append!(p, [Quadrupole(FLQ / FNQ, -k) for _ in 1:FNQ])
    append!(p, [Drift(0.9 / FN1)  for _ in 1:FN1])
    p
end

qmark(ylo, yhi) = ([-1-FLQ/2, -1-FLQ/2, nothing, -1+FLQ/2, -1+FLQ/2, nothing,
                     1-FLQ/2,  1-FLQ/2, nothing,  1+FLQ/2,  1+FLQ/2],
                   [ylo, yhi, nothing, ylo, yhi, nothing,
                    ylo, yhi, nothing, ylo, yhi])

# |Tr M| = 2 exactly is the boundary, where the motion grows linearly rather
# than staying bounded — a pure drift (k₁ = 0) sits there, so name it as its own
# case instead of rounding it into "stable".
stability(t) = t < 2 - 1e-9 ? "stable" :
               t <= 2 + 1e-9 ? "marginal (|Tr M| = 2)" : "unstable"

explorer(
    title   = "A FODO cell: QF at s = −1 m, QD at s = +1 m, equal strength",
    sliders = [Knob("β₀ [m]", [1.0, 2.9, 5.2, 10.0, 20.0];
                    fmt = b -> string(round(b, digits=1)), init = 3),
               Knob("α₀", [-2.0, -1.63, -1.12, 0.0, 1.0];
                    fmt = a -> string(round(a, digits=1)), init = 2),
               Knob("k₁ [m⁻²]", [0.0, 2.0, 4.0];
                    fmt = k -> string(round(k, digits=1)), init = 4)],
    panels  = [Panel(xlabel="s [m]", ylabel="x [mm]",
                     title="horizontal trajectories through the cell",
                     ylim=(-12.0, 12.0), height=300, legend=:bottomleft,
                     basis="100%"),
               Panel(title="(x, x′) at s = −2 m", autoscale=:frame, share="xxp",
                     height=175, basis="23%", minwidth=140, ticklabels=false),
               Panel(title="(x, x′) at QF, s = −1 m", autoscale=:frame, share="xxp",
                     height=175, basis="23%", minwidth=140, ticklabels=false),
               Panel(title="(x, x′) at QD, s = +1 m", autoscale=:frame, share="xxp",
                     height=175, basis="23%", minwidth=140, ticklabels=false),
               Panel(title="(x, x′) at s = +2 m", autoscale=:frame, share="xxp",
                     height=175, basis="23%", minwidth=140, ticklabels=false)],
    note = "The entrance size is pinned at 1.5 mm rms, so the emittance follows from "*
           "the knobs as ε = σ₀²/β₀ and β₀ trades spot size against divergence. "*
           "The four insets share one x/x′ box, wide enough for every particle at "*
           "all four stations, so the clouds can be compared directly. |Tr M| is the "*
           "trace of this 4 m cell — the stability test of the previous section, "*
           "applied as if the cell were repeated forever; it comes out the same in "*
           "both planes, which is why one number decides the FODO cell.",
) do β0, α0, k
    ε  = SIGF^2 / β0                               # holds σₓ(−2 m) = 1.5 mm
    c0 = matched_gaussian(MersenneTwister(1234), NPF, optics4DUC(β0, α0, β0, α0);
                          emitx=ε, emity=ε, emitz=1e-9, betaz=0.2)
    S, C = track_states(fodo_pieces(k), beam, c0)
    s    = S .- 2.0                                # cell centre at s = 0
    env  = [std(c[:, 1]) for c in C] .* 1e3
    sray = r4([s[j] for j in IFR])

    M   = transfer_map(Lattice(fodo_pieces(k)), beam)
    trx = abs(M[1,1] + M[2,2]); try_ = abs(M[3,3] + M[4,4])
    

    rays = [line(sray, r4([C[j][i, 1] * 1e3 for j in IFR]); panel=1,
                 color=PALETTE[1], alpha=0.35, width=1.0) for i in 1:NRF]
    insets = [points(r4([c[i, 1] * 1e3 for i in 1:NPF]), r4([c[i, 2] * 1e3 for i in 1:NPF]);
                     panel=p, color=PALETTE[p - 1], size=2.4, alpha=0.45)
              for (p, c) in zip(2:5, (C[IST[1]], C[IST[2]], C[IST[3]], C[IST[4]]))]

    (series = vcat(rays,
        [line(r4(s), r4(env);  panel=1, color=PALETTE[4], width=2.0, label="±σₓ(s)"),
         line(r4(s), r4(-env); panel=1, color=PALETTE[4], width=2.0),
         line(qmark(-12.0, 12.0)...; panel=1, color="#9aa4b2", dash=true, width=1,
              label="quadrupoles")],
        insets),
     readouts = ["ε = σ₀²/β₀"  => string(round(ε * 1e9; sigdigits=3), " nm·rad"),
                 "k₁ℓ [m⁻¹]"  => string(round(k * FLQ; sigdigits=3)),
                 "σₓ(−2 m)"   => string(round(env[IST[1]]; sigdigits=3), " mm"),
                 "σₓ at QF"   => string(round(env[IST[2]]; sigdigits=3), " mm"),
                 "σₓ at QD"   => string(round(env[IST[3]]; sigdigits=3), " mm"),
                 "σₓ(+2 m)"   => string(round(env[IST[4]]; sigdigits=3), " mm"),
                 "|Tr Mₓ|"    => string(round(trx;  digits=3)),
                 "|Tr Mᵧ|"    => string(round(try_; digits=3)),
                 "if repeated" => stability(max(trx, try_))])
end
```

In next section, we will see how to anaylize linear lattice and their properties.