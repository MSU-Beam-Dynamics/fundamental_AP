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

# The three 2×2 maps of the text, written out as plain matrices acting on (x, x′).
drift_matrix(L) = [1.0  L
                   0.0  1.0]

sector_bend_matrix(L, ρ) = [ cos(L/ρ)       ρ*sin(L/ρ)
                            -sin(L/ρ)/ρ     cos(L/ρ)]

"Thick quadrupole: trigonometric when focusing (k > 0), hyperbolic when not."
function quadrupole_matrix(L, k)
    k == 0 && return drift_matrix(L)
    q = sqrt(abs(k))
    return k > 0 ? [ cos(q*L)      sin(q*L)/q
                    -q*sin(q*L)    cos(q*L)] :
                   [ cosh(q*L)     sinh(q*L)/q
                     q*sinh(q*L)   cosh(q*L)]
end

# The unit circle in (x, x′), the set of starting conditions we map.
angles   = range(0, 2π, length=91)
circle_x = cos.(angles)
circle_y = sin.(angles)

# Each element reads the knob differently, so each says how to turn the knob
# value `u` into a matrix and into a label.
ELEMENTS = [
    (name  = "drift",
     build = u -> drift_matrix(u),
     label = u -> "drift(L=$(round(u; digits=2)))"),
    (name  = "dipole",
     build = u -> sector_bend_matrix(u, 1.5),
     label = u -> "dipole(L=$(round(u; digits=2)), ρ=1.5)"),
    (name  = "quadrupole",
     # The knob spans 0.1 … 3, remapped to k = −2 … +10 so one slider sweeps
     # from defocusing, through a drift at k = 0, to strongly focusing.
     build = u -> quadrupole_matrix(0.3, 4u - 2),
     label = u -> "quadrupole(L=0.3, k=$(round(4u - 2; digits=2)))"),
]

explorer(
    title   = "Drift, dipole and quadrupole are all just 2×2 matrices",
    sliders = [Knob("element", eachindex(ELEMENTS); fmt = i -> ELEMENTS[i].name),
               Knob("length / strength", range(0.1, 3.0, length=25);
                    fmt = u -> string(round(u; digits=2)), init = 13)],
    panels  = [Panel(xlabel="x", ylabel="x′", title="unit circle → M·(circle)",
                     xlim=(-3.2,3.2), ylim=(-3.2,3.2), equal=true, height=300)],
    statics = [line(circle_x, circle_y; color="#9aa4b2", dash=true, label="unit circle")],
    note = "det M = 1 for all three — every one of these elements is symplectic.",
) do i, u
    element = ELEMENTS[i]
    M       = element.build(u)

    image     = [M * [circle_x[j], circle_y[j]] for j in eachindex(circle_x)]
    ellipse_x = [p[1] for p in image]
    ellipse_y = [p[2] for p in image]

    (series = [line(ellipse_x, ellipse_y; color=PALETTE[i], label=element.label(u))],
     readouts = ["m₁₁"   => round(M[1,1]; digits=3),
                 "m₁₂"   => round(M[1,2]; digits=3),
                 "m₂₁"   => round(M[2,1]; digits=3),
                 "m₂₂"   => round(M[2,2]; digits=3),
                 "det M" => round(det(M); digits=6)])
end
```

### Tracking through a drift
Now reproduce the same physics with TrackPad's tracking engine: launch an ensemble of particles through a drift, watch the phase-space distribution shear:

```{code-cell} julia
:tags: [hide-input]

using StaticArrays, TrackPad, CairoMakie, Random, Statistics, TrackPadWidgets

# Reference beam for every tracking cell and widget on this page. Defined at top
# level, in the first cell that needs it, so later cells share it and the page
# keeps working whatever order the cells are executed in.
beam = Beam(3.0e9)                   # 3 GeV kinetic energy, electron by default

# Column layout of a TrackPad coordinate array: (x, pₓ, y, p_y, z, δE).
const IX, IPX, IY, IPY = 1, 2, 3, 4

"""
    track_and_record(pieces, beam, coords0)

Track the bunch `coords0` through `pieces` one element at a time, saving the
full coordinate array after every element.

Returns `(s, states)`: `s[i]` is the distance from the start of the line and
`states[i]` is a copy of the `nparticle × 6` array there. Cutting a magnet into
several short `pieces` is what gives a smooth curve *inside* the magnet —
`linepass!` on its own only reports the exit.
"""
function track_and_record(pieces, beam, coords0)
    coords = copy(coords0)
    lost   = zeros(Int, size(coords, 1))
    s      = [0.0]
    states = [copy(coords)]
    for element in pieces
        linepass!(coords, Lattice(AbstractElement[element]), beam, lost)
        push!(s, s[end] + get_length(element))
        push!(states, copy(coords))
    end
    return s, states
end

"Coordinate `col` of particle `i` at every recorded station."
trajectory(states, i, col) = [state[i, col] for state in states]

"Round plot data to 4 significant digits — finer than a screen pixel, and it
 keeps the data embedded in the page small."
plotdata(v) = round.(v; sigdigits=4)

# ---------------------------------------------------------------------------
# The shear demo: one bunch, one drift, nothing else.
# ---------------------------------------------------------------------------
const N_SHEAR = 1000                        # particles

entrance_optics = optics4DUC(1.0, 0.0, 1.0, 0.0)     # βₓ, αₓ, βᵧ, αᵧ
bunch_at_entrance = matched_gaussian(
    MersenneTwister(1234), N_SHEAR, entrance_optics;
    emitx = 20e-9,
    emity = 20e-9,
    emitz = 1e-6,
    betaz = 0.2,
)

explorer(
    title   = "A drift shears the phase-space distribution",
    sliders = [Knob("drift length [m]", 1.0:1.0:5.0;
                    fmt = L -> string(round(Int, L)), init = 2)],
    panels  = [Panel(xlabel="x [mm]", ylabel="pₓ [mrad]", title="Phase space",
                     height=380, legend=:bottomright)],
    # The entrance cloud does not depend on the knob, so it is drawn once.
    statics = [points(bunch_at_entrance[:, IX]  .* 1e3,
                      bunch_at_entrance[:, IPX] .* 1e3;
                      color="#808000", size=2.2, alpha=0.3,
                      label="initial distribution")],
    note = "x grows by x′·L while pₓ is untouched: the cloud shears, but its " *
           "area — the emittance — never changes.",
) do L
    bunch = copy(bunch_at_entrance)
    linepass!(bunch, Lattice(AbstractElement[Drift(L)]), beam, zeros(Int, N_SHEAR))

    x  = bunch[:, IX]
    px = bunch[:, IPX]

    # rms emittance ε = √(⟨x²⟩⟨pₓ²⟩ − ⟨xpₓ⟩²) — the area of the cloud, which the
    # shear cannot change.
    emittance = sqrt(var(x)*var(px) - cov(x, px)^2)

    (series = [points(x .* 1e3, px .* 1e3; color="#1e90ff", size=2.2, alpha=0.3,
                      label="after $(round(Int, L)) m")],
     readouts = ["σₓ"            => string(round(std(x)  * 1e3; sigdigits=3), " mm"),
                 "σₚₓ"           => string(round(std(px) * 1e3; sigdigits=3), " mrad"),
                 "⟨x pₓ⟩"        => string(round(cov(x, px) * 1e6; sigdigits=3), " mm·mrad"),
                 "rms emittance" => string(round(emittance * 1e9; sigdigits=4), " nm·rad")])
end
```
Now we can see how particle envelope evolves in a drift space.  You will note how different initial condition matters in their evolution:

```{code-cell} julia
:tags: [hide-input]

using Random, Statistics, StaticArrays, TrackPad, TrackPadWidgets

const DRIFT_LENGTH = 2.0                       # [m]
const N_PARTICLES  = 300                       # per frame
const SIGMA_START  = sqrt(100.0 * 20e-9)       # 1.414 mm, the entrance rms size

# β₀ spans two decades geometrically; α₀ runs from strongly converging to a waist
# at the entrance.
BETA_CHOICES  = 10.0 .^ range(2, 0, length=5)        # 100 → 1 m
ALPHA_CHOICES = collect(range(100.0, 0.0, length=5)) # 100 → 0

# Where the rms envelope is evaluated. A drift is linear, so the envelope can be
# computed analytically from the entrance cloud — no tracking needed.
envelope_s = collect(range(0.0, DRIFT_LENGTH, length=41))

explorer(
    title   = "A drift acting on a matched bunch of fixed entrance size",
    sliders = [Knob("β₀ [m]", BETA_CHOICES;  fmt = b -> string(round(b; sigdigits=3)), init = 1),
               Knob("α₀",     ALPHA_CHOICES; fmt = a -> string(round(a; digits=1)),    init = 1)],
    panels  = [Panel(xlabel="s [m]", ylabel="x [mm]", title="trajectories through the drift",
                     autoscale=:frame, height=270, legend=:bottomleft),
               Panel(xlabel="x [mm]", ylabel="pₓ [mrad]", title="phase space",
                     autoscale=:frame, height=270, legend=:bottomleft)],
    note = "The meaning of the beta and alpha functions is explained in the next chapter. "*
           "Here the entrance beam size is held fixed, so changing β₀ changes the beam "*
           "quality — the emittance — rather than the spot.",
) do β₀, α₀
    # Pinning the entrance size means the emittance is no longer free: it follows
    # from σ₀ = √(β₀ε).
    emittance = SIGMA_START^2 / β₀

    entrance = optics4DUC(β₀, α₀, β₀, α₀)      # the two planes are set up identically
    bunch_in = matched_gaussian(MersenneTwister(1234), N_PARTICLES, entrance;
                                emitx=emittance, emity=emittance,
                                emitz=1e-9, betaz=0.2)

    bunch_out = copy(bunch_in)
    linepass!(bunch_out, Lattice(AbstractElement[Drift(DRIFT_LENGTH)]),
              beam, zeros(Int, N_PARTICLES))

    # In a drift x(s) = x₀ + s·pₓ₀, so each particle's path is a straight line and
    # two points draw it exactly. `nothing` lifts the pen between particles.
    fan_s = Union{Float64,Nothing}[]
    fan_x = Union{Float64,Nothing}[]
    for i in 1:N_PARTICLES
        push!(fan_s, 0.0);           push!(fan_x, bunch_in[i,  IX] * 1e3)
        push!(fan_s, DRIFT_LENGTH);  push!(fan_x, bunch_out[i, IX] * 1e3)
        push!(fan_s, nothing);       push!(fan_x, nothing)
    end

    # The envelope, measured from the ensemble rather than from optics functions:
    # propagate every particle to s by hand and take the standard deviation.
    envelope = [std(bunch_in[:, IX] .+ s .* bunch_in[:, IPX]) * 1e3 for s in envelope_s]
    waist    = argmin(envelope)

    # Closed forms to check the measurement against, both consequences of holding
    # σ₀ fixed: the waist depth depends only on α₀, its position on both knobs.
    waist_size_formula = SIGMA_START * 1e3 / sqrt(1 + α₀^2)
    waist_s_formula    = α₀ * β₀ / (1 + α₀^2)
    divergence         = SIGMA_START * 1e3 * sqrt(1 + α₀^2) / β₀

    (series = [line(fan_s, fan_x; panel=1, color=PALETTE[1], width=0.8, alpha=0.18,
                    label="particles"),
               line(envelope_s,  envelope; panel=1, color=PALETTE[4], width=2.0,
                    label="±σₓ(s)"),
               line(envelope_s, -envelope; panel=1, color=PALETTE[4], width=2.0),
               points(bunch_in[:, IX]  .* 1e3, bunch_in[:, IPX]  .* 1e3; panel=2,
                      color=PALETTE[1], size=2.2, alpha=0.45, label="entrance"),
               points(bunch_out[:, IX] .* 1e3, bunch_out[:, IPX] .* 1e3; panel=2,
                      color=PALETTE[2], size=2.2, alpha=0.45,
                      label="after $(DRIFT_LENGTH) m")],
     readouts = ["ε = σ₀²/β₀"     => string(round(emittance; sigdigits=3), " m·rad"),
                 "σₓ at entrance" => string(round(envelope[1];     sigdigits=3), " mm"),
                 "σₓ at exit"     => string(round(envelope[end];   sigdigits=3), " mm"),
                 "waist σₓ"       => string(round(envelope[waist]; sigdigits=3), " mm"),
                 "σ₀/√(1+α₀²)"    => string(round(waist_size_formula; sigdigits=3), " mm"),
                 "waist at s"     => string(round(envelope_s[waist]; digits=2), " m"),
                 "α₀β₀/(1+α₀²)"   => string(round(waist_s_formula; digits=2), " m"),
                 "σₓ′ = √(γ₀ε)"   => string(round(divergence; sigdigits=3), " mrad")])
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

using TrackPad, TrackPadWidgets
# `beam`, `track_and_record`, `trajectory` and `IX`/`IY` come from the cell above.

const L_QUAD    = 0.3                 # quadrupole length [m]
const S_QUAD    = (0.5, 0.8)          # where it sits along the line, for the marker
const N_SLICES  = 6                   # slices, so the ray bends smoothly inside it

# Nine rays entering parallel to the axis, at the same offset in x and in y, so
# the two panels differ only through the sign of the focusing.
ray_offsets = collect(range(-4e-3, 4e-3, length=9))

"A dashed vertical pair marking where the quadrupole begins and ends."
quad_marker(ylo, yhi) = ([S_QUAD[1], S_QUAD[1], nothing, S_QUAD[2], S_QUAD[2]],
                         [ylo, yhi, nothing, ylo, yhi])

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
    beamline = AbstractElement[
        [Drift(0.1)                  for _ in 1:5]...,
        [Quadrupole(L_QUAD/N_SLICES, k) for _ in 1:N_SLICES]...,
        [Drift(0.15)                 for _ in 1:20]...,
    ]

    rays = zeros(length(ray_offsets), 6)
    rays[:, IX] .= ray_offsets
    rays[:, IY] .= ray_offsets
    s, states = track_and_record(beamline, beam, rays)

    # The whole magnet as ONE thick element, for the exact transfer matrix.
    M = transfer_map(Lattice(AbstractElement[Quadrupole(L_QUAD, k)]), beam)

    horizontal = [line(s, trajectory(states, i, IX) .* 1e3; panel=1,
                       color=PALETTE[1], alpha=0.85, width=1.4)
                  for i in eachindex(ray_offsets)]
    vertical   = [line(s, trajectory(states, i, IY) .* 1e3; panel=2,
                       color=PALETTE[2], alpha=0.85, width=1.4)
                  for i in eachindex(ray_offsets)]

    # A lens takes a parallel ray to the axis after a focal length f; for a
    # transfer matrix that is x′_out = −x_in/f, i.e. f = −1/M₂₁.
    thin_lens_f  = k == 0    ? "∞" : string(round(1/(k*L_QUAD); sigdigits=4), " m")
    thick_lens_f = M[2,1] == 0 ? "∞" : string(round(-1/M[2,1];  sigdigits=4), " m")

    (series = vcat(horizontal, vertical,
        [line(quad_marker(-9.0, 9.0)...; panel=1, color="#9aa4b2", dash=true,
              width=1, label="quad"),
         line(quad_marker(-9.0, 9.0)...; panel=2, color="#9aa4b2", dash=true,
              width=1, label="quad")]),
     readouts = ["k₁ℓ [m⁻¹]"             => round(k*L_QUAD; sigdigits=4),
                 "thin lens f = 1/(k₁ℓ)" => thin_lens_f,
                 "thick lens f = −1/M₂₁" => thick_lens_f,
                 "det M (x block)"       => round(M[1,1]*M[2,2] - M[1,2]*M[2,1]; digits=9)])
end
```

Now we can insert a short quadrupole in the middle of a 2 meter drfit space. Launch a 2 mm bunch a metre upstream and collect them at downstream: the trajectory fan and
the rms envelope give the absolute picture, and the three insets show the
phase-space cloud at entrance, at the magnet and at exit.

```{code-cell} julia
:tags: [hide-input]

using Random, Statistics, StaticArrays, TrackPad, TrackPadWidgets
# `beam`, `track_and_record`, `plotdata` and `IX`/`IPX` come from the cells above.

const L_QUAD_SHORT = 0.1          # quadrupole length [m], centred on s = 0
const L_ARM        = 0.95         # drift on each side, so s runs −1 m → +1 m
const SLICES_DRIFT = 10           # slices per drift arm, for a smooth envelope
const SLICES_QUAD  = 4            # slices inside the quadrupole
const N_CLOUD      = 500          # particles in the phase-space clouds
const N_RAYS       = 50           # trajectories drawn individually
const SIGMA_IN     = 2.0e-3       # entrance rms beam size, held at 2 mm

"The line drift – quadrupole – drift, sliced so the optics can be watched inside."
function quad_beamline(k)
    arm = [Drift(L_ARM/SLICES_DRIFT) for _ in 1:SLICES_DRIFT]
    return AbstractElement[
        arm...,
        [Quadrupole(L_QUAD_SHORT/SLICES_QUAD, k) for _ in 1:SLICES_QUAD]...,
        arm...,
    ]
end

"Index of the recorded station closest to `s_target`."
nearest_station(s, s_target) = argmin(abs.(s .- s_target))

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
) do β₀, α₀, k
    # Pinning σₓ(−1 m) = 2 mm fixes the emittance, as in the drift example above.
    emittance = SIGMA_IN^2 / β₀
    bunch_in  = matched_gaussian(MersenneTwister(1234), N_CLOUD,
                                 optics4DUC(β₀, α₀, β₀, α₀);
                                 emitx=emittance, emity=emittance,
                                 emitz=1e-9, betaz=0.2)

    s_raw, states = track_and_record(quad_beamline(k), beam, bunch_in)
    s = s_raw .- (L_ARM + L_QUAD_SHORT/2)      # put the quadrupole centre at s = 0

    # The three stations the insets show, found by position rather than by
    # counting slices — so the geometry constants above can change freely.
    i_entrance = nearest_station(s, -1.0)
    i_centre   = nearest_station(s,  0.0)
    i_exit     = nearest_station(s, +1.0)

    # A ray is straight in each drift, so it is drawn exactly by the stations
    # where the geometry changes. The envelope is a hyperbola and needs them all.
    ray_stations = unique([1, nearest_station(s, -L_QUAD_SHORT/2), i_centre,
                           nearest_station(s, +L_QUAD_SHORT/2), length(s)])
    ray_s = plotdata(s[ray_stations])
    rays  = [line(ray_s, plotdata([states[j][i, IX] * 1e3 for j in ray_stations]);
                  panel=1, color=PALETTE[1], alpha=0.4, width=1.0)
             for i in 1:N_RAYS]

    envelope = [std(state[:, IX]) for state in states] .* 1e3

    insets = [points(plotdata(states[j][:, IX]  .* 1e3),
                     plotdata(states[j][:, IPX] .* 1e3);
                     panel=panel, color=PALETTE[panel - 1], size=2.6, alpha=0.5)
              for (panel, j) in zip(2:4, (i_entrance, i_centre, i_exit))]

    quad_marker = ([-L_QUAD_SHORT/2, -L_QUAD_SHORT/2, nothing,
                     L_QUAD_SHORT/2,  L_QUAD_SHORT/2],
                   [-12.0, 12.0, nothing, -12.0, 12.0])

    (series = vcat(rays,
        [line(plotdata(s), plotdata(envelope);  panel=1, color=PALETTE[4], width=2.0,
              label="±σₓ(s)"),
         line(plotdata(s), plotdata(-envelope); panel=1, color=PALETTE[4], width=2.0),
         line(quad_marker...; panel=1, color="#9aa4b2", dash=true, width=1,
              label="quadrupole")],
        insets),
     readouts = ["ε = σ₀²/β₀"  => string(round(emittance * 1e9; sigdigits=3), " nm·rad"),
                 "σₓ′(−1 m)"   => string(round(std(states[i_entrance][:, IPX]) * 1e3;
                                               sigdigits=3), " mrad"),
                 "k₁ℓ [m⁻¹]"   => string(round(k * L_QUAD_SHORT; sigdigits=3)),
                 "f = 1/(k₁ℓ)" => k == 0 ? "∞" :
                                  string(round(1/(k * L_QUAD_SHORT); sigdigits=3), " m"),
                 "σₓ(0)"       => string(round(envelope[i_centre]; sigdigits=3), " mm"),
                 "σₓ(+1 m)"    => string(round(envelope[i_exit];   sigdigits=3), " mm"),
                 "waist σₓ"    => string(round(minimum(envelope);  sigdigits=3), " mm"),
                 "waist at s"  => string(round(s[argmin(envelope)]; digits=2), " m")])
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
# `beam`, `track_and_record`, `plotdata`, `nearest_station` and `IX`/`IPX` all
# come from the cells above.

const L_QUAD_FODO  = 0.2          # quadrupole length [m]
const S_QF         = -1.0         # focusing quadrupole centre [m]
const S_QD         = +1.0         # defocusing quadrupole centre [m]
const CELL_HALF    = 2.0          # the cell runs −2 m → +2 m
const N_CLOUD_FODO = 200          # particles in the phase-space clouds
const N_RAYS_FODO  = 10           # trajectories drawn individually
const SIGMA_IN_FODO = 1.5e-3      # entrance rms beam size, held at 1.5 mm

"""
    fodo_beamline(k)

The cell drift – QF – drift – QD – drift, spanning −2 m → +2 m with the two
quadrupoles centred at ∓1 m. Every piece is sliced so the envelope is smooth;
the slice counts only affect resolution, never the optics.
"""
function fodo_beamline(k)
    outer_drift  = (CELL_HALF + S_QF) - L_QUAD_FODO/2       # 0.9 m
    middle_drift = (S_QD - S_QF) - L_QUAD_FODO              # 1.8 m
    return AbstractElement[
        [Drift(outer_drift/6)               for _ in 1:6]...,
        [Quadrupole(L_QUAD_FODO/4, +k)      for _ in 1:4]...,
        [Drift(middle_drift/12)             for _ in 1:12]...,
        [Quadrupole(L_QUAD_FODO/4, -k)      for _ in 1:4]...,
        [Drift(outer_drift/6)               for _ in 1:6]...,
    ]
end

"Dashed verticals at the four quadrupole faces."
function fodo_markers(ylo, yhi)
    faces = [S_QF - L_QUAD_FODO/2, S_QF + L_QUAD_FODO/2,
             S_QD - L_QUAD_FODO/2, S_QD + L_QUAD_FODO/2]
    xs = Union{Float64,Nothing}[]
    ys = Union{Float64,Nothing}[]
    for (i, f) in enumerate(faces)
        i > 1 && (push!(xs, nothing); push!(ys, nothing))
        push!(xs, f, f)
        push!(ys, ylo, yhi)
    end
    return xs, ys
end

# |Tr M| = 2 exactly is the boundary, where the motion grows linearly rather than
# staying bounded — a pure drift (k₁ = 0) sits there, so it gets its own name
# instead of being rounded into "stable".
function stability_verdict(trace)
    trace < 2 - 1e-9  && return "stable"
    trace <= 2 + 1e-9 && return "marginal (|Tr M| = 2)"
    return "unstable"
end

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
) do β₀, α₀, k
    emittance = SIGMA_IN_FODO^2 / β₀          # holds σₓ(−2 m) = 1.5 mm
    bunch_in  = matched_gaussian(MersenneTwister(1234), N_CLOUD_FODO,
                                 optics4DUC(β₀, α₀, β₀, α₀);
                                 emitx=emittance, emity=emittance,
                                 emitz=1e-9, betaz=0.2)

    s_raw, states = track_and_record(fodo_beamline(k), beam, bunch_in)
    s = s_raw .- CELL_HALF                    # cell centre at s = 0

    # The four stations the insets show, found by position.
    i_start = nearest_station(s, -CELL_HALF)
    i_qf    = nearest_station(s, S_QF)
    i_qd    = nearest_station(s, S_QD)
    i_end   = nearest_station(s, +CELL_HALF)

    # A ray is straight in every drift, so the quadrupole faces plus the two ends
    # draw it exactly; the envelope, a hyperbola between the magnets, needs all.
    ray_stations = unique([1,
                           nearest_station(s, S_QF - L_QUAD_FODO/2), i_qf,
                           nearest_station(s, S_QF + L_QUAD_FODO/2),
                           nearest_station(s, S_QD - L_QUAD_FODO/2), i_qd,
                           nearest_station(s, S_QD + L_QUAD_FODO/2),
                           length(s)])
    ray_s = plotdata(s[ray_stations])
    rays  = [line(ray_s, plotdata([states[j][i, IX] * 1e3 for j in ray_stations]);
                  panel=1, color=PALETTE[1], alpha=0.35, width=1.0)
             for i in 1:N_RAYS_FODO]

    envelope = [std(state[:, IX]) for state in states] .* 1e3

    insets = [points(plotdata(states[j][:, IX]  .* 1e3),
                     plotdata(states[j][:, IPX] .* 1e3);
                     panel=panel, color=PALETTE[panel - 1], size=2.4, alpha=0.45)
              for (panel, j) in zip(2:5, (i_start, i_qf, i_qd, i_end))]

    # The one-turn trace of the whole 4 m cell, in each plane separately.
    M = transfer_map(Lattice(fodo_beamline(k)), beam)
    trace_x = abs(M[1,1] + M[2,2])
    trace_y = abs(M[3,3] + M[4,4])

    (series = vcat(rays,
        [line(plotdata(s), plotdata(envelope);  panel=1, color=PALETTE[4], width=2.0,
              label="±σₓ(s)"),
         line(plotdata(s), plotdata(-envelope); panel=1, color=PALETTE[4], width=2.0),
         line(fodo_markers(-12.0, 12.0)...; panel=1, color="#9aa4b2", dash=true,
              width=1, label="quadrupoles")],
        insets),
     readouts = ["ε = σ₀²/β₀"  => string(round(emittance * 1e9; sigdigits=3), " nm·rad"),
                 "k₁ℓ [m⁻¹]"   => string(round(k * L_QUAD_FODO; sigdigits=3)),
                 "σₓ(−2 m)"    => string(round(envelope[i_start]; sigdigits=3), " mm"),
                 "σₓ at QF"    => string(round(envelope[i_qf];    sigdigits=3), " mm"),
                 "σₓ at QD"    => string(round(envelope[i_qd];    sigdigits=3), " mm"),
                 "σₓ(+2 m)"    => string(round(envelope[i_end];   sigdigits=3), " mm"),
                 "|Tr Mₓ|"     => string(round(trace_x; digits=3)),
                 "|Tr Mᵧ|"     => string(round(trace_y; digits=3)),
                 "if repeated" => stability_verdict(max(trace_x, trace_y))])
end
```

In next section, we will see how to anaylize linear lattice and their properties.