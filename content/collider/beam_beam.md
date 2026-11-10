---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# The Beam–Beam Effect

## Two beams, two kinds of interaction

At the interaction point, two intense bunches of charged particles pass through
each other. Two kinds of interaction take place:

- the **strong and weak interactions** between individual particles, which are
  rare events and the whole purpose of the collider;
- the **electromagnetic interaction** between each particle and the collective
  field of the other bunch. It acts on every particle at every crossing. It is
  parasitic and impossible to avoid.

The force of the opposing bunch is the **beam–beam force**. The force of a
bunch on its own particles is the **space-charge force**. The beam–beam effect
is one of the most important factors that limit the luminosity of colliders.
Raising the luminosity means more particles in smaller bunches, and that
strengthens exactly this force.

## Fields of a relativistic bunch

The fields are easiest to find in the rest frame of the bunch. There it is a
static charge distribution, with $\mathbf B'=0$ and an electric field from the
Poisson equation:

$$
\nabla^2U=-\frac{\rho(x,y,z)}{\epsilon_0},\qquad \mathbf E'=-\nabla U .
$$

The laboratory fields follow from the Lorentz transformation of the fields
(see [](../math_prep/relativistic_mechanics.md)):

$$
E_\parallel=E'_\parallel,\qquad \mathbf E_\perp=\gamma\,\mathbf E'_\perp,\qquad
\mathbf B=\frac{\boldsymbol\beta\times\mathbf E}{c}.
$$

The transverse field is enhanced by $\gamma$ and the longitudinal one is not.
Equivalently, the field of a fast charge is squeezed into a thin disc
perpendicular to its motion, with an opening angle of about $1/\gamma$. For a
bunch with $\gamma$ in the hundreds or thousands, only the transverse fields
matter, so the beam–beam effect is a problem of **transverse** dynamics.

```{code-cell} julia
:tags: [hide-input]
:label: fig-moving-charge
:caption: Electric field lines of a point charge moving to the right, drawn in a plane through its path with a density proportional to the field strength. As γ grows, the field is compressed into the plane perpendicular to the motion.

using CairoMakie, TrackPadWidgets, TrackPad, Printf

"""
Angles (from the direction of motion) of n field lines in the half plane, spaced so
that their density follows the field strength E(θ) ∝ (1 - β²)/(1 - β² sin²θ)^{3/2}.
"""
function field_line_angles(γ, n)
    β2 = 1 - 1/γ^2
    θ = range(0, π, length=4001)
    E = @. (1 - β2)/(1 - β2*sin(θ)^2)^1.5
    F = cumsum(E); F ./= F[end]                          # fraction of the lines within θ
    [θ[findfirst(>=((k - 0.5)/n), F)] for k in 1:n]
end

let fig = Figure(size=(960, 330), fontsize=15)
    for (j, γ) in enumerate([1.0, 2.0, 6.0])
        ax = Axis(fig[1, j]; aspect=DataAspect(), title=γ == 1 ? "at rest (γ = 1)" : "γ = $(Int(γ))",
                  limits=(-1.1, 1.1, -1.1, 1.1))
        hidedecorations!(ax); hidespines!(ax)
        for θ in field_line_angles(γ, 14), sgn in (1, -1)
            lines!(ax, [0.1cos(θ), cos(θ)], sgn .* [0.1sin(θ), sin(θ)]; color=PALETTE[1], linewidth=1.5)
        end
        scatter!(ax, [0.0], [0.0]; color=PALETTE[4], markersize=14)
        γ > 1 && arrows2d!(ax, [0.0], [0.0], [0.45], [0.0]; color=PALETTE[4], shaftwidth=3, tipwidth=12, tiplength=12)
    end
    fig
end
```

A particle moving with velocity $\mathbf v$ feels the Lorentz force
$\mathbf F=q(\mathbf E+\mathbf v\times\mathbf B)$, and for the bunch's own field
$\mathbf v\times\mathbf B=\boldsymbol\beta_{\mathrm{particle}}\times(\boldsymbol\beta_{\mathrm{bunch}}\times\mathbf E)$.
The magnetic part either adds to the electric part or cancels it:

$$
\text{opposing beam:}\quad F_r=qE_\perp\left(1+\beta^2\right)\approx2qE_\perp,
\qquad
\text{same beam:}\quad F_r=qE_\perp\left(1-\beta^2\right)=\frac{qE_\perp}{\gamma^2}.
$$

For particles in the same bunch, the magnetic attraction of parallel currents
almost cancels the electric repulsion. That is why space charge fades at high
energy. For the opposing beam the two forces add, and the beam–beam force stays
at full strength at any energy.

## The beam–beam kick of a round Gaussian beam

Take a round Gaussian bunch with $n$ charges per unit length:

$$
\rho(x,y)=\frac{n\,e}{2\pi\sigma^2}\exp\!\left(-\frac{x^2+y^2}{2\sigma^2}\right).
$$

Gauss's law over a circle of radius $r$ gives the radial field. The force on an
opposing particle is then

$$
F_r=\frac{n\,e^2\left(1+\beta^2\right)}{2\pi\epsilon_0\,r}\left[1-\exp\!\left(-\frac{r^2}{2\sigma^2}\right)\right].
$$

A particle crossing the bunch receives a kick, the time integral of the force.
For a bunch of $N$ particles, the change of angle is

$$
\Delta r'=\frac{1}{mc\beta\gamma}\int_{-\infty}^{\infty}F_r\,dt
=\pm\frac{2Nr_0}{\gamma}\,\frac{1}{r}\left[1-\exp\!\left(-\frac{r^2}{2\sigma^2}\right)\right],
\qquad r_0=\frac{e^2}{4\pi\epsilon_0mc^2}.
$$

Here $r_0$ is the classical radius of the kicked particle
($1.535\times10^{-18}$ m for protons, $2.818\times10^{-15}$ m for electrons).
The kick points outward (defocusing) for like charges, as in proton–proton
collisions, and inward (focusing) for opposite charges, as in $e^+e^-$ or
$p\bar p$ collisions.

The left panel of the figure shows the shape of this kick:

- **inside the core** ($r\lesssim\sigma$) it is linear in $r$, so it acts like a
  thin quadrupole lens that focuses or defocuses in both planes at once;
- it **peaks** at $r\approx1.6\sigma$;
- **far outside** it falls as $1/r$, like the field of a line charge. This tail
  is what a bunch feels at the *long-range* encounters where the beams pass at a
  distance (see [](interaction_region.md)).

**Non-round beams.** Many colliders, and all $e^+e^-$ factories, use flat beams
with $\sigma_x\gg\sigma_y$. The field of an elliptical Gaussian was found by
Bassetti and Erskine in terms of the complex error function
$w(z)=e^{-z^2}\operatorname{erfc}(-iz)$:

$$
E_y+iE_x=\frac{n\,e}{2\epsilon_0\sqrt{2\pi(\sigma_x^2-\sigma_y^2)}}
\left[w\!\left(\frac{x+iy}{\sqrt{2(\sigma_x^2-\sigma_y^2)}}\right)
-\exp\!\left(-\frac{x^2}{2\sigma_x^2}-\frac{y^2}{2\sigma_y^2}\right)
w\!\left(\frac{x\frac{\sigma_y}{\sigma_x}+iy\frac{\sigma_x}{\sigma_y}}{\sqrt{2(\sigma_x^2-\sigma_y^2)}}\right)\right].
$$

This is the formula behind the beam–beam elements used in the simulations below.
The right panel shows its main consequence. The small dimension carries the
steep field: near the center, the vertical gradient is $\sigma_x/\sigma_y$ times
the horizontal one.

```{code-cell} julia
:tags: [hide-input]
:label: fig-bb-kick
:caption: Left, the beam–beam kick of a round Gaussian beam in units of 2Nr₀/(γσ), with the linear approximation of the core and the 1/r field of a line charge; the shaded curve is the beam's transverse profile. Right, the fields of a flat beam with σₓ = 4σᵧ (Bassetti–Erskine) along the two axes, against distance in units of σᵧ.

let fig = Figure(size=(980, 380), fontsize=15)
    r = range(-6, 6, length=601)
    kick = [gaussian_beam_field(x, 0.0, 1.0, 1.0)[1]/2 for x in r]      # (1/r)(1 - e^{-r²/2}) for σ = 1
    ax = Axis(fig[1, 1]; xlabel="x / σ", ylabel="Δx′  [2Nr₀ / γσ]", title="round beam",
              limits=(-6, 6, -0.75, 0.75))
    band!(ax, r, zeros(length(r)), 0.6 .* exp.(-r.^2/2); color=(PALETTE[1], 0.15))
    lines!(ax, r, r ./ 2; color=:gray55, linestyle=:dash, label="linear: x/2σ")
    rr = filter(x -> abs(x) > 1.4, r)
    lines!(ax, rr, 1 ./ rr; color=:gray55, linestyle=:dot, label="line charge: σ/x")
    lines!(ax, r, kick; color=PALETTE[2], linewidth=3, label="Gaussian beam")
    vlines!(ax, [1.585, -1.585]; color=(PALETTE[2], 0.4), linestyle=:dash)
    text!(ax, 1.65, -0.7; text="peak at 1.6σ", fontsize=13, color=:gray30)
    axislegend(ax; position=:lt, framevisible=false, labelsize=13)

    σx, σy = 4.0, 1.0
    d = range(0, 16, length=401)
    ax = Axis(fig[1, 2]; xlabel="distance from the center [σᵧ]", ylabel="field [arb. units]",
              title="flat beam, σₓ = 4σᵧ", limits=(0, 16, 0, 0.75))
    lines!(ax, d, [gaussian_beam_field(0.0, y, σx, σy)[2] for y in d]; color=PALETTE[4], linewidth=3,
           label="Eᵧ along y")
    lines!(ax, d, [gaussian_beam_field(x, 0.0, σx, σy)[1] for x in d]; color=PALETTE[1], linewidth=3,
           label="Eₓ along x")
    lines!(ax, d[d .< 1.2], 2 .* d[d .< 1.2] ./ (σy*(σx + σy)); color=PALETTE[4], linestyle=:dash)
    lines!(ax, d[d .< 4.5], 2 .* d[d .< 4.5] ./ (σx*(σx + σy)); color=PALETTE[1], linestyle=:dash)
    axislegend(ax; position=:rt, framevisible=false)
    fig
end
```

## The beam–beam parameter

A particle near the center of the opposing bunch sees only the linear part of
the kick, $\Delta x'=\pm Nr_0\,x/(\gamma\sigma^2)$: a thin lens with
$1/f=Nr_0/(\gamma\sigma^2)$ in each plane. A thin lens at a point where the Courant–Snyder function is $\beta$
shifts the tune by $\beta/(4\pi f)$. This tune shift is the **beam–beam
parameter**:

$$
\xi_{x,y}=\frac{N r_0\,\beta^*_{x,y}}{2\pi\gamma\,\sigma_{x,y}\left(\sigma_x+\sigma_y\right)}
\;\xrightarrow{\ \text{round}\ }\;
\xi=\frac{N r_0\,\beta^*}{4\pi\gamma\sigma^2}=\frac{N r_0}{4\pi\epsilon_n}.
$$

The last form, for round beams with normalized emittance
$\epsilon_n\approx\gamma\epsilon$, holds a surprise. The beam–beam parameter
depends only on the **brightness** $N/\epsilon_n$ of the opposing bunch. It does
not depend on its energy or on $\beta^*$. For RHIC's protons, $N=2\times10^{11}$
and $\epsilon_n=2.5\ \mu$m give $\xi\approx0.01$ per IP.

Written in terms of $\xi$, the luminosity of [](luminosity.md) becomes

$$
\mathcal L=\frac{fN^2}{4\pi\sigma^2}=\frac{f\,N\,\gamma\,\xi}{r_0\,\beta^*}.
$$

Once the beams are as bright as the beam–beam effect allows, with $\xi$ at its
limit, more luminosity can come only from more bunches, more particles per bunch
*at the same brightness*, higher energy, or a smaller $\beta^*$. That last point
is the reason for the low-β insertions of [](interaction_region.md).

Experience shows that $\xi$ cannot be raised beyond an empirical
**beam–beam limit**:

| | LHC | RHIC | Tevatron | B factories |
|---|---|---|---|---|
| particles | hadrons | hadrons | hadrons | leptons |
| energy (GeV) | 3500 | 250 | 980 | 3 / 9 |
| beam–beam parameter (total) | 0.03 | 0.015 | 0.03 | 0.16 |

Hadron colliders stay near 0.01 per interaction point. Lepton colliders reach
ten times more, because synchrotron radiation damps the particles back towards
the core faster than the beam–beam force can drive them out. Why the limit sits
where it does is still not fully understood.

## Linear analysis

Treat the beam–beam force as a thin lens at the IP, where $\alpha^*=0$. Writing
$\mu=2\pi\nu$ for the unperturbed phase advance and $1/f=4\pi\xi/\beta^*$
(focusing for $\xi>0$), the one-turn matrix becomes

$$
K\cdot M=
\begin{pmatrix}1&0\\-1/f&1\end{pmatrix}
\begin{pmatrix}\cos\mu&\beta^*\sin\mu\\-\sin\mu/\beta^*&\cos\mu\end{pmatrix}
=\begin{pmatrix}\cos\mu&\beta^*\sin\mu\\[0.5ex]
-\dfrac{\sin\mu}{\beta^*}-\dfrac{\cos\mu}{f}&\cos\mu-\dfrac{\beta^*\sin\mu}{f}\end{pmatrix}.
$$

Half its trace gives the new tune and its upper-right element the new $\beta^*$:

$$
\cos2\pi\nu'=\cos2\pi\nu-2\pi\xi\sin2\pi\nu,
\qquad
\frac{\beta^{*\prime}}{\beta^*}=\frac{\sin2\pi\nu}{\sin2\pi\nu'} .
$$

For small $\xi$, the tune shifts by $\Delta\nu\approx\xi$. The beam–beam
"lens" also changes $\beta^*$ itself, the **dynamic beta**. The motion is stable
as long as $|\cos2\pi\nu'|<1$, which gives

$$
\xi<\begin{cases}\dfrac{1}{2\pi}\cot(\pi\nu), & n<\nu<n+\tfrac12,\\[2ex]
-\dfrac{1}{2\pi}\tan(\pi\nu), & n+\tfrac12<\nu<n+1 .\end{cases}
$$

For like charges, replace $\xi$ by $-\xi$, which mirrors everything about the
half integer. The figure shows that this linear limit is far above the observed
values of $\xi$, except within a few hundredths of an integer or half-integer
tune. The linear picture does **not** explain the beam–beam limit.

```{code-cell} julia
:tags: [hide-input]
:label: fig-bb-linear
:caption: Linear beam–beam theory for a focusing beam–beam force. Left, the largest stable ξ against the fractional tune, with the observed beam–beam limits (gray band). Right, the dynamic beta β*′/β* for three values of ξ. Just above the half integer, the beam–beam force squeezes β* further, which e⁺e⁻ factories exploit.

let fig = Figure(size=(980, 380), fontsize=15)
    ν = range(0.001, 0.999, length=2000)
    ξmax = [v < 0.5 ? cot(π*v)/2π : -tan(π*v)/2π for v in ν]
    ax = Axis(fig[1, 1]; xlabel="fractional tune ν", ylabel="largest stable ξ", yscale=log10,
              title="linear stability limit", limits=(0, 1, 1e-3, 10))
    band!(ax, [0.0, 1.0], [0.004, 0.004], [0.15, 0.15]; color=(:gray, 0.2))
    text!(ax, 0.25, 0.006; text="observed beam–beam limits", fontsize=13, color=:gray35, align=(:center, :bottom))
    lines!(ax, ν, ξmax; color=PALETTE[1], linewidth=3)

    ax = Axis(fig[1, 2]; xlabel="fractional tune ν", ylabel="β*′ / β*", title="dynamic beta",
              limits=(0, 1, 0, 3))
    for (j, ξ) in enumerate([0.02, 0.05, 0.1])
        b = map(ν) do v
            c = cos(2π*v) - 2π*ξ*sin(2π*v)
            abs(c) < 1 || return NaN
            μ′ = acos(c); sin(2π*v) < 0 && (μ′ = 2π - μ′)
            sin(2π*v)/sin(μ′)
        end
        lines!(ax, ν, b; color=PALETTE[j + 1], linewidth=2.5, label="ξ = $ξ")
    end
    hlines!(ax, [1.0]; color=:gray55, linestyle=:dash)
    axislegend(ax; position=:lt, framevisible=false)
    fig
end
```

## Simulating the beam–beam kick

From here on, the simulations use the model collider of
[](interaction_region.md): a 519 m ring of FODO arcs with one low-β interaction
region, $\beta^*=1$ m, and 250 GeV protons with a normalized emittance of
2.5 μm. The beam–beam force is a thin element at the IP: the kick of a
Gaussian bunch of the other beam, computed with the field formulas above. This
is the **weak–strong** model. The tracked particles feel the opposing bunch,
which is held fixed. It is the standard tool for single-particle beam–beam
dynamics.

```{code-cell} julia
:tags: [hide-input]

# The model collider of the interaction-region page: lattice and matching.
using TrackPad, TrackPadWidgets, CairoMakie, StaticArrays, LinearAlgebra, Printf

beam = Beam(250.0e9; mass=TrackPad.M_PROTON, charge=1.0)    # 250 GeV protons
const EMIT = 2.5e-6/(beam.beta*beam.gamma)   # geometric rms emittance for εₙ = 2.5 μm

const LCELL = 30.0          # arc FODO cell [m]
const LQA   = 1.5           # arc quadrupole length
const NARC  = 12            # FODO cells in the arc
const LSTAR = 12.0          # IP to the first triplet quadrupole: room for the detector
const LQ1, LQ2 = 2.0, 3.5   # triplet: Q1 = Q3 short, Q2 long
const LQM   = 1.5           # matching quadrupoles Q4–Q7

arc_cell(kf, kd) = AbstractElement[
    Quadrupole(LQA/2, kf; name=:QF), Drift(LCELL/2 - LQA), Quadrupole(LQA, kd; name=:QD),
    Drift(LCELL/2 - LQA), Quadrupole(LQA/2, kf; name=:QF)]

"IP → start of the arc: the free drift L*, the triplet, and four matching quadrupoles."
ir_half(k) = AbstractElement[
    Drift(LSTAR),
    Quadrupole(LQ1, k[1]; name=:Q1), Drift(1.0), Quadrupole(LQ2, k[2]; name=:Q2), Drift(1.0),
    Quadrupole(LQ1, k[1]; name=:Q3), Drift(16.0),
    Quadrupole(LQM, k[3]; name=:Q4), Drift(9.0), Quadrupole(LQM, k[4]; name=:Q5), Drift(9.0),
    Quadrupole(LQM, k[5]; name=:Q6), Drift(9.0), Quadrupole(LQM, k[6]; name=:Q7), Drift(9.0)]

"""
The ring, starting at the IP: one half of the IR, the arc, and the mirror image
of the first half. `ip` holds whatever sits at the collision point.
"""
collider_lattice(c; ip=AbstractElement[]) = Lattice(AbstractElement[ip..., ir_half(c.k)...,
    (e for _ in 1:NARC for e in arc_cell(c.kf, c.kd))..., reverse(ir_half(c.k))...]; periodic=true)

# --- Courant–Snyder parameters from 2×2 blocks of a transfer matrix ------------
function cs_periodic(M)
    μ = acos(clamp((M[1,1] + M[2,2])/2, -1, 1)); M[1,2] < 0 && (μ = 2π - μ)
    (β = M[1,2]/sin(μ), α = (M[1,1] - M[2,2])/(2sin(μ)), μ = μ)
end
function cs_propagate(β0, α0, M)
    γ0 = (1 + α0^2)/β0
    β = M[1,1]^2*β0 - 2M[1,1]*M[1,2]*α0 + M[1,2]^2*γ0
    α = -M[1,1]*M[2,1]*β0 + (M[1,1]*M[2,2] + M[1,2]*M[2,1])*α0 - M[1,2]*M[2,2]*γ0
    (β = β, α = α, μ = mod(atan(M[1,2], M[1,1]*β0 - M[1,2]*α0), 2π))
end
blocks(M) = (M[1:2, 1:2], M[3:4, 3:4])

"Damped least squares (Levenberg–Marquardt) with a finite-difference Jacobian."
function lsq(f, p; w=1e-4, tol=1e-10)
    p0 = copy(p); g(q) = [f(q); w .* (q .- p0)]       # a weak pull keeps the solution continuous
    r = g(p); λ = 1e-3
    for _ in 1:300
        norm(f(p)) < tol && break
        J = reduce(hcat, [(g(p .+ 1e-8 .* (1:length(p) .== j)) .- r) ./ 1e-8 for j in 1:length(p)])
        while λ < 1e10
            q = p .- (J'J + λ*I) \ (J'r)
            rq = g(q)
            if norm(rq) < norm(r)
                p, r, λ = q, rq, max(λ/3, 1e-12); break
            end
            λ *= 4
        end
        λ ≥ 1e10 && break
    end
    norm(f(p)) < 1e-8 || error("matching failed")
    p
end

"Arc strengths (kf, kd) for phase advances (μx, μy) per cell."
arc_strengths(μx, μy, guess) = lsq(guess; w=0.0) do (kf, kd)
    Mx, My = blocks(transfer_map(Lattice(arc_cell(kf, kd)), beam))
    [(Mx[1,1] + Mx[2,2])/2 - cos(μx), (My[1,1] + My[2,2])/2 - cos(μy)]
end

"Match β* at the IP, a matched arc, and the tunes (Qx, Qy), starting from `seed`."
function match_collider(βstar, Qx, Qy; seed)
    k, (kf, kd), μx, μy = seed.k, (seed.kf, seed.kd), seed.μx, seed.μy
    for _ in 1:60
        kf, kd = arc_strengths(μx, μy, [kf, kd])
        cx, cy = cs_periodic.(blocks(transfer_map(Lattice(arc_cell(kf, kd)), beam)))
        ir_exit(k) = map(M -> cs_propagate(βstar, 0.0, M), blocks(transfer_map(Lattice(ir_half(k)), beam)))
        k = lsq(k) do k
            ex, ey = ir_exit(k)
            [log(ex.β/cx.β), ex.α, log(ey.β/cy.β), ey.α]
        end
        ex, ey = ir_exit(k)
        μxn, μyn = (2π*Qx - 2ex.μ)/NARC, (2π*Qy - 2ey.μ)/NARC    # trim the arcs for the tunes
        done = abs(μxn - μx) + abs(μyn - μy) < 1e-11
        μx, μy = μxn, μyn
        done && return (; βstar, k, kf, kd, μx, μy)
    end
    error("tune iteration failed")
end

"Walk β* from the matched state `c` to `target` in small steps, re-matching each time."
function squeeze(c, target, Qx, Qy)
    n = max(1, ceil(Int, abs(log(target/c.βstar))/log(1.12)))
    for b in exp.(range(log(c.βstar), log(target), length=n + 1)[2:end])
        c = match_collider(b, Qx, Qy; seed=c)
    end
    c
end

# A starting point at β* = 1 m, found once by a random search over the six strengths.
const SEED = (βstar = 1.0, k = [-0.0653, 0.0735, -0.0161, -0.0116, 0.0231, 0.0241],
              kf = 0.055, kd = -0.061, μx = deg2rad(69.6), μy = deg2rad(84.9))

const QX, QY = 3.685, 3.695           # RHIC-like working point
nominal = match_collider(1.0, QX, QY; seed=SEED)
nothing
```

```{code-cell} julia
const BSTAR = 1.0
const SIGMA = sqrt(EMIT*BSTAR)       # rms beam size at the IP (round beams)

"""
Thin weak–strong beam–beam kick at the IP, scaled to a linear tune shift ξ.
`StrongThinGaussianBeam` kicks by amplitude × field; the amplitude N r₀/γ
follows from ξ = β* amplitude/(4πσ²). Positive: like charges, defocusing.
"""
bb_kick(ξ) = StrongThinGaussianBeam(4π*SIGMA^2*ξ/BSTAR, SIGMA, SIGMA; name=:IP)

"Track `coords` turn by turn; returns x, pₓ, y, p_y histories (particles × turns)."
function track_turns(ring, coords, nturn)
    packed = GPULattice(ring, beam; dtype=Float64)
    c = copy(coords); N = size(c, 1)
    H = [zeros(N, nturn + 1) for _ in 1:4]
    for t in 0:nturn
        t > 0 && batch_ringpass!(c, packed, 1)
        for (j, col) in enumerate((1, 2, 3, 4)); H[j][:, t + 1] .= c[:, col]; end
    end
    H
end

"""
Tune from turn-by-turn data at the IP (α* = 0): the phase advance per turn in
normalized coordinates, averaged with smooth weights (a weighted Birkhoff
average), which converges much faster than a plain mean.
"""
function phase_tune(x, px)
    θ = atan.(-px .* BSTAR, x)
    n = length(θ) - 1; t = ((1:n) .- 0.5) ./ n
    w = exp.(-1 ./ (t .* (1 .- t)))
    sum(w .* mod.(diff(θ), 2π))/sum(w)/2π
end

# Check: a particle at 0.05σ should have its tune lowered by ξ.
let ξ = 0.01
    ring = collider_lattice(nominal; ip=[bb_kick(ξ)])
    X, PX, Y, PY = track_turns(ring, [0.05SIGMA 0 0.05SIGMA 0 0 0], 1024)
    @printf("beam–beam ξ = %.3f, N r₀/γ = %.3e m\n", ξ, 4π*SIGMA^2*ξ/BSTAR)
    @printf("tracked tunes (%.5f, %.5f), shifted by (%.5f, %.5f)\n", phase_tune(X[1, :], PX[1, :]),
            phase_tune(Y[1, :], PY[1, :]), phase_tune(X[1, :], PX[1, :]) - QX % 1,
            phase_tune(Y[1, :], PY[1, :]) - QY % 1)
end
```

A particle near the core sees its tune lowered by $\xi$, as the linear analysis
predicts (lowered because two proton beams repel).

## Amplitude-dependent tune shift and resonances

Only small-amplitude particles see the linear kick. A particle that oscillates
with amplitude $a$ spends part of each crossing in the weaker field outside the
core, so its tune is shifted less. For a round beam and motion in one plane,
averaging the kick over the betatron oscillation gives

$$
\frac{\Delta\nu(a)}{\xi}=\frac{4}{\alpha^2}\left[1-e^{-\alpha^2/4}\,I_0\!\left(\frac{\alpha^2}{4}\right)\right],
\qquad\alpha=\frac{a}{\sigma},
$$

where $I_0$ is a modified Bessel function. The shift is $\xi$ at the center,
about $\xi/2$ at $2\sigma$, and falls as $4\xi/\alpha^2$ far out. The
beam–beam force therefore spreads the tunes of the particles in a bunch over a
band of width about $\xi$.

That spread is what makes the beam–beam force dangerous. A particle whose tune
lands on a resonance $m\nu=p$ is driven by the nonlinear kick in step with its
own oscillation. Because the force of a centered, symmetric opposing beam is an
odd function of $x$, it drives the even resonances
$2\nu,\ 4\nu,\ 6\nu,\dots$ directly.

The widget below sets the horizontal tune at $\nu_x=0.765$, just above the
fourth-order resonance $\nu=3/4$, and raises $\xi$. The left panel shows the
phase space at the IP in normalized coordinates. The right panel shows the
tracked tune of each amplitude against the formula above.

- At $\xi=0$ the orbits are circles.
- As soon as $\xi>0.015$, the tune of the core falls below $3/4$. Particles at
  the amplitude where $\nu(a)=3/4$ are captured into **four resonance
  islands**, and their tune locks to $3/4$.

The islands are where particles can be carried to larger amplitudes; combined
with noise and with the motion in the other plane, this leads to emittance
growth and losses.

```{code-cell} julia
:tags: [hide-input]

"e^{-u} I₀(u) by quadrature of its integral representation."
eI0(u) = sum(exp(-u*(1 - cos(θ))) for θ in range(0, π, length=401)[1:end-1] .+ π/800)/400
detune(α) = α < 1e-6 ? 1.0 : 4/α^2*(1 - eI0(α^2/4))        # Δν/ξ

# Re-match the ring to νₓ = 0.765, in small steps from the nominal working point.
near_quarter = foldl((c, Q) -> match_collider(1.0, Q, QY; seed=c), range(QX, 3.765, length=5)[2:end];
                     init=nominal)
const NUX = 0.765
const XIS = [0.0, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04]
const AMP_PS = collect(0.3:0.3:6.0)             # starting amplitudes [σ] for the portraits
const AMP_Q  = collect(0.1:0.15:6.0)            # amplitudes for the tune curve

PORTRAITS = Dict(map(XIS) do ξ
    ring = collider_lattice(near_quarter; ip=[bb_kick(ξ)])
    n = length(AMP_PS)
    c = zeros(2n, 6)
    c[1:n, 1] .= AMP_PS .* SIGMA                                 # on the x̄ axis
    c[n+1:end, 1] .= AMP_PS ./ sqrt(2) .* SIGMA                  # and on the diagonal
    c[n+1:end, 2] .= AMP_PS ./ sqrt(2) .* SIGMA ./ BSTAR
    X, PX, _, _ = track_turns(ring, c, 200)
    q = zeros(length(AMP_Q), 6); q[:, 1] .= AMP_Q .* SIGMA
    QX_, QPX, _, _ = track_turns(ring, q, 1024)
    ξ => (X = round.(X ./ SIGMA; digits=2), P = round.(PX .* BSTAR ./ SIGMA; digits=2),
          ν = [phase_tune(QX_[i, :], QPX[i, :]) for i in 1:length(AMP_Q)])
end)

αs = range(0, 6, length=121)
explorer(
    title   = "Beam–beam resonance islands near ν = 3/4",
    sliders = [Knob("beam–beam parameter ξ", XIS; init = 5)],
    panels  = [Panel(xlabel="x̄ [σ]", ylabel="x̄′ [σ]", title="phase space at the IP",
                     xlim=(-6.6, 6.6), ylim=(-6.6, 6.6), equal=true, height=380, minwidth=300),
               Panel(xlabel="amplitude [σ]", ylabel="horizontal tune", title="tune against amplitude",
                     xlim=(0.0, 6.0), ylim=(0.72, 0.768), height=380, minwidth=280, legend=:bottomright)],
    statics = [line([0, 6], [0.75, 0.75]; panel=2, color=PALETTE[4], dash=true, width=1.5,
                    label="resonance ν = 3/4")],
    note    = "The islands sit at the amplitude where the tune curve crosses 3/4. Raising ξ moves "*
              "that crossing outward and widens the islands.",
) do ξ
    d = PORTRAITS[ξ]
    orbits = [points(d.X[i, :], d.P[i, :]; panel=1, color=PALETTE[mod1(i, 6)], size=1.6)
              for i in 1:size(d.X, 1)]
    (series = [orbits...,
               line(αs, NUX .- ξ .* detune.(αs); panel=2, color="#9aa4b2", width=2, label="formula"),
               points(AMP_Q, d.ν; panel=2, color=PALETTE[1], size=5, label="tracked")],
     readouts = ["core tune" => @sprintf("%.3f", NUX - ξ),
                 "tune spread" => @sprintf("%.3f", ξ)])
end
```

## The tune footprint

In two planes, every particle has its own pair of tunes $(\nu_x,\nu_y)$, set by
its two amplitudes. Mapping a grid of amplitudes from 0 to 6σ in each plane onto
the tune diagram gives the **tune footprint** of the beam. The footprint has to
fit between the resonance lines $m\nu_x+n\nu_y=p$ that the beam cannot
tolerate.

RHIC's polarized protons ran at $(\nu_x,\nu_y)\approx(0.685,0.695)$. That
working point sits between the third-order resonance $\nu=2/3$, driven by the
ring's sextupoles, and the tenth-order resonance $\nu=7/10$, both shown in the
widget below with the other resonances up to tenth order. The space between them
is only $0.033$ wide.

```{code-cell} julia
:tags: [hide-input]

const XIF = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03]
const RS = [0.1; collect(0.5:0.5:6.0)]                  # amplitude radius [σ]
const AS = range(0.02, π/2 - 0.02, length=9)            # angle between the planes
const GRID = [(r*cos(a), r*sin(a)) for a in AS, r in RS]

FOOTPRINTS = Dict(map(XIF) do ξ
    ring = collider_lattice(nominal; ip=[bb_kick(ξ)])
    c = zeros(length(GRID), 6)
    for (i, (ax, ay)) in enumerate(GRID); c[i, 1] = ax*SIGMA; c[i, 3] = ay*SIGMA; end
    X, PX, Y, PY = track_turns(ring, c, 1024)
    νx = [phase_tune(X[i, :], PX[i, :]) for i in 1:length(GRID)]
    νy = [phase_tune(Y[i, :], PY[i, :]) for i in 1:length(GRID)]
    ξ => (reshape(νx, size(GRID)), reshape(νy, size(GRID)))
end)

"Resonance lines m νx + n νy = p up to `order`, clipped to the window [lo, hi]²."
function resonance_lines(order, lo, hi)
    out = Any[]
    for m in -order:order, n in -order:order
        o = abs(m) + abs(n)
        (1 <= o <= order) || continue
        (m > 0 || (m == 0 && n > 0)) || continue            # each line once
        for p in floor(Int, min(m*lo + n*lo, m*lo + n*hi, m*hi + n*lo, m*hi + n*hi)):ceil(Int, max(m*lo + n*lo, m*lo + n*hi, m*hi + n*lo, m*hi + n*hi))
            gcd(gcd(m, n), p) == 1 || continue                # 2νx = 4/3 is the line νx = 2/3
            pts = Tuple{Float64,Float64}[]
            n != 0 && for x in (lo, hi)
                y = (p - m*x)/n; lo - 1e-12 <= y <= hi + 1e-12 && push!(pts, (x, y))
            end
            m != 0 && for y in (lo, hi)
                x = (p - n*y)/m; lo - 1e-12 <= x <= hi + 1e-12 && push!(pts, (x, y))
            end
            pts = unique(q -> round.(q; digits=9), pts)           # a corner is found twice
            length(pts) >= 2 || continue
            a, b = pts[1], pts[2]
            push!(out, line([a[1], b[1]], [a[2], b[2]]; color = o == 3 ? PALETTE[4] : o == 10 ? PALETTE[2] : "#9aa4b2",
                            width = o == 3 ? 2.2 : o == 10 ? 1.3 : max(0.6, 2.0 - 0.15o),
                            alpha = o == 3 ? 0.9 : 0.65))
        end
    end
    out
end

# Lines of constant amplitude are coloured by band, in both panels, so the two
# meshes can be matched: core (≤ 2σ), middle (2–4σ), tails (> 4σ).
band_color(r) = r <= 2 ? PALETTE[1] : r <= 4 ? PALETTE[3] : PALETTE[2]
amplitude_grid = Any[]
for i in 1:size(GRID, 1)
    push!(amplitude_grid, line(first.(GRID[i, :]), last.(GRID[i, :]); panel=2, color="#9aa4b2", width=1))
end
for (j, r) in enumerate(RS)
    push!(amplitude_grid, line(first.(GRID[:, j]), last.(GRID[:, j]); panel=2, color=band_color(r), width=1.6))
end

explorer(
    title   = "Beam–beam tune footprint at the RHIC working point",
    sliders = [Knob("beam–beam parameter ξ (total)", XIF; init = 4)],
    panels  = [Panel(xlabel="νₓ", ylabel="νᵧ", title="tune diagram", xlim=(0.65, 0.71), ylim=(0.65, 0.71),
                     equal=true, height=400, legend=:bottomright),
               Panel(xlabel="horizontal amplitude [σ]", ylabel="vertical amplitude [σ]",
                     title="the amplitudes tracked", xlim=(0.0, 6.3), ylim=(0.0, 6.3), equal=true, height=400)],
    statics = [resonance_lines(10, 0.65, 0.71)...,
               points([QX % 1], [QY % 1]; color="#111111", size=5, label="working point"),
               amplitude_grid...],
    note    = "Red: third-order resonances (2/3). Orange: tenth-order (7/10). Gray: orders 1–9, "*
              "thinner for higher order. Each mesh line of constant amplitude has the same colour "*
              "in both panels: the core of the beam (blue) is shifted most, the tails (orange) least.",
) do ξ
    νx, νy = FOOTPRINTS[ξ]
    mesh = Any[]
    for i in 1:size(νx, 1)                       # lines of constant angle
        push!(mesh, line(νx[i, :], νy[i, :]; color="#9aa4b2", width=1))
    end
    for (j, r) in enumerate(RS)                  # lines of constant amplitude
        push!(mesh, line(νx[:, j], νy[:, j]; color=band_color(r), width=1.6))
    end
    (series = mesh,
     readouts = ["core tunes" => @sprintf("(%.4f, %.4f)", νx[1, 1], νy[end, 1]),
                 "spread in νₓ" => @sprintf("%.4f", maximum(νx) - minimum(νx)),
                 "core distance to 2/3" => @sprintf("%+.4f", min(νx[1, 1], νy[end, 1]) - 2/3)])
end
```

At $\xi=0.015$ the whole footprint still fits between $2/3$ and $7/10$. Beyond
about $0.018$ the core of the beam reaches the third-order line, and the beam
loses particles and polarization. RHIC accordingly ran close to the total
beam–beam tune shift of about 0.015–0.02 that this space allows, with two
experiments colliding.

## Beyond the weak–strong picture

Several further effects shape real beam–beam limits:

- **Strong–strong interaction.** Both bunches are disturbed, and each changes
  the field the other feels. The two beams can then oscillate coherently
  against each other ($\sigma$ and $\pi$ modes), which is studied with
  self-consistent simulations of both beams.
- **Long-range encounters.** With many bunches, each bunch passes the bunches
  of the other beam at a distance several times before and after the IP, in the
  $1/r$ tail of the kick. The crossing angle of [](interaction_region.md) keeps
  these encounters ten or more beam sizes apart.
- **Compensation.** An electron beam can be shaped to exert the opposite of the
  beam–beam force on the protons. RHIC operated such **electron lenses** to
  partially compensate the head-on beam–beam effect.
- **Disruption.** In electron–ion and linear colliders, the beam–beam force can
  be strong enough to change the shape of the lighter beam during a single
  collision.

## Further reading

- *Beam Dynamics and Technologies for Future Colliders*, CERN Accelerator
  School, Zurich, 2018.
- W. Herr and T. Pieloni, "Beam–Beam Effects", arXiv:1601.05235.
