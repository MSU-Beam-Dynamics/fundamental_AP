---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Injection

## The goal

To zeroth order, injection means putting the beam into the synchrotron at its
lowest energy, the injection energy. In practice there is more to it:

- a ring with harmonic number $h\gg 1$ has many buckets to fill, often one
  bunch at a time from a smaller machine ("box-car" stacking);
- injecting the $i$-th bunch must not disturb bunches $1\ldots i-1$ already
  circulating;
- so the elements that steer the new beam onto the closed orbit have to be
  pulsed — *time-dependent* kicks, switched on for the injected beam only.


```{code-cell} julia
:tags: [hide-input]
:label: fig-bunch-filling
:caption: Box-car filling: bunches from a smaller ring A are transferred into the buckets of ring B one batch at a time, until B is full.

using CairoMakie, TrackPadWidgets

"Bunch positions on a ring of radius r centred at c, starting at the top and going clockwise."
ring_slots(c, r, n) = [(c[1] + r*sin(2π*(k - 1)/n), c[2] + r*cos(2π*(k - 1)/n)) for k in 1:n]

let cB = (0.0, 0.0), rB = 1.0, cA = (-1.05, -1.05), rA = 0.42, nA = 8, nB = 16
    stages = [("ring A filled", true, 0), ("first transfer", false, nA),
              ("A refilled", true, nA), ("B full", false, nB)]
    fig = Figure(size=(960, 280), fontsize=15)
    for (j, (title, a_full, nb)) in enumerate(stages)
        ax = Axis(fig[1, j]; aspect=DataAspect(), title=title, limits=(-1.6, 1.15, -1.6, 1.15))
        hidedecorations!(ax); hidespines!(ax)
        θ = range(0, 2π, length=200)
        for (c, r, name) in ((cB, rB, "B"), (cA, rA, "A"))
            lines!(ax, c[1] .+ r .* cos.(θ), c[2] .+ r .* sin.(θ); color=:gray60)
            text!(ax, c...; text=name, align=(:center, :center), fontsize=17)
        end
        a_full && scatter!(ax, ring_slots(cA, rA, nA); color=PALETTE[2], markersize=11)
        nb > 0 && scatter!(ax, ring_slots(cB, rB, nB)[1:nb]; color=PALETTE[1], markersize=11)
    end
    fig
end
```

## Normalized phase space

Injection and extraction are easiest to follow in **normalized phase space**.
At a point where the lattice has Courant–Snyder parameters $\beta$ and $\alpha$
(see [](../transverse_dynamics/twiss_fodo.md)), define

$$
\bar x=\frac{x}{\sqrt\beta},\qquad
\bar x'=\frac{\alpha x+\beta x'}{\sqrt\beta}.
$$

In these coordinates a matched beam is round, with rms radius $\sqrt\epsilon$, and
the linear motion from one point to another is a pure rotation by the betatron
phase advance between them. One full turn is a rotation by $2\pi Q$. All
phase-space plots on this page and the next use $\bar x$ and $\bar x'$ in units
of the rms size $\sigma=\sqrt\epsilon$ of the circulating beam.

The simulations use one small model ring throughout: four thin-lens FODO
cells, 80 m around (about the size of a medical synchrotron). For a thin-lens
FODO cell $\sin(\mu/2)=L/(2f)$, so the quadrupole strength that gives any
chosen tune follows directly. After the next figure, every plot is taken at the
point A, the centre of a focusing quadrupole, where $\alpha=0$, so
$\bar x'=\sqrt\beta\,p_x$.

The figure below shows the ring's optics, and then compares the two
descriptions of phase space at four points A–D in the first cell, one per
tab. Each tab shows a matched beam, and one particle observed at that point on
eight successive turns (numbered 0–7), in both coordinate systems.
The normalization is a linear map,

$$
\begin{pmatrix}\bar x\\ \bar x'\end{pmatrix}=
\begin{pmatrix}1/\sqrt\beta & 0\\ \alpha/\sqrt\beta & \sqrt\beta\end{pmatrix}
\begin{pmatrix}x\\ x'\end{pmatrix},
$$

with determinant 1. It therefore keeps areas, and with them the emittance. It
sends straight lines through the origin to straight lines through the origin,
so the gray rays connect the same points in both lower panels.

- **Top, the ring.** $\beta_x(s)$ (left axis) and the phase advance
  $\mu_x(s)/2\pi$ (right axis). The phase grows fastest where $\beta$ is
  small, and after one turn it reaches the tune $Q=1.27$. The point of the
  current tab is in black.
- **Bottom left, physical coordinates.** The beam fills an ellipse whose shape
  and tilt depend on the point, and the particle jumps around it in steps of
  unequal length. The axis ranges are the same in every tab.
- **Bottom right, normalized coordinates.** At every point the ellipse becomes
  the same circle, and every turn is the same clockwise rotation by $2\pi Q$.

```{code-cell} julia
:tags: [hide-input]

using StaticArrays, TrackPad, TrackPadWidgets, Random, Statistics, CairoMakie

beam = Beam(2.0e9; mass=TrackPad.M_PROTON, charge=1.0)    # 2 GeV kinetic protons

const NCELL = 4          # FODO cells around the ring
const LHALF = 10.0       # quadrupole centre to quadrupole centre [m]: an 80 m ring
const EMIT  = 1.0e-6     # rms emittance of the circulating beam [m·rad] (1 mm·mrad)

# Thin-lens FODO: sin(μ/2) = L/(2f), and each cell advances the phase by μ = 2πQ/NCELL.
inv_focal(Q) = 2sin(π*Q/NCELL)/LHALF
thin_quad(k1L, name) = ThinMultipole(0.0, nothing, SVector(0.0, k1L, 0.0, 0.0); max_order=1, name)

"""
    synchrotron(Q; b2=0, b3=0)

A ring of `NCELL` thin-lens FODO cells with tune `Q` in both planes. It starts,
and is observed, at the centre of a focusing quadrupole (α = 0), where a thin
multipole adds the sextupole and octupole kicks Δpₓ = −b2·x² − b3·x³.
"""
function synchrotron(Q; b2=0.0, b3=0.0)
    k = inv_focal(Q)
    elements = AbstractElement[ThinMultipole(0.0, nothing, SVector(0.0, 0.0, b2, b3);
                                             max_order=3, name=:NL)]
    for _ in 1:NCELL
        append!(elements, [thin_quad(k/2, :QF), Drift(LHALF), thin_quad(-k, :QD),
                           Drift(LHALF), thin_quad(k/2, :QF)])
    end
    Lattice(elements; periodic=true)
end

"β at the observation point of a periodic lattice (α = 0 there)."
beta_start(ring) = periodic_twiss(ring, beam).betax[1]

"Normalized coordinates (x̄, x̄′) of every particle, in units of σ = √EMIT."
normalized(c, β) = (c[:, 1] ./ sqrt(β*EMIT), c[:, 2] .* sqrt(β/EMIT))

"A Gaussian bunch matched to β (α = 0) with rms momentum spread σδ."
matched_beam(rng, N, β; σδ=0.0) =
    matched_gaussian(rng, N; betax=β, betay=β, emitx=EMIT, emity=EMIT,
                     longitudinal=[1.0e-6 0.0; 0.0 max(σδ, 1e-12)^2])

ring = synchrotron(1.27)

# The figure in the tabs below: the ring's optics, and phase space at A, B, C or D.

# Four observation points in the first cell, each given by the elements from the start
# of the ring (the centre of QF) up to it. A and C lie between the two halves of a thin
# quadrupole, where α = 0.
const OBS_POINTS = let k = inv_focal(1.27)
    [("A", AbstractElement[]),
     ("B", AbstractElement[thin_quad(k/2, :QF), Drift(LHALF/2)]),
     ("C", AbstractElement[thin_quad(k/2, :QF), Drift(LHALF), thin_quad(-k/2, :QD)]),
     ("D", AbstractElement[thin_quad(k/2, :QF), Drift(LHALF), thin_quad(-k, :QD), Drift(LHALF/2)])]
end

"Carry the phase-space vector `r` from the start of the ring through the elements `line`."
to_point(line, r) = isnothing(line) ? r : linepass(line, r, beam)

"""
    phase_space_figure(ip)

The ring's βₓ and phase advance (top), with observation point `ip` of `OBS_POINTS`
highlighted, and the same beam and the same particle at that point in physical
(bottom left) and normalized (bottom right) phase space. The axis ranges are the same
for every point, so the figures can be compared directly.
"""
function phase_space_figure(ip)
    β0, C, σ = beta_start(ring), NCELL*2LHALF, sqrt(EMIT)
    tw = periodic_twiss(ring, beam; sample_integrator_steps=true, max_step=0.1)
    Q  = tw.mux[end]/2π                                   # the full tune, 1.27
    ξ  = getchrom(ring, beam; centered=true)[1]           # natural chromaticity

    # β, α and s at every observation point.
    optics = map(OBS_POINTS) do (_, elems)
        isempty(elems) && return (β=β0, α=0.0, s=0.0, μ=0.0)
        t = transport_twiss(Lattice(elems), beam, optics4DUC(β0, 0.0, β0, 0.0))
        (β=t.betax[end], α=t.alphax[end] + 0.0, s=t.s[end], μ=t.mux[end])
    end
    name, elems = OBS_POINTS[ip]
    line = isempty(elems) ? nothing : Lattice(elems)
    β, α = optics[ip].β, optics[ip].α
    γ = (1 + α^2)/β

    # The same matched beam and the same particle (turns 0–7), carried to this point.
    cloud = matched_beam(MersenneTwister(3), 800, β0)
    for i in axes(cloud, 1)
        cloud[i, :] .= to_point(line, SVector{6}(cloud[i, :]))
    end
    r = SVector(2.6sqrt(β0*EMIT), 0.6sqrt(EMIT/β0), 0.0, 0.0, 0.0, 0.0)
    pts = SVector{6,Float64}[]
    for _ in 0:7
        push!(pts, to_point(line, r))
        r = linepass(ring, r, beam)                       # one more turn, back at the start
    end
    x  = [p[1] for p in pts];  xp = [p[2] for p in pts]

    norm_x(x, xp)  = x ./ sqrt(β) ./ σ                       # x̄  in units of σ
    norm_xp(x, xp) = (α .* x .+ β .* xp) ./ sqrt(β) ./ σ     # x̄′ in units of σ
    W = γ*x[1]^2 + 2α*x[1]*xp[1] + β*xp[1]^2                # the particle's invariant
    φ = range(0, 2π, length=200)
    ellipse(w) = (sqrt(w*β) .* cos.(φ), -sqrt(w/β) .* (α .* cos.(φ) .+ sin.(φ)))

    fig = Figure(size=(980, 800), fontsize=15)

    # Top: β and the phase advance around the ring, with the quadrupoles as a strip.
    top = β0/(1 - 0.17)                                   # room for the strip above the data
    axβ = Axis(fig[1, 1:2]; xlabel="s [m]", ylabel="βₓ [m]", ylabelcolor=PALETTE[1],
               title="model ring: Qₓ = Qᵧ = $(round(Q; digits=2)), natural chromaticity ξ = $(replace(string(round(ξ; digits=2)), "-" => "−"))")
    axμ = Axis(fig[1, 1:2]; ylabel="μₓ/2π", yaxisposition=:right, ylabelcolor=PALETTE[4],
               yticklabelcolor=PALETTE[4], yticks=0:0.25:1.25, ygridvisible=false)
    hidespines!(axμ); hidexdecorations!(axμ); linkxaxes!(axβ, axμ)
    base, h = top*(1 - 0.045), top*0.035
    hlines!(axβ, [base]; color=:gray50, linewidth=1)
    for n in 0:2NCELL                                    # QF above the line, QD below
        poly!(axβ, Rect(n*LHALF - 0.6, iseven(n) ? base : base - h, 1.2, h); color=(:seagreen, 0.7))
    end
    for shift in (-C, 0.0)                               # periodic: also show the end of the previous turn
        lines!(axβ, tw.s .+ shift, tw.betax; color=PALETTE[1], linewidth=2.5)
    end
    lines!(axμ, tw.s, tw.mux ./ 2π; color=PALETTE[4], linewidth=2)
    text!(axμ, C, Q; text="Q = $(round(Q; digits=2)) ", align=(:right, :bottom), color=PALETTE[4])
    for (j, (lab, _)) in enumerate(OBS_POINTS)           # this point in black, the others in gray
        c, o = j == ip ? (:black, optics[j]) : (:gray65, optics[j])
        lines!(axβ, [o.s, o.s], [0, β0]; color=c, linestyle=:dash, linewidth=1)
        scatter!(axβ, [o.s], [o.β]; color=c, markersize=j == ip ? 11 : 8)
        text!(axβ, o.s, β0; text=" " * lab, align=(:left, :bottom), color=c,
              fontsize=j == ip ? 17 : 14, font=:bold)
    end
    xlims!(axβ, -2, C + 2); ylims!(axβ, 0, top); ylims!(axμ, 0, Q/(1 - 0.17))

    # Bottom: phase space at this point, physical (left) and normalized (right).
    colors = cgrad(:viridis, 8; categorical=true)
    axL = Axis(fig[2, 1]; xlabel="x [mm]", ylabel="x′ [mrad]",
               title="physical phase space at $name (β = $(round(β; digits=1)) m, α = $(replace(string(round(α; digits=2) + 0.0), "-" => "−")))")
    axR = Axis(fig[2, 2]; xlabel="x̄ [σ]", ylabel="x̄′ [σ]", aspect=DataAspect(),
               title=ip == 1 ? "normalized phase space at A" :
                     "normalized phase space at $name ($(round(Int, rad2deg(optics[ip].μ)))° after A)")

    # The same ranges at every point: wide enough for the largest β and γ.
    lim = 1.25*sqrt(W)/σ
    xl = 1.2*sqrt(W*maximum(o.β for o in optics))*1e3
    yl = 1.2*sqrt(W*maximum((1 + o.α^2)/o.β for o in optics))*1e3
    for (ax, fx, fy, sx, sy, L) in ((axL, (a, b) -> a, (a, b) -> b, 1e3, 1e3, (xl, yl)),
                                    (axR, norm_x, norm_xp, 1.0, 1.0, (lim, lim)))
        scatter!(ax, sx .* fx(cloud[:, 1], cloud[:, 2]), sy .* fy(cloud[:, 1], cloud[:, 2]);
                 color=(:gray60, 0.35), markersize=3, strokewidth=0)
        for w in (EMIT, 4EMIT)                                   # 1σ and 2σ contours
            ex, ey = ellipse(w)
            lines!(ax, sx .* fx(ex, ey), sy .* fy(ex, ey); color=:gray50, linestyle=:dash, linewidth=1)
        end
        ex, ey = ellipse(W)                                      # the particle's own ellipse
        lines!(ax, sx .* fx(ex, ey), sy .* fy(ex, ey); color=PALETTE[1], linewidth=2)
        X, Y = sx .* fx(x, xp), sy .* fy(x, xp)
        for i in 1:8                                             # rays from the origin
            lines!(ax, [0, X[i]], [0, Y[i]]; color=(:gray40, 0.5), linewidth=1)
        end
        scatter!(ax, X, Y; color=[colors[i] for i in 1:8], markersize=14, strokecolor=:black, strokewidth=1)
        for i in 1:8                                             # labels pushed outward on screen
            d = (X[i]/L[1], Y[i]/L[2]); d = d ./ hypot(d...)
            text!(ax, X[i], Y[i]; text=string(i - 1), align=(:center, :center), fontsize=15,
                  offset=(17d[1], 17d[2]))
        end
        scatter!(ax, [0.0], [0.0]; marker=:cross, color=:black, markersize=10)
    end
    # The angle of one turn, between points 0 and 1 in normalized coordinates.
    θ0 = atan(norm_xp(x[1], xp[1]), norm_x(x[1], xp[1]))
    θ1 = θ0 - 2π*Q % (2π)
    arc = range(θ0, θ1, length=40)
    lines!(axR, 0.7 .* cos.(arc), 0.7 .* sin.(arc); color=PALETTE[4], linewidth=2.5)
    text!(axR, 1.45cos((θ0 + θ1)/2), 1.45sin((θ0 + θ1)/2); text="2πQ", color=PALETTE[4],
          align=(:center, :center), fontsize=15)
    limits!(axR, -lim, lim, -lim, lim)
    limits!(axL, -xl, xl, -yl, yl)
    rowsize!(fig.layout, 1, Relative(0.34))
    fig
end
nothing
```

```{code-cell} julia
:tags: [remove-cell]
:label: fig-phase-space-A

phase_space_figure(1)
```

```{code-cell} julia
:tags: [remove-cell]
:label: fig-phase-space-B

phase_space_figure(2)
```

```{code-cell} julia
:tags: [remove-cell]
:label: fig-phase-space-C

phase_space_figure(3)
```

```{code-cell} julia
:tags: [remove-cell]
:label: fig-phase-space-D

phase_space_figure(4)
```

::::{tab-set}
:::{tab-item} A · centre of QF
**A, the centre of a focusing quadrupole.** $\beta$ is at its largest and
$\alpha=0$. The ellipse is upright, wide in $x$ and narrow in $x'$, and
normalizing only rescales the two axes.

```{embed} #fig-phase-space-A
:remove-input: true
```
:::
:::{tab-item} B · drift after QF
**B, halfway along the drift after QF.** $\beta$ is falling steeply, so
$\alpha=-\beta'/2$ is large and positive. The beam is converging: particles
with $x>0$ have $x'<0$, which tilts the ellipse.

```{embed} #fig-phase-space-B
:remove-input: true
```
:::
:::{tab-item} C · centre of QD
**C, the centre of a defocusing quadrupole.** $\beta$ is at its smallest and
$\alpha=0$ again. The ellipse is upright once more, now narrow in $x$ and tall
in $x'$.

```{embed} #fig-phase-space-C
:remove-input: true
```
:::
:::{tab-item} D · drift before QF
**D, halfway along the drift before the next QF.** $\beta$ is growing and
$\alpha<0$. The beam is diverging, and the ellipse tilts the other way.

```{embed} #fig-phase-space-D
:remove-input: true
```
:::
::::

In physical coordinates the beam changes shape from point to point, while its
area, the emittance, stays the same. In normalized coordinates every tab shows
the same circle. Only the angle of the particle's points changes. They turn
clockwise by the phase advance from A, which is given in the panel title.

The computed optics confirm the tune of 1.27 that the quadrupoles were set
for. Without any sextupoles, the ring's natural chromaticity is about −2; this
matters below, where it is what makes an injection error grow into a larger
emittance.

## Single-turn injection

The simplest scheme uses two elements:

- a **septum magnet**: a thin wall (2–20 mm) separates the field-free
  region of the circulating beam from a high-field region that bends the
  incoming beam to run alongside it;
- a **kicker magnet**, downstream, which gives the injected beam the last
  angle it needs to land on the closed orbit.

Slow closed-orbit bumpers move the circulating orbit close to the septum, so
the angle the kicker has to supply stays small. The kicker field rises after
the last circulating bunch has passed and falls before it comes back, so only
the injected bunch is deflected. Bunches already in the ring miss the septum
because they are separated from the injected beam in *space*, and miss the kick
because they are separated from it in *time*.


```{figure} ../images/synchrotron/one_turn_injection.png
:width: 85%
:name: fig-one-turn-injection
Single-turn injection. The transfer line brings the beam through the septum magnet, closed-orbit bumpers move the circulating orbit close to the septum, and the kicker puts the injected beam onto the closed orbit. Inset: 'box-car' stacking, with the kicker field switched on for the injected bunch only.
```

```{figure} ../images/synchrotron/injection_layouts.png
:width: 90%
:name: fig-injection-layouts
Real injection layouts. Top: RHIC's six o'clock insertion, with the injection equipment of the yellow ring. Bottom: LHC injection, with the septum deflecting in the horizontal plane and the kicker in the vertical.
```

**How strong must the kick be?** Let the septum (point 0) and the kicker
(point 1) be separated by a betatron phase advance $\mu$. The transfer
matrix between two points with known Courant–Snyder parameters is

$$
M_{0\to1}=
\begin{pmatrix}
\sqrt{\dfrac{\beta_1}{\beta_0}}\,(\cos\mu+\alpha_0\sin\mu) & \sqrt{\beta_0\beta_1}\,\sin\mu\\[2ex]
-\dfrac{(1+\alpha_0\alpha_1)\sin\mu+(\alpha_1-\alpha_0)\cos\mu}{\sqrt{\beta_0\beta_1}} &
\sqrt{\dfrac{\beta_0}{\beta_1}}\,(\cos\mu-\alpha_1\sin\mu)
\end{pmatrix}.
$$

The injected beam leaves the septum at $(x_0,x_0')$ and should reach the
kicker on axis, $x_1=0$, which fixes the launch angle
$x_0'=-M_{11}x_0/M_{12}$. Its angle there is then
$x_1'=(M_{21}M_{12}-M_{22}M_{11})\,x_0/M_{12}=-x_0/M_{12}$, because
$\det M=1$. The kicker has to cancel that angle:

$$
\theta_{\mathrm{kicker}}=\frac{x_0}{\sqrt{\beta_0\beta_1}\,\sin\mu}.
$$

The kick is smallest for $\mu=\pi/2$. A quarter of a betatron oscillation turns
the septum's displacement into a pure angle at the kicker, which the kicker then
removes. In normalized phase space that is three steps: a displacement, a
rotation by 90°, and a kick back to the origin. The widget below shows those
three steps, and lets you change the phase advance of the cell between the two
magnets:

```{code-cell} julia
:tags: [hide-input]

# The septum-to-kicker section: one thin-lens FODO cell from QF centre to QF
# centre, so α = 0 and β is the same at both ends.
injection_cell(μ) = (k = 2sin(μ/2)/LHALF;
    Lattice(AbstractElement[thin_quad(k/2, :QF), Drift(LHALF), thin_quad(-k, :QD),
                            Drift(LHALF), thin_quad(k/2, :QF)]; periodic=true))

const X_SEPTUM = 5.0      # septum blade, in σ of the circulating beam
const X_INJECT = 7.0      # centre of the injected beam at the septum
const NINJ = 400

explorer(
    title   = "Single-turn injection in normalized phase space",
    sliders = [Knob("phase advance septum → kicker μ [deg]", 45.0:15.0:135.0;
                    fmt = m -> string(Int(m)), init = 4)],
    panels  = [Panel(xlabel="x̄ [σ]", ylabel="x̄′ [σ]", title="at the septum", xlim=(-12.0, 12.0),
                     ylim=(-12.0, 12.0), equal=true, height=260, legend=:bottomleft,
                     basis="30%", minwidth=200),
               Panel(xlabel="x̄ [σ]", ylabel="x̄′ [σ]", title="at the kicker, before the kick",
                     xlim=(-12.0, 12.0), ylim=(-12.0, 12.0), equal=true, height=260,
                     legend=:bottomleft, basis="30%", minwidth=200),
               Panel(xlabel="x̄ [σ]", ylabel="x̄′ [σ]", title="after the kick", xlim=(-12.0, 12.0),
                     ylim=(-12.0, 12.0), equal=true, height=260, legend=:bottomleft,
                     basis="30%", minwidth=200)],
    statics = [line([X_SEPTUM, X_SEPTUM], [-12.0, 12.0]; panel=1, color="#5b6472", width=3.0,
                    label="septum blade")],
    note = "Blue: the circulating beam. Orange: the injected beam, launched by the septum "*
           "so that it crosses the axis at the kicker. The kick it then needs, "*
           "x₀/(β sin μ), is smallest for μ = 90°, where the septum displacement has turned "*
           "into a pure angle.",
) do μdeg
    μ = deg2rad(μdeg)
    cell = injection_cell(μ)
    β = beta_start(cell)
    rng = MersenneTwister(1)
    circ = matched_beam(rng, NINJ, β)
    inj = matched_beam(rng, NINJ, β)
    x0, x0p = X_INJECT*sqrt(β*EMIT), -X_INJECT*cot(μ)*sqrt(EMIT/β)    # launch with x₁ = 0
    inj[:, 1] .+= x0; inj[:, 2] .+= x0p
    stage = [(normalized(circ, β), normalized(inj, β))]
    linepass!(circ, cell, beam, zeros(Int, NINJ)); linepass!(inj, cell, beam, zeros(Int, NINJ))
    push!(stage, (normalized(circ, β), normalized(inj, β)))
    θ = x0/(β*sin(μ))                                   # the kicker angle derived above
    linepass!(inj, Lattice(AbstractElement[HKicker(; xkick=θ)]), beam, zeros(Int, NINJ))
    push!(stage, (normalized(circ, β), normalized(inj, β)))
    series = Any[]
    for (p, ((cx, cp), (ix, ip))) in enumerate(stage)
        push!(series, points(cx, cp; panel=p, color=PALETTE[1], size=1.6, alpha=0.6,
                             label = p == 1 ? "circulating" : nothing))
        push!(series, points(ix, ip; panel=p, color=PALETTE[2], size=1.6, alpha=0.6,
                             label = p == 1 ? "injected" : nothing))
    end
    xf, pf = stage[3][2]
    (series = series,
     readouts = ["β at septum and kicker" => string(round(β; digits=1), " m"),
                 "septum launch angle x₀′" => string(round(1e3*x0p; digits=3) + 0.0, " mrad"),
                 "kicker angle θ" => string(round(1e3*θ; digits=3), " mrad"),
                 "injected centroid after the kick" => string("(", round(mean(xf); digits=2) + 0.0, ", ",
                                                              round(mean(pf); digits=2) + 0.0, ") σ")])
end
```

### Septum and kicker magnets

A **septum** is a pulsed or DC magnet whose thin septum conductor (2–20 mm)
separates a field-free region, where the circulating beam passes, from the
high-field region that deflects the injected beam. Currents of 5–25 kA flow
in the septum coil. In a **Lambertson septum** the septum is a thin steel wall
instead: the circulating beam passes through a field-free hole in the yoke, and
the deflection is perpendicular to the plane in which the two beams were
separated.


```{figure} ../images/synchrotron/septum.jpg
:width: 95%
:name: fig-septum
Septum magnets: a septum with its coil (left), the cross-section of a current septum, whose thin coil separates the field-free region from the deflecting field $B_0$ (top right), and a Lambertson septum, where a thin steel wall does the same and the circulating beam passes through the field-free hole (bottom right).
```

A **kicker** is a fast pulsed magnet, usually with a ferrite yoke. Its rise and
fall times must fit into the gaps between bunches, and its flat top must be as
long as the injected bunch train. This is why the filling pattern of a ring
leaves gaps for the kicker.


```{figure} ../images/synchrotron/injection_kicker.png
:width: 90%
:name: fig-injection-kicker
The RHIC injection kicker (left) and a measured kicker current pulse (right): the bunch is injected while the pulse is at its top.
```

## Injection errors and filamentation

If the kicker strength or the septum position is slightly wrong, the injected
beam arrives off the closed orbit and performs a coherent betatron oscillation
around it: an **injection oscillation**. If every particle had exactly the same
tune, the bunch would simply rotate in normalized phase space forever. But a
real beam has a tune spread, from

- **chromaticity** combined with the momentum spread, $\Delta Q=\xi\,\delta$, and
- **nonlinear fields**, which make the tune depend on amplitude.

Particles with different tunes rotate at different rates, so the bunch smears
around the circle: it **filaments**. The centroid oscillation decoheres, and
eventually the beam fills a circle much larger than the one it was injected
with, so the emittance has grown.

```{code-cell} julia
:tags: [hide-input]
:label: fig-injection-error
:caption: An injection error in normalized phase space, tracked through the ring of this page with an octupole added so that the tune grows with amplitude. The bunch arrives with $\bar x'$ off by 3σ. For the first few turns it rotates about the closed orbit (+) as a whole, by $2\pi Q\approx 97^\circ$ per turn. Then the particles at larger amplitude run ahead, the bunch winds into a spiral, and it ends as a centred beam filling the predicted 2σ circle (red). Grey: the 2σ circle of the matched beam.

using CairoMakie

βf = beta_start(ring)
oct_ring = synchrotron(1.27; b3 = 0.01/(βf^2*EMIT))    # Δx̄′ = −0.01 x̄³ (σ units): Q grows with amplitude
ring_circle(r) = (r .* cos.(range(0, 2π, length=121)), r .* sin.(range(0, 2π, length=121)))

"Phase-space snapshots, at the given turns, of a bunch spoiled by `spoil!` and tracked through `oct_ring`."
function filament(spoil!, turns; N=4000)
    c = matched_beam(MersenneTwister(7), N, βf)
    spoil!(c)
    packed = GPULattice(oct_ring, beam; dtype=Float64)
    snaps = Dict{Int,Any}(); t = 0
    for T in turns
        T > t && (batch_ringpass!(c, packed, T - t); t = T)
        snaps[T] = normalized(c, βf)
    end
    snaps
end

"A grid of snapshots; the last panel also shows the predicted filamented 2σ circle of radius r_final."
function snapshot_grid(snaps, turns, ncol, r_final)
    nrow = cld(length(turns), ncol)
    fig = Figure(size=(230ncol, 240nrow), fontsize=13ncol/3)    # same type size once scaled to the page
    for (k, T) in enumerate(turns)
        i, j = cld(k, ncol), mod1(k, ncol)
        ax = Axis(fig[i, j]; aspect=DataAspect(), title="turn $T", limits=(-8, 8, -8, 8),
                  xlabel="x̄ [σ]", ylabel="x̄′ [σ]")
        i < nrow && hidexdecorations!(ax; grid=false, ticks=false)
        j > 1 && hideydecorations!(ax; grid=false, ticks=false)
        X, P = snaps[T]
        scatter!(ax, X, P; color=(PALETTE[1], 0.35), markersize=2.5, strokewidth=0)
        lines!(ax, ring_circle(2.0)...; color=:gray55, linestyle=:dash)
        k == length(turns) && lines!(ax, ring_circle(r_final)...; color=PALETTE[4], linestyle=:dash, linewidth=2)
        scatter!(ax, [0.0], [0.0]; marker=:cross, color=:black, markersize=10)
    end
    fig
end

T_ERR = [0, 1, 2, 3, 4, 10, 20, 30, 45, 60, 100, 150, 200, 400, 1000, 3000]
snapshot_grid(filament(c -> (c[:, 2] .+= 3sqrt(EMIT/βf)), T_ERR), T_ERR, 4, 2sqrt(1 + 3^2/2))
```

An optical **mismatch** (a beam injected with
the wrong $\beta$ or $\alpha$; see [](../transverse_dynamics/emittance.md))
filaments the same way into a larger ellipse of the matched shape.

```{code-cell} julia
:tags: [hide-input]
:label: fig-mismatch-filamentation
:caption: An optical mismatch filaments the same way. The bunch arrives with $\beta$ three times too large, an ellipse stretched along $\bar x$ (grey: the matched 2σ circle). At first it simply rotates, by $2\pi Q$ per turn; then it shears into a spiral, and it ends as a round beam filling the predicted 2σ circle (red). Same ring and octupole as in the figure above.

T_MIS = [0, 1, 2, 3, 10, 20, 40, 70, 100, 150, 400, 3000]
snapshot_grid(filament(c -> (c[:, 1] .*= sqrt(3); c[:, 2] ./= sqrt(3)), T_MIS), T_MIS, 4,
              2sqrt((3 + 1/3)/2))
```

Both final emittances follow from averaging over the betatron phase:

$$
\epsilon_{\mathrm{dipole}}=\epsilon_0+\tfrac12\,\Delta\bar A^2,
\qquad
\epsilon_{\mathrm{mismatch}}=\tfrac12\left(\beta\gamma_0-2\alpha\alpha_0+\gamma\beta_0\right)\epsilon_0,
$$

where $\Delta\bar A$ is the normalized amplitude of the injection error, and
$(\beta_0,\alpha_0,\gamma_0)$ describe the injected beam while
$(\beta,\alpha,\gamma)$ describe the lattice. In the widget below the tune spread comes instead from
the ring's natural chromaticity of −2.1 and an rms momentum spread of
$10^{-3}$, so $\sigma_Q\approx 2\times10^{-3}$, and the beam is filamented
after a few hundred turns:

```{code-cell} julia
:tags: [hide-input]

const NFIL = 1000
const FIL_TURNS = [0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 60, 80, 100, 150, 200, 300, 400, 600]
βfil = beta_start(ring)

"rms emittance about the centroid, in units of EMIT, and the centroid x̄."
function emittance_and_centroid(c)
    X, P = normalized(c, βfil)
    X0, P0 = X .- mean(X), P .- mean(P)
    sqrt(mean(X0.^2)*mean(P0.^2) - mean(X0.*P0)^2), mean(X)
end

# (label, how the injected beam differs from the matched one, predicted ε∞/ε₀)
FIL_CASES = [("kicker error: x̄′ off by 3σ", c -> (c[:, 2] .+= 3sqrt(EMIT/βfil)), 1 + 3^2/2),
             ("β mismatch: β_beam = 3β",   c -> (c[:, 1] .*= sqrt(3); c[:, 2] ./= sqrt(3)), (3 + 1/3)/2)]

fil_runs = map(FIL_CASES) do (_, spoil!, _)
    c = matched_beam(MersenneTwister(2), NFIL, βfil; σδ=1.0e-3)
    spoil!(c)
    packed = GPULattice(ring, beam; dtype=Float64)     # packed copy: same physics, faster
    snaps = Dict(0 => normalized(c, βfil))
    hist = [emittance_and_centroid(c)]
    for t in 1:maximum(FIL_TURNS)
        batch_ringpass!(c, packed, 1)
        push!(hist, emittance_and_centroid(c))
        t in FIL_TURNS && (snaps[t] = normalized(c, βfil))
    end
    (snaps = snaps, emit = first.(hist), centroid = last.(hist))
end

circle(r) = (r .* cos.(range(0, 2π, length=121)), r .* sin.(range(0, 2π, length=121)))

explorer(
    title   = "An injection error filaments into a larger emittance",
    sliders = [Knob("turn", FIL_TURNS; fmt = string, init = 1),
               Knob("error", 1:length(FIL_CASES); fmt = i -> FIL_CASES[i][1], init = 1)],
    panels  = [Panel(xlabel="x̄ [σ]", ylabel="x̄′ [σ]", title="normalized phase space",
                     xlim=(-8.0, 8.0), ylim=(-8.0, 8.0), equal=true, height=300,
                     legend=:bottomleft, basis="40%", minwidth=240),
               Panel(xlabel="turn", ylabel="ε / ε₀", title="emittance and centroid",
                     xlim=(0.0, Float64(maximum(FIL_TURNS))), ylim=(0.0, 6.5),
                     y2label="centroid ⟨x̄⟩ [σ]", y2lim=(-4.0, 4.0), height=300,
                     legend=:right, basis="55%", minwidth=260)],
    statics = [line(circle(2.0)...; panel=1, color="#9aa4b2", dash=true, label="2σ, matched beam")],
    note = "The centroid (orange) rings at the betatron frequency and decays as the particles "*
           "dephase; the rms emittance (blue) climbs to the value predicted by averaging over "*
           "the betatron phase (dashed). Neither process can be undone once the beam has "*
           "filamented.",
) do turn, icase
    label, _, ratio = FIL_CASES[icase]
    run = fil_runs[icase]
    X, P = run.snaps[turn]
    t = 0:turn
    (series = [points(X, P; panel=1, color=PALETTE[1], size=1.6, alpha=0.6, label="particles"),
               line(circle(2sqrt(ratio))...; panel=1, color=PALETTE[4], dash=true,
                    label="2σ, filamented (predicted)"),
               line([0, maximum(FIL_TURNS)], [ratio, ratio]; panel=2, color=PALETTE[4], dash=true,
                    label="predicted ε∞/ε₀"),
               line(collect(t), run.emit[t .+ 1]; panel=2, color=PALETTE[1], width=1.6,
                    label="ε / ε₀"),
               line(collect(t), run.centroid[t .+ 1]; panel=2, color=PALETTE[2], width=1.0,
                    axis=:right, label="centroid")],
     readouts = ["turn" => turn,
                 "ε / ε₀" => round(run.emit[turn+1]; digits=3),
                 "predicted ε∞ / ε₀" => round(ratio; digits=3),
                 "centroid ⟨x̄⟩" => string(round(run.centroid[turn+1]; digits=2) + 0.0, " σ")])
end
```

To sum up: both an injection angle or position error and an optical mismatch
blow up the emittance, and the blow-up happens through filamentation, which in
turn comes from the tune spread across the beam (nonlinearity and chromatic
effects). The only cure is to correct the injection trajectory and match the
optics before the beam filaments.

## Multi-turn injection

Sometimes a single turn cannot deliver enough intensity: the source or the
previous accelerator may be too weak, or space charge may limit the density.
The beam is then injected over many turns:

- instead of one kicker, several fast dipoles make a closed-orbit **bump** at the
  septum, and its amplitude shrinks turn by turn during injection;
- each new beamlet lands just outside the septum blade, and so at a larger
  amplitude about the moving closed orbit than the one before it;
- the tune is chosen so that the earlier beamlets, as they rotate, come back to
  the septum's side of phase space only after the bump has pulled the orbit
  away.

```{figure} ../images/synchrotron/multiturn_idea.png
:width: 75%
:name: fig-multiturn-idea
Multi-turn injection: fast bumpers make a closed-orbit bump at the septum whose amplitude varies from turn to turn.
```

The beam is "painted" into transverse phase space. Liouville's theorem has a
price, though: beamlets cannot be laid on top of each other, so the stored
emittance is much larger than that of a single beamlet, and particles that rotate
back into the blade are lost. Step through the turns, and compare tunes:

```{code-cell} julia
:tags: [hide-input]

const XS = 4.0                    # septum blade [σ of the ring's reference beam]
const σB = 0.4                    # rms size of one injected beamlet [σ]
const XINJ = XS + 2σB             # beamlets arrive just outside the blade
const NB = 100                    # particles per beamlet
const NPAINT = 12                 # turns shown
const BUMP_STEP = 1.1σB           # the bump falls this much per turn

# closed orbit at the septum on turn n: it starts just inside the blade and falls to zero
bump(n) = max(XS - 2σB - BUMP_STEP*(n - 1), 0.0)
PAINT_TUNES = [0.22, 0.25, 0.28, 0.31]

"Track the painting at fractional tune q; return a snapshot after every turn."
function paint(q)
    r = synchrotron(1 + q); β = beta_start(r)
    rng = MersenneTwister(3)
    U = zeros(0, 6); tag = Int[]; lost = 0
    snaps = []
    for n in 1:NPAINT
        if bump(n) > 0                               # inject while the bump is on
            new = zeros(NB, 6)
            new[:, 1] .= ((XINJ - bump(n)) .+ σB .* randn(rng, NB)) .* sqrt(β*EMIT)
            new[:, 2] .= σB .* randn(rng, NB) .* sqrt(EMIT/β)
            U = vcat(U, new); append!(tag, fill(n, NB))
        end
        ringpass!(U, r, beam, zeros(Int, size(U, 1)), 1)
        X, P = normalized(U, β)
        hit = X .+ bump(n + 1) .> XS                 # lab position beyond the blade
        lost += count(hit)
        U, tag = U[.!hit, :], tag[.!hit]
        X, P = X[.!hit], P[.!hit]
        push!(snaps, (X = X .+ bump(n + 1), P = P, tag = copy(tag), lost = lost,
                      injected = NB*count(m -> bump(m) > 0, 1:n), orbit = bump(n + 1),
                      emit = isempty(X) ? NaN :
                             sqrt(mean((X .- mean(X)).^2)*mean((P .- mean(P)).^2) -
                                  mean((X .- mean(X)).*(P .- mean(P)))^2)))
    end
    snaps
end
paint_runs = Dict(q => paint(q) for q in PAINT_TUNES)
TURN_COLORS = ["#2563eb", "#e2711d", "#0f9d58", "#c2255c", "#7c3aed", "#0d9488",
               "#a16207", "#475569", "#db2777", "#65a30d", "#0891b2", "#9333ea"]

explorer(
    title   = "Multi-turn injection: painting phase space",
    sliders = [Knob("turn", 1:NPAINT; fmt = string, init = 1),
               Knob("fractional tune", PAINT_TUNES; fmt = q -> string(q), init = 3)],
    panels  = [Panel(xlabel="x̄ at the septum [σ]", ylabel="x̄′ [σ]", title="phase space at the septum",
                     xlim=(-7.0, 7.0), ylim=(-7.0, 7.0), equal=true, height=340,
                     legend=:bottomleft, basis="60%", minwidth=300)],
    statics = [line([XS, XS], [-7.0, 7.0]; color="#5b6472", width=3.0, label="septum blade")],
    note = "Colours mark the turn on which each beamlet was injected. The black dot is the "*
           "bumped closed orbit: it starts just inside the blade and moves away from it by "*
           "1.1 beamlet σ per turn until the bump is off. Near Q = 0.25 a beamlet comes back "*
           "to the septum side every four turns, so the tune decides how many beamlets survive.",
) do n, q
    s = paint_runs[q][n]
    series = Any[points([s.orbit], [0.0]; color="#111827", size=6.0, label="closed orbit")]
    for k in 1:n
        sel = s.tag .== k
        any(sel) || continue
        push!(series, points(s.X[sel], s.P[sel]; color=TURN_COLORS[k], size=1.8, alpha=0.75))
    end
    (series = series,
     readouts = ["bump at the septum" => string(round(s.orbit; digits=2), " σ"),
                 "particles injected" => s.injected,
                 "lost on the septum" => string(round(100*s.lost/s.injected; digits=1), " %"),
                 "stored emittance" => string(round(s.emit/σB^2; digits=1), " × one beamlet")])
end
```

### Charge-exchange (H⁻) injection

Liouville's theorem forbids stacking particles into phase space that is already
occupied — **as long as the motion is Hamiltonian**. Charge-exchange injection
gets around this by changing the particles themselves. An H⁻ beam is bent onto
the circulating proton orbit; both pass through a thin stripping foil, which
removes the two electrons, and the resulting protons simply join the circulating
beam in the same region of phase space. A bump steers the stack onto the foil during
injection and off it afterwards; any partially stripped H⁰ continue straight
to a dump.

```{figure} ../images/synchrotron/hminus_injection.png
:width: 80%
:name: fig-hminus-injection
Charge-exchange injection: the H⁻ beam and the bumped circulating protons pass through the stripping foil together; the stripped protons join the circulating beam, while partly stripped H⁰ go on to a dump.
```

Charge-exchange injection is essential for high-intensity, low-loss machines:
the best conventional multi-turn injection loses about 10 % of the beam, while
H⁻ injection can reach losses around 0.02 %. Its limits come from the foil:
heating, scattering of the circulating beam each time it crosses the foil, and
the foil's lifetime.

## Further reading

Much of the injection material here, including several of the figures, follows
the CERN Accelerator School course *Beam Injection, Extraction and Transfer*
(Erice, 2017): <https://cas.web.cern.ch/schools/erice-2017>.
