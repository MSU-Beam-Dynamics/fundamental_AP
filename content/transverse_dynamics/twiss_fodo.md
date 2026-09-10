---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Twiss Parametrization and the FODO Cell

## Twiss parametrization

To separate the description of transverse motion into a *lattice* part and a
*particle* part we adopt the **Twiss parametrization**. For a periodical
lattice the one-turn matrix is written as

$$
M=\begin{pmatrix}
\cos\Phi+\alpha\sin\Phi & \beta\sin\Phi\\
-\gamma\sin\Phi & \cos\Phi-\alpha\sin\Phi
\end{pmatrix}
$$

where $\beta(s)$, $\alpha(s)$ and $\gamma(s)$ are the **Twiss functions** that
depend only on the lattice, and $\Phi$ is the phase advance per period. They
satisfy

$$
\beta \gamma = 1 + \alpha^2,
$$

with $\beta$ measured in metres and $\alpha$ dimensionless. Turn-by-turn, the
particle coordinates evolve as

$$
\begin{pmatrix}x\\x'\end{pmatrix}_{n+1}=M\,\begin{pmatrix}x\\x'\end{pmatrix}_{n}
$$

and, if $\beta(s)$ is known everywhere,

$$
x(s)=\sqrt{2J \beta(s)}\;\cos\!\left(\psi(s)+\psi_0\right)
$$

The constant $J$ is the ***action*** of the particle, determined by its initial
condition; $\psi$ is the betatron phase advance from the origin. The pair
$(\beta, \alpha)$ fixes the shape and orientation of the ellipse traced in
phase space, while $J$ sets its area: each turn maps the point along an
ellipse without changing it.

Given a one-turn matrix at location $s_0$, the optical functions follow from

$$
\cos\Phi=\frac{m_{11}+m_{22}}{2},\qquad
\beta=\frac{m_{12}}{\sin\Phi},\qquad
\alpha=\frac{m_{11}-m_{22}}{2\sin\Phi},
$$

and optics at another location follow from
$M(s_1)=M(s_1|s_0)\,M(s_0)\,M(s_0|s_1)$, where $M(s_i|s_j)$ transfers from
$s_j$ to $s_i$. Explicitly, the transfer between two locations reads

$$
M\left(s_{1}\mid s_{0}\right)
=\begin{pmatrix}
\sqrt{\frac{\beta_{1}}{\beta_{0}}}\left(\cos\psi+\alpha_{0}\sin\psi\right) & \sqrt{\beta_{0}\beta_{1}}\sin\psi\\
-\frac{1+\alpha_{0}\alpha_{1}}{\sqrt{\beta_{0}\beta_{1}}}\sin\psi+\frac{\alpha_{0}-\alpha_{1}}{\sqrt{\beta_{0}\beta_{1}}}\cos\psi & \sqrt{\frac{\beta_{0}}{\beta_{1}}}\left(\cos\psi-\alpha_{1}\sin\psi\right)
\end{pmatrix}.
$$

### The Twiss form is the general symplectic matrix

Writing $M$ this way is not merely a change of variables. Any $2\times2$ matrix
of the Twiss form is automatically symplectic,

$$
\det M=\cos^{2}\Phi-\alpha^{2}\sin^{2}\Phi+\beta\gamma\sin^{2}\Phi=1
\qquad\text{as}\qquad \beta\gamma=1+\alpha^{2},
$$

In the following example we start from the initial condition $(1,0)$ and apply
the one-turn map repeatedly. Drag the **turn** knob to step the particle one
turn at a time, or press ▶ play to watch it hop: the start point stays marked,
the turns already taken fade into the background, and the current position is
highlighted. Every point lands on the same invariant ellipse, the one predicted
by the Twiss parameters at the observation point.

```{code-cell} julia
:tags: [hide-input]

using LinearAlgebra, TrackPadWidgets

# Every knob multiplies the frame count, and each frame carries its own copy of
# the ellipse, so the grid is deliberately small: 21 turns × 3 × 3 × 5 = 945 frames.
const NTURN = 30
const θE    = range(0, 2π, length=33)     # the invariant ellipse, drawn once per frame
r4t(v) = round.(v; sigdigits=4)           # 4 digits is finer than a pixel here

explorer(
    title   = "A Twiss-form one-turn map preserves its own ellipse",
    # `turn` is the first knob, so the ▶ play button walks the particle around
    # the ellipse; drag it by hand to step turn by turn.
    sliders = [Knob("turn", 0:NTURN; fmt = n -> string(n), init = 9),
               # β below 2 would throw the α = ±2 orbit off the ±3.2 axis, since the
               # x′ amplitude of the ellipse through (1, 0) is (1+α²)/β.
               Knob("β [m]", [2.0, 5.0, 10.0]; fmt = b -> string(round(b;digits=2)), init = 2),
               Knob("α", [-2.0, 0.0, 2.0]; fmt = a -> string(round(a;digits=2)), init = 2),
               Knob("Φ [deg]", [30, 45, 60, 71, 90, 120]; fmt = p -> string(Int(round(p))), init = 2)],
    panels  = [Panel(xlabel="x", ylabel="x′", title="turn-by-turn motion on the invariant ellipse",
                     # `equal` would stretch the wider axis to match pixel aspect —
                     # about ±7 in x on a full-width panel — so the x limit only
                     # holds with it off. x and x′ carry different units anyway.
                     xlim=(-3.5,3.5), ylim=(-3.5,3.5), height=340,
                     legend=:bottomleft)],
    note = "Matrix has det M = 1, for any (β, α): the ellipse shape change, the particle will remain on the same ellipse.",
) do n, βt, αt, Φdeg
    Φ  = deg2rad(Φdeg)
    γt = (1+αt^2)/βt
    M  = [cos(Φ)+αt*sin(Φ)  βt*sin(Φ); -γt*sin(Φ)  cos(Φ)-αt*sin(Φ)]

    # Iterate from (1, 0) and keep the whole history up to the selected turn.
    v = [1.0, 0.0]
    xs = Float64[v[1]]; ps = Float64[v[2]]
    for _ in 1:n
        v = M*v
        push!(xs, v[1]); push!(ps, v[2])
    end

    J  = (γt*xs[1]^2 + 2αt*xs[1]*ps[1] + βt*ps[1]^2)/2
    ex = @. sqrt(2J*βt)*cos(θE); ep = @. -(αt*cos(θE)+sin(θE))*sqrt(2J/βt)
    Jn = (γt*xs[end]^2 + 2αt*xs[end]*ps[end] + βt*ps[end]^2)/2

    (series = [line(r4t(ex), r4t(ep); color="#9aa4b2", dash=true, label="invariant ellipse"),
               # turns 1 … n−1: where the particle has already been
               points(r4t(xs[2:max(1, end-1)]), r4t(ps[2:max(1, end-1)]);
                      color=PALETTE[1], size=3.0, alpha=0.25, label="past turns"),
               # the start, always on the plot
               points([1.0], [0.0]; color=PALETTE[4], size=6.0, label="start (1, 0)"),
               # where the particle is now
               points(r4t([xs[end]]), r4t([ps[end]]); color=PALETTE[2], size=6.5,
                      label="turn $n")],
     readouts = ["turn"        => n,
                 "phase nΦ"    => string(round(mod(n*Φdeg, 360); digits=1), "°"),
                 ])
end
```





## Example: FODO cell

A **FODO cell** is the workhorse of accelerator lattices: focusing quad –
drift – defocusing quad – drift. For a symmetric cell starting at the middle
of the focusing quadrupole,

$$
\tfrac12\text{QF}-\text{O}-\text{QD}-\text{O}-\tfrac12\text{QF},
$$

and using thin-lens quads ($\pm f$) with equal drifts $L_1$, the cell matrix is

$$
\begin{aligned}
M_{FODO}
&=
\begin{pmatrix}1 & 0\\ -\frac{1}{2f} & 1\end{pmatrix}
\begin{pmatrix}1 & L_1\\ 0 & 1\end{pmatrix}
\begin{pmatrix}1 & 0\\ \frac{1}{f} & 1\end{pmatrix}
\begin{pmatrix}1 & L_1\\ 0 & 1\end{pmatrix}
\begin{pmatrix}1 & 0\\ -\frac{1}{2f} & 1\end{pmatrix}\\[2ex]
&=
\begin{pmatrix}
1-\frac{L_{1}^{2}}{2 f^{2}} & \frac{L_{1}}{f} \left(L_{1} + 2 f\right)\\[1ex]
\frac{L_{1}}{4 f^{3}} \left(L_{1} - 2 f\right) & 1-\frac{L_{1}^{2}}{2 f^{2}}
\end{pmatrix}.
\end{aligned}
$$

The determinant is always 1 (due to symplecticity, as expected), while stability requires
$|\mathrm{Tr}(M)| = |2 - L_1^2/f^2|\le2$.  Therefore $L_1/f < 2$ is required for stability. The phase advance per cell follows:

$$
\cos{\Phi}=1-\frac{L_1^2}{2f^2}
\quad\Leftrightarrow\quad
\sin{\frac{\Phi}{2}}=\frac{L_1}{2f},
$$

with Twiss parameters at the focusing-quad midpoint

$$
\beta_F  = \frac{M_{12}}{\sin\Phi}=\frac{2L_1\left(1+\sin(\Phi/2)\right)}{\sin\Phi},
\qquad
\alpha_F = 0 .
$$

In the following interactive example the $L_1/f$ knob scans the thin length quadrupole strength. Four views of the same cell move together: the eigenvalues of $M_{FODO}$ in the complex plane, the iterated orbit of a particle starting at $(1,0)$ to identify stable/unstable orbits, the trace of $M_{FODO}$ and the betatron tune plotted against $L_1/f$ knob across the whole scan, and the matched phase-space ellipse at the focuing quad center. 

```{code-cell} julia
:tags: [hide-input]

using LinearAlgebra, TrackPadWidgets

nper = 40
θe = range(0, 2π, length=121)

"The thin-lens cell ½QF–O–QD–O–½QF as a 2×2 map; only the ratio L₁/f matters."
function mfodo(u)
    f, L1  = 1.0, u
    halfQF = [1.0 0.0; -1/(2f) 1.0]
    QD     = [1.0 0.0;  1/f    1.0]
    O      = [1.0 L1;   0.0    1.0]
    halfQF * O * QD * O * halfQF
end

# Trace and tune across the whole scan, drawn once behind every frame. The cell
# is a palindrome, so M₁₁ = M₂₂ and the matched α is zero at this point; the
# tune is real only while |Tr M| ≤ 2, and NaN elsewhere lifts the pen.
ufine  = collect(range(0.05, 2.1, length=205))
trfine = [(M = mfodo(u); M[1,1] + M[2,2]) for u in ufine]
qfine  = [abs(t) <= 2 ? acos(t/2)/2π : NaN for t in trfine]

explorer(
    title   = "Thin-lens FODO cell ½QF–O–QD–O–½QF, iterated",
    sliders = [Knob("L₁ / f", range(0.05, 2.1, length=50);
                      fmt = u -> string(round(u; digits=2)), init = 20)],
    panels  = [Panel(xlabel="Re λ", ylabel="Im λ", title="eigenvalues of M_FODO",
                     xlim=(-3.0, 3.0), ylim=(-1.6, 1.6), equal=true, height=250),
               Panel(xlabel="period number n", ylabel="xₙ / x₀",
                     title="iterated orbit", ylim=(-6.0, 6.0), height=250),
               Panel(xlabel="L₁ / f", ylabel="Tr M", y2label="tune Q = Φ/2π",
                     title="trace and tune versus L₁/f",
                     xlim=(0.0, 2.15), ylim=(-2.8, 2.8), y2lim=(0.0, 0.56),
                     height=250, legend=:bottomleft),
               Panel(xlabel="x", ylabel="x′",
                     title="phase-space orientation at the cell centre",
                     xlim=(-1.5, 1.5), ylim=(-0.6, 0.6), height=250,
                     legend=:bottomleft)],
    statics = [line(cos.(θe), sin.(θe); panel=1, color="#9aa4b2", dash=true, label="|λ| = 1"),
               line([-3, 3], [0, 0]; panel=1, color="#9aa4b2", width=1, alpha=0.5),
               line(ufine, trfine; panel=3, color=PALETTE[1], width=2.0, label="Tr M (left)"),
               line(ufine, qfine;  panel=3, color=PALETTE[2], width=2.0, axis=:right,
                    label="Q (right)"),
               line([0.0, 2.15, NaN, 0.0, 2.15], [2.0, 2.0, NaN, -2.0, -2.0];
                    panel=3, color=PALETTE[4], dash=true, width=1.4, label="|Tr M| = 2")],
    note = "Stable ⇔ |Tr M| ≤ 2 ⇔ both eigenvalues sit on the unit circle — the same "*
           "M_FODO derived above. The cell is a palindrome, so M₁₁ = M₂₂ and the matched "*
           "α is exactly zero here: the ellipse never tilts, it only flattens as the tune "*
           "climbs towards the half-integer at L₁/f = 2.",
) do u
    M  = mfodo(u)
    tr = M[1,1] + M[2,2]
    λ  = eigvals(complex(M))

    v = [1.0, 0.0]; xs = Float64[v[1]]; ps = Float64[v[2]]
    for _ in 1:nper
        v = M * v
        push!(xs, v[1]); push!(ps, v[2])
    end

    stable = abs(tr) <= 2
    col = stable ? PALETTE[1] : PALETTE[4]
    Q   = stable ? acos(tr/2)/2π : NaN

    # Matched ellipse through the launch point (1, 0). α* = 0 by symmetry, so
    # β* = M₁₂/sin Φ sets the aspect and nothing sets a tilt.
    βs  = stable ? M[1,2]/sin(acos(tr/2)) : NaN
    ellipse = stable ?
        [line(cos.(θe), -sin.(θe) ./ βs; panel=4, color="#9aa4b2", dash=true,
              label="matched ellipse")] : []

    (series = vcat(
        [points(real.(λ), imag.(λ); panel=1, color=col, size=6.0, label="λ₊, λ₋"),
         line(0:nper, xs; panel=2, color=col, label="xₙ"),
         line([u, u], [-2.8, 2.8]; panel=3, color="#9aa4b2", width=1.4),
         points([u], [tr]; panel=3, color=PALETTE[1], size=5.5),
         points([u], [stable ? Q : NaN]; panel=3, color=PALETTE[2], size=5.5, axis=:right)],
        ellipse,
        [points(xs, ps; panel=4, color=col, size=3.0, alpha=0.55, label="turn-by-turn"),
         points([1.0], [0.0]; panel=4, color=PALETTE[4], size=6.0, label="start (1, 0)")]),
     readouts = ["Tr M"    => round(tr; digits=3),
                 "|λ|max"  => round(maximum(abs.(λ)); digits=4),
                 "Φ/cell"  => stable ? string(round(rad2deg(acos(tr/2)); digits=1), "°") : "—",
                 "tune Q"  => stable ? string(round(Q; digits=4)) : "—",
                 "β* [m]"  => stable ? string(round(βs; sigdigits=4)) : "—",
                 "α*"      => stable ? "0 (symmetry point)" : "—",
                 "verdict" => stable ? "stable" : "unstable"])
end
```

Two things are worth noticing. Inside the stable band the two eigenvalues are a
complex-conjugate pair $e^{\pm i\Phi}$ and the iterated orbit is a sampled
cosine of exactly that phase advance. At $L_1/f = 2$ the pair collides at
$\lambda=-1$ and then splits along the real axis: one eigenvalue larger than
one in modulus, and the orbit grows exponentially. The matched ellipse at the focusing quad center is always upright, due to the symmetry of the lattice.

## Simulating a FODO ring with TrackPad

Build a realistic thick-lens FODO ring and let TrackPad compute the periodic
solution of optics function:

```{code-cell} julia
:tags: [hide-input]

using StaticArrays, TrackPad, CairoMakie

# Reference beam for every tracking cell on this page. Defined at top level, in
# the first cell that needs it, so the later cells share it.
beam = Beam(3.0e9)                        # 3 GeV kinetic-energy electron

L_cell, f_eff, L_q = 4.0, 1.8, 0.4         # cell length, focal length, quad length [m]

# Thick quads approximating the thin lens f: k·L_q ≈ 1/f. The cell holds *two*
# quadrupoles, so each drift is (L_cell − 2L_q)/2 and the period is exactly L_cell.
kq = 1/(f_eff*L_q)
ring = Lattice(AbstractElement[
    Quadrupole(L_q, +kq; name=:QF),
    Drift((L_cell - 2L_q)/2),
    Quadrupole(L_q, -kq; name=:QD),
    Drift((L_cell - 2L_q)/2),
]; name=:FODO, periodic=true)

# By default the optics come back at element boundaries only — five points for
# this four-element cell, which draws β as four straight segments. The sampling
# keywords slice the lattice for optics: `sample_integrator_steps` splits thick
# multipoles at their integration steps (10 per quadrupole here) and `max_step`
# caps the length of drift and bend pieces.
tw = periodic_twiss(ring, beam; sample_integrator_steps=true, max_step=0.05)
println("tunex = ", tw.tunex, "   tuney = ", tw.tuney)
println("optics sampled at ", length(tw.s), " points over ", round(tw.s[end]; digits=2), " m")

# analytic thin-lens prediction for comparison
Φ_analytic = 2asin((L_cell/2)/(2*f_eff))
println("phase advance/cell: TrackPad = ", round(tw.tunex*2π; digits=4),
        " rad ; thin-lens model = ", round(Φ_analytic; digits=4), " rad")

# Translucent glyphs, so the beamline reads as a background band and never
# hides a curve that dips into it. Unlisted element kinds fall back to :element.
const STRIP_COLORS = Dict{Symbol,Any}(
    :quadrupole => (:seagreen,    0.55),
    :bend       => (:steelblue, 0.55),
    :sextupole  => (:tomato,  0.55),
    :rf_cavity  => (:goldenrod, 0.55),
    :element    => (:gray,      0.45),
)

"""
    twiss_figure(lat, datamax; figsize, ylabel, title, strip=0.09, headroom=0.10)

One axis carrying both the curves and the beamline. `plot_lattice!` draws the
element glyphs inside the box, in a band across the top `strip` fraction of it:
focusing quadrupoles above the band's midline, defocusing below, drifts left to
the baseline. The y range is set to `[0, a]` with

    a = datamax / (1 - strip - headroom)

so the band always has room of its own above the data, separated from it by
`headroom`. Returns `(fig, ax)`; draw into `ax` afterwards, so the curves land
on top of the glyphs.
"""
function twiss_figure(lat, datamax; figsize=(800, 380), ylabel="", title="",
                      strip=0.09, headroom=0.10)
    fig = Figure(size=figsize)
    ax  = Axis(fig[1, 1]; xlabel="s [m]", ylabel, title)
    a   = datamax / (1 - strip - headroom)
    # band spans [a(1−strip), a]: baseline on its midline, glyphs ±height about it
    plot_lattice!(ax, lat; baseline=a*(1 - strip/2), height=a*strip/2,
                  colors=STRIP_COLORS)
    # plot_lattice! sets limits to suit the strip alone; restore ours afterwards
    xlims!(ax, 0, total_length(lat))
    ylims!(ax, 0, a)
    return fig, ax
end

βmax = maximum(max.(tw.betax, tw.betay))
fig, ax = twiss_figure(ring, βmax; ylabel="β [m]", title="FODO cell Twiss functions")
lines!(ax, tw.s, tw.betax; linewidth=2, label="βₓ")
lines!(ax, tw.s, tw.betay; linewidth=2, label="βᵧ")
axislegend(ax; position=:rb)          # top of the box now belongs to the beamline
fig
```


Reading the curves, we can see that:
- $\beta_x$ peaks inside the focusing quadrupole and reaches its minimum inside the defocusing one, 
- $\beta_y$ does the opposite

This is the alternating-gradient principle, visible directly below the magnets that cause it.



```{code-cell} julia
:tags: [hide-input]

# phase advances accumulate monotonically around the ring:
μmax = max(maximum(tw.mux), maximum(tw.muy))
fig, ax = twiss_figure(ring, μmax; ylabel="betatron phase [rad]")
lines!(ax, tw.s, tw.mux; linewidth=2, label="μₓ")
lines!(ax, tw.s, tw.muy; linewidth=2, label="μᵧ")
axislegend(ax; position=:rb)
fig
```

The turn-by-turn motion of a single particle draws the matched phase-space
ellipse predicted by $(\beta, \alpha)$ at the observation point:




