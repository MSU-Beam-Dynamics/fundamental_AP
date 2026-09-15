---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Courant–Snyder Parametrization and the FODO Cell

## Courant–Snyder parametrization

To separate the description of transverse motion into a *lattice* part and a
*particle* part we adopt the **Courant–Snyder parametrization**, after Ernest Courant
and Hartland Snyder; much of the literature calls the same quantities the
*Twiss parameters*. For a periodical
lattice the one-turn matrix is written as

$$
M=\begin{pmatrix}
\cos\Phi+\alpha\sin\Phi & \beta\sin\Phi\\
-\gamma\sin\Phi & \cos\Phi-\alpha\sin\Phi
\end{pmatrix}
$$

where $\beta(s)$, $\alpha(s)$ and $\gamma(s)$ are the **Courant–Snyder functions** that
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

### The Courant–Snyder form is the general symplectic matrix

Writing $M$ this way is not merely a change of variables. Any $2\times2$ matrix
of the Courant–Snyder form is automatically symplectic,

$$
\det M=\cos^{2}\Phi-\alpha^{2}\sin^{2}\Phi+\beta\gamma\sin^{2}\Phi=1
\qquad\text{as}\qquad \beta\gamma=1+\alpha^{2},
$$

In the following example we start from the initial condition $(1,0)$ and apply
the one-turn map repeatedly. Drag the **turn** knob to step the particle one
turn at a time, or press ▶ play to watch it hop: the start point stays marked,
the turns already taken fade into the background, and the current position is
highlighted. Every point lands on the same invariant ellipse, the one predicted
by the Courant–Snyder parameters at the observation point.

```{code-cell} julia
:tags: [hide-input]

using LinearAlgebra, TrackPadWidgets

# Every knob multiplies the frame count and each frame carries its own copy of
# the ellipse, so the grid is deliberately small: 31 turns × 3 × 3 × 6 frames.
const MAX_TURN = 30
const ELLIPSE_ANGLES = range(0, 2π, length=33)

"Round plot data to 4 significant digits — finer than a screen pixel, and it
 keeps the data embedded in the page small."
plotdata(v) = round.(v; sigdigits=4)

"""
    cs_matrix(β, α, Φ)

The one-turn map written in Courant–Snyder form. Its determinant is 1 and its trace is
2cos Φ for any (β, α), which is the whole point of the parametrisation.
"""
function cs_matrix(β, α, Φ)
    γ = (1 + α^2)/β
    return [cos(Φ) + α*sin(Φ)    β*sin(Φ)
           -γ*sin(Φ)             cos(Φ) - α*sin(Φ)]
end

"The Courant–Snyder invariant 2J = γx² + 2αxx′ + βx′², halved."
action(x, xp, β, α) = ((1 + α^2)/β * x^2 + 2α*x*xp + β*xp^2) / 2

"The ellipse of constant action J for the lattice (β, α), as (x, x′) arrays."
function invariant_ellipse(β, α, J)
    x  = @. sqrt(2J*β) * cos(ELLIPSE_ANGLES)
    xp = @. -sqrt(2J/β) * (α*cos(ELLIPSE_ANGLES) + sin(ELLIPSE_ANGLES))
    return x, xp
end

explorer(
    title   = "A Courant–Snyder-form one-turn map preserves its own ellipse",
    # `turn` is the first knob, so the ▶ play button walks the particle around
    # the ellipse; drag it by hand to step turn by turn.
    sliders = [Knob("turn", 0:MAX_TURN; fmt = n -> string(n), init = 9),
               # β below 2 would throw the α = ±2 orbit off the ±3.5 axis, since
               # the x′ amplitude of the ellipse through (1, 0) is (1+α²)/β.
               Knob("β [m]", [2.0, 5.0, 10.0]; fmt = b -> string(round(b; digits=2)), init = 2),
               Knob("α", [-2.0, 0.0, 2.0];     fmt = a -> string(round(a; digits=2)), init = 2),
               Knob("Φ [deg]", [30, 45, 60, 71, 90, 120];
                    fmt = p -> string(Int(round(p))), init = 2)],
    panels  = [Panel(xlabel="x", ylabel="x′",
                     title="turn-by-turn motion on the invariant ellipse",
                     # `equal` would stretch the wider axis to match pixel aspect —
                     # about ±7 in x on a full-width panel — so the x limit only
                     # holds with it off. x and x′ carry different units anyway.
                     xlim=(-3.5,3.5), ylim=(-3.5,3.5), height=340,
                     legend=:bottomleft)],
    note = "det M = 1 for any (β, α): changing the knobs reshapes the ellipse, but the "*
           "particle never leaves the one it started on.",
) do turn, β, α, Φ_deg
    Φ = deg2rad(Φ_deg)
    M = cs_matrix(β, α, Φ)

    # Launch at (1, 0) and apply the map `turn` times, keeping the whole history.
    position = [1.0, 0.0]
    xs = [position[1]]
    xps = [position[2]]
    for _ in 1:turn
        position = M * position
        push!(xs, position[1])
        push!(xps, position[2])
    end

    # The ellipse is fixed by the action of the STARTING point; every later turn
    # must land on it, which is what the figure demonstrates.
    J = action(xs[1], xps[1], β, α)
    ellipse_x, ellipse_xp = invariant_ellipse(β, α, J)

    # Turns 1 … n−1 are drawn faint as "where the particle has already been";
    # `max(1, end-1)` keeps the range empty rather than negative at turn 0.
    past_x  = xs[2:max(1, end-1)]
    past_xp = xps[2:max(1, end-1)]

    (series = [line(plotdata(ellipse_x), plotdata(ellipse_xp); color="#9aa4b2",
                    dash=true, label="invariant ellipse"),
               points(plotdata(past_x), plotdata(past_xp); color=PALETTE[1],
                      size=3.0, alpha=0.25, label="past turns"),
               points([1.0], [0.0]; color=PALETTE[4], size=6.0, label="start (1, 0)"),
               points(plotdata([xs[end]]), plotdata([xps[end]]); color=PALETTE[2],
                      size=6.5, label="turn $turn")],
     readouts = ["turn"     => turn,
                 "phase nΦ" => string(round(mod(turn*Φ_deg, 360); digits=1), "°")])
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

with Courant–Snyder parameters at the focusing-quad midpoint

$$
\beta_F  = \frac{M_{12}}{\sin\Phi}=\frac{2L_1\left(1+\sin(\Phi/2)\right)}{\sin\Phi},
\qquad
\alpha_F = 0 .
$$

In the following interactive example the $L_1/f$ knob scans the thin length quadrupole strength. Four views of the same cell move together: the eigenvalues of $M_{FODO}$ in the complex plane, the iterated orbit of a particle starting at $(1,0)$ to identify stable/unstable orbits, the trace of $M_{FODO}$ and the betatron tune plotted against $L_1/f$ knob across the whole scan, and the matched phase-space ellipse at the focuing quad center. 

```{code-cell} julia
:tags: [hide-input]

using LinearAlgebra, TrackPadWidgets

const N_PERIODS = 40                      # how many cells the orbit is iterated through
const CIRCLE = range(0, 2π, length=121)

"""
    fodo_matrix(ratio)

The thin-lens cell ½QF – O – QD – O – ½QF as a 2×2 map, where `ratio` = L₁/f.
Only that ratio matters, so f is fixed at 1 and the drift length carries it.
Matrices multiply right to left, i.e. in the reverse of the order the particle
meets the elements — here the cell is a palindrome, so it reads the same.
"""
function fodo_matrix(ratio)
    f, L1  = 1.0, ratio
    half_QF = [1.0     0.0
              -1/(2f)  1.0]
    QD      = [1.0     0.0
               1/f     1.0]
    drift   = [1.0     L1
               0.0     1.0]
    return half_QF * drift * QD * drift * half_QF
end

trace_of(M) = M[1,1] + M[2,2]

"Phase advance per cell from the trace; real only inside the stable band."
phase_advance(trace) = abs(trace) <= 2 ? acos(trace/2) : NaN

# The trace and tune across the whole scan, drawn once behind every frame. NaN
# outside the stable band lifts the pen rather than drawing a spurious line.
scan_ratio = collect(range(0.05, 2.1, length=205))
scan_trace = [trace_of(fodo_matrix(u)) for u in scan_ratio]
scan_tune  = [phase_advance(t)/2π for t in scan_trace]

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
    statics = [line(cos.(CIRCLE), sin.(CIRCLE); panel=1, color="#9aa4b2", dash=true,
                    label="|λ| = 1"),
               line([-3, 3], [0, 0]; panel=1, color="#9aa4b2", width=1, alpha=0.5),
               line(scan_ratio, scan_trace; panel=3, color=PALETTE[1], width=2.0,
                    label="Tr M (left)"),
               line(scan_ratio, scan_tune;  panel=3, color=PALETTE[2], width=2.0,
                    axis=:right, label="Q (right)"),
               line([0.0, 2.15, NaN, 0.0, 2.15], [2.0, 2.0, NaN, -2.0, -2.0];
                    panel=3, color=PALETTE[4], dash=true, width=1.4, label="|Tr M| = 2")],
    note = "Stable ⇔ |Tr M| ≤ 2 ⇔ both eigenvalues sit on the unit circle — the same "*
           "M_FODO derived above. The cell is a palindrome, so M₁₁ = M₂₂ and the matched "*
           "α is exactly zero here: the ellipse never tilts, it only flattens as the tune "*
           "climbs towards the half-integer at L₁/f = 2.",
) do ratio
    M      = fodo_matrix(ratio)
    trace  = trace_of(M)
    λ      = eigvals(complex(M))
    stable = abs(trace) <= 2

    # Iterate the map from (1, 0) and keep the x history, one point per cell.
    position = [1.0, 0.0]
    xs  = [position[1]]
    xps = [position[2]]
    for _ in 1:N_PERIODS
        position = M * position
        push!(xs, position[1])
        push!(xps, position[2])
    end

    Φ = phase_advance(trace)
    Q = Φ/2π
    colour = stable ? PALETTE[1] : PALETTE[4]

    # Matched ellipse through the launch point (1, 0). The cell is symmetric, so
    # α* = 0 and β* = M₁₂/sin Φ sets the aspect ratio with nothing to tilt it.
    β_match = stable ? M[1,2]/sin(Φ) : NaN
    matched_ellipse = stable ?
        [line(cos.(CIRCLE), -sin.(CIRCLE) ./ β_match; panel=4, color="#9aa4b2",
              dash=true, label="matched ellipse")] : []

    (series = vcat(
        [points(real.(λ), imag.(λ); panel=1, color=colour, size=6.0, label="λ₊, λ₋"),
         line(0:N_PERIODS, xs; panel=2, color=colour, label="xₙ"),
         line([ratio, ratio], [-2.8, 2.8]; panel=3, color="#9aa4b2", width=1.4),
         points([ratio], [trace]; panel=3, color=PALETTE[1], size=5.5),
         points([ratio], [Q];     panel=3, color=PALETTE[2], size=5.5, axis=:right)],
        matched_ellipse,
        [points(xs, xps; panel=4, color=colour, size=3.0, alpha=0.55,
                label="turn-by-turn"),
         points([1.0], [0.0]; panel=4, color=PALETTE[4], size=6.0, label="start (1, 0)")]),
     readouts = ["Tr M"    => round(trace; digits=3),
                 "|λ|max"  => round(maximum(abs.(λ)); digits=4),
                 "Φ/cell"  => stable ? string(round(rad2deg(Φ); digits=1), "°") : "—",
                 "tune Q"  => stable ? string(round(Q; digits=4)) : "—",
                 "β* [m]"  => stable ? string(round(β_match; sigdigits=4)) : "—",
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
beam = Beam(3.0e9)                          # 3 GeV kinetic energy, electron

const CELL_LENGTH  = 4.0                    # one FODO period [m]
const FOCAL_LENGTH = 1.8                    # the thin lens the thick quad imitates [m]
const QUAD_LENGTH  = 0.4                    # [m]

# A thin lens of focal length f becomes a thick quadrupole of the same integrated
# strength when k·ℓ = 1/f.
quad_strength = 1/(FOCAL_LENGTH * QUAD_LENGTH)

# The cell holds TWO quadrupoles, so each drift is (cell − 2ℓ_q)/2 and the period
# comes out at exactly CELL_LENGTH.
drift_length = (CELL_LENGTH - 2*QUAD_LENGTH)/2

ring = Lattice(AbstractElement[
    Quadrupole(QUAD_LENGTH, +quad_strength; name=:QF),
    Drift(drift_length),
    Quadrupole(QUAD_LENGTH, -quad_strength; name=:QD),
    Drift(drift_length),
]; name=:FODO, periodic=true)

# By default the optics come back at element boundaries only — five points for
# this four-element cell, which draws β as four straight segments. The sampling
# keywords slice the lattice for optics: `sample_integrator_steps` splits thick
# multipoles at their integration steps (10 per quadrupole here) and `max_step`
# caps the length of drift and bend pieces.
twiss = periodic_twiss(ring, beam; sample_integrator_steps=true, max_step=0.05)

println("tunex = ", twiss.tunex, "   tuney = ", twiss.tuney)
println("optics sampled at ", length(twiss.s), " points over ",
        round(twiss.s[end]; digits=2), " m")

# The thin-lens prediction sin(Φ/2) = L₁/(2f), with L₁ the half-cell length.
phase_advance_thin_lens = 2asin((CELL_LENGTH/2)/(2*FOCAL_LENGTH))
println("phase advance/cell: TrackPad = ", round(twiss.tunex*2π; digits=4),
        " rad ; thin-lens model = ", round(phase_advance_thin_lens; digits=4), " rad")

# Translucent glyphs, so the beamline reads as a background band and never hides
# a curve that dips into it. Unlisted element kinds fall back to :element.
const STRIP_COLORS = Dict{Symbol,Any}(
    :quadrupole => (:seagreen,  0.55),
    :bend       => (:steelblue, 0.55),
    :sextupole  => (:tomato,    0.55),
    :rf_cavity  => (:goldenrod, 0.55),
    :element    => (:gray,      0.45),
)

"""
    twiss_figure(lat, datamax; figsize, ylabel, title, strip=0.09, headroom=0.10)

One axis carrying both the curves and the beamline. `plot_lattice!` draws the
element glyphs inside the box, in a band across the top `strip` fraction of it:
focusing quadrupoles above the band's midline, defocusing below, drifts left to
the baseline. The y range is set to `[0, top]` with

    top = datamax / (1 - strip - headroom)

so the band always has room of its own above the data, separated from it by
`headroom`. Returns `(fig, ax)`; draw into `ax` afterwards, so the curves land
on top of the glyphs.
"""
function twiss_figure(lat, datamax; figsize=(800, 380), ylabel="", title="",
                      strip=0.09, headroom=0.10)
    fig = Figure(size=figsize)
    ax  = Axis(fig[1, 1]; xlabel="s [m]", ylabel, title)
    top = datamax / (1 - strip - headroom)

    # The band spans [top(1−strip), top]: baseline on its midline, glyphs ±height.
    plot_lattice!(ax, lat; baseline=top*(1 - strip/2), height=top*strip/2,
                  colors=STRIP_COLORS)

    # plot_lattice! sets limits to suit the strip alone; restore ours afterwards.
    xlims!(ax, 0, total_length(lat))
    ylims!(ax, 0, top)
    return fig, ax
end

beta_max = maximum(max.(twiss.betax, twiss.betay))
fig, ax = twiss_figure(ring, beta_max; ylabel="β [m]",
                       title="FODO cell Courant–Snyder functions")
lines!(ax, twiss.s, twiss.betax; linewidth=2, label="βₓ")
lines!(ax, twiss.s, twiss.betay; linewidth=2, label="βᵧ")
axislegend(ax; position=:rb)      # the top of the box now belongs to the beamline
fig
```


Reading the curves, we can see that:
- $\beta_x$ peaks inside the focusing quadrupole and reaches its minimum inside the defocusing one, 
- $\beta_y$ does the opposite

This is the alternating-gradient principle, visible directly below the magnets that cause it.



```{code-cell} julia
:tags: [hide-input]

# The betatron phase accumulates monotonically, fastest where β is smallest —
# that is the statement dμ/ds = 1/β. Over one period it climbs to 2πQ.
phase_max = max(maximum(twiss.mux), maximum(twiss.muy))
fig, ax = twiss_figure(ring, phase_max; ylabel="betatron phase [rad]")
lines!(ax, twiss.s, twiss.mux; linewidth=2, label="μₓ")
lines!(ax, twiss.s, twiss.muy; linewidth=2, label="μᵧ")
axislegend(ax; position=:rb)
fig
```

The turn-by-turn motion of a single particle draws the matched phase-space
ellipse predicted by $(\beta, \alpha)$ at the observation point:




