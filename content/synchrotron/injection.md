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
cells, 80 m around (about the size of a medical synchrotron). Every
plot is taken at the centre of a focusing quadrupole, where $\alpha=0$, so
$\bar x'=\sqrt\beta\,p_x$. For a thin-lens FODO cell
$\sin(\mu/2)=L/(2f)$, so the quadrupole strength that gives any chosen tune
follows directly:

```{code-cell} julia
using StaticArrays, TrackPad, TrackPadWidgets, Random, Statistics

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
println("tune  = ", gettune(ring, beam))
println("β     = ", round(beta_start(ring); digits=2), " m at the observation point")
println("chromaticity ξ = ", round.(getchrom(ring, beam; centered=true); digits=3))
```

The computed optics confirm the tune of 1.27 that the quadrupoles were set for. Without
any sextupoles, the ring's natural chromaticity is about −2; this matters
below, where it is what makes an injection error grow into a larger emittance.

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
rotation by 90°, and a kick back to the origin. The widget below lets you change
the phase advance of the cell between the two magnets:

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

A **kicker** is a fast pulsed magnet, usually with a ferrite yoke. Its rise and
fall times must fit into the gaps between bunches, and its flat top must be as
long as the injected bunch train. This is why the filling pattern of a ring
leaves gaps for the kicker.

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
with, so the emittance has grown. An optical **mismatch** (a beam injected with
the wrong $\beta$ or $\alpha$; see [](../transverse_dynamics/emittance.md))
filaments the same way into a larger ellipse of the matched shape. Both final
emittances follow from averaging over the betatron phase:

$$
\epsilon_{\mathrm{dipole}}=\epsilon_0+\tfrac12\,\Delta\bar A^2,
\qquad
\epsilon_{\mathrm{mismatch}}=\tfrac12\left(\beta\gamma_0-2\alpha\alpha_0+\gamma\beta_0\right)\epsilon_0,
$$

where $\Delta\bar A$ is the normalized amplitude of the injection error, and
$(\beta_0,\alpha_0,\gamma_0)$ describe the injected beam while
$(\beta,\alpha,\gamma)$ describe the lattice. Here the tune spread comes from
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

Charge-exchange injection is essential for high-intensity, low-loss machines:
the best conventional multi-turn injection loses about 10 % of the beam, while
H⁻ injection can reach losses around 0.02 %. Its limits come from the foil:
heating, scattering of the circulating beam each time it crosses the foil, and
the foil's lifetime.
