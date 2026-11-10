---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Lattice Design and the Interaction Region

## Arcs and insertions

A collider ring is built from two kinds of sections.

- The **arcs** bend the beam around the ring. They are long chains of identical
  cells, usually the FODO cell of [](../transverse_dynamics/twiss_fodo.md), with
  dipoles between the quadrupoles and sextupoles to correct the chromaticity.
- The **insertions** are straight sections for everything else: RF cavities,
  injection and beam dump, collimation and, above all, the **interaction
  regions** (IR) where the beams collide in the middle of a detector.

The luminosity $\mathcal L\propto1/(\sigma_x\sigma_y)$ found in
[](luminosity.md) asks for the smallest possible beams at the interaction
point (IP). With a momentum spread $\delta$, the beam size at the IP is

$$
\sigma=\sqrt{\epsilon\beta^*+D^{*2}\delta^2}.
$$

An IR therefore has two jobs:

- **low β\***, to squeeze the betatron part $\sqrt{\epsilon\beta^*}$;
- **zero dispersion** $D^*$, because even a small dispersion would dominate.
  For protons at 250 GeV with a normalized emittance of 2.5 μm and
  $\beta^*=1$ m, $\sqrt{\epsilon\beta^*}\approx97\ \mu$m. A dispersion of 0.5 m
  with $\delta=5\times10^{-4}$ would add 250 μm. *Dispersion suppressors* at
  the ends of the arcs bring $D$ and $D'$ to zero before the beam enters the
  IR.

RHIC, for example, has six straight sections with six crossing points. In its
proton run, two of them housed the experiments with $\beta^*=0.65$ m, while the
others ran at a relaxed $\beta^*$ of several metres.

## A model collider

Here and on the next page we use a small model collider with one IP. It has the
features that matter for the optics:

- an **arc** of twelve 30 m FODO cells, standing in for the bending arcs
  (dipoles and dispersion are left out, since the IR is dispersion-free
  anyway);
- on each side of the IP, a free **drift** of $L^*=12$ m for the detector;
- a **final-focus triplet** of quadrupoles (Q1, Q2, Q3; Q1 and Q3 are equal);
- four **matching quadrupoles** (Q4–Q7) that connect the IR optics to the arc.

The ring is mirror-symmetric about the IP. That makes $\alpha=0$ at the IP
automatically, so only $\beta^*$ has to be set there.

```{code-cell} julia
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
```

The quadrupole strengths come from **matching**. This is the everyday task of
lattice design: adjust magnet strengths until the optics meet a list of
targets. Here there are six targets:

- $\beta_x=\beta_y=\beta^*$ at the IP;
- $(\beta,\alpha)$ at the arc entrance equal to the periodic values of the FODO
  cell in both planes, so the IR is invisible to the arcs and no
  β-beating leaks out;
- the two tunes of the ring.

The six IR quadrupoles handle the first group, with a least-squares solver. The
arc phase advance is then trimmed to put the tunes where they are wanted, and
the two steps repeat until both are satisfied.

```{code-cell} julia
:tags: [hide-input]

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
nothing
```

```{code-cell} julia
const QX, QY = 3.685, 3.695           # RHIC-like working point
nominal = match_collider(1.0, QX, QY; seed=SEED)
ring = collider_lattice(nominal)
tw = periodic_twiss(ring, beam; sample_integrator_steps=true, max_step=0.5)

println(join([
    @sprintf("circumference   %.1f m", total_length(ring)),
    @sprintf("tunes           (%.4f, %.4f)", gettune(ring, beam)...),
    @sprintf("β* at the IP    (%.3f, %.3f) m,  α* = (%.1e, %.1e)", tw.betax[1], tw.betay[1], tw.alphax[1], tw.alphay[1]),
    @sprintf("arc phase/cell  (%.1f°, %.1f°)", rad2deg(nominal.μx), rad2deg(nominal.μy)),
    @sprintf("β max           (%.0f, %.0f) m in the triplet", maximum(tw.betax), maximum(tw.betay))], "\n"))
```

The matched ring is 519 m around, with the IP at $\beta^*=1$ m and the requested
tunes. The figure shows the optics around the whole ring, with the IP in the
middle. The FODO arcs carry small, regular β functions of 10–50 m. The two
halves of the IR carry the large peaks, as in the RHIC lattice: on each side of
the IP, $\beta$ climbs to several hundred metres in the triplet.

```{code-cell} julia
:tags: [hide-input]
:label: fig-collider-optics
:caption: Courant–Snyder functions around the model collider with β* = 1 m, plotted with the IP in the middle. In the beamline strip at the top, quadrupoles above the line focus horizontally and those below defocus.

const STRIP_COLORS = Dict{Symbol,Any}(:quadrupole => (:seagreen, 0.6), :element => (:gray, 0.45))

"The same ring, but starting in the middle of the arc so that the IP sits at s = C/2."
function centered_lattice(c)
    half_arc = AbstractElement[(e for _ in 1:NARC÷2 for e in arc_cell(c.kf, c.kd))...]
    Lattice(AbstractElement[half_arc..., reverse(ir_half(c.k))..., ir_half(c.k)..., half_arc...]; periodic=true)
end

let lat = centered_lattice(nominal), C = total_length(ring)
    t = periodic_twiss(lat, beam; sample_integrator_steps=true, max_step=0.5)
    top = 470/(1 - 0.09 - 0.06)
    fig = Figure(size=(980, 400), fontsize=15)
    ax = Axis(fig[1, 1]; xlabel="distance from the IP [m]", ylabel="β [m]",
              xticks=(C/2 .+ (-250:50:250), string.(-250:50:250)))
    plot_lattice!(ax, lat; baseline=top*(1 - 0.06), height=top*0.04, colors=STRIP_COLORS)
    xlims!(ax, 0, C); ylims!(ax, 0, top)
    lines!(ax, t.s, t.betax; color=PALETTE[1], linewidth=2, label="βₓ")
    lines!(ax, t.s, t.betay; color=PALETTE[2], linewidth=2, label="βᵧ")
    vlines!(ax, [C/2]; color=:gray50, linestyle=:dash)
    text!(ax, C/2, 0.86top; text="IP", align=(:center, :bottom))
    for (x, label) in [(C/2 - 65, "IR"), (C/2 + 65, "IR"), (C/2 - 170, "arc"), (C/2 + 170, "arc")]
        text!(ax, x, 0.76top; text=label, align=(:center, :center), color=:gray35)
    end
    axislegend(ax; position=:lc)
    fig
end
```

## The free space around the IP

Between the IP and the first quadrupole the beam is in a drift. A drift
transports the Courant–Snyder parameters as

$$
\beta(s)=\beta^*-2\alpha^*s+\gamma^*s^2
=\beta^*+\frac{s^2}{\beta^*},
\qquad
\alpha(s)=-\frac{s}{\beta^*},
$$

because $\alpha^*=0$ and $\gamma^*=1/\beta^*$ at a symmetric IP. Within one
$\beta^*$ of the IP, β has only doubled. After that it grows *quadratically*, so
at the first quadrupole

$$
\beta_{Q1}\approx\frac{L^{*2}}{\beta^*}.
$$

Halving $\beta^*$ doubles the luminosity, and also doubles β in the final focus.
The phase advance over the drift is $\psi(s)=\arctan(s/\beta^*)$. It approaches
90° within a few $\beta^*$ of the IP, which is why a beam that crosses the IP at
an angle is displaced almost entirely in the triplet.

The widget below squeezes $\beta^*$ from 10 m to 0.25 m. Every step is a full
re-match of the IR quadrupoles and the tunes.

- **Top panel:** β functions from the IP through the triplet and the matching
  section into the first arc cells. β is on a logarithmic scale, so that
  $\beta^*$ stays visible next to the peaks.
- **Bottom panel:** what this means for the beams near the IP. Each beam is
  drawn as its centroid and a ±3σ envelope. The beams cross at the angle that
  keeps them ten beam sizes apart in the drift, a typical requirement to keep
  the long-range beam–beam forces small (see [](beam_beam.md)).

```{code-cell} julia
:tags: [hide-input]

const BSTARS = [10.0, 6.0, 4.0, 3.0, 2.0, 1.5, 1.0, 0.7, 0.5, 0.35, 0.25]
SQUEEZE = let out = Dict(1.0 => nominal)
    for dir in (filter(>(1.0), BSTARS), reverse(filter(<(1.0), BSTARS)))
        c = nominal
        for b in sort(dir; by = b -> abs(log(b)))
            c = squeeze(c, b, QX, QY); out[b] = c
        end
    end
    out
end

const SPAN = 140.0                   # metres shown on each side of the IP (top panel)
const NEAR = 40.0                    # bottom panel

"Mirror a glyph list to negative s, for the left half of the symmetric IR."
mirror(gs) = [(; kind=g.kind, height=g.height, centered=g.centered,
                 plot_start=-g.plot_end, plot_end=-g.plot_start) for g in gs]

"Natural chromaticity of the arc alone (NARC cells) and of the whole ring."
function chroma(c)
    arc = Lattice(AbstractElement[(e for _ in 1:NARC for e in arc_cell(c.kf, c.kd))...]; periodic=true)
    getchrom(arc, beam; centered=true), getchrom(collider_lattice(c), beam; centered=true)
end

explorer(
    title   = "Squeezing β* at the interaction point",
    sliders = [Knob("β* [m]", BSTARS; init = 7)],
    panels  = [Panel(xlabel="distance from the IP [m]", ylabel="β [m]", title="β functions, log scale",
                     xlim=(-SPAN, SPAN), ylim=(0.1, 3000.0), yscale=:log10, height=280, basis="100%",
                     legend=:bottomright),
               Panel(xlabel="distance from the IP [m]", ylabel="x [mm]",
                     title="the two beams in the crossing plane (±3σ)",
                     xlim=(-NEAR, NEAR), ylim=(-36.0, 44.0), height=280, basis="100%", legend=:bottomright)],
    note    = "The peak β in the triplet grows as L*²/β*: the beam is largest exactly where the "*
              "quadrupoles are strongest, so the triplet needs a large aperture, has to be very "*
              "accurate, and dominates the chromaticity of the ring.",
) do βs
    c = SQUEEZE[βs]
    half = Lattice(AbstractElement[ir_half(c.k)..., arc_cell(c.kf, c.kd)..., arc_cell(c.kf, c.kd)...])
    t = transport_twiss(half, beam, optics4DUC(βs, 0.0, βs, 0.0); sample_integrator_steps=true, max_step=0.5)
    s = [-reverse(t.s); t.s]; bx = [reverse(t.betax); t.betax]; by = [reverse(t.betay); t.betay]
    sd = range(-30, 30, length=121)                  # the free-drift law, near the IP only

    # Bottom panel: the crossing-angle orbit x = θc √(β β*) sin ψ and the ±3σ envelope.
    θc = 5sqrt(EMIT/βs)                           # half angle for a 10σ separation in the drift
    near = t.s .<= NEAR
    sn = [-reverse(t.s[near]); t.s[near]]
    orb = θc .* sqrt.(t.betax[near] .* βs) .* sin.(t.mux[near])
    xo = 1e3 .* [-reverse(orb); orb]              # beam 1; beam 2 is its mirror image
    env = 1e3 .* 3 .* sqrt.(EMIT .* [reverse(t.betax[near]); t.betax[near]])
    glyphs = lattice_plot_data(Lattice(ir_half(c.k)))
    strip = lattice_strip([glyphs; mirror(glyphs)], (-NEAR, NEAR); panel=2, ylim=(-36.0, 44.0), strip=0.1)

    (ξarc, ξring) = chroma(c)
    iq = findall(t.s .<= LSTAR + 2LQ1 + LQ2 + 2.0)
    (series = [line(s, bx; panel=1, color=PALETTE[1], width=2, label="βₓ"),
               line(s, by; panel=1, color=PALETTE[2], width=2, label="βᵧ"),
               line(sd, βs .+ sd.^2 ./ βs; panel=1, color="#9aa4b2", dash=true, width=1.3, label="β* + s²/β*"),
               strip...,
               line(sn, xo; panel=2, color=PALETTE[1], width=2, label="beam 1"),
               line(sn, xo .+ env; panel=2, color=PALETTE[1], width=1, alpha=0.6),
               line(sn, xo .- env; panel=2, color=PALETTE[1], width=1, alpha=0.6),
               line(sn, -xo; panel=2, color=PALETTE[4], width=2, label="beam 2"),
               line(sn, -xo .+ env; panel=2, color=PALETTE[4], width=1, alpha=0.6),
               line(sn, -xo .- env; panel=2, color=PALETTE[4], width=1, alpha=0.6)],
     readouts = ["σ* at the IP" => @sprintf("%.0f μm", 1e6*sqrt(EMIT*βs)),
                 "β max (triplet)" => @sprintf("%.0f m", maximum(max.(t.betax[iq], t.betay[iq]))),
                 "L*²/β*" => @sprintf("%.0f m", LSTAR^2/βs),
                 "half crossing angle" => @sprintf("%.0f μrad", 1e6θc),
                 "chromaticity of the arcs" => @sprintf("(%.1f, %.1f)", ξarc...),
                 "chromaticity of the ring" => @sprintf("(%.1f, %.1f)", ξring...),
                 "luminosity ∝ 1/β*" => @sprintf("×%.1f vs β* = 1 m", 1/βs)])
end
```

Three things grow as $\beta^*$ shrinks, and together they set how far an IR can
be pushed:

- **Aperture.** The beam size in the triplet scales as
  $\sqrt{\beta_{\max}}\propto L^*/\sqrt{\beta^*}$. The crossing-angle orbit grows
  the same way, because the crossing angle itself must grow as
  $1/\sqrt{\beta^*}$ to keep the beams separated.
- **Chromaticity.** A quadrupole of strength $k$ and length $\ell$ contributes
  $-\beta k\ell/4\pi$ to the chromaticity. In the triplet β is huge, so the IR
  soon contributes more chromaticity than all the arc cells together, and the
  sextupoles have to correct it.
- **Field quality.** Magnet errors kick the beam in proportion to the beam size
  there. Errors in the triplet are amplified accordingly.

## Interaction region facts

- **Functions:** squeeze the beams to a small $\beta^*$, match the optics and
  the dispersion to the arcs, and leave room for the detector.
- **Detector:** the space around the IP is a drift of length $2L^*$. Over it,
  $\beta(s)=\beta^*+s^2/\beta^*$.
- **Final focus:** the triplet sits right outside the detector. Its quadrupoles
  are the strongest in the ring, and β is highest inside them.
- **Shared aperture:** both beams pass through the triplet, so it needs aperture
  for two beams separated by the crossing angle. Separation dipoles after the
  triplet guide the beams into their own vacuum chambers.
- **Matching:** matching quadrupoles and a dispersion suppressor connect the IR
  to the arcs.
- **Cost:** a rule of thumb is that ten metres of interaction region cost as
  much as a kilometre of regular arc.

The newest designs push all of this further. The Electron–Ion Collider crosses
its electron and hadron beams at 25 mrad, restores head-on collisions with crab
cavities, and integrates its first magnets into the detector. Belle II at
SuperKEKB uses $\beta_y^*=0.3$ mm with the nano-beam scheme of
[](luminosity.md).

```{figure} ../images/collider/detectors.jpg
:width: 100%
:name: fig-detectors
Detectors around interaction points. Top: PHENIX (left) and STAR (right) at RHIC. Bottom: ATLAS at the LHC (left) and a cut-away drawing of Belle II at SuperKEKB (right).
```
