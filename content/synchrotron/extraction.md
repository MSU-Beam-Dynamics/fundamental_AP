---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Extraction

## Schemes

Extraction can be treated as injection run backwards, but nonlinear resonances
open up more possibilities:

- **fast extraction**, the reverse of injection:
  - single-turn extraction,
  - non-resonant multi-turn extraction;
- **resonant extraction**, which uses a nonlinear resonance:
  - slow extraction over thousands of turns,
  - low-loss multi-turn extraction.

The simulations below use the same model ring as the injection page: four
thin-lens FODO cells, observed at the centre of a focusing quadrupole where
$\alpha=0$, with a thin multipole there for the sextupole and octupole kicks.
Phase space is again normalized, $\bar x=x/\sqrt\beta$ and
$\bar x'=\sqrt\beta\,x'$, in units of the rms beam size $\sigma=\sqrt\epsilon$:

```{code-cell} julia
using StaticArrays, TrackPad, TrackPadWidgets, Random, Statistics

beam = Beam(2.0e9; mass=TrackPad.M_PROTON, charge=1.0)    # 2 GeV kinetic protons

const NCELL = 4          # FODO cells around the ring
const LHALF = 10.0       # quadrupole centre to quadrupole centre [m]: an 80 m ring
const EMIT  = 1.0e-6     # rms emittance of the circulating beam [m·rad] (1 mm·mrad)
const σ     = sqrt(EMIT) # the unit of normalized phase space [√m]

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

"Normalized coordinates (x̄, x̄′) of every particle, in units of σ."
normalized(c, β) = (c[:, 1] ./ (sqrt(β)*σ), c[:, 2] .* sqrt(β) ./ σ)

"A Gaussian beam matched to β (α = 0)."
matched_beam(rng, N, β) = matched_gaussian(rng, N; betax=β, betay=β, emitx=EMIT, emity=EMIT,
                                           longitudinal=[1.0e-6 0.0; 0.0 1.0e-24])

β = beta_start(synchrotron(4/3))
println("β at the observation point: ", round(β; digits=2), " m, so σ = ",
        round(1e3*sqrt(β*EMIT); digits=2), " mm")
```

## Fast single-turn extraction

Fast extraction removes the whole beam in one turn:

- slow closed-orbit bumpers (rising over milliseconds) bring the circulating
  beam close to the septum;
- a fast kicker, fired in a gap in the bunch train, deflects the beam;
- a betatron phase advance $\mu$ downstream, the deflected beam crosses the septum
  blade and the septum bends it into the extraction line.

```{figure} ../images/synchrotron/single_turn_extraction.png
:width: 95%
:name: fig-single-turn-extraction
Fast single-turn extraction: closed-orbit bumpers move the circulating beam close to the septum, the kicker deflects it into the septum, and a defocusing quadrupole in between enhances the kick. After M. Fraser, CERN Accelerator School, 2017.
```

The displacement a kick $\theta$ produces at the septum follows from the same
transfer matrix as for injection,

$$
\Delta x_{\mathrm{septum}}=\theta\,\sqrt{\beta_{\mathrm{kicker}}\beta_{\mathrm{septum}}}\,\sin\mu ,
$$

so the most efficient layout again has $\mu\approx\pi/2$ between kicker and
septum. A defocusing quadrupole in between, where the kicked beam is far off
axis, adds to the deflection: a "kick-enhancing" quadrupole. With
$\beta\approx 40$ m at both ends, a 1 mrad kick moves the beam by 40 mm at the
septum.

## Non-resonant multi-turn extraction

Some machines downstream want the beam spread over a few turns rather than
delivered in one.

```{figure} ../images/synchrotron/multiturn_extraction.png
:width: 80%
:name: fig-multiturn-extraction
Non-resonant multi-turn extraction: closed-orbit bumpers push the circulating beam across the septum, which peels part of it off into the extraction channel.
```

CERN's PS fed the SPS this way for decades with its
*continuous transfer*:

- slow bumpers move the beam so that the septum blade cuts into it;
- the horizontal tune is set to $Q_h=0.25$, so the beam turns by 90° in
  normalized phase space every turn and a different slice faces the blade each
  turn;
- after four slices have been shaved off, the remaining core is kicked out on
  the fifth turn.

The widget below colours each particle of the circulating beam by the turn on
which it leaves:

```{code-cell} julia
:tags: [hide-input]

ring_ct = synchrotron(1.25)
β_ct = beta_start(ring_ct)
const NCT = 3000
const BLADE = 0.05      # septum blade thickness [σ]: particles landing in it are lost
SLICE_COLORS = (PALETTE[1], PALETTE[2], PALETTE[3], PALETTE[4], PALETTE[5])

"Turn on which each particle leaves (1–4 at the septum, 5 = kicked core, 0 = lost on the blade)."
function slice_beam(xs)
    c = matched_beam(MersenneTwister(4), NCT, β_ct)
    X0, P0 = normalized(c, β_ct)
    turn = zeros(Int, NCT); alive = trues(NCT)
    for n in 1:4
        X, _ = normalized(c, β_ct)
        on_blade = alive .& (xs .< X .< xs + BLADE)
        out = alive .& (X .>= xs + BLADE)
        turn[out] .= n; alive .&= .!(out .| on_blade)
        ringpass!(c, ring_ct, beam, zeros(Int, NCT), 1)
    end
    turn[alive] .= 5
    X0, P0, turn
end

explorer(
    title   = "Five-turn continuous transfer at Q = 0.25",
    sliders = [Knob("septum blade position [σ]", [0.4, 0.6, 0.8, 1.0, 1.2]; fmt = string, init = 2)],
    panels  = [Panel(xlabel="x̄ [σ]", ylabel="x̄′ [σ]", title="the beam before extraction, by exit turn",
                     xlim=(-4.0, 4.0), ylim=(-4.0, 4.0), equal=true, height=340, legend=:right,
                     basis="60%", minwidth=300)],
    note = "Slices 1–4 are shaved off by the septum on successive turns while the beam rotates "*
           "by 90° per turn; the core (5) is kicked out last. No septum position makes all "*
           "five slices equal, and the slices also differ in shape, so the extracted turns "*
           "have different emittances.",
) do xs
    X0, P0, turn = slice_beam(xs)
    series = Any[line([xs, xs], [-4.0, 4.0]; color="#5b6472", width=3.0, label="septum blade")]
    for n in 1:5
        sel = turn .== n
        push!(series, points(X0[sel], P0[sel]; color=SLICE_COLORS[n], size=1.8, alpha=0.75,
                             label = n < 5 ? "turn $n" : "turn 5 (kicked)"))
    end
    pct(n) = string(round(100*count(==(n), turn)/NCT; digits=1), " %")
    (series = series,
     readouts = vcat(["slice $n" => pct(n) for n in 1:5], ["lost on the blade" => pct(0)]))
end
```

This scheme extracts the beam within a few turns, but it has drawbacks:

- particles hitting the septum blade are lost, so the blade has to be very thin;
- the slices do not share the same intensity, emittance or centroid.

A nonlinear resonance can do the job better.

## Third-integer resonant slow extraction

Fixed-target experiments and medical therapy want a slow, steady beam spread
over thousands of turns (a "spill"), not a few-turn pulse. Resonant slow
extraction does this with a sextupole and a tune close to a third integer,
$Q\approx 1/3$ or $2/3$. The sextupole kick $\Delta\bar x'=-S\bar x^2$,
with normalized strength $S=\tfrac12\beta^{3/2}k_2L$, drives the third-order
resonance. Near it, phase space splits into a **stable triangle**, with three
unstable fixed points at its corners, and three separatrix arms leading out of
those corners (Kobayashi's theory). The triangle's inscribed radius shrinks in
proportion to the distance from the resonance, $\delta Q=Q-Q_r$:

$$
h=\frac{4\pi\,\delta Q}{S},
\qquad
A_{\mathrm{stable}}=3\sqrt3\,h^2=\frac{48\sqrt3\,\pi^2}{S^2}\,(\delta Q)^2 .
$$

Particles outside the triangle stream out along the arms, gaining amplitude
every third turn, until they jump over an electrostatic septum wire. Ramping the
tune slowly toward the resonance shrinks the triangle, which peels the beam off
from the outside in and sets the extraction rate. Particles that land on the
wire itself are lost; that fraction is roughly the wire thickness divided by
the step a particle makes in three turns.

```{code-cell} julia
:tags: [hide-input]
:label: fig-slow-extraction-portrait
:caption: Phase portrait near the third-integer resonance, at $Q=1+1/3+0.0017$ with a sextupole that makes the triangle's inscribed radius $h=2.5\sigma$, as at the start of the simulation below. Particles inside the triangle (blue) circulate on closed curves. Particles started just outside it (orange) hop from arm to arm every turn, a little farther out each time they return, and leave along the arms; one arm crosses the septum wire.

using CairoMakie

let δq = 0.335 - 1/3, S = 4π*δq/2.5            # S: normalized sextupole strength, h = 4πδq/S
    r3 = synchrotron(1 + 1/3 + δq; b2 = -S/(σ*β^1.5))
    β3 = beta_start(r3)
    x0 = [0.5, 1.0, 1.5, 2.0, 2.4, 2.8, -2.55, -2.65, -2.8, -3.0]   # start at x̄ = x0 σ, x̄′ = 0
    c  = zeros(length(x0), 6); c[:, 1] .= x0 .* (sqrt(β3)*σ)
    X = [Float64[] for _ in x0]; P = [Float64[] for _ in x0]
    for _ in 1:900
        ringpass!(c, r3, beam, zeros(Int, length(x0)), 1)
        x, p = normalized(c, β3)
        for i in eachindex(x0)
            (isfinite(x[i]) && abs(x[i]) < 12 && abs(p[i]) < 12) || continue
            push!(X[i], x[i]); push!(P[i], p[i])
        end
    end
    fig = Figure(size=(540, 500), fontsize=15)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="x̄ [σ]", ylabel="x̄′ [σ]", limits=(-10, 10, -10, 10))
    for i in eachindex(x0)
        bound = length(X[i]) == 900 && maximum(abs, X[i]) < 6
        scatter!(ax, X[i], P[i]; color=(bound ? PALETTE[1] : PALETTE[2], 0.6), markersize=2.5, strokewidth=0)
    end
    vlines!(ax, [8.0]; color=:gray30, linewidth=3)
    text!(ax, 7.7, -9.4; text="septum wire", align=(:right, :bottom), fontsize=13)
    fig
end
```

The simulation below tracks this process: the tune starts at $1/3+0.0017$,
where the triangle has a 2.5σ radius, and is ramped onto the resonance. Compare a **linear** ramp with one **shaped** so that
the fraction of the beam outside the triangle grows linearly in time:

```{code-cell} julia
:tags: [hide-input]

const NSE      = 2500
const RAMP     = 7000           # turns to reach the resonance
const SE_TURNS = 8000           # then 1000 turns on resonance
const SE_SNAP  = 0:400:SE_TURNS
const δQ0      = 0.335 - 1/3    # initial distance from the resonance
const S_SE     = 4π*δQ0/2.5     # normalized sextupole strength [1/σ]: h = 2.5σ at the start
const X_SEP    = 8.0            # electrostatic septum wire [σ]
const WIRE     = 0.02           # wire thickness [σ]

# h = 4πδQ/S is the inscribed radius; a circle of equal area has radius 1.286 h.
r_equal_area(δq) = 1.286*4π*δq/S_SE
function δQ_shaped(t)
    u = clamp(t/RAMP, 0.0, 1.0)
    F0 = exp(-r_equal_area(δQ0)^2/2)         # Gaussian fraction outside at the start
    F  = F0 + (1 - F0)*u                     # grows linearly to 1
    sqrt(-2log(F))*S_SE/(1.286*4π)
end
δQ_linear(t) = δQ0*clamp(1 - t/RAMP, 0.0, 1.0)
SE_RAMPS = [("linear tune ramp", δQ_linear), ("shaped tune ramp", δQ_shaped)]

"Track the slow extraction for one tune ramp, checking the septum every turn."
function slow_extraction(δQ)
    c = matched_beam(MersenneTwister(3), NSE, β)
    b2 = -S_SE/(σ*β^1.5)                   # gives Δx̄′ = −S x̄² in units of σ
    extracted = zeros(Int, SE_TURNS); lost = zeros(Int, SE_TURNS)
    snaps = Dict(0 => normalized(c, β))
    for t0 in 0:20:SE_TURNS-1               # the tune is stepped every 20 turns
        packed = GPULattice(synchrotron(1 + 1/3 + δQ(t0); b2), beam; dtype=Float64)
        for t in t0+1:t0+20
            batch_ringpass!(c, packed, 1)
            X, _ = normalized(c, β)
            out = .!isfinite.(X) .| (X .> X_SEP)
            on_wire = out .& (X .< X_SEP + WIRE)
            lost[t] = count(on_wire); extracted[t] = count(out) - lost[t]
            c = c[.!out, :]
            t in SE_SNAP && (snaps[t] = normalized(c, β))
        end
    end
    (snaps = snaps, extracted = extracted, lost = lost)
end
se_runs = [slow_extraction(f) for (_, f) in SE_RAMPS]

spill(v) = [100*sum(v[k:k+99])/NSE for k in 1:100:length(v)-99]   # % of the beam per 100 turns
SPILL_T = collect(50:100:SE_TURNS)

explorer(
    title   = "Third-integer slow extraction",
    sliders = [Knob("turn", collect(SE_SNAP); fmt = string, init = 1),
               Knob("tune ramp", 1:2; fmt = i -> SE_RAMPS[i][1], init = 1)],
    panels  = [Panel(xlabel="x̄ [σ]", ylabel="x̄′ [σ]", title="circulating beam",
                     xlim=(-12.0, 12.0), ylim=(-12.0, 12.0), equal=true, height=320,
                     legend=:bottomleft, basis="45%", minwidth=260),
               Panel(xlabel="turn", ylabel="spill [% of the beam per 100 turns]", title="spill",
                     xlim=(0.0, Float64(SE_TURNS)), ylim=(0.0, 3.5),
                     y2label="(Q − 1/3) × 10³", y2lim=(0.0, 1.8), height=320, legend=:topright,
                     basis="50%", minwidth=260)],
    statics = [line([X_SEP, X_SEP], [-12.0, 12.0]; panel=1, color="#5b6472", width=2.0,
                    label="septum wire")],
    note = "The stable triangle shrinks as the tune approaches 1/3; particles outside it leave "*
           "along the separatrix arms and jump the septum wire. A linear ramp extracts the dense "*
           "core late, so the spill grows; slowing the ramp near the resonance flattens it.",
) do turn, iramp
    run = se_runs[iramp]
    δQ = SE_RAMPS[iramp][2]
    X, P = run.snaps[turn]
    nbin = turn ÷ 100
    tq = collect(0:100:SE_TURNS)
    series = Any[points(X, P; panel=1, color=PALETTE[1], size=1.6, alpha=0.6, label="particles"),
                 line(tq, 1e3 .* δQ.(tq); panel=2, color=PALETTE[2], dash=true, axis=:right,
                      label="tune − 1/3"),
                 points([turn], [1e3*δQ(turn)]; panel=2, color=PALETTE[2], size=5.0, axis=:right)]
    nbin > 0 && push!(series, line(SPILL_T[1:nbin], spill(run.extracted)[1:nbin]; panel=2,
                                   color=PALETTE[1], width=1.6, label="extracted"))
    (series = series,
     readouts = ["tune Q" => round(1 + 1/3 + δQ(turn); digits=5),
                 "triangle radius h = 4πδQ/S" => string(round(4π*δQ(turn)/S_SE; digits=2) + 0.0, " σ"),
                 "extracted" => string(round(100*sum(run.extracted[1:turn])/NSE; digits=1), " %"),
                 "lost on the wire" => string(round(100*sum(run.lost[1:turn])/NSE; digits=1), " %"),
                 "circulating" => string(round(100*length(X)/NSE; digits=1), " %")])
end
```

## Resonant low-loss multi-turn extraction

Resonant extraction still loses beam on the septum, because particles keep
arriving at the wire everywhere along the arm. The low-loss alternative first
**splits the beam in phase space** and only then extracts it:

- with a sextupole and an octupole, the fourth-order resonance $Q=1/4$ develops
  four stable **islands**; the octupole makes the tune fall with amplitude, so
  once the bare tune is above $1/4$, some amplitude is exactly on resonance;
- the tune is ramped slowly across $1/4$: the islands appear at the centre, move
  outward, and carry particles with them, leaving four equally populated
  beamlets around a core;
- because $Q\approx1/4$, the islands swap places every turn. A bump and a kicker
  then push them across a septum placed in the *empty* space between the
  beamlets, one island per turn, with almost no particles hitting the blade.

This is the principle of CERN's Multi-Turn Extraction (MTE), which replaced the
continuous transfer between the PS and the SPS. The simulation below
ramps the tune from 0.249 to 0.258 over 10 000 turns:

```{code-cell} julia
:tags: [hide-input]

const NIS       = 2000
const IS_TURNS  = 11000
const IS_SNAP   = 0:500:IS_TURNS
const S_IS      = 0.075         # sextupole: Δx̄′ = +S x̄² (units of σ)
const O_IS      = 2.5e-4        # octupole:  Δx̄′ = +O x̄³
const X_SEP_IS  = 8.0           # where a septum would go, between the core and the outer island
q_island(t) = 0.249 + 1.0e-5*min(t, 10000)/11
β_is = beta_start(synchrotron(1.25))             # β near the quarter-integer tune

island_run = let
    c = matched_beam(MersenneTwister(5), NIS, β_is)
    b2 = -S_IS/(σ*β_is^1.5); b3 = -O_IS/(σ^2*β_is^2)   # the thin multipole kicks Δpₓ = −b2 x² − b3 x³
    snaps = Dict(0 => normalized(c, β_is))
    for t0 in 0:11:IS_TURNS-1                    # the tune is stepped every 11 turns
        packed = GPULattice(synchrotron(1 + q_island(t0); b2, b3), beam; dtype=Float64)
        batch_ringpass!(c, packed, 11)
        for t in t0+1:t0+11
            t in IS_SNAP && (snaps[t] = normalized(c, β_is))
        end
    end
    snaps
end

explorer(
    title   = "Capturing the beam in fourth-order islands",
    sliders = [Knob("turn", collect(IS_SNAP); fmt = string, init = 1)],
    panels  = [Panel(xlabel="x̄ [σ]", ylabel="x̄′ [σ]", title="normalized phase space",
                     xlim=(-15.0, 15.0), ylim=(-15.0, 15.0), equal=true, height=360,
                     legend=:bottomleft, basis="60%", minwidth=300)],
    statics = [line([X_SEP_IS, X_SEP_IS], [-15.0, 15.0]; color="#5b6472", width=2.0, dash=true,
                    label="where the septum goes")],
    note = "As the tune crosses 1/4 the islands emerge from the centre and drift outward, each "*
           "carrying a beamlet. Snapshots are 500 turns apart, so the four islands appear in "*
           "fixed places even though they trade places every turn. The septum then sits where "*
           "there is no beam.",
) do turn
    X, P = island_run[turn]
    r = hypot.(X, P)
    core = count(<(3.0), r)
    (series = [points(X, P; color=PALETTE[1], size=1.6, alpha=0.6, label="particles")],
     readouts = ["tune Q" => round(1 + q_island(turn); digits=5),
                 "core (r < 3σ)" => string(round(100*core/NIS; digits=1), " %"),
                 "outside the core" => string(round(100*(NIS - core)/NIS; digits=1), " %"),
                 "per island, if shared equally" => string(round(25*(NIS - core)/NIS; digits=1), " %")])
end
```

## Further reading

The injection and extraction material in these notes largely follows the CERN
Accelerator School course *Beam Injection, Extraction and Transfer* (Erice,
2017): <https://cas.web.cern.ch/schools/erice-2017>.
