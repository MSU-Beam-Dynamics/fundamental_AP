---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Dispersion and Chromaticity

## Dispersion

So far, we have considered the motion of particles at the reference momentum $P_0$. What will be the dynamical effect for an off-momentum particle?

Obviously, an off-momentum particle ($\delta=(P-P_0)/P_0 \ne 0$) bends differently in a dipole:

```{figure} ../images/dispersive_orbit_in_dipole.png
:width: 250px
:name: fig:dispersiveorbit
Dispersive orbit of an off-momentum particle inside a dipole.
```

Hill's equation gains a source term,

$$
x''+\left(\frac{1}{\rho^2}+k(s)\right)x=\frac{\delta}{\rho},
$$

and writing 

$$
x=x_\beta+ D \delta
$$

defines the **dispersion function**

$$
D(s)=\frac{\partial x(s)}{\partial \delta},
\qquad
D''+\left(\frac{1}{\rho^2}+k(s)\right)D=\frac{1}{\rho}.
$$

Dispersion obeys the same equation as $x$ but is driven only by dipoles —
quadrupoles modify it but do not create it. 

Below is an example of a simple FODO lattice with two dipoles per cell, with the twiss parameters and the dispersion function plotted along the cell. 

```{code-cell} julia
:tags: [hide-input]

using TrackPad, TrackPadWidgets

beam = Beam(3.0e9)

# One FODO cell serves every example on this page:
#     QF – SF – O – B – O – QD – SD – O – B – O
# 0.5 m quadrupoles, 0.2 m sextupoles, 2 m dipoles, 0.5 m of drift after each
# magnet group, so the period is 7 m. Each sextupole is carved out of the drift
# that follows its quadrupole, which leaves the dipole positions — and therefore
# every result below — exactly where they were before the sextupoles existed.
const KQ    = 0.9                     # m⁻², stable for every angle on the knob
const LQ_D  = 0.5
const LS_D  = 0.2
const LB_D  = 2.0
const LD_D  = 0.5
const CELL  = 2*(LQ_D + LB_D + 2*LD_D)

fodo_bend(θ; k = KQ, k2f = 0.0, k2d = 0.0) = AbstractElement[
    Quadrupole(LQ_D, +k; name=:QF), Sextupole(LS_D, k2f; name=:SF), Drift(LD_D - LS_D),
    SBend(LB_D, θ; name=:B1),       Drift(LD_D),
    Quadrupole(LQ_D, -k; name=:QD), Sextupole(LS_D, k2d; name=:SD), Drift(LD_D - LS_D),
    SBend(LB_D, θ; name=:B2),       Drift(LD_D)]

# Element spans, used both for the glyph band and for the α_C integral.
const QF_S, QD_S = (0.0, LQ_D), (LQ_D + LB_D + 2LD_D, LQ_D + LB_D + 2LD_D + LQ_D)
const B1_S = (LQ_D + LD_D, LQ_D + LD_D + LB_D)
const B2_S = (QD_S[2] + LD_D, QD_S[2] + LD_D + LB_D)
inbend(s) = (B1_S[1] < s < B1_S[2]) || (B2_S[1] < s < B2_S[2])

# The beamline band across the top of the β axis. Geometry does not depend on the
# bend angle, so it is built once from TrackPad's backend-independent glyph data
# and handed to the widget as statics.
const BMAX  = 15.0                    # left-axis range, leaving the band its own room
BEAMLINE = lattice_strip(lattice_plot_data(Lattice(fodo_bend(deg2rad(10.0)); periodic=true)),
                         CELL; ymax=BMAX)

ANGLES = collect(range(0.0, 20.0, length=9))     # bending angle per dipole [deg]

explorer(
    title   = "A FODO cell with two bends: dipoles create the dispersion",
    sliders = [Knob("bend angle per dipole [deg]", ANGLES;
                    fmt = a -> string(round(a; digits=1), "°"), init = 5)],
    panels  = [Panel(xlabel="s [m]", ylabel="β [m]", y2label="Dₓ [m]",
                     title="Twiss functions and dispersion along the cell",
                     ylim=(0.0, BMAX), y2lim=(0.0, 3.2), height=320,
                     legend=:bottomleft, basis="100%")],
    statics = BEAMLINE,
    note = "Set the angle to zero and Dₓ vanishes identically: quadrupoles shape the "*
           "dispersion but only dipoles create it. Opening the bend leaves βᵧ and Qᵧ "*
           "untouched — a sector bend has no vertical focusing — while the 1/ρ² term "*
           "adds horizontal focusing, so βₓ shrinks and Qₓ creeps up.",
) do adeg
    θ   = deg2rad(adeg)
    lat = Lattice(fodo_bend(θ); periodic=true)
    tw  = periodic_twiss(lat, beam;      sample_integrator_steps=true, max_step=0.05)
    dp  = periodic_dispersion(lat, beam; sample_integrator_steps=true, max_step=0.05)

    # α_C = (1/C)∮ D/ρ ds, and 1/ρ = θ/L_B is nonzero only inside the two dipoles.
    # dp.dx is d/dδ_E; beam.beta ≈ 1 here, so this is the momentum convention too.
    αC = 0.0
    for i in 1:length(dp.s)-1
        inbend((dp.s[i] + dp.s[i+1])/2) || continue
        αC += (dp.dx[i] + dp.dx[i+1])/2 * (θ/LB_D) * (dp.s[i+1] - dp.s[i])
    end
    αC /= CELL

    (series = [line(tw.s, tw.betax; panel=1, color=PALETTE[1], width=2.0, label="βₓ"),
               line(tw.s, tw.betay; panel=1, color=PALETTE[2], width=2.0, label="βᵧ"),
               line(dp.s, dp.dx;    panel=1, color=PALETTE[4], width=2.0, axis=:right,
                    label="Dₓ (right)")],
     readouts = ["bend angle"    => string(round(adeg; digits=1), "° per dipole"),
                 "ρ = L_B/θ"     => θ == 0 ? "∞" : string(round(LB_D/θ; sigdigits=4), " m"),
                 "Qₓ"            => round(tw.tunex; digits=4),
                 "Qᵧ"            => round(tw.tuney; digits=4),
                 "βₓ max"        => string(round(maximum(tw.betax); digits=3), " m"),
                 "βᵧ max"        => string(round(maximum(tw.betay); digits=3), " m"),
                 "Dₓ min / max"  => string(round(minimum(dp.dx); digits=3), " / ",
                                           round(maximum(dp.dx); digits=3), " m"),
                 "α_C"           => round(αC; sigdigits=4),
                 "cells for 2π"  => adeg == 0 ? "—" : string(round(360/(2adeg); digits=1))])
end
```

$D_x$ is exactly zero at zero bend angle.  Bending manget creates dispersion, while disperson function transport in other components like particles does, because disperson function satisfy the same Hill's equation outside bending magnets.






A direct consequence of dispersive effect is that the path length depends on energy:

$$
\frac{\Delta C}{C}=\frac{1}{C}\int\frac{D(s)\,\delta}{\rho}ds\equiv\alpha_C \delta,
$$

with $\alpha_C$ the **momentum compaction factor**. The revolution-time
dependence combines path-length and speed changes:

$$
\frac{\Delta T}{T}
=\left(\alpha_C-\frac{1}{\gamma_0^2}\right)\delta
\equiv \eta\,\delta,
$$



defining the **phase-slip factor** $\eta$.

That integral runs over the whole cell, but $1/\rho$ is zero everywhere except
inside the two dipoles, so only the dispersion *inside the bends* contributes.
Putting the two factors on one set of axes makes that plain, and the readouts
carry the integral itself.

```{code-cell} julia
:tags: [hide-input]

# Reuses the cell, the constants and the geometry of the dispersion example above.
INVRHO(θ) = θ/LB_D                       # 1/ρ inside a sector bend of length L_B

"1/ρ(s) as an exact step: zero outside the bends, θ/L_B inside."
function invrho_step(θ)
    k = INVRHO(θ)
    ([0.0, B1_S[1], B1_S[1], B1_S[2], B1_S[2], B2_S[1], B2_S[1], B2_S[2], B2_S[2], CELL],
     [0.0, 0.0,     k,       k,       0.0,     0.0,     k,       k,       0.0,     0.0])
end

explorer(
    title   = "Only the dispersion inside the bends contributes to α_C",  
    sliders = [Knob("bend angle per dipole [deg]", ANGLES;
                    fmt = a -> string(round(a; digits=1), "°"), init = 5)],
    panels  = [Panel(xlabel="s [m]", ylabel="Dₓ [m]", y2label="1/ρ [m⁻¹]",
                     title="the two factors of the integrand",
                     xlim=(0.0, CELL), ylim=(0.0, 3.2), y2lim=(0.0, 0.25), height=300,
                     legend=:bottomleft, basis="100%")],
    note = "Dₓ is nonzero all the way round — quadrupoles and drifts transport it — but "*
           "1/ρ switches it off outside the dipoles, so the product Dₓ/ρ is nonzero only "*
           "there. Dividing its integral by the 7 m cell length gives α_C.",
) do adeg
    θ   = deg2rad(adeg)
    lat = Lattice(fodo_bend(θ); periodic=true)
    dp  = periodic_dispersion(lat, beam; sample_integrator_steps=true, max_step=0.02)

    sr, kr = invrho_step(θ)

    # ∫Dₓ/ρ ds by trapezoid on the sampled grid, over the bends alone.
    area = 0.0
    for i in 1:length(dp.s)-1
        inbend((dp.s[i] + dp.s[i+1])/2) || continue
        area += (dp.dx[i] + dp.dx[i+1])/2 * INVRHO(θ) * (dp.s[i+1] - dp.s[i])
    end

    (series = [line(dp.s, dp.dx; panel=1, color=PALETTE[4], width=2.0, label="Dₓ (left)"),
               line(sr, kr;      panel=1, color=PALETTE[3], width=2.0, axis=:right,
                    label="1/ρ (right)")],
     readouts = ["bend angle"        => string(round(adeg; digits=1), "° per dipole"),
                 "1/ρ in the bends"  => string(round(INVRHO(θ); sigdigits=4), " m⁻¹"),
                 "Dₓ in bend 1"      => string(round(minimum(dp.dx[inbend.(dp.s)]); digits=3), " … ",
                                               round(maximum(dp.dx[inbend.(dp.s)]); digits=3), " m"),
                 "∫Dₓ/ρ ds"          => string(round(area; sigdigits=4), " m"),
                 "cell length C"     => string(CELL, " m"),
                 "α_C = area / C"    => round(area/CELL; sigdigits=4),
                 "1/γ_t² = α_C"      => adeg == 0 ? "—" : string("γ_t = ", round(1/sqrt(area/CELL); digits=3))])
end
```






## Chromaticity

Off-momentum particles also see a different quadrupole focal length
($k=G/(B\rho)$ shrinks with momentum). Consequently the phase advance varies
with energy — the linear ***chromaticity*** is defined as

$$
\xi=\frac{d\nu}{d\delta}.
$$

Before any of that, it is worth seeing the effect on single trajectories. Send
a nearly parallel bunch into one short focusing quadrupole and colour every
particle by its momentum offset: the stiffer particles (red) are bent less and
cross the axis late, the softer ones (blue) are bent more and cross early. The
focal point is smeared along $s$, and that smear *is* the chromaticity.

```{code-cell} julia
:tags: [hide-input]

using TrackPad, TrackPadWidgets, Random, Statistics

const NCH   = 245                # particles, 35 per colour band
const NBAND = 7                  # δ bands, blue (low momentum) → red (high)
const LQCH  = 0.1                # quadrupole length [m], centred on s = 0
const ARM   = 1.0 - LQCH/2       # drift on each side, so s runs −1 m → +1 m
# β = 20 m with ε = 32 nm·rad puts 0.8 mm and 0.04 mrad at the entrance: a bunch
# wide enough to see and divergent enough to ignore, i.e. essentially parallel.
const ENTCH = optics4DUC(20.0, 0.0, 20.0, 0.0)
const EMCH  = 32e-9

# Blue → grey → red, chosen to stay legible on both book themes (a bwr ramp's
# white middle disappears on the light background).
const BANDS = ["#2166ac", "#4393c3", "#92c5de", "#9aa4b2", "#f4a582", "#d6604d", "#b2182b"]

"The line as a TrackPad lattice, split so s = ±1 and the quadrupole centre are boundaries."
ch_lattice(k1) = Lattice(AbstractElement[
    Drift(ARM), Quadrupole(LQCH/2, k1; name=:Q), Quadrupole(LQCH/2, k1; name=:Q),
    Drift(ARM)]; name=:chromdemo)

"A Gaussian bunch matched to the entrance optics, with a chosen rms momentum spread."
ch_beam(σδ) = matched_gaussian(MersenneTwister(7), NCH, ENTCH;
                               emitx=EMCH, emity=EMCH,
                               longitudinal=[0.01^2 0.0; 0.0 σδ^2])

"Track through the lattice, keeping every coordinate at each element boundary."
function ch_track(k1, σδ)
    c = ch_beam(σδ); flags = zeros(Int, NCH)
    S = Float64[-1.0]; C = [copy(c)]
    for e in ch_lattice(k1)
        linepass!(c, Lattice(AbstractElement[e]), beam, flags)
        push!(S, S[end] + get_length(e)); push!(C, copy(c))
    end
    S, C
end

"Split the bunch into equal-population bands ordered by momentum."
band_index(c) = (p = sortperm(c[:, 6]); b = zeros(Int, length(p));
                 for (rank, i) in enumerate(p); b[i] = min(NBAND, 1 + (rank-1)*NBAND ÷ length(p)); end; b)

const CHYL   = (-4.5, 4.5)      # σ_δ is a Gaussian rms, so the tails reach past 3 mm
BEAMLINE_CH = lattice_strip(lattice_plot_data(ch_lattice(14.3)), (-1.0, 1.0);
                            ylim=CHYL, offset=-1.0)

SPREAD = [0.0, 0.05, 0.10, 0.15, 0.20]                # rms δ
FOCAL  = [0.5, 0.7, 1.0]                              # nominal focal length [m]

explorer(
    title   = "One quadrupole, one bunch: the focal point smears with momentum",
    sliders = [Knob("σ_δ", SPREAD; fmt = d -> string(round(Int, d*100), " %"), init = 5),
               Knob("focal length [m]", FOCAL;
                    fmt = f -> string(round(f; digits=2)), init = 2)],
    panels  = [Panel(xlabel="s [m]", ylabel="x [mm]",
                     title="trajectories, coloured by momentum offset",
                     xlim=(-1.0, 1.0), ylim=CHYL, height=320,
                     legend=:bottomleft, basis="100%"),
               Panel(title="(x, pₓ) at s = −1 m", autoscale=:frame, share="chr",
                     height=175, basis="31%", minwidth=150, ticklabels=false),
               Panel(title="(x, pₓ) at the quadrupole", autoscale=:frame, share="chr",
                     height=175, basis="31%", minwidth=150, ticklabels=false),
               Panel(title="(x, pₓ) at s = +1 m", autoscale=:frame, share="chr",
                     height=175, basis="31%", minwidth=150, ticklabels=false)],
    statics = BEAMLINE_CH,
    note = "The bunch is a Gaussian matched to β = 20 m at the entrance, so it arrives "*
           "0.8 mm wide and almost parallel. TrackPad gives every particle the same pₓ "*
           "kick — it is the drift that divides by (1+δ), so stiffer particles turn "*
           "through a smaller angle and focus further downstream. Set σ_δ to zero and "*
           "the fan collapses to one focal point; open it and the point becomes a "*
           "segment. The three insets share one x/pₓ box.",
) do σδ, f
    k1 = 1/(f*LQCH)
    S, C = ch_track(k1, σδ)
    bands = band_index(C[1])

    traj = map(1:NBAND) do b
        xs = Union{Float64,Nothing}[]; ys = Union{Float64,Nothing}[]
        for i in findall(==(b), bands)
            for j in eachindex(S); push!(xs, S[j]); push!(ys, C[j][i, 1]*1e3); end
            push!(xs, nothing); push!(ys, nothing)
        end
        line(xs, ys; panel=1, color=BANDS[b], width=0.9, alpha=0.55,
             label = b == 1 ? "δ < 0" : b == NBAND ? "δ > 0" : nothing)
    end

    insets = [points(round.(C[j][idx, 1] .* 1e3; sigdigits=4),
                     round.(C[j][idx, 2] .* 1e3; sigdigits=4);
                     panel=p, color=BANDS[b], size=2.6, alpha=0.55)
              for (p, j) in zip(2:4, (1, 3, 5))
              for (b, idx) in enumerate([findall(==(bb), bands) for bb in 1:NBAND])]

    (series = vcat(traj, insets),
     readouts = ["nominal f"        => string(round(f; digits=2), " m"),
                 "k₁"               => string(round(k1; digits=3), " m⁻²"),
                 "σ_δ"              => string(round(Int, σδ*100), " %"),
                 "σₓ at entrance"   => string(round(std(C[1][:,1])*1e3; sigdigits=3), " mm"),
                 "σₓ′ at entrance"  => string(round(std(C[1][:,2])*1e3; sigdigits=3), " mrad"),
                 "focus of δ = −σ_δ" => string(round(f*(1-σδ); sigdigits=3), " m"),
                 "focus of δ = +σ_δ" => string(round(f*(1+σδ); sigdigits=3), " m")])
end
```




Both dispersion and chromaticity are inevitable in lattices containing dipoles
and quadrupoles, and many beam-dynamics problems stem from them, so their
control is an essential part of accelerator design. Sextupoles correct the
chromaticity by producing another energy dependent focusing effect, but introduce nonlinearity, requiring attention to the dynamic aperture of accelerator.



```{code-cell} julia
:tags: [hide-input]

using TrackPad, TrackPadWidgets, StaticArrays

# The page's cell, with the bends fixed at 5° and k₁ on the knob. The sextupoles
# are present but unpowered, so this really is the *natural* chromaticity.
const ΘNAT = deg2rad(5.0)
fodo_k(k) = fodo_bend(ΘNAT; k = k)

δs = collect(range(-1e-2, 1e-2, length=13))

"Tune at each δ, tracked from the off-momentum closed orbit."
function tunecurve(lat)
    qx = Float64[]; qy = Float64[]
    for d in δs
        o = find_closed_orbit_4d(lat, beam; dp=d)
        t = gettune(lat, beam; reference=SVector{6,Float64}(o[1],o[2],o[3],o[4],0.0,d))
        push!(qx, t[1]); push!(qy, t[2])
    end
    qx, qy
end

"Least-squares slope dQ/dδ of a tune curve — the chromaticity the plot shows."
slope(q) = (δ̄ = sum(δs)/length(δs); q̄ = sum(q)/length(q);
            sum((δs .- δ̄) .* (q .- q̄)) / sum((δs .- δ̄).^2))

explorer(
    title   = "Natural chromaticity of the FODO cell (dipoles at 5°, no sextupoles)",
    sliders = [Knob("k₁ [m⁻²]", range(0.3, 1.1, length=17);
                    fmt = k -> string(round(k;digits=2)), init = 13)],
    # The tune moves by only a few parts in a thousand across the energy scan,
    # so the axis has to follow the frame or the slope — the whole point — is flat.
    panels  = [Panel(xlabel="δ_E", ylabel="tune", title="tune versus energy",
                     autoscale=:frame, height=300, legend=:bottomleft)],
    note = "ξ = −(1/4π)∮β(s)k(s)ds: stronger quadrupoles focus more but also chromatically "*
           "detune more, and the cell runs out of stability just above k₁ = 1.197 m⁻². Real rings "*
           "correct this to small positive values with sextupole families — the next "*
           "figure. The readouts check TrackPad's getchrom against the slope of the "*
           "curve actually plotted.",
) do k
    lat = Lattice(fodo_k(k); periodic=true)
    qx, qy = tunecurve(lat)
    ξ = getchrom(lat, beam; centered=true)
    (series = [line(δs, qx; color=PALETTE[1], label="Qₓ"),
               line(δs, qy; color=PALETTE[2], label="Qᵧ")],
     readouts = ["Qₓ(0)"            => round(qx[7]; digits=4),
                 "Qᵧ(0)"            => round(qy[7]; digits=4),
                 "ξₓ (getchrom)"    => round(ξ[1]; digits=4),
                 "ξᵧ (getchrom)"    => round(ξ[2]; digits=4),
                 "dQₓ/dδ from curve" => round(slope(qx); digits=4),
                 "dQᵧ/dδ from curve" => round(slope(qy); digits=4),
                 "ΔQₓ over the scan" => round(maximum(qx) - minimum(qx); sigdigits=3)])
end
```

### A sextupole is needed to fix it

A sextupole kicks as $\Delta p_x \propto -(x^2-y^2)$, which does not directly introduce energy dependent kicks.  However, when placed where the beam is **dispersive**, so that $x=x_\beta+D\delta$, and
the square produces a cross term $\propto D\,\delta\,x_\beta$ — a kick
proportional to the betatron amplitude *and* to the momentum offset. That is a
quadrupole whose strength tracks $\delta$, which is exactly the error the
chromatic quadrupole made. Launch the same bunch with an initial dispersion and
watch the fan of focal points close up:

```{code-cell} julia
:tags: [hide-input]

using TrackPad, TrackPadWidgets, Random, Statistics

const LSX  = 0.1                 # sextupole length [m]
const SGDL = 0.02                # rms momentum spread for this demo
const K1SX = 1/(0.7*LQCH)        # the same 0.7 m quadrupole as above

"Drift, sextupole, quadrupole, drift — s = 0 at the quadrupole centre."
sx_lattice(k2) = Lattice(AbstractElement[
    Drift(ARM - LSX), Sextupole(LSX, k2; name=:S),
    Quadrupole(LQCH/2, K1SX; name=:Q), Quadrupole(LQCH/2, K1SX; name=:Q),
    Drift(ARM)]; name=:sextdemo)

"The same matched bunch, but launched with dispersion D: x = x_β + Dδ."
sx_beam(D) = matched_gaussian(MersenneTwister(11), NCH, ENTCH;
                              emitx=EMCH, emity=EMCH,
                              longitudinal=[0.01^2 0.0; 0.0 SGDL^2],
                              dispersion=[D, 0.0, 0.0, 0.0])

function sx_track(k2, D)
    c = sx_beam(D); flags = zeros(Int, NCH)
    S = Float64[-1.0]; C = [copy(c)]
    for e in sx_lattice(k2)
        linepass!(c, Lattice(AbstractElement[e]), beam, flags)
        push!(S, S[end] + get_length(e)); push!(C, copy(c))
    end
    S, C
end

"""
Where each momentum band comes to a betatron waist, drifting on from `c`.

Each band still spans a little δ, and with dispersion that alone smears x by
several mm — far more than the 0.8 mm betatron beam — so the band's own linear
x(δ) and pₓ(δ) trend is removed first. Without that the metric would report a
difference between D = 0 and D = 0.5 at k₂ = 0, where the lattices are identical.
"""
function band_waists(c, bands)
    sg = range(0.0, 1.6, length=161)
    map(1:NBAND) do b
        idx = findall(==(b), bands)
        x, px, d = c[idx,1], c[idx,2], c[idx,6]
        dd = d .- mean(d); den = sum(dd.^2)
        if den > 0
            x  = x  .- (sum(dd .* (x  .- mean(x)))/den) .* dd
            px = px .- (sum(dd .* (px .- mean(px)))/den) .* dd
        end
        env = [std(x .+ s .* px ./ (1 .+ d)) for s in sg]
        sg[argmin(env)]
    end
end

# The dispersive beam is 30 mm across at D = 0.5 m and 4 mm at D = 0, so a fixed
# axis would either clip one or flatten the other into a thread. The panel scales
# to the frame instead, and the beamline band is rebuilt to match — it has to be,
# since lattice_strip places the glyphs relative to the axis range it is given.
const SXHEAD = 0.70          # fraction of the half-range the data may occupy

K2S = collect(range(0.0, 100.0, length=11))
DS  = [0.0, 0.25, 0.5]

explorer(
    title   = "A sextupole in a dispersive spot undoes the chromatic focusing error",
    sliders = [Knob("k₂ [m⁻³]", K2S; fmt = k -> string(round(Int, k)), init = 1),
               Knob("initial dispersion D [m]", DS;
                    fmt = d -> string(round(d; digits=2)), init = 3)],
    panels  = [Panel(xlabel="s [m]", ylabel="x [mm]",
                     title="trajectories, coloured by momentum offset",
                     xlim=(-1.0, 1.0), autoscale=:frame, height=300,
                     legend=:bottomleft, basis="100%"),
               Panel(xlabel="band mean δ", ylabel="waist position [m]",
                     title="where each momentum band focuses",
                     # fixed, unlike the panel above: comparing the slope across
                     # settings is the point, and the waists span 0.59-0.72 m
                     xlim=(-0.05, 0.05), ylim=(0.56, 0.76), height=260,
                     legend=:bottomleft, basis="100%")],
    statics = [line([-0.05, 0.05], [0.7, 0.7]; panel=2, color="#9aa4b2",
                    dash=true, width=1.2, label="on-momentum focus")],
    note = "With D = 0 the sextupole sits on the beam axis and cannot help: the lower "*
           "curve keeps its slope whatever k₂ does. Give the bunch dispersion and the "*
           "slope — which is the chromaticity — can be driven to zero, then past it into "*
           "over-correction. Only the product k₂·D matters, so halving D doubles the k₂ "*
           "needed. The waist is measured on each band's betatron residual, with its own "*
           "x(δ) trend removed, so the lower panel shows focusing and not the dispersion "*
           "itself. The sextupole is still nonlinear: it flattens the slope, not the "*
           "curvature.",
) do k2, D
    S, C = sx_track(k2, D)
    bands = band_index(C[1])
    δmean = [mean(C[1][findall(==(b), bands), 6]) for b in 1:NBAND]
    waist = band_waists(C[end-1], bands)          # state at the quadrupole exit

    traj = map(1:NBAND) do b
        xs = Union{Float64,Nothing}[]; ys = Union{Float64,Nothing}[]
        for i in findall(==(b), bands)
            for j in eachindex(S); push!(xs, S[j]); push!(ys, C[j][i, 1]*1e3); end
            push!(xs, nothing); push!(ys, nothing)
        end
        line(xs, ys; panel=1, color=BANDS[b], width=0.9, alpha=0.5,
             label = b == 1 ? "δ < 0" : b == NBAND ? "δ > 0" : nothing)
    end

    # slope of waist position against δ — the chromaticity of this little line
    δ̄, w̄ = mean(δmean), mean(waist)
    slope = sum((δmean .- δ̄) .* (waist .- w̄)) / sum((δmean .- δ̄).^2)

    # Zoom to this frame, then lay the beamline band in the space above the data.
    xall  = maximum(maximum(abs, C[j][:, 1]) for j in eachindex(C)) * 1e3
    ytop  = xall / SXHEAD
    strip = lattice_strip(lattice_plot_data(sx_lattice(k2)), (-1.0, 1.0);
                          ylim=(-ytop, ytop), offset=-1.0)

    (series = vcat(traj, strip,
        [line(δmean, waist; panel=2, color="#5b6472", width=1.4)],
        [points([δmean[b]], [waist[b]]; panel=2, color=BANDS[b], size=7.0)
         for b in 1:NBAND]),
     readouts = ["k₂"              => string(round(Int, k2), " m⁻³"),
                 "D"               => string(round(D; digits=2), " m"),
                 "k₂·D"            => string(round(k2*D; digits=1), " m⁻²"),
                 "σₓ at entrance"  => string(round(std(C[1][:,1])*1e3; sigdigits=3), " mm"),
                 "waist spread"    => string(round(maximum(waist)-minimum(waist); digits=3), " m"),
                 "d(waist)/dδ"     => string(round(slope; digits=2), " m per unit δ"),
                 "verdict"         => abs(slope) < 0.2 ? "corrected" :
                                      slope > 0 ? "under-corrected" : "over-corrected"])
end
```






## Correcting chromaticity with two sextupole families

A single sextupole moves both $\xi_x$ and $\xi_y$ at once, so canceling both
chromaticities needs (at least) two independently powered families — typically
one near the focusing quadrupoles, where $\beta_x$ and the dispersion are large,
and one near the defocusing quadrupoles, where $\beta_y$ is large instead.
Watch the two knobs pull the operating point around the $(\xi_x,\xi_y)$ plane:

```{code-cell} julia
:tags: [hide-input]

using TrackPad, TrackPadWidgets, StaticArrays, LinearAlgebra

# The page's cell once more — now with its two sextupole families powered.
chromcell(k2f, k2d) = Lattice(fodo_bend(ΘNAT; k2f, k2d); periodic=true)

δs = collect(range(-4e-3, 4e-3, length=13))
function tunecurve(lat)
    qx = Float64[]; qy = Float64[]
    for d in δs
        try
            o = find_closed_orbit_4d(lat, beam; dp=d)
            t = gettune(lat, beam; reference=SVector{6,Float64}(o[1],o[2],o[3],o[4],0.0,d))
            push!(qx, t[1]); push!(qy, t[2])
        catch
            push!(qx, NaN); push!(qy, NaN)
        end
    end
    qx, qy
end

chrom(k2f, k2d) = getchrom(chromcell(k2f, k2d), beam; centered=true, closed_orbit=true)

# Two families, two chromaticities: the response is a 2×2 matrix, so solve it
# rather than hunting by hand. ξ₀ is the natural chromaticity of the bare cell.
const ξ0 = collect(chrom(0.0, 0.0))
const RESP = let h = 5.0, a = chrom(h, 0.0), b = chrom(0.0, h)
    [ (a[1]-ξ0[1])/h  (b[1]-ξ0[1])/h ; (a[2]-ξ0[2])/h  (b[2]-ξ0[2])/h ]
end
const KSOL = RESP \ (-ξ0)                       # strengths that zero both planes

# Grids built around the solution, so the corrected point is always the middle
# knob position and stays reachable if the cell is ever retuned.
SF_GRID = collect(range(0.0, 2KSOL[1], length=11))
SD_GRID = collect(range(2KSOL[2], 0.0, length=11))

explorer(
    title   = "Two sextupole families steer the two chromaticities",
    sliders = [Knob("k₂ focusing family, SF [m⁻³]", SF_GRID;
                      fmt = v -> string(round(v; digits=1)), init = 1),
               Knob("k₂ defocusing family, SD [m⁻³]", SD_GRID;
                      fmt = v -> string(round(v; digits=1)), init = 11)],
    panels  = [Panel(xlabel="δ_E = (E−E₀)/(P₀c)", ylabel="tune", title="tune versus energy",
                     ylim=(0.265, 0.280), height=260, legend=:bottomleft),
               Panel(xlabel="ξₓ", ylabel="ξᵧ", title="chromaticity plane",
                     xlim=(-0.7, 0.7), ylim=(-0.7, 0.7), height=260, legend=:topleft)],
    statics = [line([0.0, 0.0], [-0.7, 0.7]; panel=2, color="#9aa4b2", width=1),
               line([-0.7, 0.7], [0.0, 0.0]; panel=2, color="#9aa4b2", width=1),
               points([ξ0[1]], [ξ0[2]]; panel=2, color="#9aa4b2", size=5.0,
                      label="no sextupoles")],
    note = "The sextupoles sit where the dispersion is large, so an off-momentum particle "*
           "passes off-centre and feels an extra quadrupole kick ∝ k₂Dδ. One family alone "*
           "cannot zero both chromaticities — the two rows of the response matrix are not "*
           "parallel, which is exactly why two families are needed. Note that the "*
           "on-momentum tunes never move: a sextupole is a nonlinear element and does not "*
           "touch the linear optics.",
) do k2f, k2d
    lat = chromcell(k2f, k2d)
    qx, qy = tunecurve(lat)
    ξ = try chrom(k2f, k2d) catch; (NaN, NaN) end
    (series = [line(δs, qx; panel=1, color=PALETTE[1], label="Qₓ"),
               line(δs, qy; panel=1, color=PALETTE[2], label="Qᵧ"),
               points([ξ[1]], [ξ[2]]; panel=2, color=PALETTE[4], size=6.5, label="(ξₓ, ξᵧ)")],
     readouts = ["ξₓ"          => round(ξ[1]; digits=4),
                 "ξᵧ"          => round(ξ[2]; digits=4),
                 "Qₓ(0)"       => round(qx[7]; digits=5),
                 "Qᵧ(0)"       => round(qy[7]; digits=5),
                 "natural ξ"   => string("(", round(ξ0[1]; digits=3), ", ",
                                         round(ξ0[2]; digits=3), ")"),
                 "both zero at" => string("SF = ", round(KSOL[1]; digits=2),
                                          ", SD = ", round(KSOL[2]; digits=2))])
end
```
