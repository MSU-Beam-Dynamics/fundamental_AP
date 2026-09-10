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




:::{note}
TrackPad's sixth coordinate is $\delta_E=(E-E_0)/(P_0c)$, not
$\delta_P=(P-P_0)/P_0$; they coincide only to first order at ultrarelativistic
energies. `getchrom`'s legacy keyword names (`dp`, `dpp`) also refer to
$\delta_E$.
:::

## Chromaticity

Off-momentum particles also see a different quadrupole focal length
($k=G/(B\rho)$ shrinks with momentum). Consequently the phase advance varies
with energy — the linear ***chromaticity*** is defined as

$$
\xi=\frac{d\nu}{d\delta}.
$$

Both dispersion and chromaticity are inevitable in lattices containing dipoles
and quadrupoles, and many beam-dynamics problems stem from them, so their
control is an essential part of accelerator design. Sextupoles correct the
chromaticity but introduce nonlinearity, requiring attention to the dynamic
aperture.

TrackPad computes the chromaticity by finite differences of the tune versus
energy. Reusing the same 7 m cell as above — 0.5 m quadrupoles, 2 m dipoles,
0.5 m drifts — with each dipole fixed at $5°$, the knob is now the quadrupole
strength. It is the *natural* chromaticity in the sense that nothing corrects
it: the cell has no sextupoles. Both planes come out negative at every setting,
which is the generic result — a particle with $\delta>0$ is stiffer, the
quadrupoles focus it less, and its tune falls.

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
