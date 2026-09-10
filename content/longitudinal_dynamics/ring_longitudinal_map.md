---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Longitudinal Dynamics in a Ring

## Harmonic number and the turn-by-turn map

In a ring, RF cavities are installed at isolated locations. To synchronize
with the beam's arrival each revolution, the RF frequency must be an integer
multiple of the revolution frequency:

$$
f_{\mathrm{RF}}=hf_0=\frac{h}{T_0},
$$

where $h$ is the **harmonic number**.

```{figure} ../images/ring_harmonic.png
:width: 300px
:name: fig:harmonic
RF buckets around the ring for harmonic number h.
```

Describing the dynamics turn by turn with subscript $n$ for the value before
entering the cavity on the $(n+1)^{\text{th}}$ pass:

$$
\delta_{n+1}=
\delta_{n} + \frac{eV}{\beta^2 E_0} \left(\sin \phi_n-\sin \phi_s\right),
\qquad
\phi_{n+1}= \phi_{n} + 2\pi h \eta  \,\delta_{n+1},
$$

or, with $\Delta E$ instead of $\delta$ (convenient when the reference energy
ramps, as in a synchrotron):

$$
\Delta E_{n+1}
=\Delta E_{n} + eV \left(\sin \phi_n-\sin \phi_s\right),
\qquad
\phi_{n+1} = \phi_{n} + \frac{2\pi h \eta}{\beta^2 E_0}\, \Delta E_{n+1}.
$$

Note that the phase update uses the *new* energy — this makes the Jacobian of
the map exactly 1, so the phase-space area in $(\delta,\phi)$ or $(\Delta
E,\phi)$ is conserved while parameters stay constant.

## Fixed points and phase stability

The map has fixed points where $(\delta,\phi)\rightarrow(\delta,\phi)$:
$(0,\phi_s)$ and $(0,\pi-\phi_s)$. One is stable, the other unstable,
depending on the sign of the slip factor $\eta$. Consider a proton synchrotron
based on RHIC parameters:

| Parameter | Value |
| :--- | ------: |
| Injection energy | 20 GeV (kinetic) |
| Top energy | 250 GeV |
| Mass | 0.938 GeV |
| RF voltage | 5 MV |
| RF frequency | 28.15 MHz |
| Harmonic | 360 |
| $\alpha_c$ | 0.0018 |

Assume $\sin\phi_s=0.5$, i.e. $\phi_s=\pi/6$ or $\phi_s=5\pi/6$. Track two
families of particles at both energies:

```{code-cell} julia
using CairoMakie

"""
    longitudinal_map(turns, φ₀, ΔE₀; ...)

Turn-by-turn longitudinal map in (φ, ΔE) with slip factor held constant.
Returns (Φ_hist, ΔE_hist), each (turns+1) × nparticles.
"""
function longitudinal_map(turns, φ₀, ΔE₀; sin_φs=0.5, αc=0.0018,
                          E0_tot=250e9, mass=0.938e9, V=5e6, h=360)
    γ0  = E0_tot/mass
    β²  = 1 - 1/γ0^2
    η   = αc - 1/γ0^2                       # slip factor
    φ, ΔE = copy(φ₀), copy(ΔE₀)
    Φh, Eh = [copy(φ)], [copy(ΔE)]
    for _ in 1:turns
        ΔE .+= V .* (sin.(φ) .- sin_φs)
        φ  .+= (2π*h*η/(β²*E0_tot)) .* ΔE   # uses UPDATED energy ⇒ area preserving
        push!(Φh, copy(φ)); push!(Eh, copy(ΔE))
    end
    return reduce(hcat, Φh)', reduce(hcat, Eh)'
end

mass, αc, V, h = 0.938e9, 0.0018, 5e6, 360
sin_φs = 0.5; φs = asin(sin_φs)
turns = 2100; npar = 5

fig = Figure(size=(900, 380))
for (panel, (label, E0_tot)) in enumerate([("250 GeV (above transition)", 250e9), ("20 GeV (below transition)", 20e9)])
    γ = E0_tot/mass; β² = 1-1/γ^2; η = αc-1/γ^2
    # bucket half-height estimate for axis scaling
    yf = √(cos(φs) - (π-2φs)*sin(φs)/2)
    height = 2*√(V/(2π*β²*E0_tot*h*abs(η)))*yf

    # family near the UNSTABLE fixed point π−φs (drifts out or gets extracted)
    Φu = fill(π-φs, npar); Eu = collect(range(height/npar, 0.99height; length=npar))
    Hu = longitudinal_map(turns, Φu, Eu; sin_φs=sin_φs, αc=αc, E0_tot=E0_tot, mass=mass, V=V, h=h)
    # family near the STABLE fixed point φs
    Φs = fill(φs, npar); Es = collect(range(0, 0.99height; length=npar))
    Hs = longitudinal_map(turns, Φs, Es; sin_φs=sin_φs, αc=αc, E0_tot=E0_tot, mass=mass, V=V, h=h)

    ax = Axis(fig[1, panel]; limits=((-π, 2π), (-2height, 2height)),
              xlabel="Phase [rad]", title=label)
    for i in 1:npar
        lines!(ax, Hu[1][:, i], Hu[2][:, i] ./ β²; color=(:darkorange, 0.55))
        lines!(ax, Hs[1][:, i], Hs[2][:, i] ./ β²; color=(:steelblue, 0.55))
    end
end
fig
```

:::{note}
`δ = ΔE/(β²E₀)` converts between the two vertical axes above; near
transition the difference matters, which is why $(\phi,\Delta E)$ is preferred
during acceleration.
:::

Summarizing the observed stability pattern:

| | $\phi_s < \pi/2$ | $\phi_s > \pi/2$ | Phase space rotation | |
|:--|:--:|:--:|:---|:---|
| High energy ($\eta>0$) | Unstable | Stable | Clockwise | **Above transition** |
| Low energy ($\eta<0$) | Stable | Unstable | Counter-clockwise | **Below transition** |

## Acceleration and bucket distortion

When the beam is actually accelerated, the phase-space area of $(\phi,\delta)$
is no longer conserved because $\beta^2E_0$ changes; the area of $(\phi,\Delta
E)$ stays nearly conserved provided $2\pi h\eta/(\beta^2E_0)\approx$
constant. At high acceleration rates even that fails — relevant e.g. to medical
therapy synchrotrons:

```{code-cell} julia
# fast acceleration: (φ, ΔE) orbits spiral outward (medical-synchrotron-like)
mass, αc, V, h = 0.938e9, 0.05, 0.1e6, 1
E0_tot = 45e6 + mass                      # low-energy proton
γ = E0_tot/mass; β² = 1-1/γ^2; η = αc-1/γ^2
sin_φs = 0.5; φs = asin(sin_φs)
yf = √(cos(φs) - (π-2φs)*sin(φs)/2)
height = 2*√(V/(2π*β²*E0_tot*h*abs(η)))*yf

npar = 4
Φ0 = fill(π-φs, npar)
E0 = range(0, 0.3e6; length=npar) .|> Float64
H = longitudinal_map(400, Φ0, E0; sin_φs=sin_φs, αc=αc, E0_tot=E0_tot, mass=mass, V=V, h=h)

fig = Figure()
ax = Axis(fig[1, 1]; xlabel="Phase [rad]", ylabel="ΔE [MeV]",
          title="Low-energy ring: distorted (non-closed) buckets")
for i in axes(H[1], 2)
    lines!(ax, H[1][:, i], H[2][:, i] ./ 1e6)
end
fig
```

## From map to differential equations

Over many turns the discrete changes can be approximated as continuous:

$$
\dot\delta=\frac{eV\omega_0}{2\pi\beta^2E_0}\left(\sin\phi-\sin\phi_s\right),
\qquad
\dot\phi=h\omega_0\eta\,\delta ,
$$

which combine into the pendulum-like second-order equation

$$
\ddot\phi=\frac{eVh\eta\omega_0^2}{2\pi\beta^2E_0}\left(\sin\phi-\sin\phi_s\right).
$$

```{figure} ../images/phase_stability_region.png
:width: 500px
:name: fig:phasestability
Stable phase regions versus slip factor sign.
```

## Small-amplitude approximation and synchrotron tune

Linearizing around $\Delta\phi=\phi-\phi_s$,

$$
\ddot{\Delta\phi}=\frac{eVh\eta\omega_0^2}{2\pi\beta^2E_0}\cos\phi_s\;\Delta\phi,
$$

so stable oscillation requires

$$
\eta\cos\phi_s < 0,
$$

with small-amplitude synchrotron tune

$$
Q_s=\sqrt{\frac{-eVh\eta\cos\phi_s}{2\pi\beta^2E_0}}
\equiv\nu_s\sqrt{\left|\cos\phi_s\right|}.
$$

Let us verify $Q_s$ against the simulated map via spectral analysis — and see
how it changes away from the small-amplitude limit the formula assumes:

```{code-cell} julia
:tags: [hide-input]

using TrackPadWidgets

# RHIC-like top energy: same longitudinal_map used throughout this chapter
mass, αc, V, h = 0.938e9, 0.0018, 5e6, 360
E0_tot = 250e9
γ = E0_tot/mass; β² = 1-1/γ^2; η = αc-1/γ^2
sin_φs = 0.5
φs_stable = π - asin(sin_φs)      # above transition (η>0): cos φs < 0 required!
νs = √(V*h*abs(η)/(2π*β²*E0_tot))
Qs_small = νs * √abs(cos(φs_stable))
turns = 4096

explorer(
    title   = "Synchrotron tune softens with amplitude",
    sliders = [Knob("amplitude Δφ₀ [rad]", range(0.05, 1.0, length=20);
                    fmt = a -> string(round(a; digits=2)), init = 1)],
    panels  = [Panel(xlabel="phase φ [rad]", ylabel="ΔE [MeV]", title="tracked orbit",
                     height=280, legend=:bottomleft)],
    note = "Qs measured from zero-crossings of the tracked orbit versus the small-amplitude "*
           "formula Qs = νs√|cos φs|: real buckets run slower far from the centre, and Qs → 0 "*
           "at the separatrix.",
) do amp
    H = longitudinal_map(turns, [φs_stable+amp], [0.0];
                         sin_φs=sin_φs, αc=αc, E0_tot=E0_tot, mass=mass, V=V, h=h)
    signal = H[1][:,1] .- φs_stable
    ncross = sum(sign.(signal[1:end-1]) .* sign.(signal[2:end]) .< 0)
    Qs_meas = ncross/(2*length(signal))
    nshow = 2400
    (series = [line(H[1][1:nshow,1], H[2][1:nshow,1]./1e6; color=PALETTE[1], label="(φ, ΔE)"),
               points([φs_stable], [0.0]; color=PALETTE[3], size=5.0, label="stable fixed point")],
     readouts = ["Δφ₀" => round(amp; digits=3),
                 "Qs measured" => round(Qs_meas; sigdigits=4),
                 "Qs small-amplitude" => round(Qs_small; sigdigits=4),
                 "ratio" => round(Qs_meas/Qs_small; digits=3)])
end
```

## The pendulum analogy

For $\sin\phi_s=0$ ($\phi_s=0$ or $\pi$) the longitudinal Hamiltonian maps
one-to-one onto the physical pendulum:

| Longitudinal motion | | Pendulum motion |
|:---|:---:|:---|
| $H=\frac{1}{2}h\omega_0\eta\delta^2+\frac{eV\omega_0}{2\pi\beta^2E_0}\left(\cos\phi\pm 1\right)$ | Hamiltonian | $H=\frac{1}{2}m v^2+mgl\left(1-\cos\phi\right)$ |
| $h\omega_0\eta$ | 'Mass' | $m$ |
| $0$ ($\eta<0$), $\pi$ ($\eta>0$) | Stable phase | $0$ |
| $\pm\sqrt{\dfrac{2eV}{\pi\beta^2E_0h\lvert\eta\rvert}}$ | Bucket height | $2\sqrt{gl}$ |
| $\omega_0\sqrt{\frac{-eVh\eta\cos\phi_s}{2\pi\beta^2E_0}}$ | $\omega$ (small amplitude) | $\sqrt{g/l}$ |

## Longitudinal phase space area

The longitudinal emittance mirrors its transverse sibling; the common variable
pair is $(\Delta t, \Delta E)$, whose area is invariant under multiple RF
frequencies and slow acceleration:

$$
\mathcal{A}_{rms}=\pi\sqrt{\langle\Delta E\rangle^2\langle\Delta t\rangle^2-
\overline{\Delta E \Delta t}^2}
\approx \pi\langle\Delta E\rangle\langle\Delta t\rangle,
\qquad
\mathcal{A}_{95\%}\sim 6\mathcal{A}_{rms}.
$$

Its usual unit is eV·s. Mismatched injection and filamentation quickly wash
out correlations, leaving the rms estimate.

## The bucket, tracked turn by turn

The table above summarized which fixed point is stable in which regime; here
is that statement produced by actually tracking particles with TrackPad's
`AccelCavity` and `LongitudinalRFMap` — the same two elements, turn after turn.
The dashed curve is the analytic separatrix through the *unstable* fixed point;
whichever of $\phi_s$ or $\pi-\phi_s$ is **not** on that curve is the stable one:

```{code-cell} julia
:tags: [hide-input]

using TrackPad, TrackPadWidgets, StaticArrays

mass_p, αc, Vrf, hrf, frf = TrackPad.M_PROTON, 0.0018, 5.0e6, 360.0, 28.15e6
krf = 2π*frf/2.99792458e8
REGIME = [("20 GeV,  φs = φ",    20.0e9, false),
          ("20 GeV,  φs = π−φ",  20.0e9, true),
          ("250 GeV, φs = φ",   250.0e9, false),
          ("250 GeV, φs = π−φ", 250.0e9, true)]
NSAMP = 140

explorer(
    title   = "RF bucket: which fixed point is phase-stable depends on the sign of η",
    sliders = [Knob("|φ| [deg]", range(0.0, 70.0, length=8);
                      fmt = p -> string(Int(round(p))), init = 4),
               Knob("energy and phase choice", 1:4; fmt = i -> REGIME[i][1], init = 1)],
    panels  = [Panel(xlabel="RF phase φ [rad]", ylabel="ΔE [MeV]",
                     title="longitudinal phase space", xlim=(-π, 2π), ylim=(-400.0, 400.0),
                     height=300, legend=:bottomleft)],
    note = "Orbits come from TrackPad's AccelCavity + LongitudinalRFMap — exactly the "*
           "turn-by-turn map written above. The dashed curve is the analytic separatrix "*
           "through the unstable fixed point.",
) do φ0deg, ireg
    lbl, Etot, flip = REGIME[ireg]
    K0 = Etot - mass_p
    bm = Beam(K0; mass=mass_p, charge=1.0)
    β, γ = bm.beta, bm.gamma
    η   = αc - 1/γ^2
    P0c = β*Etot
    φs  = flip ? π - deg2rad(φ0deg) : deg2rad(φ0deg)
    cav  = AccelCavity(0.0; volt=Vrf, freq=frf, h=hrf, phis=φs, energy=K0, charge=1.0)
    ring = Lattice(AbstractElement[cav, LongitudinalRFMap(αc, cav)]; periodic=true)

    stable = η*cos(φs) < 0
    φ_u  = stable ? π - φs : φs                     # the other fixed point
    νs   = sqrt(Vrf*hrf*abs(η)/(2π*β^2*Etot))
    Qs   = νs*sqrt(abs(cos(φs)))

    pref = -Vrf/(π*hrf*η*β^2*Etot)
    φg = collect(range(-π, 2π, length=181))
    up = Union{Float64,Nothing}[]; dn = Union{Float64,Nothing}[]
    for φ in φg
        v = pref*(cos(φ) - cos(φ_u) + (φ - φ_u)*sin(φs))
        if v >= 0
            d = sqrt(v)*β^2*Etot/1e6
            push!(up, d); push!(dn, -d)
        else
            push!(up, nothing); push!(dn, nothing)
        end
    end

    np = 3
    stride = max(1, round(Int, 1.8/max(Qs, 1e-5)/NSAMP))
    c = zeros(np, 6)
    c[:, 6] .= collect(range(60.0e6, 240.0e6, length=np)) ./ P0c
    flags = zeros(Int, np)
    Φh = [Float64[] for _ in 1:np]; Eh = [Float64[] for _ in 1:np]
    for _ in 1:NSAMP
        ringpass!(c, ring, bm, flags, stride)
        for j in 1:np
            push!(Φh[j], φs - krf*c[j,5]); push!(Eh[j], c[j,6]*P0c/1e6)
        end
    end

    (series = vcat(
        [line(φg, up; color=PALETTE[4], dash=true, width=1.6, label="separatrix"),
         line(φg, dn; color=PALETTE[4], dash=true, width=1.6)],
        [line(Φh[j], Eh[j]; color=PALETTE[1], alpha=0.85, width=1.2,
              label = j == 1 ? "tracked orbits" : nothing) for j in 1:np],
        [points([φs, φ_u], [0.0, 0.0]; color=PALETTE[3], size=5.5, label="fixed points")]),
     readouts = ["η" => round(η; sigdigits=4),
                 "γ / γ_tr" => string(round(γ; digits=1), " / ", round(1/sqrt(αc); digits=1)),
                 "φs" => string(Int(round(rad2deg(φs))), "°"),
                 "η cos φs" => round(η*cos(φs); sigdigits=3),
                 "verdict" => stable ? "φs is phase-stable" : "φs is unstable",
                 "Qs" => stable ? string(round(Qs; sigdigits=4)) : "—",
                 "turns shown" => string(NSAMP*stride)])
end
```
