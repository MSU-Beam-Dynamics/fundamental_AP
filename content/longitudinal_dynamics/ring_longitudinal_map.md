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
| Injection energy | 20 GeV (total) |
| Top energy | 250 GeV (total) |
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
    longitudinal_map(turns, φ₀, ΔE₀; sin_φs, αc, E0_tot, mass, V, h, ramp=false)

Turn-by-turn longitudinal map in (φ, ΔE), with ΔE in eV. With `ramp=true` the
reference energy grows by eV sin φs every turn and β², η follow it; otherwise
they stay constant. Returns `(φ, ΔE, δ)`, each (turns+1) × nparticles, where
δ = ΔE/(β²E₀) is taken against the reference energy of that turn.
"""
function longitudinal_map(turns, φ₀, ΔE₀; sin_φs=0.5, αc=0.0018,
                          E0_tot=250e9, mass=0.938e9, V=5e6, h=360, ramp=false)
    E  = E0_tot
    β² = 1 - (mass/E)^2
    φ, ΔE = float.(φ₀), float.(ΔE₀)            # work on copies
    Φh, Eh, δh = [copy(φ)], [copy(ΔE)], [ΔE ./ (β²*E)]
    for _ in 1:turns
        ΔE .+= V .* (sin.(φ) .- sin_φs)         # RF kick
        ramp && (E += V*sin_φs)                 # the reference gains eV sin φs
        γ  = E/mass
        β² = 1 - 1/γ^2
        η  = αc - 1/γ^2                         # slip factor
        φ .+= (2π*h*η/(β²*E)) .* ΔE             # uses UPDATED energy ⇒ area preserving
        push!(Φh, copy(φ)); push!(Eh, copy(ΔE)); push!(δh, ΔE ./ (β²*E))
    end
    rows(hist) = permutedims(reduce(hcat, hist))
    return rows(Φh), rows(Eh), rows(δh)
end

mass, αc, V, h = 0.938e9, 0.0018, 5e6, 360
sin_φs = 0.5; φs = asin(sin_φs)
turns = 2100; npar = 5

fig = Figure(size=(900, 420))
for (panel, (label, E0_tot)) in enumerate([("250 GeV (above transition)", 250e9), ("20 GeV (below transition)", 20e9)])
    γ = E0_tot/mass; β² = 1-1/γ^2; η = αc-1/γ^2
    # bucket half-height, in δ = ΔE/(β²E₀): sets the starting points and the axes
    yf = √(cos(φs) - (π-2φs)*sin(φs)/2)
    height = 2*√(V/(2π*β²*E0_tot*h*abs(η)))*yf
    toΔE = β²*E0_tot                        # the map runs in ΔE [eV]
    pars = (; sin_φs, αc, E0_tot, mass, V, h)

    # one family starts at π−φs, the other at φs; the sign of η decides which stays bounded
    Hu = longitudinal_map(turns, fill(π-φs, npar),
                          collect(range(height/npar, 0.99height; length=npar)) .* toΔE; pars...)
    Hs = longitudinal_map(turns, fill(φs, npar),
                          collect(range(0, 0.99height; length=npar)) .* toΔE; pars...)

    ax = Axis(fig[1, panel]; limits=((-π, 2π), (-2height, 2height)), title=label,
              xlabel="Phase φ [rad]", ylabel = panel == 1 ? "δ = ΔE/(β²E₀)" : "")
    for i in 1:npar
        lines!(ax, Hu[1][:, i], Hu[3][:, i]; color=(:darkorange, 0.55))
        lines!(ax, Hs[1][:, i], Hs[3][:, i]; color=(:steelblue, 0.55))
    end
end
Legend(fig[2, :], [LineElement(color=:steelblue), LineElement(color=:darkorange)],
       ["start at φs", "start at π − φs"]; orientation=:horizontal, framevisible=false)
fig
```

:::{note}
The map runs in $\Delta E$; the figure plots $\delta=\Delta E/(\beta^2E_0)$.
While $E_0$ is constant the two differ only by a scale factor, but once the
reference energy ramps they do not, which is why $(\phi,\Delta E)$ is preferred
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
# fast acceleration in a low-energy proton ring (medical-synchrotron-like)
mass, αc, V, h = 0.938e9, 0.05, 0.1e6, 1
E0_tot = 45e6 + mass                      # 45 MeV kinetic
γ = E0_tot/mass; β² = 1-1/γ^2; η = αc-1/γ^2
sin_φs = 0.5; φs = asin(sin_φs)           # η < 0 here, so φs is the stable phase
yf = √(cos(φs) - (π-2φs)*sin(φs)/2)
Ĥ  = 2*√(V/(2π*β²*E0_tot*h*abs(η)))*yf * β²*E0_tot   # bucket half-height in ΔE [eV]

turns = 600                               # the reference gains eV sin φs = 50 keV per turn
ΔE₀   = [0.3, 0.6, 0.9] .* Ĥ              # three particles inside the bucket
pars  = (; sin_φs, αc, E0_tot, mass, V, h)
Hflat = longitudinal_map(turns, fill(φs, 3), ΔE₀; pars..., ramp=false)
Hramp = longitudinal_map(turns, fill(φs, 3), ΔE₀; pars..., ramp=true)

fig = Figure(size=(900, 420))
a1 = Axis(fig[1, 1]; xlabel="Phase φ [rad]", ylabel="δ = ΔE/(β²E₀)",
          title="(φ, δ): orbits spiral inward")
a2 = Axis(fig[1, 2]; xlabel="Phase φ [rad]", ylabel="ΔE [MeV]",
          title="(φ, ΔE): orbits do not close")
for (i, c) in enumerate((:steelblue, :darkorange, :seagreen))
    lines!(a1, Hflat[1][:, i], Hflat[3][:, i]; color=(:gray, 0.7), linestyle=:dash)
    lines!(a1, Hramp[1][:, i], Hramp[3][:, i]; color=c)
    lines!(a2, Hflat[1][:, i], Hflat[2][:, i] ./ 1e6; color=(:gray, 0.7), linestyle=:dash)
    lines!(a2, Hramp[1][:, i], Hramp[2][:, i] ./ 1e6; color=c)
end
Legend(fig[2, :], [LineElement(color=:gray, linestyle=:dash), LineElement(color=:steelblue)],
       ["no ramp (closed orbits)", "reference ramped by 50 keV per turn"];
       orientation=:horizontal, framevisible=false)
fig
```

Over these 600 turns the reference climbs from 45 to 75 MeV, $\beta^2E_0$ grows
by about 60 % and $2\pi h\eta/(\beta^2E_0)$ falls by about 40 % — too fast for
the orbits to follow. In $(\phi,\delta)$ every orbit spirals inward; in
$(\phi,\Delta E)$ they narrow in phase while growing in energy, so neither
picture closes on itself.

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

Let us verify $Q_s$ against the simulated map by timing the oscillation of the
tracked phase — and see how it changes away from the small-amplitude limit the
formula assumes:

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

"""
    oscillation_tune(signal)

Oscillations per turn of `signal`, from the turn numbers of its upward zero
crossings, each interpolated between the two turns that straddle it. Timing
whole periods, rather than counting crossings, resolves Qs to well under 1 %.
"""
function oscillation_tune(signal)
    t = Float64[]
    for n in 1:length(signal)-1
        a, b = signal[n], signal[n+1]
        a < 0 <= b && push!(t, n + a/(a - b))
    end
    length(t) < 2 ? NaN : (length(t) - 1)/(t[end] - t[1])
end

# The separatrix through the unstable fixed point π − φs turns back at the phase
# beyond φs where cos φ + φ sin φs takes the same value as at π − φs.
φu = π - φs_stable
on_separatrix(φ) = cos(φ) + φ*sin_φs - cos(φu) - φu*sin_φs
Δφ_sep = let a = φs_stable, b = φs_stable + 2.0          # bisect for that root
    for _ in 1:80
        m = (a + b)/2
        on_separatrix(a)*on_separatrix(m) <= 0 ? (b = m) : (a = m)
    end
    (a + b)/2 - φs_stable
end

# Qs over a fine grid of amplitudes, packed close to the separatrix where it falls steeply
amps_curve = [collect(0.02:0.02:1.18); Δφ_sep .- 10.0 .^ -(2:0.5:5)]
Hc = longitudinal_map(20000, φs_stable .+ amps_curve, zeros(length(amps_curve));
                      sin_φs=sin_φs, αc=αc, E0_tot=E0_tot, mass=mass, V=V, h=h)
Qs_curve = [oscillation_tune(Hc[1][:, j] .- φs_stable) for j in eachindex(amps_curve)]

explorer(
    title   = "Synchrotron tune softens with amplitude",
    sliders = [Knob("amplitude Δφ₀ [rad]", [collect(0.05:0.05:1.15); 1.19];
                    fmt = a -> string(round(a; digits=2)), init = 1)],
    panels  = [Panel(xlabel="phase φ [rad]", ylabel="ΔE [MeV]", title="tracked orbit",
                     height=280, legend=:bottomleft),
               Panel(xlabel="amplitude Δφ₀ [rad]", ylabel="Qs × 10³", title="tune versus amplitude",
                     xlim=(0.0, 1.3), ylim=(0.0, 1.5), height=280, legend=:bottomleft)],
    statics = [line(amps_curve, 1e3 .* Qs_curve; panel=2, color=PALETTE[1], label="tracked"),
               line([0.0, 1.3], fill(1e3*Qs_small, 2); panel=2, color=PALETTE[2], dash=true,
                    label="small-amplitude formula"),
               line([Δφ_sep, Δφ_sep], [0.0, 1.5]; panel=2, color=PALETTE[4], dash=true,
                    label="separatrix")],
    note = "Qs is timed from the zero crossings of the tracked phase and compared with the "*
           "small-amplitude formula Qs = νs√|cos φs|. Real buckets run slower far from the "*
           "centre, and Qs reaches 0 at the separatrix (Δφ₀ = $(round(Δφ_sep; digits=3)) rad "*
           "here) — but only logarithmically: an orbit 10⁻⁵ rad inside it still has about a "*
           "third of the small-amplitude tune.",
) do amp
    H = longitudinal_map(turns, [φs_stable+amp], [0.0];
                         sin_φs=sin_φs, αc=αc, E0_tot=E0_tot, mass=mass, V=V, h=h)
    Qs_meas = oscillation_tune(H[1][:,1] .- φs_stable)
    nshow = 2400
    (series = [line(H[1][1:nshow,1], H[2][1:nshow,1]./1e6; panel=1, color=PALETTE[1], label="(φ, ΔE)"),
               points([φs_stable], [0.0]; panel=1, color=PALETTE[3], size=5.0, label="stable fixed point"),
               points([amp], [1e3*Qs_meas]; panel=2, color=PALETTE[1], size=6.0, label="this orbit")],
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

The figure below tracks Gaussian bunches through the turn-by-turn map of this
chapter, in the RHIC-like ring with $\phi_s=150^\circ$ of the tune figure above.
Its time coordinate is $\Delta t=(\phi-\phi_s)/\omega_{\mathrm{RF}}$, positive
for a particle that arrives late. Four bunches are tracked at 250 GeV:

- **matched**: the ratio $\sigma_{\Delta E}/\sigma_{\Delta t}$ fits the
  small-amplitude bucket, so the bunch rotates in place and keeps its area;
- **half the matched length**, with twice the energy spread and the same area:
  the bunch tumbles at twice the synchrotron frequency. While it is tilted, the
  correlation term matters and $\pi\sigma_{\Delta t}\sigma_{\Delta E}$
  overestimates the area; as the bunch filaments, its true rms area grows;
- **half the matched length with a tenth of the area**: the same mismatch in a
  smaller bunch, which stays near the bucket centre where the motion is nearly
  linear. All its particles share almost the same synchrotron tune, so it keeps
  tumbling as one body: $\pi\sigma_{\Delta t}\sigma_{\Delta E}$ swings back and
  forth while the rms area holds still. The tune spread grows with the square of
  the amplitude, so this bunch filaments too, but several times later;
- **injected 0.3 rad late**, with the matched area: the whole bunch oscillates
  about the bucket centre and smears out around it.

What happens to the area during acceleration, including across transition, is
the subject of the last section.

```{code-cell} julia
:tags: [hide-input]

using Random

# The RHIC-like ring of the tune figure at 250 GeV: φs = 150°, above transition
mass, αc, V, h = 0.938e9, 0.0018, 5e6, 360
sin_φs = 0.5; φs = π - asin(sin_φs)
ω_rf = 2π*28.15e6                       # Δt = (φ − φs)/ω_rf, positive = arrives late
const E_LE = 250e9                      # total energy of the bunches [eV]
slip(E)  = αc - (mass/E)^2              # η at total energy E
beta2(E) = 1 - (mass/E)^2

"σ_ΔE/σ_φ of a bunch matched to the small-amplitude bucket at total energy E [eV/rad]."
function matched_ratio(E)
    Qs = sqrt(V*h*abs(slip(E))/(2π*beta2(E)*E))*sqrt(abs(cos(φs)))
    Qs*beta2(E)*E/(h*abs(slip(E)))
end

"""
    bunch_moments(Δt, ΔE)

Centroid, second moments, the rms area π√(σt²σE² − ⟨ΔtΔE⟩²), its uncorrelated
estimate π σt σE, and r95: the ellipse shaped like the rms one that holds 95 %
of the particles has area r95 × A_rms.
"""
function bunch_moments(Δt, ΔE)
    N  = length(Δt)
    mt, mE = sum(Δt)/N, sum(ΔE)/N
    dt, dE = Δt .- mt, ΔE .- mE
    vt, vE, c = sum(abs2, dt)/N, sum(abs2, dE)/N, sum(dt .* dE)/N
    d  = vt*vE - c^2
    r2 = sort((dE.^2 .* vt .- 2c .* dt .* dE .+ dt.^2 .* vE) ./ d)
    (; mt, mE, vt, vE, c, rms = π*sqrt(d), uncorr = π*sqrt(vt*vE), r95 = r2[ceil(Int, 0.95N)])
end

"Ellipse of the bunch's rms shape scaled to area k × A_rms."
function moment_ellipse(m, k; n=90)
    L11 = sqrt(m.vt); L21 = m.c/L11; L22 = sqrt(m.vE - L21^2)
    θ = range(0, 2π, length=n)
    m.mt .+ sqrt(k)*L11 .* cos.(θ), m.mE .+ sqrt(k) .* (L21 .* cos.(θ) .+ L22 .* sin.(θ))
end

# The bucket runs from the unstable fixed point φu to the far tip φs + Δφ_sep (tune figure).
φ_lo, φ_hi = π - φs, φs + Δφ_sep
"Upper half of the separatrix at total energy E, as (Δt [ns], ΔE [MeV])."
function separatrix(E; n=80)
    pref = -V/(π*h*slip(E)*beta2(E)*E)
    φg = collect(range(φ_lo, φ_hi, length=n))
    ΔEsep = [beta2(E)*E*sqrt(max(pref*(cos(φ) - cos(φ_lo) + (φ - φ_lo)*sin_φs), 0.0))/1e6
             for φ in φg]
    round.((φg .- φs) ./ ω_rf .* 1e9; digits=2), round.(ΔEsep; digits=1)
end

const NLE = 1500
const NPLOT = 800                       # particles drawn; the statistics use all NLE
const σφ0 = 0.15                        # rms length of the matched bunch [rad of RF phase]
LE_SNAP  = [0, 50, 100, 200, 300, 400, 800, 1000, 2000, 3000, 5000, 8000, 12000,
            20000, 30000, 60000]
LE_TRACE = [10:10:1000; 1050:50:6000; 6100:100:60000]        # turns where the areas are logged

# (label, length factor r, area factor a, phase error [rad]): the bunch is r times
# its matched length, and its rms area a times that of the matched bunch
LE_CASES = [("matched",                       1.0, 1.0, 0.0),
            ("half the matched length",       0.5, 1.0, 0.0),
            ("half length, 1/10 of the area", 0.5, 0.1, 0.0),
            ("injected 0.3 rad late",         1.0, 1.0, 0.3)]
# fully filamented rms area in linear theory, in units of the starting area
filamented(r, a, Δφ) = (r^2 + 1/r^2)/2 + (Δφ/(σφ0*sqrt(a)))^2/2

function emittance_run(r, a, Δφ)
    rng = MersenneTwister(7)
    σφ = σφ0*sqrt(a)                       # matched length for this area
    φ  = φs + Δφ .+ (r*σφ) .* randn(rng, NLE)
    ΔE = (σφ*matched_ratio(E_LE)/r) .* randn(rng, NLE)
    t = 0
    alive = trues(NLE)
    A0 = bunch_moments((φ .- φs) ./ ω_rf, ΔE).rms
    snaps = Dict{Int,Any}(); trace = (t = Int[], rms = Float64[], uncorr = Float64[])
    for target in sort(unique([LE_SNAP; LE_TRACE]))
        if target > t
            Φh, Eh, _ = longitudinal_map(target - t, φ, ΔE; sin_φs, αc, E0_tot=E_LE, mass, V, h)
            φ, ΔE = Φh[end, :], Eh[end, :]
            t = target
        end
        alive .&= φ_lo .< φ .< φ_hi             # a particle that leaves the bucket is lost
        Δt = (φ[alive] .- φs) ./ ω_rf
        m  = bunch_moments(Δt, ΔE[alive])
        if t in LE_TRACE
            push!(trace.t, t); push!(trace.rms, round(m.rms/A0; digits=3))
            push!(trace.uncorr, round(m.uncorr/A0; digits=3))
        end
        shown = 1:min(NPLOT, length(Δt))
        t in LE_SNAP && (snaps[t] = (Δt = round.(Δt[shown] .* 1e9; digits=2),
                                     ΔE = round.(Int, ΔE[alive][shown] ./ 1e6),
                                     m = m, A0 = A0, lost = NLE - count(alive)))
    end
    (snaps = snaps, trace = trace)
end
le_runs = [emittance_run(c[2:end]...) for c in LE_CASES]
sep_t, sep_E = separatrix(E_LE)

explorer(
    title   = "Longitudinal emittance: matched and mismatched bunches",
    sliders = [Knob("turn", LE_SNAP; fmt = string, init = 1),
               Knob("bunch", 1:length(LE_CASES); fmt = i -> LE_CASES[i][1], init = 2,
                    newline = true)],
    panels  = [Panel(xlabel="Δt [ns]", ylabel="ΔE [MeV]", title="the bunch in (Δt, ΔE)",
                     xlim=(-12.5, 7.5), ylim=(-750.0, 750.0), height=320, legend=:bottomleft,
                     basis="48%", minwidth=280),
               Panel(xlabel="turn", ylabel="area / starting rms area", title="area versus turn",
                     xscale=:log10, xlim=(10.0, 6.0e4), ylim=(0.0, 4.0), height=320,
                     legend=:topleft, basis="48%", minwidth=280)],
    statics = [line(sep_t, sep_E; panel=1, color=PALETTE[4], dash=true, width=1.6, label="separatrix"),
               line(sep_t, -sep_E; panel=1, color=PALETTE[4], dash=true, width=1.6)],
    note = "The dashed line is the linear-theory area after full filamentation, "*
           "(r² + 1/r²)/2 + (Δφ/σφ)²/2 for a bunch r times the matched length injected Δφ "*
           "off phase; the bucket's nonlinearity makes the real growth a little larger. "*
           "A₉₅ is the area of the ellipse, shaped like the rms ellipse, that holds 95 % of "*
           "the particles: about 6 A_rms while the bunch is Gaussian, but not once it has "*
           "filamented.",
) do turn, icase
    label, r, a, Δφ = LE_CASES[icase]
    run = le_runs[icase]
    s = run.snaps[turn]
    ex, ey = moment_ellipse(s.m, s.m.r95)
    k = searchsortedlast(run.trace.t, turn)
    series = Any[points(s.Δt, s.ΔE; panel=1, color=PALETTE[1], size=1.6, alpha=0.55, label="protons"),
                 line(round.(ex .* 1e9; digits=3), round.(ey ./ 1e6; digits=1); panel=1, color=PALETTE[2], width=2.0, label="95 % ellipse"),
                 line([10.0, 6.0e4], fill(filamented(r, a, Δφ), 2); panel=2, color="#5b6472",
                      dash=true, label="filamented, linear theory")]
    if k > 0
        push!(series, line(run.trace.t[1:k], run.trace.uncorr[1:k]; panel=2, color=PALETTE[2],
                           width=1.4, label="π σt σE"))
        push!(series, line(run.trace.t[1:k], run.trace.rms[1:k]; panel=2, color=PALETTE[1],
                           width=2.0, label="rms area"))
        push!(series, points([Float64(turn)], [run.trace.rms[k]]; panel=2, color=PALETTE[1], size=5.0))
    end
    eVs(A) = string(round(A; sigdigits=3), " eV·s")
    (series = series,
     readouts = ["rms area" => eVs(s.m.rms),
                 "π σt σE" => eVs(s.m.uncorr),
                 "rms area / start" => round(s.m.rms/s.A0; digits=3),
                 "A₉₅ / A_rms" => round(s.m.r95; digits=2),
                 "σt, σE" => string(round(sqrt(s.m.vt)*1e9; digits=2), " ns, ",
                                    round(sqrt(s.m.vE)/1e6; digits=1), " MeV"),
                 "lost" => string(s.lost, " of ", NLE)])
end
```

## The bucket, tracked turn by turn

The table above summarized which fixed point is stable in which regime; here
is that statement produced by actually tracking particles through the RF cavity
and the rest of the ring, turn after turn.
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
                     title="longitudinal phase space", xlim=(-π, 2π), autoscale=:frame,
                     height=300, legend=:bottomleft)],
    note = "The orbits are tracked turn by turn with exactly the map written above. "*
           "The dashed curve is the analytic separatrix "*
           "through the unstable fixed point. The ΔE axis follows the bucket height, "*
           "which grows with |η|/(β²E₀) and shrinks as |φs| moves away from 0 or π.",
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
    φ_u  = stable ? π - φs : φs                     # the unstable fixed point
    φ_c  = π - φ_u                                  # the stable one, centre of the bucket
    νs   = sqrt(Vrf*hrf*abs(η)/(2π*β^2*Etot))
    Qs   = νs*sqrt(abs(cos(φs)))

    # separatrix through φ_u: ΔE(φ) = ±β²E₀·√sep(φ); its value at φ_c is the bucket height
    pref = -Vrf/(π*hrf*η*β^2*Etot)
    sep(φ) = pref*(cos(φ) - cos(φ_u) + (φ - φ_u)*sin(φs))
    Ĥ    = sqrt(sep(φ_c))*β^2*Etot                  # bucket half-height [eV]
    ymax = 1.6*Ĥ/1e6                                # plotted ΔE range [MeV]
    clip(v) = abs(v) <= ymax ? v : nothing

    φg = collect(range(-π, 2π, length=181))
    tips = Float64[]                                # where the separatrix closes: sep = 0
    for i in 1:length(φg)-1
        a, b = φg[i], φg[i+1]
        sep(a)*sep(b) < 0 || continue
        for _ in 1:60
            m = (a + b)/2
            sep(a)*sep(m) <= 0 ? (b = m) : (a = m)
        end
        push!(tips, (a + b)/2)
    end
    pts = sort!(vcat([(φ, sep(φ)) for φ in φg], [(φ, 0.0) for φ in tips]))
    φsep = first.(pts)
    up = [v >= 0 ? clip( sqrt(v)*β^2*Etot/1e6) : nothing for (_, v) in pts]
    dn = [v >= 0 ? clip(-sqrt(v)*β^2*Etot/1e6) : nothing for (_, v) in pts]

    np = 3
    stride = max(1, round(Int, 1.8/max(Qs, 1e-5)/NSAMP))
    c = zeros(np, 6)
    c[:, 6] .= [0.25, 0.55, 0.85] .* Ĥ ./ P0c       # start at φs, scaled to the bucket
    flags = zeros(Int, np)
    Φh = [Float64[] for _ in 1:np]; Eh = [Union{Float64,Nothing}[] for _ in 1:np]
    for _ in 1:NSAMP
        ringpass!(c, ring, bm, flags, stride)
        for j in 1:np
            push!(Φh[j], φs - krf*c[j,5]); push!(Eh[j], clip(c[j,6]*P0c/1e6))
        end
    end

    (series = vcat(
        [line(φsep, up; color=PALETTE[4], dash=true, width=1.6, label="separatrix"),
         line(φsep, dn; color=PALETTE[4], dash=true, width=1.6)],
        [line(Φh[j], Eh[j]; color=PALETTE[1], alpha=0.85, width=1.2,
              label = j == 1 ? "tracked orbits" : nothing) for j in 1:np],
        [points([φs, φ_u], [0.0, 0.0]; color=PALETTE[3], size=5.5, label="fixed points")]),
     readouts = ["η" => round(η; sigdigits=4),
                 "γ / γ_tr" => string(round(γ; digits=1), " / ", round(1/sqrt(αc); digits=1)),
                 "φs" => string(Int(round(rad2deg(φs))), "°"),
                 "η cos φs" => round(η*cos(φs); sigdigits=3),
                 "verdict" => stable ? "φs is phase-stable" : "φs is unstable",
                 "Qs" => stable ? string(round(Qs; sigdigits=4)) : "—",
                 "bucket half-height" => string(round(Ĥ/1e6; sigdigits=3), " MeV"),
                 "turns shown" => string(NSAMP*stride)])
end
```

## Crossing transition

The RHIC-like ring of this chapter injects protons at 20 GeV, below the
transition energy $\gamma_t=1/\sqrt{\alpha_c}\approx 23.6$ (22.1 GeV total), and
accelerates them to 250 GeV, so its ramp crosses transition. (In the real RHIC
it is the heavy ions that are injected below transition; protons arrive above
it.[^gtjump]) There $\eta$ changes sign, and by
the stability condition $\eta\cos\phi_s<0$ the stable phase moves from $\phi_s$
to $\pi-\phi_s$: with $\sin\phi_s=0.5$, from $30^\circ$ to $150^\circ$. At the
moment of crossing the RF system has to shift its phase by
$\pi-2\phi_s=120^\circ$, and it has to do so abruptly. On its way from
$30^\circ$ to $150^\circ$ the phase passes the crest, so during the switch the
synchronous particle gains more than the $eV\sin\phi_s$ that the magnet ramp
asks for, and the whole bunch ends up off energy. (Centring the switch on
transition at least keeps the bunch on the focusing slope throughout: the phase
is below $90^\circ$ while $\eta<0$ and above $90^\circ$ once $\eta>0$.)

Transition is hard on the bunch even with a perfect jump. As $\eta\to 0$ the
synchrotron motion freezes, so the motion stops being adiabatic, and the bunch
becomes short with a large energy spread.

The figure below ramps a matched bunch from 20 GeV with the map of this chapter,
2.5 MeV per turn, so transition comes at turn 844. The phase is switched either
over 10 turns or over 100 turns, centred on the crossing. $\Delta t$ is measured
from the synchronous phase that the RF program sets on that turn.

```{code-cell} julia
:tags: [hide-input]

# RHIC-like ramp through transition: the reference gains eV sin φs = 2.5 MeV per turn
E_inj = 20e9                                   # injection, total energy [eV]
E_tr  = mass/sqrt(αc)                          # transition: γ = γt = 1/√αc
n_tr  = (E_tr - E_inj)/(V*sin_φs)              # turn on which the bunch reaches transition
φs_below, φs_above = asin(sin_φs), π - asin(sin_φs)

const NTR = 1500
TR_TURNS  = 20000
TR_SNAP   = [0, 200, 400, 600, 700, 800, 850, 900, 1000, 1200, 1500, 2000, 3000,
             4000, 6000, 8000, 12000, 16000, 20000]
TR_TRACE  = [100:10:3000; 3050:50:TR_TURNS]    # turns where the area and centroid are logged
TR_SWITCH = [10, 100]                          # turns taken by the phase switch

"Bucket edges in φ for stable phase φs: from the unstable fixed point π − φs to the far tip."
bucket_edges(φs) = φs < π/2 ? (φs - Δφ_sep, π - φs) : (π - φs, φs + Δφ_sep)

"Upper half of the separatrix at total energy E around stable phase φs, as (Δt [ns], ΔE [MeV])."
function bucket(E, φs; n=80)
    φu   = π - φs
    pref = -V/(π*h*slip(E)*beta2(E)*E)
    φg   = collect(range(bucket_edges(φs)..., length=n))
    ΔEs  = [beta2(E)*E*sqrt(max(pref*(cos(φ) - cos(φu) + (φ - φu)*sin_φs), 0.0))/1e6 for φ in φg]
    round.((φg .- φs) ./ ω_rf .* 1e9; digits=2), round.(ΔEs; digits=1)
end

function cross_transition(nswitch)
    rng = MersenneTwister(11)
    φ   = φs_below .+ σφ0 .* randn(rng, NTR)
    ΔE  = (σφ0*matched_ratio(E_inj)) .* randn(rng, NTR)
    E, ψ = E_inj, 0.0                          # ψ: how far the RF phase has been shifted
    t_on  = round(Int, n_tr - nswitch/2)       # the switch is centred on transition
    alive = trues(NTR)
    snaps = Dict{Int,Any}(); trace = (t = Int[], A = Float64[], mE = Float64[])
    for t in 0:TR_TURNS
        if t > 0
            Φh, Eh, _ = longitudinal_map(1, φ, ΔE; sin_φs, αc, E0_tot=E, mass, V, h, ramp=true)
            φ, ΔE = Φh[end, :], Eh[end, :]
            E += V*sin_φs
            if t_on < t <= t_on + nswitch      # the RF wave moves by (π − 2φs)/nswitch per turn
                φ .+= (φs_above - φs_below)/nswitch
                ψ += (φs_above - φs_below)/nswitch
            end
        end
        φsync = φs_below + ψ                   # the synchronous phase the RF program sets
        switching = t_on < t < t_on + nswitch
        if !switching                          # a particle outside the bucket is lost
            lo, hi = bucket_edges(φsync)
            alive .&= lo .< φ .< hi
        end
        Δt = (φ[alive] .- φsync) ./ ω_rf
        m  = bunch_moments(Δt, ΔE[alive])
        if t in TR_TRACE
            push!(trace.t, t); push!(trace.A, round(m.rms; sigdigits=4))
            push!(trace.mE, round(m.mE/1e6; digits=1))
        end
        shown = 1:min(NPLOT, length(Δt))
        t in TR_SNAP && (snaps[t] = (Δt = round.(Δt[shown] .* 1e9; digits=2),
                                     ΔE = round.(Int, ΔE[alive][shown] ./ 1e6),
                                     m = m, E = E, φsync = φsync, switching = switching,
                                     lost = NTR - count(alive)))
    end
    (snaps = snaps, trace = trace)
end
tr_runs = [cross_transition(n) for n in TR_SWITCH]

TR_COLORS = (PALETTE[1], PALETTE[2])
tr_statics = Any[]
for (j, (n, run)) in enumerate(zip(TR_SWITCH, tr_runs))
    push!(tr_statics, line(run.trace.t, run.trace.A; panel=2, color=TR_COLORS[j], label="switch in $n turns"))
    push!(tr_statics, line(run.trace.t, run.trace.mE; panel=3, color=TR_COLORS[j], label="switch in $n turns"))
end
for p in (2, 3)
    lim = p == 2 ? (0.0, 0.45) : (-200.0, 200.0)
    push!(tr_statics, line([n_tr, n_tr], collect(lim); panel=p, color="#5b6472", dash=true,
                           label = p == 2 ? "transition" : nothing))
end

explorer(
    title   = "Crossing transition: the RF phase has to jump",
    sliders = [Knob("turn", TR_SNAP; fmt = string, init = 1),
               Knob("phase switch", 1:length(TR_SWITCH); fmt = i -> "in $(TR_SWITCH[i]) turns",
                    init = 1, newline = true)],
    panels  = [Panel(xlabel="Δt [ns]", ylabel="ΔE [MeV]", title="the bunch in (Δt, ΔE)",
                     xlim=(-12.5, 12.5), ylim=(-600.0, 600.0), height=300, legend=:bottomleft,
                     basis="100%"),
               Panel(xlabel="turn", ylabel="rms area [eV·s]", title="area",
                     xscale=:log10, xlim=(100.0, 2.0e4), ylim=(0.0, 0.45), height=220,
                     legend=:topleft, basis="48%", minwidth=260),
               Panel(xlabel="turn", ylabel="⟨ΔE⟩ [MeV]", title="bunch centroid",
                     xscale=:log10, xlim=(100.0, 2.0e4), ylim=(-200.0, 200.0), height=220,
                     legend=:bottomleft, basis="48%", minwidth=260)],
    statics = tr_statics,
    note = "The RF phase is switched from φs = 30° to 150°, centred on the transition turn "*
           "($(round(Int, n_tr))). The dashed bucket is drawn when it fits the plot: near "*
           "transition η → 0 and its height diverges, and during the switch there is no "*
           "stable phase to draw it around.",
) do turn, isw
    run = tr_runs[isw]
    s = run.snaps[turn]
    series = Any[points(s.Δt, s.ΔE; panel=1, color=TR_COLORS[isw], size=1.6, alpha=0.55,
                        label="protons")]
    if !s.switching
        bt, bE = bucket(s.E, s.φsync)
        if maximum(bE) < 600
            push!(series, line(bt, bE; panel=1, color=PALETTE[4], dash=true, width=1.6, label="bucket"))
            push!(series, line(bt, -bE; panel=1, color=PALETTE[4], dash=true, width=1.6))
        end
    end
    k = searchsortedlast(run.trace.t, turn)
    if k > 0
        push!(series, points([Float64(turn)], [run.trace.A[k]]; panel=2, color=TR_COLORS[isw], size=5.0))
        push!(series, points([Float64(turn)], [run.trace.mE[k]]; panel=3, color=TR_COLORS[isw], size=5.0))
    end
    γ = s.E/mass
    (series = series,
     readouts = ["energy" => string(round(s.E/1e9; digits=2), " GeV"),
                 "γ / γt" => round(γ*sqrt(αc); digits=3),
                 "η" => round(slip(s.E); sigdigits=3),
                 "programmed φs" => string(round(Int, rad2deg(s.φsync)), "°"),
                 "rms area" => string(round(s.m.rms; sigdigits=3), " eV·s"),
                 "σt, σE" => string(round(sqrt(s.m.vt)*1e9; digits=2), " ns, ",
                                    round(sqrt(s.m.vE)/1e6; digits=1), " MeV"),
                 "⟨ΔE⟩" => string(round(s.m.mE/1e6; digits=1), " MeV"),
                 "lost" => string(s.lost, " of ", NTR)])
end
```

With the 10-turn switch the bunch crosses transition with its rms area intact,
0.119 → 0.121 eV·s, even though it is squeezed to about half its length on the
way. With the 100-turn switch it comes out about 160 MeV off energy, oscillates
about the bucket centre, and filaments to three times its area by 70 GeV.

Real machines go further and also make transition itself pass quickly, with a
$\gamma_t$ jump: a set of quadrupoles whose current is reversed just as the beam
arrives, changing $\alpha_c$ and so moving $\gamma_t$ past the beam's $\gamma$.
RHIC has eight such quadrupoles in each sextant, and reverses their current in
40–60 ms; that is short enough for its ions to cross transition with no
particle loss and negligible growth of the bunch area.[^gtjump]

[^gtjump]: J. Mi *et al.*, "RHIC Gamma Transition Jump Power Supply Prototype
    Test", PAC 2001, BNL-68035, <https://www.osti.gov/servlets/purl/783176>.
