---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# The Synchrotron

## Why a synchrotron?

A charged particle in a uniform magnetic field circulates with angular frequency

$$
\omega=\frac{qB}{\gamma m}.
$$

A cyclotron holds $B$ and the RF frequency fixed and lets the orbit spiral
outward as the particle gains energy. That works only while $\gamma\approx 1$:
once the particle becomes relativistic its revolution frequency drops and it
slips out of step with the RF. A synchrotron turns the problem around and keeps
the **orbit** fixed instead:

| | Magnetic field | Revolution frequency | Orbit |
|:--|:--:|:--:|:--:|
| Cyclotron | fixed | changes for a relativistic particle | grows with energy |
| Synchrotron | ramped | kept in step by the RF | kept constant |

With a fixed orbit the magnets only have to cover a thin ring rather than a
whole disk, which is what makes very high energies affordable. A synchrotron
works in cycles: inject at low field, ramp the field and the RF together while
accelerating, hold at flat top while the beam is extracted or stored, then ramp
back down for the next cycle.


```{figure} ../images/synchrotron/synchrotron_cycle.png
:width: 70%
:name: fig-synchrotron-cycle
A synchrotron cycle. The beam is injected at low field, accelerated while the bending field $B$ and the momentum $pc$ ramp together, held at flat top while it is ejected, and the magnets then recover for the next cycle. The orbit stays fixed throughout.
```

## Early synchrotrons

- **Cosmotron**, Brookhaven National Laboratory — the first proton synchrotron
  to pass 1 GeV, reaching 3.3 GeV in the 1950s. It was a *weak-focusing*
  machine, so its large beam needed large magnet apertures.
- **ADA** (*Anello Di Accumulazione*), Frascati, 1961–1964 — a 250 MeV
  electron–positron storage ring, the first machine to store and collide beams
  of matter and antimatter.


```{figure} ../images/synchrotron/first_synchrotron.jpg
:width: 55%
:name: fig-first-synchrotron
The first synchrotron.
```

```{figure} ../images/synchrotron/cosmotron.jpg
:width: 80%
:name: fig-cosmotron
The Cosmotron at Brookhaven National Laboratory, a weak-focusing proton synchrotron that reached 3.3 GeV.
```

```{figure} ../images/synchrotron/ada.jpg
:width: 60%
:name: fig-ada
ADA (*Anello Di Accumulazione*) at Frascati, the 250 MeV ring in which electrons and positrons were first stored and collided.
```

## Operating principle

As the momentum grows, the bending field is ramped so that the design radius
$\rho$ never changes:

$$
P(t)=qB(t)\rho
\quad\Longleftrightarrow\quad
\gamma(t)\beta(t)mc=qB(t)\rho .
$$

The angular frequency of the motion through the bends is then

$$
\omega(t)=\frac{qB(t)}{m\gamma(t)}=\frac{\beta(t)c}{\rho},
$$

and over the whole circumference $C$ the revolution frequency is
$f_0=\beta c/C$. The RF has to run at a harmonic $f_{\mathrm{RF}}=hf_0$ of it,
so it must follow $\beta(t)$ through the ramp. Only in the ultra-relativistic
limit, $\beta\to 1$, do the orbit and the frequency stay synchronized with a
fixed RF frequency. Expressed in practical units, the magnetic rigidity is

$$
B\rho\,[\mathrm{T\,m}] = \frac{P\,[\mathrm{GeV}/c]}{0.2998}.
$$


```{code-cell} julia
:tags: [hide-input]
:label: fig-synchronization
:caption: Synchronization in a proton synchrotron of fixed radius. The revolution frequency, relative to that of a particle moving at the speed of light, is $\beta$ (left axis); the magnetic rigidity $B\rho$ (right axis) grows with the momentum. Above a few GeV the frequency hardly changes, while $B\rho$ keeps growing.

using CairoMakie, TrackPadWidgets

mp = 0.93827                                    # proton rest energy [GeV]
K  = exp10.(range(-1, 2, length=300))           # kinetic energy [GeV]
β  = @. sqrt(1 - (mp/(K + mp))^2)               # f / f(v = c)
Bρ = @. sqrt((K + mp)^2 - mp^2)/0.299792458     # P/q in T·m, P in GeV/c

fig = Figure(size=(720, 360), fontsize=15)
ax1 = Axis(fig[1, 1]; xscale=log10, xlabel="kinetic energy [GeV]",
           ylabel="revolution frequency  f / f(v = c)", limits=(nothing, (0.3, 1.05)))
ax2 = Axis(fig[1, 1]; xscale=log10, yscale=log10, yaxisposition=:right, ylabel="Bρ [T·m]",
           limits=(nothing, (0.1, 1000)))
hidespines!(ax2); hidexdecorations!(ax2); linkxaxes!(ax1, ax2)
l1 = lines!(ax1, K, β; color=PALETTE[1], linewidth=2.5)
l2 = lines!(ax2, K, Bρ; color=PALETTE[2], linewidth=2.5)
axislegend(ax1, [l1, l2], ["frequency (left)", "Bρ (right)"]; position=:rc, framevisible=false)
fig
```


## Energy dynamic range and accelerator complexes

No single accelerator takes particles from the source all the way to a high
top energy:

- **The frequency changes a lot** at low energy, where $\beta$ is far from 1.
- **Space charge** is strongest at low energy.
- **The magnets have a limited dynamic range**: remanent fields spoil the field
  quality at very low field and the iron saturates at high field.

So a linac takes the beam to tens or hundreds of MeV, and each synchrotron after
it boosts the energy by a factor of roughly 10–100. Reaching high energy takes a
*chain* of machines, and moving the beam into and out of each one — injection
and extraction, the subject of the next two pages — becomes essential. The same
chain also re-fills storage rings such as the LHC and feeds external
experiments.

The two complexes below are typical:

- **CERN (LHC)**: Linac → Proton Synchrotron Booster (PSB, 50 MeV → 1.4 GeV)
  → Proton Synchrotron (PS, 1.4 → 25 GeV, 628 m) → Super Proton Synchrotron
  (SPS, 25 → 450 GeV, 7 km) → LHC (450 GeV → 7 TeV). The complex also serves
  the Low Energy Ion Ring (LEIR), the Antiproton Decelerator (AD), ISOLDE and
  fixed-target beams.[^liu]
- **Brookhaven (RHIC)**: EBIS ion source / linac → Booster (200 MeV protons, or
  2 MeV/nucleon ions → 100 MeV/nucleon) → AGS (Alternating Gradient Synchrotron,
  ions to about 9 GeV/nucleon) → AGS-to-RHIC transfer line → RHIC.


```{figure} ../images/synchrotron/lhc_complex.png
:width: 90%
:name: fig-lhc-complex
The CERN accelerator complex that feeds the LHC.
```

```{figure} ../images/synchrotron/cern_synchrotrons.jpg
:width: 100%
:name: fig-cern-synchrotrons
The synchrotrons of the LHC injector chain: the Proton Synchrotron Booster (50 MeV → 1.4 GeV), the Low Energy Ion Ring, the Proton Synchrotron (1.4 → 25 GeV, 628 m) and the Super Proton Synchrotron (25 → 450 GeV, 7 km).
```

```{figure} ../images/synchrotron/rhic_cern_complexes.png
:width: 100%
:name: fig-rhic-cern-complexes
The RHIC complex (left): 1 EBIS ion source, 2 linac, 3 Booster, 4 AGS, 5 AGS-to-RHIC line, 6 RHIC. The CERN complex is shown on the right for comparison.
```

```{figure} ../images/synchrotron/rhic_booster_ags.jpg
:width: 100%
:name: fig-rhic-booster-ags
Synchrotrons in the RHIC complex: the Booster, from 200 MeV protons or 2 MeV/nucleon ions to 100 MeV/nucleon (left), and the Alternating Gradient Synchrotron, to about 9 GeV/nucleon (right).
```

[^liu]: These are the energies of the LHC injector chain before the 2021
    upgrade; since then Linac4 injects H⁻ into the PSB at 160 MeV and the PSB
    extracts at 2 GeV.

The widget below follows a proton through the CERN chain. For each ring it marks
the stretch of $\beta$ (which sets the revolution and RF frequency) and of
$B\rho$ (which sets the bending field) that the ring has to cover:

```{code-cell} julia
:tags: [hide-input]

using TrackPad, TrackPadWidgets

const CLIGHT = 2.99792458e8
proton(K) = Beam(K; mass=TrackPad.M_PROTON, charge=1.0)   # K = kinetic energy [eV]
rigidity(K) = (b = proton(K); b.beta*(K + b.mass)/CLIGHT)    # Bρ = P/q  [T·m]

# (name, injection and top kinetic energy [eV], circumference [m])
RINGS = [("PSB",  50.0e6,  1.4e9,    157.08),
         ("PS",    1.4e9, 25.0e9,    628.32),
         ("SPS",  25.0e9, 450.0e9,  6911.5),
         ("LHC", 450.0e9,  7.0e12, 26658.9)]
RING_COLORS = (PALETTE[2], PALETTE[3], PALETTE[4], PALETTE[5])

Kgrid(a, b; n=60) = exp10.(range(log10(a), log10(b), length=n))
energy_label(K) = K < 1e9 ? string(round(K/1e6; sigdigits=3), " MeV") :
                  K < 1e12 ? string(round(K/1e9; sigdigits=3), " GeV") :
                             string(round(K/1e12; sigdigits=3), " TeV")
freq_label(f) = f < 1e6 ? string(round(f/1e3; sigdigits=3), " kHz") :
                          string(round(f/1e6; sigdigits=4), " MHz")

Kall = Kgrid(10e6, 10e12; n=200)
statics = Any[line(Kall ./ 1e9, [proton(K).beta for K in Kall]; panel=1, color="#9aa4b2", width=1.2),
              line(Kall ./ 1e9, rigidity.(Kall); panel=2, color="#9aa4b2", width=1.2)]
for ((name, Ki, Kt, _), col) in zip(RINGS, RING_COLORS)
    Ks = Kgrid(Ki, Kt)
    push!(statics, line(Ks ./ 1e9, [proton(K).beta for K in Ks]; panel=1, color=col, width=2.0,
                        alpha=0.45, label=name))
    push!(statics, line(Ks ./ 1e9, rigidity.(Ks); panel=2, color=col, width=2.0, alpha=0.45))
end

explorer(
    title   = "One ring, one energy range: the LHC injector chain",
    sliders = [Knob("ring", 1:length(RINGS); fmt = i -> RINGS[i][1], init = 1)],
    panels  = [Panel(xlabel="kinetic energy [GeV]", ylabel="β = v/c  (∝ revolution frequency)",
                     title="velocity", xscale=:log10, xlim=(0.01, 1.0e4), ylim=(0.0, 1.05),
                     height=260, legend=:bottomright),
               Panel(xlabel="kinetic energy [GeV]", ylabel="Bρ [T·m]", title="magnetic rigidity",
                     xscale=:log10, yscale=:log10, xlim=(0.01, 1.0e4), ylim=(0.1, 5.0e4),
                     height=260, legend=:bottomright)],
    statics = statics,
    note = "Each ring covers a factor of 10–20 in Bρ. Only the first one sees a large change "*
           "in β: the PSB's RF has to sweep by almost a factor of three, while from the SPS "*
           "upward the revolution frequency hardly moves.",
) do i
    name, Ki, Kt, C = RINGS[i]
    col = RING_COLORS[i]
    Ks = Kgrid(Ki, Kt)
    bi, bt = proton(Ki).beta, proton(Kt).beta
    (series = [line(Ks ./ 1e9, [proton(K).beta for K in Ks]; panel=1, color=col, width=4.0),
               points([Ki, Kt] ./ 1e9, [bi, bt]; panel=1, color=col, size=5.0),
               line(Ks ./ 1e9, rigidity.(Ks); panel=2, color=col, width=4.0),
               points([Ki, Kt] ./ 1e9, rigidity.([Ki, Kt]); panel=2, color=col, size=5.0)],
     readouts = ["kinetic energy" => string(energy_label(Ki), " → ", energy_label(Kt)),
                 "β" => string(round(bi; sigdigits=4), " → ", round(bt; sigdigits=4)),
                 "revolution frequency" => string(freq_label(bi*CLIGHT/C), " → ", freq_label(bt*CLIGHT/C)),
                 "RF frequency swing" => string(round(100*(bt/bi - 1); sigdigits=3), " %"),
                 "Bρ" => string(round(rigidity(Ki); sigdigits=3), " → ", round(rigidity(Kt); sigdigits=3),
                                " T·m  (×", round(rigidity(Kt)/rigidity(Ki); sigdigits=3), ")")])
end
```
