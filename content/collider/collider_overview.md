---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Colliders and Center-of-Mass Energy

## Why collide two beams

A **collider** brings two beams of particles into head-on collision, instead of
sending one beam onto a stationary target. The aim is the same as in every
high-energy experiment. The more energy a collision makes available, the heavier
the particles it can create and the smaller the structures it can resolve. Over
a century this has taken physics from atoms to nuclei, then to quarks, and to
the particles of the Standard Model.

Colliders can be linear accelerators or rings. This part of the notes deals with
**ring colliders**. These are synchrotrons that store two beams for hours and
collide the same bunches millions of times. Three numbers characterize such a
machine, and each has a section of its own:

- the **center-of-mass energy** $E_{cm}$, which decides what can be produced
  (this page);
- the **luminosity** $\mathcal L$, which decides how often it is produced
  ([](luminosity.md));
- the **beam–beam parameter** $\xi$, which measures how strongly the two beams
  disturb each other and in practice caps the luminosity
  ([](beam_beam.md)).

The optics that squeeze the beams at the collision point are the subject of
[](interaction_region.md).

## The first colliders

The first collider was AdA (*Anello di Accumulazione*), built at Frascati in
1961 on Bruno Touschek's proposal. Electrons and positrons have opposite charge,
so they circulate in opposite directions in the same magnetic field. AdA
therefore stored both beams in a single ring of magnets, at up to 250 MeV each.
Moved to Orsay to use a stronger injector, it recorded the first
electron–positron collisions.

```{figure} ../images/synchrotron/ada.jpg
:width: 55%
:name: fig-ada-collider
AdA, the first electron–positron storage ring and the first collider (Frascati, 1961).
```

Protons of equal charge need two rings with opposite fields. CERN's
**Intersecting Storage Rings** (ISR, 1971–1984) were two interleaved rings of
about 300 m diameter. They crossed at eight points and collided protons of up to
31.4 GeV per beam. The ISR was the first hadron collider, and the place where
much of the beam physics of colliders was first learned.

```{figure} ../images/collider/isr.jpg
:width: 50%
:name: fig-isr
The Intersecting Storage Rings at CERN, the first hadron collider.
```

## Colliders today

The largest collider is the **LHC** at CERN (26.7 km, protons at 6.8 TeV per
beam, also lead ions), fed by the injector chain shown in
[](../synchrotron/synchrotron_principle.md).

**RHIC** at Brookhaven is 3.8 km around. It collided heavy ions (up to
100 GeV per nucleon) and polarized protons (up to 255 GeV), with the STAR and
PHENIX experiments. Its hadron ring is becoming part of the **Electron–Ion
Collider** (EIC), which adds an electron ring in the same tunnel.

**SuperKEKB** at KEK collides 7 GeV electrons with 4 GeV positrons for the
Belle II experiment, at the highest luminosity yet reached.

```{figure} ../images/collider/colliders_today.jpg
:width: 100%
:name: fig-colliders-today
Left: RHIC and its injectors at Brookhaven, with the experiments PHENIX and STAR. Right: the SuperKEKB electron–positron collider at KEK, with its injector linac and positron damping ring.
```

## Fixed target or collider?

Particle 1 has energy $E_1$, mass $m_1$ and Lorentz factor
$\gamma_1=E_1/(m_1c^2)$. Particle 2 has $E_2$, $m_2$ and $\gamma_2$. The energy
available for producing new particles is the energy in the center-of-mass frame.
That energy is a Lorentz invariant, so it can be evaluated in the laboratory
(see [](../math_prep/relativistic_mechanics.md)):

$$
E_{cm}^2=\left(E_1+E_2\right)^2-c^2\left(\mathbf P_1+\mathbf P_2\right)^2 .
$$

For a head-on collision, $\mathbf P_1=\hat{\mathbf z}\,\gamma_1 m_1\beta_1 c$ and
$\mathbf P_2=-\hat{\mathbf z}\,\gamma_2 m_2\beta_2 c$. Then

$$
E_{cm}^2=\left(m_1\gamma_1+m_2\gamma_2\right)^2c^4-\left(m_1\gamma_1\beta_1-m_2\gamma_2\beta_2\right)^2c^4 .
$$

Using $\gamma^2-\gamma^2\beta^2=1$ for each particle gives

$$
E_{cm}=c^2\sqrt{m_1^2+m_2^2+2m_1m_2\gamma_1\gamma_2\left(1+\beta_1\beta_2\right)}
\;\approx\;2\sqrt{E_1E_2}\qquad(\gamma_1,\gamma_2\gg1).
$$

For a **fixed target**, $\gamma_2=1$ and $\beta_2=0$:

$$
E_{cm}=c^2\sqrt{m_1^2+m_2^2+2m_1m_2\gamma_1}\;\approx\;\sqrt{2E_1\,m_2c^2}.
$$

A collider with equal beams turns all of its energy into $E_{cm}=2E$. On a
fixed target, $E_{cm}$ grows only as the square root of the beam energy. The
rest of the energy goes into the motion of the center of mass, which carries the
collision products forward. A car hitting a parked car does far less damage
than two cars meeting head-on at the same speed.

```{code-cell} julia
:tags: [hide-input]
:label: fig-ecm
:caption: Center-of-mass energy of proton–proton collisions against the energy of one beam, for a collider (two equal beams) and for a fixed hydrogen target. Matching the LHC with a fixed target would need a beam of about 10¹⁷ eV, the energy of rare cosmic rays.

using CairoMakie, TrackPadWidgets

const MP = 0.93827      # proton rest energy [GeV]

"E_cm [GeV] for two beams of total energy E1, E2 [GeV] colliding head-on (E2 = MP: fixed target)."
function ecm(E1, E2; m1=MP, m2=MP)
    p1, p2 = sqrt(E1^2 - m1^2), sqrt(max(E2^2 - m2^2, 0.0))
    sqrt((E1 + E2)^2 - (p1 - p2)^2)
end

let E = 10.0 .^ range(0, 8, length=400)
    fig = Figure(size=(760, 460), fontsize=15)
    ax = Axis(fig[1, 1]; xscale=log10, yscale=log10, xlabel="beam energy [GeV]",
              ylabel="center-of-mass energy [GeV]", limits=(1, 4e8, 1, 4e4))
    lines!(ax, E, ecm.(E, E); color=PALETTE[1], linewidth=3, label="collider, two equal beams")
    lines!(ax, E, ecm.(E, MP); color=PALETTE[2], linewidth=3, label="fixed target")

    colliders = [("ISR", 31.4), ("RHIC", 250.0), ("Tevatron", 980.0), ("LHC", 7000.0)]
    for (name, Eb) in colliders
        scatter!(ax, [Eb], [ecm(Eb, Eb)]; color=PALETTE[1], markersize=11)
        text!(ax, Eb, ecm(Eb, Eb); text=" " * name, align=(:right, :bottom), offset=(-4, 4))
    end
    for (name, Eb, al) in [("SPS", 450.0, (:right, :bottom)), ("CEBAF", 12.0, (:left, :top))]
        scatter!(ax, [Eb], [ecm(Eb, MP)]; color=PALETTE[2], markersize=11)
        text!(ax, Eb, ecm(Eb, MP); text=name, align=al, offset=al[1] == :right ? (-4, 4) : (4, -4))
    end
    # the fixed-target beam that would match the LHC
    Elhc = ecm(7000.0, 7000.0); Eft = (Elhc^2 - 2MP^2)/(2MP)
    lines!(ax, [7000.0, Eft], [Elhc, Elhc]; color=:gray50, linestyle=:dash)
    scatter!(ax, [Eft], [Elhc]; color=PALETTE[2], marker=:star5, markersize=16)
    text!(ax, Eft, Elhc; text="$(round(Eft/1e8; digits=1))×10⁸ GeV", align=(:right, :bottom),
          offset=(-8, 6), color=:gray30)
    axislegend(ax; position=:lt)
    fig
end
```

The two-cars picture also explains **asymmetric colliders**, whose beams have
different energies. Then the center of mass moves in the laboratory with
$\beta_{cm}\gamma_{cm}=c\,(P_1-P_2)/E_{cm}$. HERA collided 27.5 GeV electrons with
920 GeV protons, for $E_{cm}=2\sqrt{27.5\times920}\approx318$ GeV. SuperKEKB
collides 7 GeV electrons with 4 GeV positrons. That gives
$E_{cm}=2\sqrt{7\times4}=10.58$ GeV, the mass of the $\Upsilon(4S)$ resonance,
which decays into pairs of B mesons. The boost $\beta\gamma\approx0.28$ carries
each B meson about 0.1 mm before it decays. That separation is what lets Belle II
measure the decay times.

## Types of colliders

| | hadron colliders | lepton colliders ($e^+e^-$) | lepton–hadron colliders |
|---|---|---|---|
| particles | composite: quarks and gluons | point-like | an electron probes a proton or nucleus |
| collision energy | spread out: each collision involves one quark or gluon from each proton, carrying a variable share of its energy | known precisely, set by the beams | known for the lepton, spread for the hadron |
| typical use | discovery: W and Z (CERN SppS), top quark (Tevatron), Higgs boson (LHC) | precision and "factories" tuned to a resonance: Z at LEP, $\Upsilon(4S)$ at the B factories | deep inelastic scattering, the inner structure of protons and nuclei: HERA, EIC |

## Ring colliders and their limits

A ring is economical. The RF system gives each particle only a small energy gain
per turn, so the same cavities are used again and again (the synchrotron
principle). Once at top energy, the bunches pass through each other once per
turn, typically for about ten hours.

**Hadron rings are limited by the magnetic field.** The bending radius follows
from $p=eB\rho$, or $p\,[\mathrm{GeV}/c]=0.2998\,B\,[\mathrm T]\;\rho\,[\mathrm m]$.
At 7 TeV and the LHC's 8.3 T, $\rho\approx2.8$ km, so the dipoles alone are
17.7 km long. Doubling the energy means doubling either the field or the
circumference. That is why collider magnets went from iron-dominated
normal-conducting dipoles (below 2 T) to superconducting Nb–Ti (RHIC 3.5 T,
LHC 8.3 T), and why studies of future machines rely on Nb₃Sn at 14–16 T.

**Electron rings are limited by synchrotron radiation.** A particle on a circle
of radius $\rho$ radiates, each turn, an energy

$$
U_0=C_\gamma\frac{E^4}{\rho},\qquad C_\gamma=8.85\times10^{-5}\ \mathrm{m/GeV^3}\ \text{for electrons},
$$

and $C_\gamma$ scales with the particle mass as $1/m^4$. LEP, in the tunnel that
now holds the LHC, reached 104.5 GeV per beam. There each electron lost about
3.4 GeV per turn, and superconducting RF had to restore it on every turn. A
proton loses only a few keV per turn at 7 TeV in the same tunnel.

```{code-cell} julia
:tags: [hide-input]
:label: fig-ring-limits
:caption: Left, the total length of the bending dipoles, 2πρ, needed to hold a proton beam of a given energy, for several dipole fields, with four colliders (the ring itself is longer still, since straight sections and quadrupoles fill a third of the LHC and more than half of RHIC). Right, the energy radiated per turn on a bending radius of 3 km (the LEP/LHC tunnel) by electrons, muons and protons.

let fig = Figure(size=(980, 420), fontsize=15)
    # --- dipole field sets the size of a hadron ring --------------------------
    ax = Axis(fig[1, 1]; xscale=log10, yscale=log10, xlabel="proton energy [TeV]",
              ylabel="length of dipoles, 2πρ [km]", limits=(0.1, 100, 0.2, 300),
              title="hadron rings: magnetic field")
    E = 10.0 .^ range(-1, 2, length=100)                 # TeV
    for (j, (B, tech)) in enumerate([(1.8, "1.8 T normal conducting"), (3.5, "3.5 T Nb–Ti (RHIC)"),
                                     (8.3, "8.3 T Nb–Ti (LHC)"), (14.0, "14 T Nb₃Sn (FCC-hh)")])
        lines!(ax, E, 2π .* (E .* 1e3 ./ (0.2998B)) ./ 1e3; color=PALETTE[j], linewidth=2.5, label=tech)
    end
    # bending radius of each machine at its top energy
    for (name, E0, ρ) in [("RHIC", 0.255, 242.8), ("Tevatron", 0.98, 754.0), ("LHC", 7.0, 2804.0),
                          ("FCC-hh", 50.0, 11900.0)]
        scatter!(ax, [E0], [2π*ρ/1e3]; color=:black, markersize=10)
        text!(ax, E0, 2π*ρ/1e3; text=name, align=(:left, :top), offset=(5, -3), fontsize=13)
    end
    axislegend(ax; position=:lt, framevisible=false, labelsize=12)

    # --- synchrotron radiation sets the energy of a lepton ring ---------------
    ax = Axis(fig[1, 2]; xscale=log10, yscale=log10, xlabel="beam energy [GeV]",
              ylabel="energy lost per turn [GeV]", limits=(10, 1e4, 1e-9, 1e3),
              title="synchrotron radiation, ρ = 3 km")
    E = 10.0 .^ range(1, 4, length=100)
    me, mμ, mp = 0.000511, 0.10566, 0.93827
    for (j, (name, m)) in enumerate([("electron", me), ("muon", mμ), ("proton", mp)])
        lines!(ax, E, 8.85e-5 .* (me/m)^4 .* E .^ 4 ./ 3000; color=PALETTE[j], linewidth=2.5, label=name)
    end
    scatter!(ax, [104.5], [8.85e-5*104.5^4/3026]; color=:black, markersize=10)
    text!(ax, 104.5, 8.85e-5*104.5^4/3026; text="LEP2", align=(:right, :bottom), offset=(-5, 3), fontsize=13)
    scatter!(ax, [7000.0], [8.85e-5*(me/mp)^4*7000^4/2804]; color=:black, markersize=10)
    text!(ax, 7000.0, 8.85e-5*(me/mp)^4*7000^4/2804; text="LHC", align=(:right, :bottom), offset=(-5, 3), fontsize=13)
    axislegend(ax; position=:lt, framevisible=false)
    fig
end
```

An $e^+e^-$ ring above LEP's energy therefore needs a much larger radius. That
is the route taken by the 91 km FCC-ee and CEPC designs. The alternative is a
**linear collider** (ILC, CLIC), which radiates no synchrotron light but must
reach full energy in a single pass and uses each bunch only once.

## Muon colliders

A muon is a lepton like the electron, but 207 times heavier. At the same energy
and bending radius it radiates $207^4\approx1.8\times10^9$ times less. A muon
ring can therefore reach multi-TeV energies in a compact footprint. The higher
mass also raises the rate of producing a Higgs boson directly
($\mu^+\mu^-\to H$) by a factor $207^2\approx4\times10^4$ over $e^+e^-$.

The difficulties come from how the muons are made and from how long they live.
A proton beam on a target produces pions, which decay into muons and neutrinos.
The muons emerge in a large, hot cloud of phase space. They have to be
*cooled* (by ionization cooling, the only method fast enough), accelerated and
collided before they decay. Their lifetime is $\tau=2.2\ \mu$s at rest and
$\gamma\tau$ in the laboratory. The number of turns a muon survives is
$\gamma\tau$ divided by the revolution time $2\pi\rho/(\beta c)$, with
$\rho=p/(eB)$:

$$
N_{\mathrm{turns}}=\frac{eB\,\tau}{2\pi m_\mu}\approx 300\times B\,[\mathrm T],
$$

independent of energy, where $B$ is the field averaged around the ring. Even with
10 T of average field, a muon makes only about 3000 turns. A hadron collider
fill lasts some $10^8$–$10^9$ turns.

```{code-cell} julia
:tags: [hide-input]
:label: fig-muon-chain
:caption: The stages of a muon collider. Every step after the target has to be completed within a few muon lifetimes.

let fig = Figure(size=(980, 210), fontsize=15)
    ax = Axis(fig[1, 1]; limits=(0, 100, 0, 20)); hidedecorations!(ax); hidespines!(ax)
    stages = [("proton\ndriver", PALETTE[1]), ("target\np → π", PALETTE[2]), ("decay\nπ → μ ν", PALETTE[3]),
              ("ionization\ncooling", PALETTE[4]), ("fast\nacceleration", PALETTE[5]), ("collider\nring", PALETTE[6])]
    w, gap = 13.0, 4.0
    for (i, (label, col)) in enumerate(stages)
        x0 = 1 + (i - 1)*(w + gap)
        poly!(ax, Rect(x0, 5, w, 10); color=(col, 0.18), strokecolor=col, strokewidth=2)
        text!(ax, x0 + w/2, 10; text=label, align=(:center, :center), fontsize=14)
        i < length(stages) && arrows2d!(ax, [x0 + w + 0.3], [10.0], [gap - 0.6], [0.0]; color=:gray30,
                                        shaftwidth=2, tipwidth=9, tiplength=9)
    end
    text!(ax, 1 + 2.5(w + gap), 1.5; text="muon lifetime: 2.2 μs at rest, γ × 2.2 μs in the laboratory",
          align=(:center, :center), fontsize=13, color=:gray30)
    fig
end
```
