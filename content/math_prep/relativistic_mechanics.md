---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Relativistic Mechanics

When an object's speed is comparable to the speed of light, relativistic
mechanics is required. The relativistic factors are defined as

$$
\boldsymbol{\beta}=\frac{\mathbf{v}}{c},\qquad
\gamma=\frac{1}{\sqrt{1-\beta^2}}
$$

The momentum of the particle is no longer $\mathbf{P}=m\mathbf{v}$ but

$$
\mathbf{P}=\gamma m\mathbf{v},
$$

which can be viewed as the object being "heavier" by a factor of $\gamma$. The
energy of the object reads

$$
E=\sqrt{P^2c^2+(mc^2)^2}=\gamma mc^2
$$

## Lorentz transformation

```{figure} ../images/lorentz.png
:width: 200px
:name: fig:lorentz
Two inertial frames $K$ and $K'$ in relative motion along $z$.
```

If the coordinate system $K'$ moves away from $K$ with velocity
$\boldsymbol{\beta}c$ directed along $\hat z$, the coordinates of both systems
are linked by

$$
x'=x,\qquad y'=y,\qquad
z'=\gamma(z-\beta ct),\qquad
t'=\gamma\left(t-\frac{\beta z}{c}\right)
$$

or, in matrix form,

$$
\begin{pmatrix}x'\\y'\\z'\\ct'\end{pmatrix}
=
\begin{pmatrix}
1 & 0 & 0 & 0\\
0 & 1 & 0 & 0\\
0 & 0 & \gamma & -\gamma\beta\\
0 & 0 & -\gamma\beta & \gamma
\end{pmatrix}
\begin{pmatrix}x\\y\\z\\ct\end{pmatrix}
$$

## Time dilation

Two events happen a time $\Delta t$ apart at the same location in frame $K$.
For an observer in the boosted frame $K'$, the time difference becomes

$$
\Delta t'=\gamma \Delta t
$$

## Length contraction

If an object is at rest in frame $K$ with length $L$ along $z$, an observer in
$K'$ measures both ends simultaneously ($\Delta t'=0$). Using the inverse
Lorentz transformation,

$$
L=\Delta z=\gamma \Delta z' + \gamma\beta c\,\Delta t'
\quad\Rightarrow\quad
\Delta z' = \frac{L}{\gamma}
$$

## 4-Vectors and Lorentz invariants

The Lorentz transformation links a space-like variable and a time-like
variable; grouping them into 4-vectors is very useful. Common examples:

- time–space four-vector: $(ct, \vec{r})=(ct,(x,y,z))$
- energy–momentum four-vector: $(E/c, \vec{P})=(E/c,(P_x,P_y,P_z))$
- electromagnetic four-potential: $(\phi/c,\vec A)=(\phi/c,(A_x,A_y,A_z))$

The dot product of two four-vectors $A=(a_0,a_1,a_2,a_3)$ and
$B=(b_0,b_1,b_2,b_3)$ uses the metric signature $(+,-,-,-)$:

$$
A\cdot B=
\begin{pmatrix}a_{0}\\a_{1}\\a_{2}\\a_{3}\end{pmatrix}^{T}
\begin{pmatrix}
1 & 0 & 0 & 0\\
0 & -1 & 0 & 0\\
0 & 0 & -1 & 0\\
0 & 0 & 0 & -1
\end{pmatrix}
\begin{pmatrix}b_{0}\\b_{1}\\b_{2}\\b_{3}\end{pmatrix}
$$

This dot product is invariant under Lorentz transformations.

## Numerical examples in Julia

Relativistic kinematics appears everywhere in accelerator physics — for
instance inside TrackPad's `Beam`, which stores kinetic energy $K_0$, rest mass
$m_0c^2$, and derives $\beta_0$ and $\gamma_0$. However you compute them, $\beta$
and $\gamma$ are locked to the same curve $\gamma=1/\sqrt{1-\beta^2}$; slide the
kinetic energy and watch TrackPad's own $(\beta_0,\gamma_0)$ ride along it:

```{code-cell} julia
:tags: [hide-input]

using TrackPad, TrackPadWidgets

const M_E = 0.510998950e6            # electron rest energy [eV]

# The analytic curve γ(β), drawn once as the backdrop.
γ_curve = range(1.0001, 60.0, length=200)
β_curve = @. sqrt(1 - 1/γ_curve^2)

# Kinetic energies to step through: 10 keV … 100 GeV, evenly spaced in decades.
kinetic_energies = 10.0 .^ range(4.0, 11.0, length=36)

"Format an energy as e.g. \"3.2e7 eV\" — compact enough for a slider label."
function energy_label(K)
    decade  = floor(Int, log10(K))
    mantissa = round(K / 10.0^decade; digits=1)
    return string(mantissa, "e", decade, " eV")
end

explorer(
    title   = "TrackPad's Beam(K): γ and β are locked to the same curve",
    sliders = [Knob("kinetic energy K", eachindex(kinetic_energies);
                    fmt = i -> energy_label(kinetic_energies[i]), init = 20)],
    panels  = [Panel(xlabel="β", ylabel="γ", title="γ = 1/√(1−β²)",
                     xlim=(0.0,1.02), ylim=(0.0,20.0), height=300, legend=:bottomleft)],
    statics = [line(β_curve, γ_curve; color="#9aa4b2", label="analytic γ(β)")],
    note = "The marker is TrackPad's own (β₀, γ₀) for Beam(K); it lands on the analytic curve "*
           "by construction — γ = 1+K/mc² and β = √(1−1/γ²) always agree.",
) do i
    K    = kinetic_energies[i]
    beam = Beam(K; mass=M_E)

    E_total   = beam.energy + beam.mass      # beam.energy is KINETIC energy
    momentum  = beam.beta * E_total          # P₀c, in eV

    # γ runs off the top of the panel at high K, so the marker is clamped to the
    # axis; the readouts still report the true value.
    γ_plotted = min(beam.gamma, 20.0)

    (series = [points([beam.beta], [γ_plotted]; color=PALETTE[1], size=7.0,
                      label="Beam(K)")],
     readouts = ["γ"          => round(beam.gamma; sigdigits=6),
                 "1+K/mc²"    => round(1 + beam.energy/M_E; sigdigits=6),
                 "β"          => round(beam.beta; sigdigits=6),
                 "√(1−1/γ²)"  => round(sqrt(1 - 1/beam.gamma^2); sigdigits=6),
                 "E₀ [GeV]"   => round(E_total/1e9; sigdigits=5),
                 "P₀c [GeV]"  => round(momentum/1e9; sigdigits=5)])
end
```

The four-vector invariant for the energy-momentum vector, $E^2-P^2c^2=(mc^2)^2$,
is preserved under boosts: a boost slides a particle's point along the
$(P_zc,E)$ mass-shell hyperbola but never off it.

```{code-cell} julia
:tags: [hide-input]

using StaticArrays, LinearAlgebra, TrackPadWidgets

# Four-vectors here are ordered (x, y, z, ct), so a boost along z mixes only the
# last two components. The metric has signature (−,−,−,+) to match that ordering.
lorentz_gamma(β) = 1/sqrt(1 - β^2)

function boost_along_z(β)
    γ = lorentz_gamma(β)
    return @SMatrix [1  0   0     0
                     0  1   0     0
                     0  0   γ    -γ*β
                     0  0  -γ*β   γ]
end

const METRIC = Diagonal([-1.0, -1.0, -1.0, 1.0])
minkowski_dot(A, B) = dot(A, METRIC * B)

# A test particle, as (Pₓ, P_y, P_z, E)/c in GeV.
P_lab = @SVector [0.0, 0.0, 2.5, 3.0]
mass_squared = minkowski_dot(P_lab, P_lab)          # (mc²)², the invariant

# The mass shell E² − (P_zc)² = (mc²)², i.e. every state this particle could be
# boosted into.
Pz_curve = range(-6.0, 6.0, length=161)
E_curve  = @. sqrt(Pz_curve^2 + mass_squared)

explorer(
    title   = "The energy-momentum invariant survives a boost",
    sliders = [Knob("boost β", range(-0.95, 0.95, length=39);
                    fmt = b -> string(round(b; digits=2)), init = 20)],
    panels  = [Panel(xlabel="Pz c [GeV]", ylabel="E [GeV]",
                     title="mass-shell E² − (Pzc)² = (mc²)²",
                     xlim=(-6.2,6.2), ylim=(0.0,7.0), height=300, legend=:bottomleft)],
    statics = [line(Pz_curve, E_curve; color="#9aa4b2", label="mass shell"),
               points([P_lab[3]], [P_lab[4]]; color=PALETTE[3], size=6.0,
                      label="P (rest frame value)")],
    note = "Boosting slides the point along the SAME hyperbola — that is exactly the "*
           "statement P·P is Lorentz-invariant.",
) do β
    P_boosted = boost_along_z(β) * P_lab

    (series = [points([P_boosted[3]], [P_boosted[4]]; color=PALETTE[1], size=7.0,
                      label="Λ(β)·P")],
     readouts = ["β"            => round(β; digits=3),
                 "Pz′c [GeV]"   => round(P_boosted[3]; sigdigits=5),
                 "E′ [GeV]"     => round(P_boosted[4]; sigdigits=5),
                 "P·P"          => round(mass_squared; sigdigits=6),
                 "(ΛP)·(ΛP)"    => round(minkowski_dot(P_boosted, P_boosted); sigdigits=6)])
end
```

## How relativistic is my beam?

One number decides more of the design of an accelerator than any other: the
velocity $\beta$, which fixes how long a particle spends inside each element and
therefore how the accelerating structures must be spaced. Its behaviour is
entirely set by the ratio of kinetic energy to rest mass, so the *same* curve
serves every species — only the position along it changes. Slide through nine
decades of kinetic energy for an electron and a proton:

```{code-cell} julia
:tags: [hide-input]

using TrackPad, TrackPadWidgets
# `energy_label` was defined in the first code cell of this page.

# 10 keV … 10 TeV, evenly spaced in decades.
kinetic_scan = 10.0 .^ range(4, 13, length=46)

SPECIES = [(name="electron", mass=TrackPad.M_ELECTRON, charge=-1.0),
           (name="proton",   mass=TrackPad.M_PROTON,   charge=+1.0)]

# One β(K) curve per species, drawn behind every frame.
beta_curves = [line(kinetic_scan,
                    [Beam(K; mass=s.mass, charge=s.charge).beta for K in kinetic_scan];
                    color=PALETTE[i], label=s.name)
               for (i, s) in enumerate(SPECIES)]

explorer(
    title   = "Velocity saturates long before the energy does",
    sliders = [Knob("kinetic energy", eachindex(kinetic_scan);
                    fmt = i -> energy_label(kinetic_scan[i]), init = 30),
               Knob("particle", eachindex(SPECIES); fmt = j -> SPECIES[j].name)],
    panels  = [Panel(xlabel="kinetic energy K [eV]", ylabel="β = v/c",
                     xscale=:log10, xlim=(8.0e3, 1.3e13), ylim=(0.0, 1.05),
                     height=300, legend=:bottomright)],
    statics = beta_curves,
    note = "The horizontal axis is logarithmic in K itself — nine decades of energy, and "*
           "each species climbs the same curve, displaced by the ratio of the rest masses.",
) do i, j
    species = SPECIES[j]
    K       = kinetic_scan[i]
    beam    = Beam(K; mass=species.mass, charge=species.charge)

    (series = [points([K], [beam.beta]; color=PALETTE[j], size=6.5, label=species.name)],
     readouts = ["particle" => species.name,
                 "K"     => energy_label(K),
                 "γ"     => round(beam.gamma; sigdigits=5),
                 "β"     => round(beam.beta;  sigdigits=6),
                 "1 − β" => string(round(1 - beam.beta; sigdigits=3)),
                 "P₀c"   => string(round(beam.beta*(beam.energy + beam.mass)/1e9;
                                         sigdigits=4), " GeV")])
end
```

At 10 MeV an electron is already at $\beta = 0.9987$ — a linac designed for it
can use identical cells from there on, because the particle never gets
appreciably faster. A proton reaches the same $\beta$ only near 10 GeV, which is
why proton linacs must stretch their cell lengths as the beam accelerates.

The second number that governs every design is the magnetic **rigidity**
$B\rho = P/q$, which behaves in exactly the opposite way: it never saturates. We
meet it in the Transverse Dynamics chapter, where it becomes the natural
normalization of the whole subject.
