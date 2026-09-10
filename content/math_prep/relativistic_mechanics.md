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

using StaticArrays, TrackPad, TrackPadWidgets

mₑ = 0.510998950e6          # electron rest mass [eV]
γcurve = range(1.0001, 60.0, length=200)
βcurve = @. sqrt(1 - 1/γcurve^2)
Ks = 10.0 .^ range(4.0, 11.0, length=36)

explorer(
    title   = "TrackPad's Beam(K): γ and β are locked to the same curve",
    sliders = [Knob("kinetic energy K", eachindex(Ks);
                    fmt = i -> let K = Ks[i]
                        string(round(K/10.0^floor(log10(K)); digits=1), "e", Int(floor(log10(K))), " eV")
                    end, init = 20)],
    panels  = [Panel(xlabel="β", ylabel="γ", title="γ = 1/√(1−β²)",
                     xlim=(0.0,1.02), ylim=(0.0,20.0), height=300, legend=:bottomleft)],
    statics = [line(βcurve, γcurve; color="#9aa4b2", label="analytic γ(β)")],
    note = "The marker is TrackPad's own (β₀, γ₀) for Beam(K); it lands on the analytic curve "*
           "by construction — γ = 1+K/mc² and β = √(1−1/γ²) always agree.",
) do i
    K = Ks[i]
    beam = Beam(K; mass=mₑ)
    E0 = beam.energy + beam.mass
    P0c = beam.beta*E0
    (series = [points([beam.beta], [min(beam.gamma,20.0)]; color=PALETTE[1], size=7.0, label="Beam(K)")],
     readouts = ["γ" => round(beam.gamma; sigdigits=6),
                 "1+K/mc²" => round(1 + beam.energy/mₑ; sigdigits=6),
                 "β" => round(beam.beta; sigdigits=6),
                 "√(1−1/γ²)" => round(sqrt(1-1/beam.gamma^2); sigdigits=6),
                 "E₀ [GeV]" => round(E0/1e9; sigdigits=5),
                 "P₀c [GeV]" => round(P0c/1e9; sigdigits=5)])
end
```

The four-vector invariant for the energy-momentum vector, $E^2-P^2c^2=(mc^2)^2$,
is preserved under boosts: a boost slides a particle's point along the
$(P_zc,E)$ mass-shell hyperbola but never off it.

```{code-cell} julia
:tags: [hide-input]

using StaticArrays, LinearAlgebra, TrackPadWidgets

g(β) = 1/sqrt(1-β^2)
Λ(β) = @SMatrix [1 0 0 0;
                 0 1 0 0;
                 0 0 g(β) -g(β)*β;
                 0 0 -g(β)*β g(β)]        # boost along z; order (x,y,z,ct)
ηmetric = Diagonal([-1.0,-1.0,-1.0,1.0])
P0 = @SVector [0.0, 0.0, 2.5, 3.0]        # (Px,Py,Pz,E)/c, GeV-ish units
m2 = dot(P0, ηmetric*P0)
Pzc = range(-6.0, 6.0, length=161)
Ecurve = @. sqrt(Pzc^2 + m2)

explorer(
    title   = "The energy-momentum invariant survives a boost",
    sliders = [Knob("boost β", range(-0.95, 0.95, length=39);
                    fmt = b -> string(round(b; digits=2)), init = 20)],
    panels  = [Panel(xlabel="Pz c [GeV]", ylabel="E [GeV]", title="mass-shell E² − (Pzc)² = (mc²)²",
                     xlim=(-6.2,6.2), ylim=(0.0,7.0), height=300, legend=:bottomleft)],
    statics = [line(Pzc, Ecurve; color="#9aa4b2", label="mass shell"),
               points([P0[3]], [P0[4]]; color=PALETTE[3], size=6.0, label="P (rest frame value)")],
    note = "Boosting slides the point along the SAME hyperbola — that is exactly the "*
           "statement P·P is Lorentz-invariant.",
) do β
    Pb = Λ(β)*P0
    m2b = dot(Pb, ηmetric*Pb)
    (series = [points([Pb[3]], [Pb[4]]; color=PALETTE[1], size=7.0, label="Λ(β)·P")],
     readouts = ["β" => round(β; digits=3),
                 "Pz′c [GeV]" => round(Pb[3]; sigdigits=5),
                 "E′ [GeV]" => round(Pb[4]; sigdigits=5),
                 "P·P" => round(m2; sigdigits=6),
                 "(ΛP)·(ΛP)" => round(m2b; sigdigits=6)])
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

Ks   = 10.0 .^ range(4, 13, length=46)               # 10 keV … 10 TeV
SPEC = [("electron", TrackPad.M_ELECTRON, -1.0), ("proton", TrackPad.M_PROTON, 1.0)]
klabel(K) = string(round(K/10.0^floor(log10(K)); digits=1), "e", Int(floor(log10(K))), " eV")

background = Any[]
for (i, (nm, m, q)) in enumerate(SPEC)
    push!(background, line(Ks, [Beam(K; mass=m, charge=q).beta for K in Ks];
                           color=PALETTE[i], label=nm))
end

explorer(
    title   = "Velocity saturates long before the energy does",
    sliders = [Knob("kinetic energy", eachindex(Ks); fmt = i -> klabel(Ks[i]), init = 30),
               Knob("particle", 1:2; fmt = j -> SPEC[j][1])],
    panels  = [Panel(xlabel="kinetic energy K [eV]", ylabel="β = v/c",
                     xscale=:log10, xlim=(8.0e3, 1.3e13), ylim=(0.0, 1.05),
                     height=300, legend=:bottomright)],
    statics = background,
    note = "The horizontal axis is logarithmic in K itself — nine decades of energy, and "*
           "each species climbs the same curve, displaced by the ratio of the rest masses.",
) do i, j
    nm, m, q = SPEC[j]
    b = Beam(Ks[i]; mass=m, charge=q)
    (series = [points([Ks[i]], [b.beta]; color=PALETTE[j], size=6.5, label=nm)],
     readouts = ["particle" => nm,
                 "K"     => klabel(Ks[i]),
                 "γ"     => round(b.gamma; sigdigits=5),
                 "β"     => round(b.beta;  sigdigits=6),
                 "1 − β" => string(round(1 - b.beta; sigdigits=3)),
                 "P₀c"   => string(round(b.beta*(b.energy+b.mass)/1e9; sigdigits=4), " GeV")])
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
