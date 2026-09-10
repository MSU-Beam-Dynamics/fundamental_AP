---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Transverse Dynamics: Fields and Magnets

Transverse dynamics concerns transporting charged particles to a specific
location (a **transport line**) or moving them with periodic boundary
conditions (an **accelerator ring**).

```{figure} ../images/FRIB_FS2.png
:width: 400px
:name: fig:frib
Transport line example: the FRIB folding segment 2.
```

```{figure} ../images/ags_ring.png
:width: 400px
:name: fig:ags
Ring example: the Alternating Gradient Synchrotron (AGS).
```

To affect the motion of a charged particle we need an external magnetic or
electric field acting through the Lorentz force:

$$
\mathbf{F}=q\left(\mathbf{E}+\mathbf{v}\times\mathbf{B}\right)
$$

## Magnetic or electric fields?

For equal Lorentz force, the required magnetic and electric field strengths are
related by the particle velocity. The table lists the electric field needed to
match the force of a 1 T magnetic field:

| $\gamma-1=\frac{E_k}{E_0}$ | $\beta c$ | $\mathbf{B}$ (T) | $\mathbf{E}$ (MV/m) |
|------|------|------|------|
| $10^{-5}$ | $0.0045\,c$ | 1 | 1.34 |
| $10^{-4}$ | $0.0141\,c$ | 1 | 4.24 |
| $10^{-3}$ | $0.0447\,c$ | 1 | 13.4 |
| $10^{-2}$ | $0.140\,c$ | 1 | 42.1 |
| $10^{-1}$ | $0.417\,c$ | 1 | 125 |

Considering the electric break-down limit of materials, only non-relativistic
particles can be effectively manipulated by electric fields. Magnetic fields
are therefore used to bend/steer/focus relativistic particles.

## Design orbit and coordinate system

The geometry of the accelerator is determined by dipole magnets. Assume a
magnetic field constant in some region and along the vertical direction,
$\mathbf{B}=B_0\hat{y}$, with the particle initially traveling in the $\hat s$
direction in the $x$-$s$ plane:

```{figure} ../images/coordinate_system.png
:width: 400px
:name: fig:frenet
Design orbit and the Frenet–Serret coordinate system.
```

Note that the coordinate system moves along the ideal trajectory: $\hat{s}$ is
always tangent to it and $\hat{x}$ points along the normal direction. This is
the **Frenet–Serret coordinate system**.

The equation of motion is

$$
\frac{d\mathbf{P}}{dt}=q\mathbf{v}\times\mathbf{B}
$$

The momentum change lies purely in the $x$-$s$ plane; the velocity stays
perpendicular to $\mathbf B$, so the energy does not change:

$$
m\gamma \frac{v^2}{\rho}=qvB_0
$$

giving the bending radius

$$
\rho=\frac{\left|\mathbf{P}\right|}{qB_0}
$$

We define the **rigidity** $B_0\rho=P/q$ to quantify "how hard to bend the
beam". It depends only on the particle momentum amplitude and is the natural
normalization of transverse motion — after normalization most calculations
become independent of particle energy. The dipoles defining the geometry are
sometimes called *main dipoles*: a charged particle with the correct energy,
injected at the correct position and angle, exactly follows the ideal
trajectory.

In the following examples, six cases of different particles with two energies are 
considered. Energies are quoted **per nucleon**, as is conventional for ions, 
and the gold ion is taken fully stripped ($^{197}\mathrm{Au}^{79+}$), so $A=197$, 
$Z=79$. For a charge state $Z$ the rigidity is

$$
B\rho=\frac{P}{q}=\frac{\gamma\beta\,A\,m_{u}c^{2}}{Z\,e\,c},
$$

with $m_{u}c^{2}$ the rest energy per nucleon. The factor $A/Z$ is what makes
ions different: it is 1 for a proton (and, trivially, for an electron), so those
two need no such correction, while for stripped gold it is 2.494.

| $K$ / nucleon | Particle | $q/e$ | $A/Z$ | $mc^2/A$ | $\gamma$ | $B\rho$ [T·m] |
|:---|:---|--:|--:|--:|--:|--:|
| 1 TeV | proton | 1 | 1 | 938.272 MeV | 1066.8 | 3339 |
| 1 TeV | $^{197}$Au$^{79+}$ | 79 | 2.494 | 931.131 MeV | 1075.0 | 8326 |
| 1 TeV | electron | -1 |  | 0.511 MeV | $1.957\times10^{6}$ | 3336 |
| 100 MeV | proton | 1 | 1 | 938.272 MeV | 1.1066 | 1.483 |
| 100 MeV | $^{197}$Au$^{79+}$ | 79 | 2.494 | 931.131 MeV | 1.1074 | 3.685 |
| 100 MeV | electron | -1 |  | 0.511 MeV | 196.70 | 0.3353 |

Three things are worth reading off it. **Equal energy per nucleon does not mean
equal rigidity**: gold needs exactly $A/Z = 2.494$ times the proton's $B\rho$ in
both rows. In a acclerator that is limited by magnet technology, that factor is 
why a ring built for 250 GeV protons stops out near 100 GeV/nucleon for gold.


Although the rigidity is linearly dependent on momentum, at high energy, the calculation
can be approximated by using beam energy. In the following figure, it shows the rigidity 
as a function of beam kinetic energy for different particles. 

```{code-cell} julia
:tags: [hide-input]

using TrackPad, TrackPadWidgets

Brho(b) = b.beta * (b.energy + b.mass) / (abs(b.charge) * 2.99792458e8)  # [T·m]
Ks   = 10.0 .^ range(4, 13, length=46)                  # 10 keV … 10 TeV
SPEC = [("electron", TrackPad.M_ELECTRON, -1.0), ("proton", TrackPad.M_PROTON, 1.0)]
klabel(K) = string(round(K/10.0^floor(log10(K)); digits=1), "e", Int(floor(log10(K))), " eV")

background = Any[]
for (i, (nm, m, q)) in enumerate(SPEC)
    push!(background, line(Ks, [Brho(Beam(K; mass=m, charge=q)) for K in Ks];
                           color=PALETTE[i], label=nm))
end

explorer(
    title   = "Rigidity does not saturate",
    sliders = [Knob("kinetic energy", eachindex(Ks); fmt = i -> klabel(Ks[i]), init = 30),
               Knob("particle", 1:2; fmt = j -> SPEC[j][1])],
    panels  = [Panel(xlabel="kinetic energy K [eV]", ylabel="Bρ [T·m]",
                     xscale=:log10, yscale=:log10, xlim=(8.0e3, 1.3e13),
                     height=300, legend=:bottomright)],
    statics = background,
    note = "Bρ = P/q is the natural normalisation of transverse motion: TrackPad's k₁ is "*
           "G/(Bρ). At low energy the two species differ by their rest masses; once both "*
           "are ultra-relativistic the curves merge, because Bρ then depends only on energy.",
) do i, j
    nm, m, q = SPEC[j]
    b = Beam(Ks[i]; mass=m, charge=q)
    (series = [points([Ks[i]], [Brho(b)]; color=PALETTE[j], size=6.5, label=nm)],
     readouts = ["particle" => nm,
                 "K"   => klabel(Ks[i]),
                 "β"   => round(b.beta; sigdigits=6),
                 "P₀c" => string(round(b.beta*(b.energy+b.mass)/1e9; sigdigits=4), " GeV"),
                 "Bρ"  => string(round(Brho(b); sigdigits=4), " T·m"),
                 "B for ρ = 10 m" => string(round(Brho(b)/10; sigdigits=4), " T")])
end
```

This one plot explains a lot of accelerator engineering. A 3 GeV electron needs
only $B\rho \approx 10\,\mathrm{T\,m}$, so a light source can bend its beam with
modest 1–1.5 T iron dipoles. A 7 TeV proton needs
$B\rho \approx 2.3\times10^{4}\,\mathrm{T\,m}$, and since superconducting NbTi
dipoles stop at about 8.3 T the only remaining free parameter is the bending
radius — which is why the LHC is 27 km around.

## Magnetic multipoles

Accelerators use many magnet types, grouped by the order of the multipole
expansion around the center of the beam pipe:

| Magnet type | Usage |
|:-----------|-------|
| Dipoles | Form the geometry of the accelerator |
| Quadrupoles | Focus/defocus the beam |
| Sextupoles | Chromatic correction, nonlinear effects |
| Octupoles | Damping of collective effects |
| Solenoids | Focusing for low-energy beams |
| Correctors | Steer the beam |

In the vacuum chamber the region is source-free (no charge or current), and
the field can be expanded in the **Beth representation**:

$$
B_y+iB_x = B_0 \sum_{n=0}^\infty \left(b_n + i a_n\right) \left(x+iy\right)^{n}
\quad \text{(U.S. convention)}
$$

Here $B_0$ is the main dipole field strength. The coefficients $b_n$ and
$a_n$ are the $2(n+1)^{\text{th}}$ multipole coefficients; the $b$ set are the
*normal* components, the $a$ set the *skew* components.

To reach high fields, accelerator magnets use high-permeability iron to boost
the flux density for the same coil current,

$$
\mathbf{B}=\mu_r \mu_0 \mathbf{H}.
$$

```{figure} ../images/Permeability_of_ferromagnet_by_Zureks.svg
:name: fig:permeability
Permeability $\mu_r$ of ferromagnetic material versus applied field $H$.
```

Iron (99.8 % pure) has an initial $\mu_r\approx150$ reaching $5000$;
99.95 % pure iron starts near $10000$ and peaks around $200000$. At such high
permeability, magnetic field lines meet the iron surface perpendicularly —
the shape of the pole faces then dictates the field pattern inside the gap,
which is why each multipole has its characteristic pole geometry.

### Dipole

Keeping only the lowest-order coefficient ($b_0=1$, $a_0=0$ for main dipoles):

$$
B_y=B_0 b_0, \qquad B_x=B_0 a_0 .
$$

```{figure} ../images/H_type_and_C_Type_dipole.png
:width: 400px
:name: fig:hctype
H-type and C-type dipoles.
```

Since dipoles bend the ideal trajectory, their edges also affect the motion.
Useful types include sector dipoles (zero edge angles on both sides) and
rectangular dipoles (both edge angles equal to half the bending angle).

```{figure} ../images/dipoles_3.png
:width: 600px
:name: fig:dipolestypes
General, sector and rectangular dipole geometries.
```

### Quadrupole

A quadrupole has $b_1$ (normal) and/or $a_1$ (skew) nonzero. The normal
quadrupole field is

$$
\mathbf{B}=B_0 b_1\left(y\,\hat{x}+x\,\hat{y}\right),
$$

and the skew quadrupole field is

$$
\mathbf{B}=B_0 a_1\left(x\,\hat{x}-y\,\hat{y}\right).
$$

```{code-cell} julia
:tags: [hide-input]
# Field-line patterns of normal and skew quadrupoles: B ∝ ∇Ψ with Ψ = B₀b₁xy
using CairoMakie

xs = range(-5, 5, length=41)
X = [x for x in xs, _ in xs]
Y = [y for _ in xs, y in xs]

sel(a) = vec(a[1:3:end, 1:3:end])
fig = Figure(size=(800, 380))
ax1 = Axis(fig[1, 1]; aspect=DataAspect(), title="Normal quad")
contour!(ax1, xs, xs, X .* Y; levels=[-6,-4,-2,2,4,6], colormap = :coolwarm, linewidth=2)
arrows2d!(ax1, sel(X), sel(Y), sel(Y), sel(X); color=:gray, lengthscale=0.08)

ax2 = Axis(fig[1, 2]; aspect=DataAspect(), title="Skew quad")
contour!(ax2, xs, xs, (X.^2 .- Y.^2)./2; levels=[-6,-4,-2,2,4,6], colormap = :coolwarm, linewidth=2)
arrows2d!(ax2, sel(X), sel(Y), sel(X), -sel(Y); color=:gray, lengthscale=0.08)
fig
```

### Sextupole

A sextupole has $b_2$ (normal) and/or $a_2$ (skew) nonzero:

$$
\mathbf{B}=B_0 b_2\left[2xy\,\hat{x}+\left(x^2-y^2\right)\hat{y}\right],
\qquad
\mathbf{B}=B_0 a_2\left[\left(x^2-y^2\right)\hat{x}- 2xy\,\hat{y}\right].
$$

```{code-cell} julia
:tags: [hide-input]

using CairoMakie

xs = range(-5, 5, length=41)
X = [x for x in xs, _ in xs]
Y = [y for _ in xs, y in xs]

sel(a) = vec(a[1:3:end, 1:3:end])

# Normal sextupole: Ψ ∝ (3x²y − y³)/3 ; skew sextupole: Ψ ∝ (x³ − 3xy²)/3
Zn = (3 .* X.^2 .* Y .- Y.^3) ./ 3
Zs = (X.^3 .- 3 .* X .* Y.^2) ./ 3

fig = Figure(size=(800, 380))
ax1 = Axis(fig[1, 1]; aspect=DataAspect(), title="Normal sextupole")
contour!(ax1, xs, xs, Zn; levels= [-6,-4,-2,2,4,6].* 3.0 ,colormap = :coolwarm, linewidth=2)
uₙ, vₙ = 2 .* X .* Y, X.^2 .- Y.^2          # B ∝ ∇Ψ components
arrows2d!(ax1, sel(X), sel(Y), sel(uₙ), sel(vₙ); color=:gray, lengthscale=0.01)

ax2 = Axis(fig[1, 2]; aspect=DataAspect(), title="Skew sextupole")
contour!(ax2, xs, xs, Zs; levels=[-6,-4,-2,2,4,6].*3.0, colormap = :coolwarm, linewidth=2)
uₛ, vₛ =X.^2 .- Y.^2, (-2) .* X .* Y         # B ∝ ∇Ψ components
arrows2d!(ax2, sel(X), sel(Y), sel(uₛ), sel(vₛ); color=:gray, lengthscale=0.01)
fig
```

In TrackPad these elements appear as `Drift(L)`, `Quadrupole(L, k1)`,
`Sextupole(L, k2)`, `Octupole(L, k3)`, `SBend(L, angle, e1, e2)` and friends,
with strengths normalized by rigidity ($k_1=G/(B\rho)$ in
$\mathrm{m}^{-2}$, $k_2$ not pre-divided by factorials).

## What each multipole order does to a beam

The field pictures above become dynamics as soon as a particle flies through
them. Because the field of the $n^{\text{th}}$ multipole grows like $x^{n}$, the
angular kick it delivers does too:

$$
\Delta p_x \propto -k_n\,\ell\,x^{n},
$$

which is why only the quadrupole ($n=1$) acts as a *lens*: it is the one order
whose kick is proportional to the offset, so every particle in the bunch sees
the same focal length. Launch a fan of parallel rays through a single magnet and
compare:

```{code-cell} julia
:tags: [hide-input]

using TrackPad, TrackPadWidgets

beam = Beam(3.0e9)

"Track a bunch element by element, recording x and y at every boundary."
function track_s(pieces, beam, coords0)
    c = copy(coords0); flags = zeros(Int, size(c, 1))
    S = Float64[0.0]; X = [copy(c[:, 1])]; Y = [copy(c[:, 3])]
    for e in pieces
        linepass!(c, Lattice(AbstractElement[e]), beam, flags)
        push!(S, S[end] + get_length(e))
        push!(X, copy(c[:, 1])); push!(Y, copy(c[:, 3]))
    end
    S, reduce(hcat, X), reduce(hcat, Y)      # X[particle, step]
end

KINDS = [("quadrupole", (L, g) -> Quadrupole(L, g), "k₁ [m⁻²]", 4.0,   "Δpₓ ∝ x"),
         ("sextupole",  (L, g) -> Sextupole(L, g),  "k₂ [m⁻³]", 400.0, "Δpₓ ∝ x²"),
         ("octupole",   (L, g) -> Octupole(L, g),   "k₃ [m⁻⁴]", 4.0e4, "Δpₓ ∝ x³")]
Lm   = 0.2
xin  = collect(range(-15e-3, 15e-3, length=41))
xfan = collect(range(-15e-3, 15e-3, length=11))

explorer(
    title   = "One multipole, one ray fan",
    sliders = [Knob("magnet", 1:3; fmt = i -> KINDS[i][1]),
               Knob("relative strength", range(0.0, 1.0, length=11);
                      fmt = f -> string(Int(round(100f)), " %"), init = 7)],
    panels  = [Panel(xlabel="entrance x [mm]", ylabel="Δpₓ at exit [mrad]",
                     title="kick versus offset", height=250, legend=:bottomright),
               Panel(xlabel="s [m]", ylabel="x [mm]", title="ray fan",
                     ylim=(-20.0, 20.0), height=250, legend=:bottomleft)],
    note = "Each magnet's strength scale is chosen so that the three orders give a comparable "*
           "kick at x = 15 mm — the difference you see is the shape of the kick, not its size.",
) do i, f
    name, make, unit, gmax, law = KINDS[i]
    g   = f * gmax
    mag = make(Lm, g)

    c0 = zeros(length(xin), 6); c0[:, 1] .= xin
    linepass!(c0, Lattice(AbstractElement[mag]), beam, zeros(Int, length(xin)))

    pieces = AbstractElement[Drift(0.1) for _ in 1:5]
    append!(pieces, [make(Lm/4, g) for _ in 1:4])
    append!(pieces, [Drift(0.15)    for _ in 1:14])
    cf = zeros(length(xfan), 6); cf[:, 1] .= xfan
    S, X, _ = track_s(pieces, beam, cf)

    (series = vcat(
        [line(xin .* 1e3, c0[:, 2] .* 1e3; panel=1, color=PALETTE[i], label=name)],
        [line(S, X[j, :] .* 1e3; panel=2, color=PALETTE[i], alpha=0.85, width=1.4)
         for j in eachindex(xfan)],
        [line([0.5, 0.5, nothing, 0.7, 0.7], [-20.0, 20.0, nothing, -20.0, 20.0];
              panel=2, color="#9aa4b2", dash=true, width=1, label="magnet")]),
     readouts = [unit => round(g; sigdigits=4),
                 "scaling" => law,
                 "Δpₓ at x = 5 mm"  => string(round(1e3*c0[findmin(abs.(xin .- 5e-3))[2], 2];
                                                    sigdigits=3), " mrad"),
                 "Δpₓ at x = 15 mm" => string(round(1e3*c0[end, 2]; sigdigits=3), " mrad")])
end
```

The quadrupole sends every ray to one point — a lens. The sextupole leaves the
core of the beam untouched and bends the two tails to the *same* side, because
$x^{2}$ is even; the octupole bends them to opposite sides but with a strength
that grows as $x^{3}$. Both are indispensable (sextupoles correct chromaticity,
octupoles damp collective instabilities) and both are dangerous: a particle at
twice the amplitude gets four or eight times the kick, which is the origin of
the dynamic-aperture problem.
