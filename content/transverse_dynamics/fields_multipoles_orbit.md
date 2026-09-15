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

const C_LIGHT = 2.99792458e8         # [m/s]

"""
    rigidity(beam)

Bρ = P/q in T·m. `beam.energy` is KINETIC energy, so the total energy is
`beam.energy + beam.mass` and P₀c = β·E₀. Dividing by the charge is what makes
this correct for an ion and not just for a proton or electron.
"""
rigidity(beam) = beam.beta * (beam.energy + beam.mass) / (abs(beam.charge) * C_LIGHT)

# 10 keV … 10 TeV, evenly spaced in decades.
kinetic_scan = 10.0 .^ range(4, 13, length=46)

SPECIES = [(name="electron", mass=TrackPad.M_ELECTRON, charge=-1.0),
           (name="proton",   mass=TrackPad.M_PROTON,   charge=+1.0)]

"Format an energy as e.g. \"3.2e7 eV\" — compact enough for a slider label."
function energy_label(K)
    decade   = floor(Int, log10(K))
    mantissa = round(K / 10.0^decade; digits=1)
    return string(mantissa, "e", decade, " eV")
end

# One Bρ(K) curve per species, drawn behind every frame.
rigidity_curves = [line(kinetic_scan,
                        [rigidity(Beam(K; mass=s.mass, charge=s.charge)) for K in kinetic_scan];
                        color=PALETTE[i], label=s.name)
                   for (i, s) in enumerate(SPECIES)]

explorer(
    title   = "Rigidity does not saturate",
    sliders = [Knob("kinetic energy", eachindex(kinetic_scan);
                    fmt = i -> energy_label(kinetic_scan[i]), init = 30),
               Knob("particle", eachindex(SPECIES); fmt = j -> SPECIES[j].name)],
    panels  = [Panel(xlabel="kinetic energy K [eV]", ylabel="Bρ [T·m]",
                     xscale=:log10, yscale=:log10, xlim=(8.0e3, 1.3e13),
                     height=300, legend=:bottomright)],
    statics = rigidity_curves,
    note = "Bρ = P/q is the natural normalisation of transverse motion: TrackPad's k₁ is "*
           "G/(Bρ). At low energy the two species differ by their rest masses; once both "*
           "are ultra-relativistic the curves merge, because Bρ then depends only on energy.",
) do i, j
    species = SPECIES[j]
    K       = kinetic_scan[i]
    beam    = Beam(K; mass=species.mass, charge=species.charge)
    Bρ      = rigidity(beam)

    (series = [points([K], [Bρ]; color=PALETTE[j], size=6.5, label=species.name)],
     readouts = ["particle" => species.name,
                 "K"   => energy_label(K),
                 "β"   => round(beam.beta; sigdigits=6),
                 "P₀c" => string(round(beam.beta*(beam.energy + beam.mass)/1e9;
                                       sigdigits=4), " GeV"),
                 "Bρ"  => string(round(Bρ; sigdigits=4), " T·m"),
                 # The field a 10 m bending radius would need — the number that
                 # decides whether the magnet can be iron or has to be supercon.
                 "B for ρ = 10 m" => string(round(Bρ/10; sigdigits=4), " T")])
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

using CairoMakie

# A multipole field in the gap is source-free, so it derives from a scalar
# potential Ψ: B ∝ ∇Ψ. Contours of Ψ are therefore the field LINES, and the pole
# faces of a real magnet are cut along them. Each magnet below is specified by
# its Ψ and by the ∇Ψ we draw as arrows.

xs = range(-5, 5, length=41)

# Grid convention on this page: index 1 runs along x, index 2 along y, which is
# the order Makie's contour!(x, y, z) expects.
grid_x = [x for x in xs, _ in xs]
grid_y = [y for _ in xs, y in xs]

# One arrow every third grid point, so the quiver stays readable.
every = 3
sample(a) = vec(a[1:every:end, 1:every:end])

"""
    multipole_panel!(ax, title, Ψ, Bx, By; levels, arrowscale)

Draw one magnet: contours of the potential `Ψ` (the field lines) with the field
`(Bx, By) = ∇Ψ` on top. All three are matrices on the (grid_x, grid_y) grid.
"""
function multipole_panel!(ax, title, Ψ, Bx, By; levels, arrowscale)
    ax.title = title
    contour!(ax, xs, xs, Ψ; levels, colormap=:coolwarm, linewidth=2)
    arrows2d!(ax, sample(grid_x), sample(grid_y), sample(Bx), sample(By);
              color=:gray, lengthscale=arrowscale)
end

# Normal quadrupole: Ψ = xy, so B ∝ (y, x) — zero on the axis, growing linearly.
# Skew quadrupole: the same field rotated by 45°, Ψ = (x² − y²)/2, B ∝ (x, −y).
fig = Figure(size=(800, 380))
quad_levels = [-6, -4, -2, 2, 4, 6]

multipole_panel!(Axis(fig[1, 1]; aspect=DataAspect()), "Normal quad",
                 grid_x .* grid_y, grid_y, grid_x;
                 levels=quad_levels, arrowscale=0.08)

multipole_panel!(Axis(fig[1, 2]; aspect=DataAspect()), "Skew quad",
                 (grid_x.^2 .- grid_y.^2) ./ 2, grid_x, -grid_y;
                 levels=quad_levels, arrowscale=0.08)

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

# Same construction as the quadrupole figure, one order higher: the potential is
# now cubic in the coordinates, so the field grows like x² instead of like x.
# `multipole_panel!`, `grid_x`, `grid_y` and `sample` come from the cell above.

# Normal sextupole: Ψ = (3x²y − y³)/3, giving B ∝ (2xy, x² − y²).
# Skew sextupole:   Ψ = (x³ − 3xy²)/3, giving B ∝ (x² − y², −2xy).
Ψ_normal = (3 .* grid_x.^2 .* grid_y .- grid_y.^3) ./ 3
Ψ_skew   = (grid_x.^3 .- 3 .* grid_x .* grid_y.^2) ./ 3

sext_levels = [-6, -4, -2, 2, 4, 6] .* 3.0

fig = Figure(size=(800, 380))

multipole_panel!(Axis(fig[1, 1]; aspect=DataAspect()), "Normal sextupole",
                 Ψ_normal, 2 .* grid_x .* grid_y, grid_x.^2 .- grid_y.^2;
                 levels=sext_levels, arrowscale=0.01)

multipole_panel!(Axis(fig[1, 2]; aspect=DataAspect()), "Skew sextupole",
                 Ψ_skew, grid_x.^2 .- grid_y.^2, -2 .* grid_x .* grid_y;
                 levels=sext_levels, arrowscale=0.01)

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

beam = Beam(3.0e9)                   # 3 GeV electron, used by every tracking cell here

# Column layout of a TrackPad coordinate array: (x, pₓ, y, p_y, z, δE).
const IX, IPX = 1, 2

"""
    track_and_record(pieces, beam, coords0)

Track the bunch `coords0` through `pieces` one element at a time, saving the
full coordinate array after every element.

Returns `(s, states)`: `s[i]` is the distance from the start of the line and
`states[i]` is a copy of the `nparticle × 6` array there. Cutting a magnet into
several short `pieces` is what gives a smooth curve *inside* the magnet —
`linepass!` on its own only reports the exit.
"""
function track_and_record(pieces, beam, coords0)
    coords = copy(coords0)
    lost   = zeros(Int, size(coords, 1))
    s      = [0.0]
    states = [copy(coords)]
    for element in pieces
        linepass!(coords, Lattice(AbstractElement[element]), beam, lost)
        push!(s, s[end] + get_length(element))
        push!(states, copy(coords))
    end
    return s, states
end

"Coordinate `col` of particle `i` at every recorded station."
trajectory(states, i, col) = [state[i, col] for state in states]

# The three magnets to compare. `strength_max` is chosen per order so that all
# three deliver a similar kick at the edge of the fan (x = 15 mm); what differs
# is the SHAPE of the kick, which is the point of the figure.
MAGNETS = [(name="quadrupole", build=Quadrupole, unit="k₁ [m⁻²]",
            strength_max=4.0,   law="Δpₓ ∝ x"),
           (name="sextupole",  build=Sextupole,  unit="k₂ [m⁻³]",
            strength_max=400.0, law="Δpₓ ∝ x²"),
           (name="octupole",   build=Octupole,   unit="k₃ [m⁻⁴]",
            strength_max=4.0e4, law="Δpₓ ∝ x³")]

# Beamline geometry: 0.5 m of drift, the magnet, then 2.1 m to watch the fan.
const L_MAGNET  = 0.2
const S_MAGNET  = (0.5, 0.7)                   # where the magnet sits, for the marker
const N_SLICES  = 4                            # magnet slices, for a smooth ray inside it

# Two sets of starting offsets: a fine one to measure the kick-versus-offset
# curve, and a coarse one to draw as a ray fan.
offsets_fine = collect(range(-15e-3, 15e-3, length=41))
offsets_fan  = collect(range(-15e-3, 15e-3, length=11))

"A dashed vertical pair marking where the magnet begins and ends."
magnet_marker(ylo, yhi) = ([S_MAGNET[1], S_MAGNET[1], nothing, S_MAGNET[2], S_MAGNET[2]],
                           [ylo, yhi, nothing, ylo, yhi])

explorer(
    title   = "One multipole, one ray fan",
    sliders = [Knob("magnet", eachindex(MAGNETS); fmt = i -> MAGNETS[i].name),
               Knob("relative strength", range(0.0, 1.0, length=11);
                      fmt = f -> string(Int(round(100f)), " %"), init = 7)],
    panels  = [Panel(xlabel="entrance x [mm]", ylabel="Δpₓ at exit [mrad]",
                     title="kick versus offset", height=250, legend=:bottomright),
               Panel(xlabel="s [m]", ylabel="x [mm]", title="ray fan",
                     ylim=(-20.0, 20.0), height=250, legend=:bottomleft)],
    note = "Each magnet's strength scale is chosen so that the three orders give a comparable "*
           "kick at x = 15 mm — the difference you see is the shape of the kick, not its size.",
) do i, fraction
    magnet   = MAGNETS[i]
    strength = fraction * magnet.strength_max

    # Panel 1: send a fine comb of offsets through the magnet alone and read off
    # the exit angle. Every particle starts parallel to the axis (pₓ = 0), so the
    # exit pₓ IS the kick.
    comb = zeros(length(offsets_fine), 6)
    comb[:, IX] .= offsets_fine
    linepass!(comb, Lattice(AbstractElement[magnet.build(L_MAGNET, strength)]),
              beam, zeros(Int, length(offsets_fine)))
    kick = comb[:, IPX]

    # Panel 2: the same magnet inside a beamline, so the rays can be watched.
    drift_before = [Drift(0.1)  for _ in 1:5]
    magnet_parts = [magnet.build(L_MAGNET/N_SLICES, strength) for _ in 1:N_SLICES]
    drift_after  = [Drift(0.15) for _ in 1:14]
    beamline = AbstractElement[drift_before..., magnet_parts..., drift_after...]

    fan = zeros(length(offsets_fan), 6)
    fan[:, IX] .= offsets_fan
    s, states = track_and_record(beamline, beam, fan)

    rays = [line(s, trajectory(states, j, IX) .* 1e3; panel=2,
                 color=PALETTE[i], alpha=0.85, width=1.4)
            for j in eachindex(offsets_fan)]

    # Report the kick at two offsets; a factor-of-3 step in x shows the power law.
    at_5mm  = argmin(abs.(offsets_fine .- 5e-3))
    at_15mm = argmin(abs.(offsets_fine .- 15e-3))

    (series = vcat(
        [line(offsets_fine .* 1e3, kick .* 1e3; panel=1, color=PALETTE[i],
              label=magnet.name)],
        rays,
        [line(magnet_marker(-20.0, 20.0)...; panel=2, color="#9aa4b2", dash=true,
              width=1, label="magnet")]),
     readouts = [magnet.unit => round(strength; sigdigits=4),
                 "scaling"   => magnet.law,
                 "Δpₓ at x = 5 mm"  => string(round(kick[at_5mm]  * 1e3; sigdigits=3), " mrad"),
                 "Δpₓ at x = 15 mm" => string(round(kick[at_15mm] * 1e3; sigdigits=3), " mrad")])
end
```

The quadrupole sends every ray to one point — a lens. The sextupole leaves the
core of the beam untouched and bends the two tails to the *same* side, because
$x^{2}$ is even; the octupole bends them to opposite sides but with a strength
that grows as $x^{3}$. Both are indispensable (sextupoles correct chromaticity,
octupoles damp collective instabilities) and both are dangerous: a particle at
twice the amplitude gets four or eight times the kick, which is the origin of
the dynamic-aperture problem.
