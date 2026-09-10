---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# RF Acceleration

## Accelerating structures

Modern accelerators use RF/SRF cavities as accelerating devices:

```{figure} ../images/rf_cavity.png
:width: 400px
:name: fig:rfcavity
An RF cavity (typically operating in its TM₀₁₀ mode).
```

The accelerating field of the TM010-like mode is written as

$$
E(s,t)=E_0(s)\sin\left(\omega t +\phi_0\right),
$$

where $E_0(s)$ is the field envelope, $t$ the time and $\phi_0$ the initial
phase at $t=0$. The envelope potential is $\int_{-\infty}^{\infty}E_0(s)\,ds$.

## Energy gain of an off-phase particle

This is not yet the potential a particle experiences. Consider a particle with
6D coordinates $(x, p_x, y, p_y, z, \delta)$ passing through the cavity,
assuming its velocity does not change appreciably inside. This holds when

- the particle is already ultra-relativistic ($\gamma \gg 1$), or
- the energy gain in the cavity is much smaller than the particle's momentum.

If the reference particle arrives at $s=0$ at $t_0$, our particle is at
$s=\beta c (t-t_0)+z$ when the reference is at $s$. Its energy gain through
the field is

$$
\Delta U=e\int_{-\infty}^{\infty} E_0(s)
\sin\left(\frac{\omega (s-z)}{\beta c}+\phi_s\right)ds ,
$$

where $\phi_s$ absorbs the $t_0$ terms and becomes the **synchronous phase**.
For an even envelope $E_0(s)=E_0(-s)$ this simplifies to

$$
\Delta U = eV\sin \phi,
\qquad
T=\frac{\int_{-\infty}^{\infty} E_0(s) \cos\frac{\omega s}{\beta c}\, ds}
        {\int_{-\infty}^{\infty} E_0(s)\, ds},
$$

defining the equivalent cavity **voltage** $V$ and the **transit-time factor**
$T$. For the reference particle,

$$
\Delta U_0 = eV\sin \phi_s .
$$

Therefore after one cavity passage the energy deviation changes by

$$
\Delta E|_{1}=\Delta E|_{0}+ eV \left(\sin \phi-\sin \phi_s\right).
$$

## The RF kick in TrackPad

TrackPad implements exactly this thin-lens physics in `RFCavity`. Its phase
convention is

$$
\phi=-2\pi f\,\frac{z+\ell_{\mathrm{lag}}}{c}-\phi_{\mathrm{lag}},
\qquad
\delta_E^{+}=\delta_E^{-}-\frac{\widehat q V}{P_0c}\,\sin\phi,
$$

with $\widehat q=q/e$ the *signed* charge.

:::{caution}
A directly constructed `RFCavity` defaults to `energy=0`, which **disables**
its kick. When building cavities by hand, always pass
`energy=beam.energy` and `charge=beam.charge` (file readers such as PALS/MAD-X
set both automatically).
:::

## The kick as a function of arrival phase

Let TrackPad trace out the whole curve: $\delta_E$ after the cavity as a
function of where in the RF wavelength the particle arrives.
$\phi_{\mathrm{lag}}$ slides that curve left and right — it sets which $z$ sees
zero crossing, i.e. which particle is synchronous:

```{code-cell} julia
:tags: [hide-input]

using StaticArrays, TrackPad, TrackPadWidgets

const CL = 2.99792458e8
beam10 = Beam(3.0e9)
f_rf   = 100.0e6
λ_rf   = CL / f_rf
zs = collect(range(-λ_rf/2, λ_rf/2, length=121))
zn = collect(range(-0.15, 0.15, length=61))

kick_curve(cav, zgrid) =
    [linepass(cav, SVector{6}(0.0,0.0,0.0,0.0,z,0.0), beam10)[6] for z in zgrid]

explorer(
    title   = "The RF kick seen by a particle arriving at longitudinal position z",
    sliders = [Knob("φ_lag [deg]", range(0.0, 2π, length=25);
                      fmt = p -> string(Int(round(rad2deg(p)))), init = 7),
               Knob("V [MV]", (1.0, 2.0, 4.0); fmt = v -> string(v), init = 2)],
    panels  = [Panel(xlabel="z [m]", ylabel="δ_E after the cavity ×10⁻³",
                     title="one full RF wavelength", height=250, legend=:bottomleft),
               Panel(xlabel="z [cm]", ylabel="δ_E − δ_E(0)  ×10⁻³",
                     title="zoom on a bunch (±15 cm)", height=250, legend=:bottomleft)],
    note = "z > 0 means the particle arrives early. A negative slope at z = 0 restores an "*
           "early particle towards the synchronous one — that is longitudinal focusing.",
) do φlag, VMV
    V = VMV * 1e6
    cav = Lattice(AbstractElement[RFCavity(0.0, V, f_rf; philag=φlag,
                                           energy=beam10.energy, charge=beam10.charge)])
    dfull = kick_curve(cav, zs)
    dnear = kick_curve(cav, zn)
    d0 = linepass(cav, SVector{6}(zeros(6)...), beam10)[6]
    slope = (dnear[end] - dnear[1]) / (zn[end] - zn[1])
    φ0 = -φlag
    (series = [line(zs, dfull .* 1e3; panel=1, color=PALETTE[1], label="δ_E(z)"),
               points([0.0], [d0*1e3]; panel=1, color=PALETTE[4], size=5.0, label="reference"),
               line(zn .* 100, (dnear .- d0) .* 1e3; panel=2, color=PALETTE[1],
                    label="δ_E − δ_E(0)")],
     readouts = ["φ at z = 0" => string(Int(round(rad2deg(mod(φ0, 2π)))), "°"),
                 "ΔE of the reference" => string(round(d0*beam10.beta*(beam10.energy+beam10.mass)/1e3;
                                                       sigdigits=4), " keV"),
                 "dδ_E/dz" => string(round(slope; sigdigits=3), " m⁻¹"),
                 "longitudinally" => slope < 0 ? "focusing" : (slope > 0 ? "defocusing" : "neutral")])
end
```

This is the same physics that the next chapter turns into a turn-by-turn map: a
particle sitting where the slope is negative gets nudged back every time it
passes the cavity, and repeated small nudges are what build the RF bucket.
