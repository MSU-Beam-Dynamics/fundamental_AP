---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Dynamics in a Linac Accelerator

## The gap-by-gap map

Consider a standing-wave linac: a series of short accelerating gaps. Between
gaps $n-1$ and $n$ there is a drift of length $l_{n-1}$. In that drift the
ideal particle has velocity $\beta_{n-1}$ and energy $E_{n-1}$; its arrival
time (phase) at gap $n$ is $t_n$ ($\phi_n$). Subscript $s$ marks reference
quantities.

The energy change across gap $n$ is

$$
\Delta E_n=E_{n}-E_{s,n}
=\Delta E_{n-1}+eV\left(\sin\left(\phi_{s,n}+\Delta\phi_n\right)-\sin\phi_{s,n}\right),
$$

and for small phase deviations,

$$
\Delta E_n \approx \Delta E_{n-1}+eV \cos(\phi_{s,n})\,\Delta\phi_n .
$$

The phase change accumulates from the velocity difference in the drift:

$$
\Delta \phi_n=\omega_0\left(t_{n}-t_{s,n}\right)
=\Delta \phi_{n-1}+\frac{\omega}{c}\left(\frac{l_{n-1}}{\beta_{n-1}}-\frac{l_{n-1}}{\beta_{s,n-1}}\right)
\approx
\Delta \phi_{n-1}-\frac{\omega}{mc^3}\,\frac{l_{n-1}}{\beta_{s,n-1}^3 \gamma_{s,n-1}^3}\,\Delta E_{n-1}.
$$

:::{note}
A linac has no bends, so its momentum compaction factor vanishes,
$\alpha_c=0$, and the slip factor $\eta=-1/\gamma_s^2<0$: a linac is always
"below transition".
:::

## Numerical gap tracking

```{code-cell} julia
using CairoMakie

"""
    track_linac(; ngap, φs, V, f_rf, K0, m_p, dE0)

Gap-by-gap linearized longitudinal map of a standing-wave linac.
Returns (Δφ history, ΔE history).
"""
function track_linac(; ngap=120, φs=deg2rad(25.0), V=2.0e6, f_rf=800.0e6,
                     K0=5.0e6, m_p=0.938e9, dE0=1.0e6)
    λ = 2.99792458e8 / f_rf            # RF wavelength [m]
    β_of(γ) = sqrt(1 - 1/γ^2)
    # φs > 0 with cos φs > 0 gives both acceleration and phase stability
    E, dE, dφ = K0 + m_p, dE0, 0.0     # reference energy & injection errors
    hist_E, hist_φ = Float64[dE], Float64[dφ]
    for _ in 1:ngap
        γs = E / m_p
        l  = β_of(γs) * λ / 2                              # drift ≈ half RF wavelength
        dφ -= (2π/λ) * l / (β_of(γs)^3 * γs^3) * dE / m_p  # phase slip through drift
        dE += V * cos(φs) * dφ                             # gap kick (linearized)
        E  += V * sin(φs)                                  # reference ramps by eV·sin φs
        push!(hist_E, dE); push!(hist_φ, dφ)
    end
    return hist_φ, hist_E
end

hist_dφ, hist_E = track_linac()
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="Δφ [rad]", ylabel="ΔE [MeV]",
          title="Linac longitudinal motion about the synchronous phase")
scatterlines!(ax, hist_dφ, hist_E ./ 1e6; markersize=5)
fig
```

## Differential approximation and small-amplitude motion

If acceleration is slow and gaps are short, the discrete map becomes the ODE

$$
\frac{d\Delta E}{ds}=e\mathcal{E}\left(\sin\left(\phi_{s}+\Delta\phi\right)-\sin\phi_{s}\right),
\qquad
\frac{d\Delta \phi}{ds}=-\frac{\omega}{mc^3\beta_{s}^3 \gamma_{s}^3}\,\Delta E,
$$

with $\mathcal{E}$ the accelerating gradient. This is the ring equation again,
with the replacement $\eta \rightarrow -1/\gamma_s^2$. Linearizing:

$$
\frac{d^2\Delta E}{ds^2}= -\frac{e\mathcal{E}\omega\cos\phi_{s}}{mc^3\beta_{s}^3 \gamma_{s}^3}\,\Delta E ,
$$

with synchrotron oscillation wave number

$$
k_s^{2}= \frac{e\mathcal{E}\omega\cos\phi_{s}}{mc^3\beta_{s}^3 \gamma_{s}^3}.
$$

Because $\eta=-1/\gamma_s^2<0$, stability requires $\cos\phi_{s}>0$: a
particle arriving late gains more voltage, climbs toward higher velocity, and
therefore arrives earlier at the next gap — a restoring motion around the
synchronous phase. As $\gamma_s$ grows along the linac, the slip weakens like
$\gamma_s^{-3}$ and the oscillation adiabatically damps.

## Scanning the synchronous phase

The same `track_linac` map used above, now swept across $\phi_s$: stable and
bunching for $\phi_s<90^\circ$, and no longer confined once $\cos\phi_s$
changes sign, even though the reference particle keeps gaining energy for any
$\phi_s<180^\circ$.

```{code-cell} julia
:tags: [hide-input]

using TrackPadWidgets

explorer(
    title   = "Linac: the synchronous phase decides whether the bunch holds together",
    sliders = [Knob("φs [deg]", range(5.0, 175.0, length=35);
                      fmt = p -> string(Int(round(p))), init = 5),
               Knob("gap voltage [MV]", (1.0, 2.0, 4.0); fmt = v -> string(v), init = 2)],
    panels  = [Panel(xlabel="Δφ [rad]", ylabel="ΔE [MeV]", title="longitudinal phase space",
                     xlim=(-1.6, 1.6), ylim=(-4.0, 4.0), height=250, legend=:bottomleft),
               Panel(xlabel="gap number", ylabel="ΔE [MeV]", title="energy error versus gap",
                     ylim=(-4.0, 4.0), height=250, legend=:bottomleft)],
    note = "A linac has no bends, so αc = 0 and η = −1/γs² < 0 always: stability needs "*
           "cos φs > 0, i.e. φs < 90°. The reference particle gains eV sin φs per gap either "*
           "way, so past 90° the beam still accelerates while the bunch falls apart.",
) do pdeg, VMV
    φs = deg2rad(pdeg)
    hφ, hE = track_linac(; φs=φs, V=VMV*1e6)
    ok = cos(φs) > 0
    col = ok ? PALETTE[1] : PALETTE[4]
    (series = [line(hφ, hE ./ 1e6; panel=1, color=col, width=1.4, label="trajectory"),
               points([hφ[1]], [hE[1]/1e6]; panel=1, color=PALETTE[3], size=5.0, label="injection"),
               line(0:length(hE)-1, hE ./ 1e6; panel=2, color=col, width=1.4, label="ΔE")],
     readouts = ["cos φs" => round(cos(φs); digits=3),
                 "sin φs" => round(sin(φs); digits=3),
                 "verdict" => ok ? "phase-stable" : "phase-unstable",
                 "|ΔE| after 120 gaps" => string(round(abs(hE[end])/1e6; sigdigits=3), " MeV")])
end
```
