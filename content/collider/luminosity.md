---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Luminosity

## Event rate and cross section

The center-of-mass energy decides *which* processes a collider can reach. The
**luminosity** $\mathcal L$ decides *how often* they happen. For a process with
cross section $\sigma_p$, the event rate is

$$
\frac{dR}{dt}=\sigma_p\,\mathcal L .
$$

The cross section measures the probability of the process. It is set by
physics, not by the machine, and is quoted in barns
($1\ \mathrm b=10^{-24}\ \mathrm{cm}^2$). The luminosity is the machine's share:
the number of collision opportunities per unit area and time, in
$\mathrm{cm^{-2}s^{-1}}$. Its time integral $\int\mathcal L\,dt$, quoted in
inverse femtobarns ($1\ \mathrm{fb^{-1}}=10^{39}\ \mathrm{cm^{-2}}$), counts the
events an experiment collects.

The numbers span many orders of magnitude. At the LHC, the total inelastic
proton–proton cross section is about 80 mb, and producing a Higgs boson has a
cross section of about 60 pb, a billion times smaller. At
$\mathcal L=2\times10^{34}\ \mathrm{cm^{-2}s^{-1}}$ the detectors therefore see
about $1.6\times10^9$ inelastic collisions per second, and about one Higgs boson.

## Fixed target and collider

For a beam on a **fixed target**, the luminosity is the flux of incident
particles $\mathcal F=dN/dt$ times the areal density of target particles, a
number density $n$ times a length $\ell$:

$$
\mathcal L=\mathcal F\,n\,\ell .
$$

For example, $10^{12}$ protons per second on one metre of liquid hydrogen
($n=4.2\times10^{22}\ \mathrm{cm^{-3}}$) give
$\mathcal L\approx4\times10^{36}\ \mathrm{cm^{-2}s^{-1}}$. That is more than any
collider reaches, but at a far lower center-of-mass energy.

In a **collider**, two bunches with $N_1$ and $N_2$ particles and normalized
densities $\rho_1$ and $\rho_2$ pass through each other $f$ times per second
($f$ is the number of colliding bunch pairs times the revolution frequency).
The luminosity is their overlap:

$$
\mathcal L=f\,K\,N_1N_2\int\!\!\int\!\!\int\!\!\int
\rho_1(x,y,s,t)\,\rho_2(x,y,s,t)\;dx\,dy\,ds\,dt ,
\qquad
K=\sqrt{(\mathbf v_1-\mathbf v_2)^2-\frac{(\mathbf v_1\times\mathbf v_2)^2}{c^2}} .
$$

The kinematic factor $K$ equals $2c$ for a head-on collision. The interaction
point (IP) is at $s=0$.

Take Gaussian bunches moving in opposite directions,

$$
\rho_{1,2}(x,y,s,t)=\frac{1}{(2\pi)^{3/2}\sigma_x\sigma_y\sigma_s}
\exp\!\left(-\frac{x^2}{2\sigma_x^2}-\frac{y^2}{2\sigma_y^2}-\frac{(s\mp ct)^2}{2\sigma_s^2}\right).
$$

Each transverse integral gives $\int\rho_x^2\,dx=1/(2\sqrt\pi\,\sigma_x)$. The
longitudinal part, with $K=2c$, integrates to one: the variables
$u=s-ct$, $v=s+ct$ have $ds\,dt=du\,dv/2c$. Together,

$$
\boxed{\;\mathcal L=\frac{f\,N_1N_2}{4\pi\,\sigma_x\sigma_y}\;}
$$

for identical bunches. If the two beams have different sizes, each $\sigma$ is
replaced by the combined width, which gives
$\mathcal L=f N_1N_2/\big(2\pi\sqrt{\sigma_{1x}^2+\sigma_{2x}^2}\sqrt{\sigma_{1y}^2+\sigma_{2y}^2}\big)$.
The bunch length has dropped out entirely, but only because the transverse
sizes were taken constant along the bunch. The next sections show when that
fails.

The formula, evaluated for RHIC and for the LHC:

```{code-cell} julia
using TrackPad, TrackPadWidgets, CairoMakie, Printf, Random

"Head-on luminosity of round Gaussian beams [cm⁻²s⁻¹]."
luminosity(; N, nbunch, frev, σ) = nbunch*frev*N^2/(4π*σ^2) * 1e-4      # m⁻² → cm⁻²

"rms beam size at the IP from the normalized emittance and β*."
function ip_size(beam, εn, βstar)
    ε = εn/(beam.beta*beam.gamma)                   # geometric emittance
    sqrt(ε*βstar)
end

rhic = Beam(255e9 - TrackPad.M_PROTON; mass=TrackPad.M_PROTON, charge=1.0)    # 255 GeV total
lhc  = Beam(7e12  - TrackPad.M_PROTON; mass=TrackPad.M_PROTON, charge=1.0)    # 7 TeV design

σ_rhic = ip_size(rhic, 2.5e-6, 0.65)
σ_lhc  = ip_size(lhc, 3.75e-6, 0.55)
println(@sprintf("RHIC: σ* = %.0f μm, L = %.1e cm⁻²s⁻¹", σ_rhic*1e6,
                 luminosity(N=2.0e11, nbunch=111, frev=78.2e3, σ=σ_rhic)), "\n",
        @sprintf("LHC:  σ* = %.1f μm, L = %.1e cm⁻²s⁻¹", σ_lhc*1e6,
                 luminosity(N=1.15e11, nbunch=2808, frev=11245, σ=σ_lhc)))
```

The LHC value is 20% above the design luminosity of
$10^{34}\ \mathrm{cm^{-2}s^{-1}}$. The difference is the crossing angle
discussed below. For RHIC the formula gives a few $10^{32}$, the scale it
actually reached. Some numbers from colliders past and present:

| collider | energy (GeV) | luminosity ($10^{30}\,\mathrm{cm^{-2}s^{-1}}$) | $\sigma_x/\sigma_y$ (μm) | $N$ per bunch ($10^{10}$) |
|---|---|---|---|---|
| SppS ($p\bar p$) | 315 × 315 | 6 | 60 / 30 | 10 |
| Tevatron ($p\bar p$) | 1000 × 1000 | 100 | 30 / 30 | 30 / 8 |
| HERA ($e p$) | 30 × 920 | 40 | 250 / 50 | 3 / 7 |
| LHC ($pp$) | 7000 × 7000 | 10000 | 17 / 17 | 11 |
| RHIC ($pp$) | 250 × 250 | 100 | 200 / 200 | 20 |
| PEP-II ($e^+e^-$) | 9 × 3 | 8000 | 150 / 5 | 2 / 6 |

Two routes lead to high luminosity. Hadron colliders use many protons per
bunch; $e^+e^-$ factories use many bunches and very flat beams. In both, the
transverse size at the IP is squeezed as far as the optics allow, which is the
subject of [](interaction_region.md).

## Watching two bunches collide

The widget below follows two bunches through the IP, from $ct=-3\sigma_s$ to
$+3\sigma_s$. The left panel shows the particles of each bunch in the
$(s,x)$ plane. The axes are scaled to $\sigma_s$ and to the transverse size
$\sigma_x^*$ at the IP; a real bunch is a thousand times longer than it is
wide. The right panel shows the rate of collisions at each instant, normalized
to the peak of an ideal head-on crossing (gray). The shaded area so far is the
luminosity accumulated up to that time.

The *configuration* slider steps through the effects of the next three
sections: the hourglass effect, a crossing angle, and crab crossing. Press
▶ to watch the bunches pass through each other.

```{code-cell} julia
:tags: [hide-input]

# Units: s in σ_s, x in σ_x*, time as ct in σ_s.  A configuration is
# (name, β*/σ_s, Piwinski angle Φ = θ_c σ_s/σ_x*, crab cavities on?).
const CONFIGS = [("head-on", Inf, 0.0, false),
                 ("hourglass, β* = σₛ/2", 0.5, 0.0, false),
                 ("crossing angle, Φ = 1", Inf, 1.0, false),
                 ("crossing angle, Φ = 1, crab cavities", Inf, 1.0, true)]
const TIMES = collect(range(-3, 3, length=25))
const NDRAW = 160

# Particles at the IP: longitudinal position z, transverse position and angle.
rng = MersenneTwister(11)
const BUNCH = [(z = randn(rng, NDRAW), x = randn(rng, NDRAW), xp = randn(rng, NDRAW)) for _ in 1:2]

hourglass(s, βs) = isinf(βs) ? 1.0 : sqrt(1 + (s/βs)^2)       # σ(s)/σ*

"Particle positions (s, x) of bunch `b` (1 moves to +s, 2 to -s) at time t."
function bunch_at(b, t, βs, Φ, crab)
    B = BUNCH[b]; sgn = b == 1 ? 1 : -1
    path = B.z .+ t                                # distance travelled past the IP
    s = sgn .* path
    x = B.x .+ (isinf(βs) ? 0.0 : 1/βs) .* B.xp .* path    # drift: the hourglass
    x = x .+ sgn .* Φ .* s                         # each beam is tilted by ±θ_c
    crab && (x = x .- Φ .* B.z)                    # the crab kick rotates the bunch
    round.(s; digits=2), round.(x; digits=2)
end

"Collision rate at time t (Gaussian overlap integrated over x, y and s)."
function rate(t, βs, Φ, crab)
    s = range(-6, 6, length=1201); ds = step(s)
    g(u) = exp(-u^2/2)/sqrt(2π)
    sum(s) do si
        σ = hourglass(si, βs)
        d = crab ? 0.0 : 2Φ*si                      # separation of the two bunch axes
        2sqrt(π) * g(si - t)*g(si + t) * exp(-d^2/(4σ^2))/σ^2 * ds
    end
end

const TFINE = collect(range(-4, 4, length=97))
RATES = Dict(c[1] => [rate(t, c[2], c[3], c[4]) for t in TFINE] for c in CONFIGS)
REF = RATES["head-on"]
total(r) = sum(r)*(TFINE[2] - TFINE[1])

explorer(
    title   = "Two bunches passing through each other at the interaction point",
    sliders = [Knob("ct / σₛ", TIMES; fmt = t -> @sprintf("%+.2f", t), init = 9),
               Knob("configuration", 1:length(CONFIGS); fmt = i -> CONFIGS[i][1], init = 1, newline = true)],
    panels  = [Panel(xlabel="s / σₛ", ylabel="x / σₓ*", title="the two bunches",
                     xlim=(-5.0, 5.0), ylim=(-6.0, 6.0), height=300, legend=:bottomright),
               Panel(xlabel="ct / σₛ", ylabel="collision rate", title="collisions per unit time",
                     xlim=(-3.0, 3.0), ylim=(0.0, 1.08), height=300, legend=:bottomright)],
    statics = [line(TFINE, REF; panel=2, color="#9aa4b2", dash=true, width=1.5, label="head-on, no hourglass")],
    note    = "Hourglass: particles drift away from the axis on either side of the IP, so the beam "*
              "is small only near s = 0, and collisions far from it are diluted. Crossing angle: "*
              "the bunches overlap only near the IP, over a length σₓ*/θ_c. Crab cavities tilt each "*
              "bunch by θ_c so they meet head-on and recover the full luminosity.",
) do t, ic
    name, βs, Φ, crab = CONFIGS[ic]
    s1, x1 = bunch_at(1, t, βs, Φ, crab)
    s2, x2 = bunch_at(2, t, βs, Φ, crab)
    r = RATES[name]
    upto = TFINE .<= t + 1e-9
    series = Any[points(s1, x1; panel=1, color=PALETTE[1], size=3, alpha=0.6, label="bunch 1 →"),
                 points(s2, x2; panel=1, color=PALETTE[4], size=3, alpha=0.6, label="← bunch 2"),
                 line(TFINE, r; panel=2, color=PALETTE[3], width=1.5),
                 line(TFINE[upto], r[upto]; panel=2, color=PALETTE[3], width=2.5, fill=true,
                      fillalpha=0.3, label="this configuration")]
    if !isinf(βs)                                   # the hourglass envelope ±σ(s)
        sg = range(-5, 5, length=81)
        push!(series, line(sg, hourglass.(sg, βs); panel=1, color="#9aa4b2", width=1.2),
                      line(sg, -hourglass.(sg, βs); panel=1, color="#9aa4b2", width=1.2))
    end
    if Φ > 0 && !crab                                # the two beam axes
        push!(series, line([-5, 5], [-5Φ, 5Φ]; panel=1, color="#9aa4b2", dash=true, width=1.2),
                      line([-5, 5], [5Φ, -5Φ]; panel=1, color="#9aa4b2", dash=true, width=1.2))
    end
    (series = series,
     readouts = ["β*/σₛ" => isinf(βs) ? "∞" : string(βs),
                 "Piwinski angle Φ" => string(Φ),
                 "luminosity / head-on" => @sprintf("%.3f", total(r)/total(REF)),
                 "collected so far" => @sprintf("%.0f%%", 100*total(r .* upto)/total(r))])
end
```

## The hourglass effect

Around the IP the beam travels through a drift, and there the Courant–Snyder
function grows as (see [](interaction_region.md))

$$
\beta(s)=\beta^*+\frac{s^2}{\beta^*},
\qquad
\sigma(s)=\sigma^*\sqrt{1+s^2/\beta^{*2}} .
$$

The beam is narrowest at the IP and widens on both sides. Seen from the side,
the envelope has the shape of an hourglass. A bunch of length $\sigma_s$
collides over $|s|\lesssim\sigma_s$. If $\sigma_s$ is comparable to $\beta^*$,
many collisions happen where the beams are already wider, and the luminosity
drops. Putting $\sigma(s)$ into the overlap integral for round beams, with
$u=s/\sigma_s$, gives

$$
\frac{\mathcal L}{\mathcal L_0}=\int_{-\infty}^{\infty}\frac{1}{\sqrt\pi}\,
\frac{e^{-u^2}}{1+u^2\sigma_s^2/\beta^{*2}}\,du
=\sqrt\pi\,t\,e^{t^2}\operatorname{erfc}(t),\qquad t=\frac{\beta^*}{\sigma_s}.
$$

Neglecting the hourglass effect requires $\sigma_s\ll\beta^*$. This is why
squeezing $\beta^*$ is useful only down to about the bunch length. Making
$\beta^*$ smaller than the bunch length gains little luminosity.

## Crossing angle

With many bunches, the two beams have to be separated outside the IP. Otherwise
every bunch would meet the bunches of the other beam at several parasitic
crossings. The standard solution is to let the beams cross at a small full angle
$\phi=2\theta_c$. The bunches then overlap only along a length of about
$\sigma_x^*/\theta_c$ instead of along their full length. The loss is governed by
the **Piwinski angle**

$$
\Phi=\frac{\sigma_s}{\sigma_x^*}\tan\frac{\phi}{2},
\qquad
\frac{\mathcal L}{\mathcal L_0}=\frac{1}{\sqrt{1+\Phi^2}} .
$$

In the LHC, $\phi=285\ \mu$rad, $\sigma_s=7.55$ cm and $\sigma_x^*=16.7\ \mu$m
give $\Phi=0.64$ and a factor 0.84: the 20% found above. SuperKEKB goes to the
other extreme. With $\Phi\approx25$ it gives up a factor 25 in geometric overlap
on purpose. The overlap region, only $\sigma_x^*/\theta_c\approx0.25$ mm long, is
shorter than the bunch, so $\beta_y^*$ can be squeezed to 0.3 mm without an
hourglass penalty. This is the "nano-beam" scheme.

```{code-cell} julia
:tags: [hide-input]
:label: fig-geometric-factors
:caption: Geometric luminosity reduction factors. Left, the hourglass effect against β*/σₛ. Right, a crossing angle against the Piwinski angle Φ.

let fig = Figure(size=(900, 360), fontsize=15)
    t = 10.0 .^ range(-1.3, 1, length=200)
    u = range(-10, 10, length=8001); du = step(u)
    R = [sum(@. exp(-u^2)/(sqrt(π)*(1 + u^2/ti^2)))*du for ti in t]
    ax = Axis(fig[1, 1]; xscale=log10, xlabel="β* / σₛ", ylabel="ℒ / ℒ₀", title="hourglass effect",
              limits=(0.05, 10, 0, 1.05), xticks=([0.1, 0.3, 1, 3, 10], ["0.1", "0.3", "1", "3", "10"]))
    lines!(ax, t, R; color=PALETTE[1], linewidth=3)
    for (tv, name) in [(0.5, "β* = σₛ/2"), (4.0, "LHC, β* = 30 cm")]
        Rv = sum(@. exp(-u^2)/(sqrt(π)*(1 + u^2/tv^2)))*du
        scatter!(ax, [tv], [Rv]; color=:black, markersize=10)
        text!(ax, tv, Rv; text=name, align=tv > 1 ? (:right, :top) : (:left, :top),
              offset=tv > 1 ? (-6, -6) : (6, -4))
    end
    Φ = range(0, 4, length=200)
    ax = Axis(fig[1, 2]; xlabel="Piwinski angle Φ", ylabel="ℒ / ℒ₀", title="crossing angle",
              limits=(0, 4, 0, 1.05))
    lines!(ax, Φ, @. 1/sqrt(1 + Φ^2); color=PALETTE[2], linewidth=3)
    for (Φv, name) in [(0.64, "LHC"), (1.0, "Φ = 1")]
        scatter!(ax, [Φv], [1/sqrt(1 + Φv^2)]; color=:black, markersize=10)
        text!(ax, Φv, 1/sqrt(1 + Φv^2); text=name, align=(:left, :bottom), offset=(6, 4))
    end
    fig
end
```

## Crab crossing

A **crab cavity** is an RF cavity operated to give a transverse kick that
depends on the longitudinal position, $\Delta x'\propto\sin(k z)\approx kz$. The
head of the bunch is kicked one way, the tail the other, and the center not at
all. Placed a quarter of a betatron oscillation upstream of the IP, the kick
becomes a transverse displacement proportional to $z$. Each bunch then arrives
tilted by $\theta_c$, like a crab walking sideways, and the two bunches meet
head-on even though their centers cross at an angle. A second set of cavities
after the IP removes the tilt again. Crab crossing was first used at KEKB, and is
part of the High-Luminosity LHC and the EIC designs. Set the last configuration
in the widget to see it.

## Integrated luminosity

During a store the luminosity decays. Particles are lost in collisions and to
the beam–beam effect, and the emittances grow from intrabeam scattering. A
simple model is an exponential decay with a time constant $\tau$:

$$
\mathcal L(t)=\mathcal L_0\,e^{-t/\tau}.
$$

After a run of length $t_r$, the beams are dumped. Preparing the next store
(ramping down, injecting, accelerating, squeezing) takes a time $t_p$, during
which nothing is collected. The quantity that counts is the luminosity averaged
over the whole cycle:

$$
\langle\mathcal L\rangle=\frac{\mathcal L_0\,\tau\left(1-e^{-t_r/\tau}\right)}{t_r+t_p}.
$$

Setting $d\langle\mathcal L\rangle/dt_r=0$ gives the condition
$e^{t_r/\tau}=1+(t_r+t_p)/\tau$, whose approximate solution is

$$
t_r\approx\tau\ln\!\left(1+\sqrt{2t_p/\tau}+t_p/\tau\right).
$$

Long preparation times favor long stores. A collider that refills quickly can
dump earlier and run nearer to its peak luminosity.

```{code-cell} julia
:tags: [hide-input]

avgL(tr, τ, tp) = τ*(1 - exp(-tr/τ))/(tr + tp)               # in units of ℒ₀

"Optimal run time from e^{t/τ} = 1 + (t + t_p)/τ, solved by bisection."
function optimal_run(τ, tp)
    f(t) = exp(t/τ) - 1 - (t + tp)/τ
    lo, hi = 1e-6, 50τ
    for _ in 1:100
        mid = (lo + hi)/2
        f(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    (lo + hi)/2
end
approx_run(τ, tp) = τ*log(1 + sqrt(2tp/τ) + tp/τ)

const TAUS = [3.0, 5.0, 8.0, 12.0, 20.0, 30.0]       # luminosity lifetime [h]
const TPREP = [0.5, 1.0, 2.0, 3.0, 5.0, 8.0]         # preparation time [h]

explorer(
    title   = "Choosing the length of a store",
    sliders = [Knob("luminosity lifetime τ [h]", TAUS; init = 3),
               Knob("preparation time t_p [h]", TPREP; init = 3)],
    panels  = [Panel(xlabel="time [h]", ylabel="ℒ / ℒ₀", title="three days of operation",
                     xlim=(0.0, 72.0), ylim=(0.0, 1.05), height=260, basis="100%", legend=:bottomright),
               Panel(xlabel="run time t_r [h]", ylabel="⟨ℒ⟩ / ℒ₀", title="average over a cycle",
                     xlim=(0.0, 40.0), ylim=(0.0, 1.0), height=240, basis="100%", legend=:right)],
    note    = "Each store is ended at the optimum t_r. The dashed line is the average over a "*
              "cycle, which is what the integrated luminosity grows with.",
) do τ, tp
    tr = optimal_run(τ, tp)
    T, L = Float64[], Float64[]
    t0 = 0.0
    while t0 < 72
        ts = range(t0, t0 + tr, length=16)
        append!(T, [t0; ts; t0 + tr; t0 + tr + tp]); append!(L, [0.0; exp.(-(ts .- t0)/τ); 0.0; 0.0])
        t0 += tr + tp
    end
    trs = range(0.05, 40, length=120)
    av = avgL(tr, τ, tp)
    (series = [line(T, L; panel=1, color=PALETTE[1], width=1.8, fill=true, fillalpha=0.2, label="ℒ(t)"),
               line([0, 72], [av, av]; panel=1, color=PALETTE[4], dash=true, width=2, label="cycle average"),
               line(trs, avgL.(trs, τ, tp); panel=2, color=PALETTE[2], width=2.5, label="⟨ℒ⟩(t_r)"),
               points([tr], [av]; panel=2, color=PALETTE[4], size=9, label="optimum"),
               points([approx_run(τ, tp)], [avgL(approx_run(τ, tp), τ, tp)]; panel=2, color=PALETTE[3],
                      size=6, label="approximation")],
     readouts = ["optimal t_r" => @sprintf("%.1f h", tr),
                 "approximation" => @sprintf("%.1f h", approx_run(τ, tp)),
                 "⟨ℒ⟩ / ℒ₀" => @sprintf("%.2f", av),
                 "stores per day" => @sprintf("%.1f", 24/(tr + tp))])
end
```
