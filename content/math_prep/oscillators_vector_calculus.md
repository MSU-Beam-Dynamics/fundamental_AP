---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Differential Equations and Vector Calculus

## Ordinary differential equations

Here we limit the discussion to first- and second-order ordinary differential
equations. Typical examples are:

- First-order homogeneous equation: $\dfrac{dx}{dt}+b(x)=0$
- First-order inhomogeneous equation: $\dfrac{dx}{dt}+b(x)=f(t)$
- Second-order homogeneous equation:
  $\dfrac{d^2x}{dt^2}+a(x)\dfrac{dx}{dt}+b(x)=0$
- Second-order inhomogeneous equation:
  $\dfrac{d^2x}{dt^2}+a(x)\dfrac{dx}{dt}+b(x)=f(t)$

### Harmonic oscillator

The harmonic oscillator satisfies a second-order differential equation,

$$
\frac{d^2x}{dt^2}+\omega^2 x=0,
$$

with solution

$$
x=Ae^{i\omega t}+Be^{-i\omega t}
= A\cos(\omega t+\phi_0),
$$

where the constants are fixed by initial conditions such as $x(0)$ and
$x'(0)$. Every particle executing betatron motion is (to first order) a
harmonic oscillator with $s$ playing the role of time.

```{code-cell} julia
using CairoMakie

ω = 5.0
t = range(0, 10, length=1000)
x = cos.(ω .* t)
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="Time", ylabel="Oscillation amplitude")
lines!(ax, t, x; label="x(t) = cos(ωt)")
axislegend(ax; position=:rt)
fig
```

### Damped harmonic oscillator

If damping exists, the motion obeys

$$
\frac{d^2x}{dt^2}+2\alpha\frac{dx}{dt}+\omega^2 x=0,
\qquad
x=Ae^{\lambda_1 t}+Be^{\lambda_2 t},
\qquad
\lambda_{1/2}=-\alpha\pm i\sqrt{\omega^2-\alpha^2}.
$$

In accelerators we almost always have weak damping only, i.e.
$\omega \gg \alpha$: an oscillation whose amplitude decays like $e^{-\alpha t}$.

```{code-cell} julia
ω, α = 5.0, 0.2
t = range(0, 100, length=10000)
x = @. exp(-α*t) * cos(sqrt(ω^2 - α^2)*t)
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="Time", ylabel="Oscillation amplitude")
lines!(ax, t, x; label="damped oscillator")
lines!(ax, t, exp.(-α .* t); linestyle=:dash, color=(:black, 0.5), label="e^{-αt} envelope")
axislegend(ax; position=:rt)
fig
```

Damping is a genuine one-parameter family, and the qualitative change happens
exactly at $\alpha = \omega$, where the two roots $\lambda_\pm$ collide on the
real axis. Sweep $\alpha$ and follow the roots and the motion together:

```{code-cell} julia
:tags: [hide-input]

using TrackPadWidgets

ω = 5.0                                  # undamped angular frequency [1/s]
t = range(0, 12, length=241)

"""
    root_locus(ω; αmax = 12)

The path both roots λ± = −α ± √(α²−ω²) sweep out as the damping α runs from 0
upwards: a vertical line at Re λ = −α while α < ω (underdamped), then the real
axis once the square root turns real.

Returned as a single polyline with a `nothing` between the two branches, which
lifts the pen so they are not joined by a spurious connecting segment.
"""
function root_locus(ω; αmax = 12)
    re = Union{Float64,Nothing}[]
    im = Union{Float64,Nothing}[]
    for branch in (+1, -1)
        isempty(re) || (push!(re, nothing); push!(im, nothing))
        for α in range(0, αmax, length=200)
            discriminant = α^2 - ω^2
            push!(re, -α)
            push!(im, discriminant < 0 ? branch*sqrt(-discriminant) : 0.0)
        end
    end
    return re, im
end

locus_re, locus_im = root_locus(ω)

"""
    displacement(α, t)

x(t) for ẍ + 2αẋ + ω²x = 0 started from x(0) = 1, ẋ(0) = 0. Underdamped and
overdamped need different closed forms, so they are written separately rather
than carried through complex arithmetic.
"""
function displacement(α, t)
    discriminant = α^2 - ω^2
    if discriminant < 0                              # underdamped: decaying cosine
        ω_damped = sqrt(-discriminant)
        return [exp(-α*τ) * cos(ω_damped*τ) for τ in t]
    else                                             # overdamped: two real decays
        fast = -α + sqrt(discriminant)
        slow = -α - sqrt(discriminant)
        weight = fast == slow ? 0.0 : -slow/(fast - slow)   # fixes x(0)=1, ẋ(0)=0
        return [weight*exp(fast*τ) + (1 - weight)*exp(slow*τ) for τ in t]
    end
end

explorer(
    title   = "Damped oscillator ẍ + 2αẋ + ω²x = 0   (ω = 5)",
    sliders = [Knob("α  [1/s]", range(0.0, 9.0, length=37);
                      fmt = a -> string(round(a; digits=2)), init = 3)],
    panels  = [Panel(xlabel="t", ylabel="x(t)", title="displacement",
                     ylim=(-1.15, 1.15), height=250),
               Panel(xlabel="Re λ", ylabel="Im λ", title="roots λ± = −α ± √(α²−ω²)",
                     xlim=(-12.5, 1.0), ylim=(-6.0, 6.0), height=250,
                     legend=:bottomleft)],
    statics = [line(locus_re, locus_im; panel=2, color="#9aa4b2", dash=true, label="λ(α)"),
               line([-12.5, 1.0], [0, 0]; panel=2, color="#9aa4b2", width=1, alpha=0.5)],
    note = "Accelerators live at the far left of the underdamped branch: α/ω is typically 10⁻⁴.",
) do α
    discriminant = α^2 - ω^2
    underdamped  = discriminant < 0

    # The two roots, as a pair of points to mark on the locus.
    root_re = underdamped ? [-α, -α] :
              [-α + sqrt(discriminant), -α - sqrt(discriminant)]
    root_im = underdamped ? [sqrt(-discriminant), -sqrt(-discriminant)] : [0.0, 0.0]

    regime = discriminant < -1e-9 ? "underdamped" :
             discriminant >  1e-9 ? "overdamped"  : "critical"

    (series = [line(t, displacement(α, t); panel=1, color=PALETTE[1], label="x(t)"),
               line(t, exp.(-α .* t); panel=1, color=PALETTE[2], dash=true,
                    label="e^(−αt)"),
               points(root_re, root_im; panel=2, color=PALETTE[1], size=6.0)],
     readouts = ["regime"    => regime,
                 "α/ω"       => round(α/ω; digits=3),
                 "damped ω′" => underdamped ? string(round(sqrt(-discriminant); digits=3)) : "0",
                 "1/e time"  => α > 0 ? string(round(1/α; digits=3)) : "∞"])
end
```

The root locus on the right is the whole story in one picture: the roots slide
down the vertical line $\mathrm{Re}\,\lambda = -\alpha$ until they meet at
$-\omega$, then separate along the real axis. Betatron and synchrotron motion
sit at the very top of that vertical branch, where the imaginary part — the
oscillation — dominates and the real part only slowly bleeds away the amplitude.

## Vector operations

Accelerator science deals with charged-particle motion in electromagnetic
fields, so Maxwell's equations are the starting point of many problems:

| | Differential form | Integral form |
|:--|:--|:--|
| Gauss's law | $\nabla\cdot\mathbf{E} = \rho/\varepsilon_0$ | $\oint \mathbf{E}\cdot d\mathbf{S}=\frac{1}{\varepsilon_0}\iiint \rho\, dV$ |
| no monopoles | $\nabla\cdot\mathbf{B} = 0$ | $\oint \mathbf{B}\cdot d\mathbf{S}=0$ |
| Faraday's law | $\nabla\times\mathbf{E} = -\partial\mathbf{B}/\partial t$ | $\oint \mathbf{E}\cdot d\mathbf{l}=-\frac{d}{dt}\int \mathbf{B}\cdot d\mathbf{S}$ |
| Ampère's law | $\nabla \times \mathbf{B} = \mu_0\left(\mathbf{J} + \varepsilon_0 \frac{\partial \mathbf{E}} {\partial t}\right)$ | $\oint \mathbf{B}\cdot d\mathbf{l}=\mu_0\left(\int\mathbf{J}\cdot d\mathbf{S}+\varepsilon_0\frac{d}{dt}\int \mathbf{E}\cdot d\mathbf{S}\right)$ |

Therefore vector operations are used constantly.

### Gradient, divergence and curl

**Gradient** of a scalar function is a vector:

$$
\nabla\psi=\left(\frac{\partial}{\partial x}, \frac{\partial}{\partial y}, \frac{\partial}{\partial z}\right)\psi
$$

**Divergence** of a vector is a scalar function:

$$
\nabla\cdot\mathbf{f}=
\frac{\partial f_x}{\partial x}+ \frac{\partial f_y}{\partial y}+ \frac{\partial f_z}{\partial z}
$$

**Curl** of a vector is another vector function:

$$
\nabla\times\mathbf{f} =
\begin{vmatrix} \hat{x} & \hat{y} & \hat{z}\\
{\partial/\partial x} & {\partial/\partial y} & {\partial/\partial z} \\
  f_x & f_y & f_z \end{vmatrix}
$$

### Gauss's divergence theorem

$$
\int_S \mathbf{f}\cdot d\mathbf{S}=\int_V \nabla\cdot\mathbf{f}\,dV
$$

Useful for finding the electric field of a given charge distribution.

### Stokes' theorem

$$
\int_{S} (\nabla \times \mathbf{f})\cdot d\mathbf{S} = \oint_{\partial S} \mathbf{f}\cdot d\mathbf{r}
$$

Useful for finding the magnetic field of a given current distribution.

### Example: gradient of a 2D potential

The electrostatic potential of three parallel line charges is a classic
exercise: compute the field as (minus) the gradient and visualize both.

```{code-cell} julia
using CairoMakie

# Three parallel line charges, given as (x, y, charge per unit length).
charges = [(3.45, -3.38, -0.7),
           (-2.48,  2.96,  1.2),
           (4.47,   3.70, -0.5)]

"Potential of a set of parallel line charges: each contributes −q·ln(r²)/2."
potential(x, y) = sum(-q * log((x - cx)^2 + (y - cy)^2) / 2 for (cx, cy, q) in charges)

xs = range(-5, 5, length=500)
ys = range(-5, 5, length=500)

# Makie's heatmap/contour convention is z[row, column] = z[y, x], so the grids
# are built with y varying down the rows and x across the columns.
grid_x = [x for _ in ys, x in xs]
grid_y = [y for y in ys, _ in xs]
Φ      = potential.(grid_x, grid_y)

"""
    central_difference(Φ, h; dims)

∂Φ/∂x (`dims=2`) or ∂Φ/∂y (`dims=1`) by the central difference
(Φ[i+1] − Φ[i−1]) / 2h. The first and last row/column are left at zero; the
strided sampling used for the arrows never reaches them.
"""
function central_difference(Φ, h; dims)
    d = zero(Φ)
    if dims == 2
        d[:, 2:end-1] .= ((@view Φ[:, 3:end]) .- (@view Φ[:, 1:end-2])) ./ 2h
    else
        d[2:end-1, :] .= ((@view Φ[3:end, :]) .- (@view Φ[1:end-2, :])) ./ 2h
    end
    return d
end

dΦ_dx = central_difference(Φ, step(xs); dims=2)
dΦ_dy = central_difference(Φ, step(ys); dims=1)

# One arrow every 24th grid point; 500² arrows would be a solid black square.
every = 24
sample(a) = vec(a[every:every:end, every:every:end])

fig = Figure(size=(850, 400))

ax_potential = Axis(fig[1, 1]; aspect=DataAspect(), title="Potential Φ")
contour!(ax_potential, xs, ys, Φ'; levels=[-5,-3,-1.5,-1,-0.5,0.5,1,1.5,3,5])

ax_field = Axis(fig[1, 2]; aspect=DataAspect(), title="Field E = −∇Φ")
arrows2d!(ax_field, sample(grid_x), sample(grid_y),
          -sample(dΦ_dx), -sample(dΦ_dy); lengthscale=0.28)

fig
```

The arrows of $\mathbf E=-\nabla\Phi$ cross the equipotential contours at right
angles — the same geometry intuition that later helps reading quadrupole and
sextupole field maps.
