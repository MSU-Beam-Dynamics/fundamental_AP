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

ω = 5.0
t = range(0, 12, length=241)

locus = let                       # where λ±(α) can sit, for all α ≥ 0
    xs = Union{Float64,Nothing}[]; ys = Union{Float64,Nothing}[]
    for branch in (1, -1)
        # a `nothing` lifts the pen, so the two branches are not joined by a
        # spurious segment running back from (−12, 0) to (0, −ω)
        isempty(xs) || (push!(xs, nothing); push!(ys, nothing))
        for a in range(0, 12, length=200)
            d = a^2 - ω^2
            push!(xs, -a); push!(ys, d < 0 ? branch*sqrt(-d) : 0.0)
        end
    end
    line(xs, ys; panel=2, color="#9aa4b2", dash=true, label="λ(α)")
end

explorer(
    title   = "Damped oscillator ẍ + 2αẋ + ω²x = 0   (ω = 5)",
    sliders = [Knob("α  [1/s]", range(0.0, 9.0, length=37);
                      fmt = a -> string(round(a; digits=2)), init = 3)],
    panels  = [Panel(xlabel="t", ylabel="x(t)", title="displacement",
                     ylim=(-1.15, 1.15), height=250),
               Panel(xlabel="Re λ", ylabel="Im λ", title="roots λ± = −α ± √(α²−ω²)",
                     xlim=(-12.5, 1.0), ylim=(-6.0, 6.0), height=250, legend=:bottomleft)],
    statics = [locus, line([-12.5, 1.0], [0, 0]; panel=2, color="#9aa4b2", width=1, alpha=0.5)],
    note = "Accelerators live at the far left of the underdamped branch: α/ω is typically 10⁻⁴.",
) do α
    d = α^2 - ω^2
    x = if d < 0
        [exp(-α*τ)*cos(sqrt(-d)*τ) for τ in t]
    else
        r1, r2 = -α + sqrt(d), -α - sqrt(d)
        A = r1 == r2 ? 0.0 : -r2/(r1 - r2)
        [A*exp(r1*τ) + (1 - A)*exp(r2*τ) for τ in t]
    end
    (series = [line(t, x; panel=1, color=PALETTE[1], label="x(t)"),
               line(t, exp.(-α .* t); panel=1, color=PALETTE[2], dash=true, label="e^(−αt)"),
               points(d < 0 ? [-α, -α] : [-α + sqrt(d), -α - sqrt(d)],
                       d < 0 ? [sqrt(-d), -sqrt(-d)] : [0.0, 0.0];
                       panel=2, color=PALETTE[1], size=6.0)],
     readouts = ["regime"    => d < -1e-9 ? "underdamped" : (d > 1e-9 ? "overdamped" : "critical"),
                 "α/ω"       => round(α/ω; digits=3),
                 "damped ω′" => d < 0 ? string(round(sqrt(-d); digits=3)) : "0",
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

# parallel line charges (x, y, strength)
charges = [(3.45, -3.38, -0.7), (-2.48, 2.96, 1.2), (4.47, 3.70, -0.5)]
potential(x, y) = sum(-c * log((x-cx)^2 + (y-cy)^2) / 2 for (cx, cy, c) in charges)

xs = range(-5, 5, length=500)
ys = range(-5, 5, length=500)
X = [x for _ in ys, x in xs]        # grid matrices use the z[y, x] convention
Y = [y for y in ys, _ in xs]
Φ  = potential.(X, Y)
hₓ, h_y = step(xs), step(ys)

∂Φ∂x = zero(Φ); ∂Φ∂y = zero(Φ)
∂Φ∂x[:, 2:end-1] .= (@view Φ[:, 3:end]) .- (@view Φ[:, 1:end-2])
∂Φ∂y[2:end-1, :] .= (@view Φ[3:end, :]) .- (@view Φ[1:end-2, :])
∂Φ∂x ./= 2hₓ; ∂Φ∂y ./= 2h_y
# boundary entries stay zero; the strided sampling below never touches them

skip = 24
sel(a) = vec(a[skip:skip:end, skip:skip:end])

fig = Figure(size=(850, 400))
ax1 = Axis(fig[1, 1]; aspect=DataAspect(), title="Potential Φ")
contour!(ax1, xs, ys, Φ'; levels=[-5,-3,-1.5,-1,-0.5,0.5,1,1.5,3,5])
ax2 = Axis(fig[1, 2]; aspect=DataAspect(), title="Field E = −∇Φ")
arrows2d!(ax2, sel(X), sel(Y), -sel(∂Φ∂x), -sel(∂Φ∂y); lengthscale=0.28)
fig
```

The arrows of $\mathbf E=-\nabla\Phi$ cross the equipotential contours at right
angles — the same geometry intuition that later helps reading quadrupole and
sextupole field maps.
