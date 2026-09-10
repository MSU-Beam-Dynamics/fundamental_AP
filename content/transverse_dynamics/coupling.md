---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Transverse Coupling

## The solenoid

So far the two transverse planes evolved independently. Magnets that mix them
exist — the most common is the **solenoid**, whose longitudinal field
$B_\parallel$ focuses *both* planes while rotating the transverse phase space.

```{figure} ../images/Solenoid_field.png
:width: 400px
:name: fig:solenoid
Solenoid field lines.
```

For a hard-edge solenoid of length $L_0$, define the rotation angle and the
relative solenoid strength:

$$
\theta_s=\frac{qB_{\parallel}L_0}{2P_0}\equiv g_s L_0,
\qquad
g_s\equiv \frac{qB_{\parallel}}{2P_0}.
$$

The linear map couples all four coordinates $(x, x', y, y')$ is given by:

$$
\begin{align}
\begin{pmatrix}
x \\
x' \\
y \\
y'
\end{pmatrix}
&=M_\text{sol}
\begin{pmatrix}
x \\
x' \\
y \\
y'
\end{pmatrix} \\
&=
\left[\begin{matrix}\cos^{2}{\left(\theta_{s}\right)} & \frac{\sin{\left(2\theta_{s}\right)}}{2g_{s}} & -\frac{\sin{\left(2\theta_{s}\right)}}{2} & -\frac{\sin^{2}{\left(\theta_{s}\right)}}{g_{s}}\\
-\frac{g_{s}\sin{\left(2\theta_{s}\right)}}{2} & \cos^{2}{\left(\theta_{s}\right)} & g_{s}\sin^{2}{\left(\theta_{s}\right)} & -\frac{\sin{\left(2\theta_{s}\right)}}{2}\\
\frac{\sin{\left(2\theta_{s}\right)}}{2} & \frac{\sin^{2}{\left(\theta_{s}\right)}}{g_{s}} & \cos^{2}{\left(\theta_{s}\right)} & \frac{\sin{\left(2\theta_{s}\right)}}{2g_{s}}\\
-g_{s}\sin^{2}{\left(\theta_{s}\right)} & \frac{\sin{\left(2\theta_{s}\right)}}{2} & -\frac{g_{s}\sin{\left(2\theta_{s}\right)}}{2} & \cos^{2}{\left(\theta_{s}\right)}
\end{matrix}\right]
\begin{pmatrix}
x \\
x' \\
y \\
y'
\end{pmatrix}
\end{align}
$$

The solenoid map can be factorizes as
$M_{sol}=M_{\text{focus}}\,M_{\text{rotation}}\,M_{\text{focus}}$, i.e. two
equal-strength focusing lenses acting on both planes separated by a rotation:

$$
\begin{equation}
M_\text{focus}=\left[\begin{matrix}\cos{\left(\frac{\theta_{s}}{2}\right)} & \frac{1}{g_{s}}\sin{\left(\frac{\theta_{s}}{2}\right)} & 0 & 0\\
-g_{s}\sin{\left(\frac{\theta_{s}}{2}\right)} & \cos{\left(\frac{\theta_{s}}{2}\right)} & 0 & 0\\
0 & 0 & \cos{\left(\frac{\theta_{s}}{2}\right)} & \frac{1}{g_{s}}\sin{\left(\frac{\theta_{s}}{2}\right)}\\
0 & 0 & -g_{s}\sin{\left(\frac{\theta_{s}}{2}\right)} & \cos{\left(\frac{\theta_{s}}{2}\right)}
\end{matrix}\right]
\end{equation}
$$

, and

$$
\begin{equation}
M_\text{rotation}=\left[\begin{matrix}\cos{\left(\theta_s \right)} & 0 & - \sin{\left(\theta_s \right)} & 0\\0 & \cos{\left(\theta_s \right)} & 0 & - \sin{\left(\theta_s \right)}\\\sin{\left(\theta_s \right)} & 0 & \cos{\left(\theta_s \right)} & 0\\0 & \sin{\left(\theta_s \right)} & 0 & \cos{\left(\theta_s \right)}\end{matrix}\right]
\end{equation}
$$

The following example illustrates the transverse beam trajactory in a solenoid of strength $g_s=2$. The slider knob controls the solenoid length. The transverse axes are exaggerated to show the helical motion of the particles.

```{code-cell} julia
:tags: [hide-input]

using StaticArrays, TrackPad, TrackPadWidgets, LinearAlgebra

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
    S, reduce(hcat, X), reduce(hcat, Y)
end

const GS    = 2.0                        # g_s = qB∥/(2P₀) [m⁻¹], held fixed
const L0MAX = 2.0                       # longest solenoid on the knob [m]
const NSL   = 24                         # slices, so we see the orbit inside the magnet
L0S = collect(range(0.1, L0MAX, length=30))
x0  = collect(range(0.5e-3, 2.5e-3, length=5))

"The hard-edge M_sol of the text, with K = g_s — what the tracking must reproduce."
function Msol(L, K)
    C, S = cos(K*L), sin(K*L)
    [ C^2      S*C/K    -S*C     -S^2/K ;
     -K*S*C    C^2     K*S^2   -S*C   ;
     S*C     S^2/K    C^2     S*C/K ;
      -K*S^2   S*C     -K*S*C   C^2   ]
end

# ---- a hand-rolled 3-D view -------------------------------------------------
# The widget draws 2-D polylines, so the (s, x, y) curves are projected here and
# handed over as ordinary lines. World axes are (s, x, y); the camera is an
# orthonormal (right, up) pair at azimuth `a` and elevation `e`.
const SSPAN = 6.0            # the LONGEST solenoid is drawn SSPAN units long
const TBOX  = 3.0            # half-width of the reference box, in mm

function project(s, x, y, a, e)
    u =    -sin(a)*(s * SSPAN/L0MAX) +  cos(a)*(x * 1e3)
    v = -cos(a)*sin(e)*(s * SSPAN/L0MAX) - sin(a)*sin(e)*(x * 1e3) + cos(e)*(y * 1e3)
    u, v
end

"The wireframe box around a solenoid of length L, as one polyline with pen lifts."
function box_lines(L, a, e)
    c = [(-TBOX, -TBOX), (TBOX, -TBOX), (TBOX, TBOX), (-TBOX, TBOX)]
    us = Union{Float64,Nothing}[]; vs = Union{Float64,Nothing}[]
    add!(s, x, y) = (p = project(s, x*1e-3, y*1e-3, a, e); push!(us, p[1]); push!(vs, p[2]))
    lift!()       = (push!(us, nothing); push!(vs, nothing))
    for sface in (0.0, L)                        # the two end frames
        for (x, y) in c; add!(sface, x, y); end
        add!(sface, c[1][1], c[1][2]); lift!()
    end
    for (x, y) in c                              # the four long edges
        add!(0.0, x, y); add!(L, x, y); lift!()
    end
    us, vs
end

# Azimuth matters: at −90° the camera looks straight down the x axis, so x and y
# land on the same screen direction and the plot collapses to a flat ribbon. At
# −55° the three world axes are 46–81° apart on screen, which is what reads as 3-D.
const AZ0  = -55.0
const ELEV = deg2rad(25.0)
VIEWS = [-30.0, -15.0, 0.0, 15.0]                       # azimuth offsets from AZ0

explorer(
    title   = "Particles traveling through a solenoid ",
    sliders = [Knob("L₀ [m]", L0S; fmt = L -> string(round(L; digits=2)), init = 11),
               Knob("3-D view [deg]", VIEWS; fmt = v -> string(round(Int, v)), init = 2)],
    panels  = [Panel(xlabel="x [mm]", ylabel="y [mm]", title="looking down the beamline",
                     equal=true, height=330, basis="33%", minwidth=230),
               Panel(title="the same rays in (s, x, y) — transverse scale exaggerated",
                     equal=true, ticklabels=false, height=330, basis="63%",
                     minwidth=330, legend=:bottomleft)],
    note = "The knob decides how far traveling in the solenoid. Particles enter with only x offset and initial bended by the edge field. ",
) do L0, view
    ks     = -2GS                                  # TrackPad's ks is 2gₛ
    pieces = AbstractElement[Solenoid(L0/NSL, ks) for _ in 1:NSL]
    c0 = zeros(length(x0), 6); c0[:, 1] .= x0
    S, X, Y = track_s(pieces, beam, c0)
    θs = GS*L0

    Mtr = Matrix(transfer_map(Lattice(pieces), beam))[1:4, 1:4]
    err = opnorm(Mtr - Msol(L0, GS))

    a = deg2rad(AZ0 + view)
    bx, by = box_lines(L0, a, ELEV)
    proj_ray(j) = begin
        p = [project(S[i], X[j, i], Y[j, i], a, ELEV) for i in eachindex(S)]
        [q[1] for q in p], [q[2] for q in p]
    end
    axis3 = let p0 = project(0.0, 0.0, 0.0, a, ELEV), p1 = project(L0, 0.0, 0.0, a, ELEV)
        ([p0[1], p1[1]], [p0[2], p1[2]])
    end

    (series = vcat(
        [line(X[j, :].*1e3, Y[j, :].*1e3; panel=1, color=PALETTE[1], alpha=0.9, width=1.5)
         for j in eachindex(x0)],
        [points(X[:, end].*1e3, Y[:, end].*1e3; panel=1, color=PALETTE[4], size=4.0, label="exit")],
        [line(bx, by; panel=2, color="#9aa4b2", width=1.0, alpha=0.55, label="solenoid, 6 mm bore"),
         line(axis3...; panel=2, color="#9aa4b2", width=1.0, dash=true, label="s axis")],
        [let (u, v) = proj_ray(j)
             line(u, v; panel=2, color=PALETTE[1], width=1.6, alpha=0.95,
                  label = j == 1 ? "trajectories" : nothing)
         end for j in eachindex(x0)],
        [let pe = [project(L0, X[j, end], Y[j, end], a, ELEV) for j in eachindex(x0)]
             points([q[1] for q in pe], [q[2] for q in pe]; panel=2, color=PALETTE[4],
                    size=4.0, label="exit")
         end]),
     readouts = ["gₛ (fixed)"        => string(GS, " m⁻¹"),
                 "L₀"                => string(round(L0; digits=3), " m"),
                 "θₛ = gₛL₀"         => string(round(rad2deg(θs); digits=1), "°"),
                 # Rotating the exit plane back by θₛ must leave the motion planar.
                 # atan(y, x) would flip by 180° once a ray focuses through the axis,
                 # so this is the check that survives the whole knob.
                 "max |y| in Larmor frame" =>
                     string(round(maximum(abs, @. sin(θs)*X[:,end] + cos(θs)*Y[:,end]) * 1e3;
                                  sigdigits=3), " mm"),
                 "‖M_track − M_sol‖" => string(round(err; sigdigits=3)),
                 ])
end
```
