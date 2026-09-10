---
kernelspec:
  name: julia-_trackpad-_-fundamental-ap_-1.12
  display_name: Julia (TrackPad — Fundamental AP)
---

# Emittance

## Defining the emittance

Plotting every particle's location and momentum $(x,x')$ fills an
ellipse-like region of phase space. The area occupied by the ensemble is
measured by the ***emittance***. The ***rms emittance*** is defined as

$$
\epsilon_{\text{rms}}=\sqrt{\sigma_x^2\sigma_{x'}^2-\sigma_{xx'}^2},
$$

where $\sigma_{xx'}$ is the correlation between $x$ and $x'$. For a Gaussian
beam the fraction of particles inside ellipses scaled from $\epsilon_{\text{rms}}$:

| ratio $\epsilon/\epsilon_{\text{rms}}$ | Percentage |
|---|---|
| 2 | 63 % |
| 4 | 86 % |
| 6 | 95 % |

## The $\sigma-$matrix

If we defind the covariance matrix as

$$
\sigma = \begin{bmatrix}
\sigma_x^2 & \sigma_{xx'} \\
\sigma_{xx'} & \sigma_{x'}^2
\end{bmatrix},
$$

then the emittance is simply the square root of the determinant of this matrix.

$$
\epsilon_\text{rms} = \sqrt{\det(\sigma)} = \sqrt{\sigma_x^2\sigma_{x'}^2-\sigma_{xx'}^2}.
$$

If the beam is Gaussian and matched to the lattice, then the $\sigma$-matrix is proportional to the Twiss matrix:

$$
\sigma = \epsilon_\text{rms} \begin{bmatrix}
\beta & -\alpha \\
-\alpha & \gamma
\end{bmatrix}.
$$

However, it is not necessary that the beam is matched to the lattice.  The calculation of the $\sigma$-matrix and the emittance is valid for any distribution of particles in phase space. 

One can prove that the emittance is an invariant of linear transport, meaning that it does not change as the beam propagates through a linear lattice.

However, the if the beam is not matched the lattice, the effective emittance will vary along the lattice.


## Emittance through a FODO cell

The definition earns its keep because of what it does *not* do: a linear,
uncoupled lattice moves particles around phase space but never changes the area
they occupy. Take the same 4 m FODO cell of the previous chapter and watch both
emittances — the rms one read off $\Sigma$, and the 100 % one the lattice has to
hold — run flat along the cell while $\sigma_x$, in the readouts, breathes by
more than a factor of two.

```{code-cell} julia
:tags: [hide-input]

using StaticArrays, TrackPad, TrackPadWidgets, Random, Statistics, LinearAlgebra

beam = Beam(3.0e9)

# The same FODO cell as the previous chapter: 4 m period, ℓq = 0.5 m, f ≈ 1.8 m.
const L_cell, f_eff, ℓq = 4.0, 1.8, 0.5
const kq   = 1/(f_eff*ℓq)
const NPE  = 400            # particles in the cloud
const EPS0 = 200e-9         # geometric emittance at the entrance [m·rad]

# Slice the cell so the optics are sampled inside the magnets too: 2 slices per
# quadrupole and 8 per drift give 20 pieces, hence 21 observation stations.
function fodo_slices()
    p = AbstractElement[]
    append!(p, [Quadrupole(ℓq/2, +kq) for _ in 1:2])
    append!(p, [Drift((L_cell - 2ℓq)/2/8) for _ in 1:8])
    append!(p, [Quadrupole(ℓq/2, -kq) for _ in 1:2])
    append!(p, [Drift((L_cell - 2ℓq)/2/8) for _ in 1:8])
    p
end

ringE = Lattice(fodo_slices(); periodic=true)
twE   = periodic_twiss(ringE, beam)              # β, α at all 21 station boundaries
const βM, αM = twE.betax[1], twE.alphax[1]       # matched values at the cell start

"Track the bunch piece by piece, keeping the full coordinate array at each boundary."
function track_cloud(pieces, coords0)
    c = copy(coords0); flags = zeros(Int, size(c, 1))
    S = Float64[0.0]; C = [copy(c)]
    for e in pieces
        linepass!(c, Lattice(AbstractElement[e]), beam, flags)
        push!(S, S[end] + get_length(e)); push!(C, copy(c))
    end
    S, C
end

"Beam matrix Σ of the cloud and its rms emittance √det Σ."
function sigma_matrix(c)
    x, xp = c[:, 1], c[:, 2]
    Σ = [var(x) cov(x, xp); cov(x, xp) var(xp)]
    Σ, sqrt(max(det(Σ), 0.0))
end

"Single-particle emittances εᵢ = γx² + 2αxx′ + βx′², measured against the LATTICE."
function lattice_emittances(c, β, α)
    γ = (1 + α^2)/β
    [γ*c[i,1]^2 + 2α*c[i,1]*c[i,2] + β*c[i,2]^2 for i in axes(c, 1)]
end

"The 1σ ellipse of the distribution: uᵀΣ⁻¹u = 1, area π·ε_rms."
function sigma_ellipse(Σ; np = 97)
    Lc = cholesky(Symmetric(Σ)).L
    θ  = range(0, 2π, length=np)
    p  = [Lc * [cos(t), sin(t)] for t in θ]
    [q[1] for q in p], [q[2] for q in p]
end

"The lattice ellipse γx² + 2αxx′ + βx′² = ε, oriented by the lattice, not the beam."
function twiss_ellipse(β, α, ε; np = 97)
    θ = range(0, 2π, length=np)
    (@. sqrt(ε*β)*cos(θ)), (@. -sqrt(ε/β)*(α*cos(θ) + sin(θ)))
end

r4e(v) = round.(v; sigdigits=4)
STATIONS = 0:20
# s of each station, so the knob can be labelled in metres rather than by index.
const SGRID = vcat(0.0, cumsum([get_length(e) for e in fodo_slices()]))

# The beamline band across the top of the emittance panel, built from TrackPad's
# backend-independent glyph data. Drawn from the UNSLICED cell so each magnet is
# one box rather than one per integration slice.
const YMAX  = 42000.0        # top-panel range, leaving the band its own room
BEAMLINE = lattice_strip(
    lattice_plot_data(Lattice(AbstractElement[
        Quadrupole(ℓq, +kq), Drift((L_cell - 2ℓq)/2),
        Quadrupole(ℓq, -kq), Drift((L_cell - 2ℓq)/2)]; periodic=true)),
    L_cell; ymax=YMAX)

explorer(
    title   = "Two ellipses: the beam's own, and the one the lattice must hold",
    sliders = [Knob("position in FODO cell", STATIONS;
                    fmt = n -> string(round(SGRID[n+1]; digits=3), " m"), init = 1),
               Knob("β₀ / β_match", [0.5, 1.0, 1.5];
                    fmt = r -> string(round(r; digits=2)), init = 2),
               Knob("α₀ − α_match", [-2.0, -1.0, 0.0, 1.0];
                    fmt = d -> string(round(d; digits=2)), init = 2)],
    # Two axes because ε₁₀₀ runs 12–40× ε_rms; on a shared scale the rms curve
    # would be pinned to the baseline and its flatness impossible to read.
    panels  = [Panel(xlabel="s [m]", ylabel="ε₁₀₀% from the lattice [nm·rad]",
                     y2label="ε_rms from Σ [nm·rad]",
                     title="both emittances along the cell",
                     ylim=(0.0, YMAX), y2lim=(0.0, 400.0), height=270,
                     legend=:bottomleft, basis="100%"),
               Panel(xlabel="x [mm]", ylabel="x′ [mrad]",
                     title="phase space at the selected station",
                     # wide enough for the 100 % ellipse of the most mismatched setting
                     xlim=(-16.0, 16.0), ylim=(-5.0, 5.0), height=320,
                     legend=:bottomleft)],
    statics = BEAMLINE,
    note = "The dashed ellipse is the beam's own 1σ ellipse, read off Σ; the solid one "*
           "is the lattice ellipse inflated to the largest single-particle action in the "*
           "bunch, so it holds every particle. Matched, the two are concentric and "*
           "similar; mismatched, they tilt against each other — and it is the lattice "*
           "ellipse, not the rms one, that has to fit through the aperture. The band "*
           "along the top of the first panel is the cell itself: QF above the line, QD "*
           "below, drift on it.",
) do n, ratio, dα
    c0 = matched_gaussian(MersenneTwister(2024), NPE,
                          optics4DUC(βM*ratio, αM + dα, βM*ratio, αM + dα);
                          emitx=EPS0, emity=EPS0, emitz=1e-9, betaz=0.2)
    S, C = track_cloud(fodo_slices(), c0)

    σx    = [std(c[:, 1]) for c in C] .* 1e3
    εrms  = [sigma_matrix(c)[2] for c in C] .* 1e9
    ε100  = [maximum(lattice_emittances(C[j], twE.betax[j], twE.alphax[j]))
             for j in eachindex(C)] .* 1e9

    j     = n + 1
    Σ, εr = sigma_matrix(C[j])
    βl, αl = twE.betax[j], twE.alphax[j]
    εi    = lattice_emittances(C[j], βl, αl)
    εmax  = maximum(εi)

    ex, ep = sigma_ellipse(Σ)
    tx, tp = twiss_ellipse(βl, αl, εmax)

    (series = [
        line(r4e(S), r4e(ε100); panel=1, color=PALETTE[4], width=2.0,
             label="ε₁₀₀% (left)"),
        line(r4e(S), r4e(εrms); panel=1, color=PALETTE[2], width=2.0, axis=:right,
             label="ε_rms (right)"),
        line([S[j], S[j]], [0.0, YMAX*0.91]; panel=1, color="#9aa4b2", width=1.4),
        points(r4e(C[j][:,1] .* 1e3), r4e(C[j][:,2] .* 1e3); panel=2, color=PALETTE[1],
               size=2.4, alpha=0.4, label="particles"),
        line(r4e(ex .* 1e3), r4e(ep .* 1e3); panel=2, color=PALETTE[2], dash=true,
             width=2.0, label="1σ ellipse from Σ (ε_rms)"),
        line(r4e(tx .* 1e3), r4e(tp .* 1e3); panel=2, color=PALETTE[4],
             width=2.0, label="lattice ellipse at max action (100 %)")],
     readouts = ["s"              => string(round(S[j]; digits=2), " m"),
                 "σₓ"             => string(round(σx[j]; sigdigits=4), " mm"),
                 "ε_rms"          => string(round(εr*1e9; sigdigits=5), " nm·rad"),
                 "ε_100%"         => string(round(εmax*1e9; sigdigits=5), " nm·rad"),
                 "ε_100% / ε_rms" => round(εmax/εr; digits=2),
                 "ε_rms / ε_rms(0)"   => round(εr*1e9/εrms[1]; digits=6),
                 "ε_100% / ε_100%(0)" => round(εmax*1e9/ε100[1]; digits=6),
                 "lattice β, α"   => string(round(βl; digits=3), " m, ", round(αl; digits=3)),
                 "particles held" => string(count(<=(εmax + 1e-18), εi), "/", NPE)])
end
```

Three readings are worth taking.

**Both emittances are invariant.** The two ratios $\epsilon/\epsilon(0)$ hold at
1 to about a part in $10^{4}$ from station to station, and the residue is the
thick-quadrupole integrator rather than physics. That is true of
$\epsilon_{\text{rms}}$, which is a property of the distribution, *and* of
$\epsilon_{100\%}$, which is the largest single-particle action
$\gamma x^2+2\alpha xx'+\beta x'^2$ measured against the lattice — each
particle's action is conserved, so the biggest one is too.

**Mismatch moves the size, not the area.** Turning either knob changes
$\sigma_x$ everywhere — read it off the readouts — and leaves
$\epsilon_{\text{rms}}$ at 200.5 nm·rad in every one of the twelve settings: a
badly matched beam is not a lower-quality beam, it is the same beam in the wrong
shape. What mismatch *does* cost is the 100 % ellipse, whose curve lifts from
2313 to 33538 nm·rad as the beam is thrown further from the periodic solution.

**The two ellipses are answering different questions.** Matched, they are
concentric and similar, and $\epsilon_{100\%}/\epsilon_{\text{rms}}\approx12$ —
just the tail of a Gaussian, since $\epsilon_i/2\epsilon_{\text{rms}}$ is
exponentially distributed and the maximum of $N$ draws grows like $\ln N$.
Mismatched, the beam ellipse tilts away from the lattice ellipse, and the
lattice ellipse has to inflate — up to 167 times the rms area at the far corner
of the two knobs — to enclose the same particles. The rms number is what you quote; the lattice
ellipse is what has to fit through the aperture.
