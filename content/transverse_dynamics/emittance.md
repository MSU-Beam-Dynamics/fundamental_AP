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

If the beam is Gaussian and matched to the lattice, then the $\sigma$-matrix is proportional to the Courant–Snyder matrix:

$$
\sigma = \epsilon_\text{rms} \begin{bmatrix}
\beta & -\alpha \\
-\alpha & \gamma
\end{bmatrix}.
$$

However, it is not necessary that the beam is matched to the lattice.  In this case, the Courant–Snyder parameters only reflect the shape/orientation of the beam distribution in the phase space.  The emittance calculated is usually refered as the geometric emittance. 

However, the lattice design may have an alternative shape and orientation, charecterized by the lattice Courant–Snyder parameters, as they are only determined by the lattice structure and not by the beam properties. One may calculate another beam emittance, sometimes referred as the effective emittance, to characterize the phase space area needed to accomodate the mismatched beam.

One can prove that both the geometric emittance and the effective emittance are invariants of linear transport, meaning that they do not change as the beam propagates through a linear lattice.




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

# Column layout of a TrackPad coordinate array: (x, pₓ, y, p_y, z, δE).
const IX, IPX = 1, 2

# The same FODO cell as the previous chapter: 4 m period, ℓq = 0.5 m, f ≈ 1.8 m.
const CELL_LENGTH  = 4.0
const FOCAL_LENGTH = 1.8
const QUAD_LENGTH  = 0.5
const QUAD_STRENGTH = 1/(FOCAL_LENGTH * QUAD_LENGTH)

const N_PARTICLES    = 400        # particles in the cloud
const EMITTANCE_IN   = 200e-9     # geometric emittance at the entrance [m·rad]

"""
    fodo_slices()

The cell cut into 20 pieces — 2 per quadrupole, 8 per drift — so the optics and
the cloud are recorded at 21 stations, inside the magnets as well as between
them. Slicing changes only where things are observed, never the optics.
"""
function fodo_slices()
    drift_slice = (CELL_LENGTH - 2*QUAD_LENGTH)/2/8
    return AbstractElement[
        [Quadrupole(QUAD_LENGTH/2, +QUAD_STRENGTH) for _ in 1:2]...,
        [Drift(drift_slice)                        for _ in 1:8]...,
        [Quadrupole(QUAD_LENGTH/2, -QUAD_STRENGTH) for _ in 1:2]...,
        [Drift(drift_slice)                        for _ in 1:8]...,
    ]
end

lattice = Lattice(fodo_slices(); periodic=true)
optics  = periodic_twiss(lattice, beam)        # β, α at all 21 station boundaries

# The matched values at the cell start, which the knobs deliberately depart from.
const BETA_MATCHED  = optics.betax[1]
const ALPHA_MATCHED = optics.alphax[1]

"""
    track_cloud(pieces, coords0)

Track the bunch piece by piece, keeping the full coordinate array at every
boundary. Returns `(s, states)`.
"""
function track_cloud(pieces, coords0)
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

"""
    beam_matrix(state)

The 2×2 beam matrix Σ = [⟨x²⟩ ⟨xx′⟩; ⟨xx′⟩ ⟨x′²⟩] of the cloud, and the rms
emittance √det Σ that goes with it. This is the beam describing ITSELF — no
lattice enters.
"""
function beam_matrix(state)
    x  = state[:, IX]
    xp = state[:, IPX]
    Σ  = [var(x)      cov(x, xp)
          cov(x, xp)  var(xp)]
    return Σ, sqrt(max(det(Σ), 0.0))
end

"""
    single_particle_actions(state, β, α)

εᵢ = γx² + 2αxx′ + βx′² for every particle, measured against the LATTICE
ellipse rather than against the beam's own. The largest of these is the
emittance the aperture actually has to accommodate.
"""
function single_particle_actions(state, β, α)
    γ = (1 + α^2)/β
    return [γ*state[i, IX]^2 + 2α*state[i, IX]*state[i, IPX] + β*state[i, IPX]^2
            for i in axes(state, 1)]
end

"The 1σ ellipse of the distribution itself: the curve uᵀΣ⁻¹u = 1."
function beam_ellipse(Σ; npoints = 97)
    L = cholesky(Symmetric(Σ)).L       # maps the unit circle onto the 1σ ellipse
    θ = range(0, 2π, length=npoints)
    pts = [L * [cos(t), sin(t)] for t in θ]
    return [p[1] for p in pts], [p[2] for p in pts]
end

"The lattice ellipse γx² + 2αxx′ + βx′² = ε, oriented by the lattice, not the beam."
function lattice_ellipse(β, α, ε; npoints = 97)
    θ = range(0, 2π, length=npoints)
    x  = @. sqrt(ε*β) * cos(θ)
    xp = @. -sqrt(ε/β) * (α*cos(θ) + sin(θ))
    return x, xp
end

"Round plot data to 4 significant digits — finer than a screen pixel."
plotdata(v) = round.(v; sigdigits=4)

# s of each station, so the knob can be labelled in metres rather than by index.
const STATION_S = vcat(0.0, cumsum([get_length(e) for e in fodo_slices()]))

# The beamline band across the top of the emittance panel, built from TrackPad's
# backend-independent glyph data. Drawn from the UNSLICED cell, so each magnet is
# one box rather than one box per integration slice.
const PANEL_TOP = 42000.0        # top-panel range, leaving the band its own room
const BEAMLINE = lattice_strip(
    lattice_plot_data(Lattice(AbstractElement[
        Quadrupole(QUAD_LENGTH, +QUAD_STRENGTH), Drift((CELL_LENGTH - 2*QUAD_LENGTH)/2),
        Quadrupole(QUAD_LENGTH, -QUAD_STRENGTH), Drift((CELL_LENGTH - 2*QUAD_LENGTH)/2)];
        periodic=true)),
    CELL_LENGTH; ymax=PANEL_TOP)

explorer(
    title   = "Two ellipses: the beam's own, and the one the lattice must hold",
    sliders = [Knob("position in FODO cell", 0:20;
                    fmt = n -> string(round(STATION_S[n+1]; digits=3), " m"), init = 1),
               Knob("β₀ / β_match", [0.5, 1.0, 1.5];
                    fmt = r -> string(round(r; digits=2)), init = 2),
               Knob("α₀ − α_match", [-2.0, -1.0, 0.0, 1.0];
                    fmt = d -> string(round(d; digits=2)), init = 2)],
    # Two axes because ε₁₀₀ runs 12–40× ε_rms; on a shared scale the rms curve
    # would be pinned to the baseline and its flatness impossible to read.
    panels  = [Panel(xlabel="s [m]", ylabel="ε₁₀₀% from the lattice [nm·rad]",
                     y2label="ε_rms from Σ [nm·rad]",
                     title="both emittances along the cell",
                     ylim=(0.0, PANEL_TOP), y2lim=(0.0, 400.0), height=270,
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
) do station, beta_ratio, alpha_offset
    entrance = optics4DUC(BETA_MATCHED*beta_ratio, ALPHA_MATCHED + alpha_offset,
                          BETA_MATCHED*beta_ratio, ALPHA_MATCHED + alpha_offset)
    bunch = matched_gaussian(MersenneTwister(2024), N_PARTICLES, entrance;
                             emitx=EMITTANCE_IN, emity=EMITTANCE_IN,
                             emitz=1e-9, betaz=0.2)

    s, states = track_cloud(fodo_slices(), bunch)

    # The two emittances along the whole cell. ε_rms is the beam's own area and
    # is invariant; ε₁₀₀ is measured against the lattice and is invariant too,
    # but only because the lattice ellipse rotates with the beam.
    rms_curve      = [beam_matrix(state)[2] for state in states] .* 1e9
    envelope_curve = [maximum(single_particle_actions(states[j],
                                                      optics.betax[j], optics.alphax[j]))
                      for j in eachindex(states)] .* 1e9

    # Everything below refers to the one station the knob selects.
    j = station + 1
    Σ, rms_emittance = beam_matrix(states[j])
    β_here, α_here   = optics.betax[j], optics.alphax[j]
    actions          = single_particle_actions(states[j], β_here, α_here)
    largest_action   = maximum(actions)

    beam_x,    beam_xp    = beam_ellipse(Σ)
    lattice_x, lattice_xp = lattice_ellipse(β_here, α_here, largest_action)

    (series = [
        line(plotdata(s), plotdata(envelope_curve); panel=1, color=PALETTE[4],
             width=2.0, label="ε₁₀₀% (left)"),
        line(plotdata(s), plotdata(rms_curve); panel=1, color=PALETTE[2],
             width=2.0, axis=:right, label="ε_rms (right)"),
        line([s[j], s[j]], [0.0, PANEL_TOP*0.91]; panel=1, color="#9aa4b2", width=1.4),
        points(plotdata(states[j][:, IX]  .* 1e3),
               plotdata(states[j][:, IPX] .* 1e3); panel=2, color=PALETTE[1],
               size=2.4, alpha=0.4, label="particles"),
        line(plotdata(beam_x .* 1e3), plotdata(beam_xp .* 1e3); panel=2,
             color=PALETTE[2], dash=true, width=2.0,
             label="1σ ellipse from Σ (ε_rms)"),
        line(plotdata(lattice_x .* 1e3), plotdata(lattice_xp .* 1e3); panel=2,
             color=PALETTE[4], width=2.0,
             label="lattice ellipse at max action (100 %)")],
     readouts = ["s"              => string(round(s[j]; digits=2), " m"),
                 "σₓ"             => string(round(std(states[j][:, IX])*1e3; sigdigits=4), " mm"),
                 "ε_rms"          => string(round(rms_emittance*1e9; sigdigits=5), " nm·rad"),
                 "ε_100%"         => string(round(largest_action*1e9; sigdigits=5), " nm·rad"),
                 "ε_100% / ε_rms" => round(largest_action/rms_emittance; digits=2),
                 # Both ratios stay at 1: each emittance is separately conserved.
                 "ε_rms / ε_rms(0)"   => round(rms_emittance*1e9/rms_curve[1]; digits=6),
                 "ε_100% / ε_100%(0)" => round(largest_action*1e9/envelope_curve[1]; digits=6),
                 "lattice β, α"   => string(round(β_here; digits=3), " m, ",
                                            round(α_here; digits=3)),
                 "particles held" => string(count(<=(largest_action + 1e-18), actions),
                                            "/", N_PARTICLES)])
end
```

From the above example, we can see that:

**Both emittances are invariant.**, which mismatch presents, the effective emittance is larger than the geometric emittance, and the lattice ellipse has to be inflated to hold the same particles.


**\%100 emittance** only make sense when the beam distribution is truncated.  For a real Gaussian beam, the effective emittance is infinite.  In practice, one may define a certain percentage of the beam distribution to be covered in the lattice ellipse, e.g. 90\%, 95\%, 99\%, etc.  
