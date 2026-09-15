# Fundamental Accelerator Physics

These lecture notes introduce the fundamental physics of particle accelerators, which are part of the PHY 862 course.  The notes are built using
[Jupyter Book 2](https://jupyterbook.org) site.  The simulations are powered by the
[TrackPad.jl](https://github.com/MSU-Beam-Dynamics/TrackPad.jl), an accelerator tracking package that supports auto-differentiation, power series analysis and GPU accelerations.

The notes are organized in three parts:

1. **Math Preparation** — linear algebra, relativistic mechanics, differential
   equations, and vector calculus, with numerical examples in Julia.
2. **Transverse Dynamics** — magnets and multipoles, Hill's equation, transfer
   matrices, Courant–Snyder functions, FODO cells, emittance, coupling, dispersion and
   chromaticity. Simulations are included for better illustration using
   TrackPad.jl.
3. **Longitudinal Dynamics in a Ring** — RF acceleration, the longitudinal turn-by-turn
   map, phase stability, synchrotron motion and bucket dynamics.

## Interactive figures

Several sections carry a figure with sliders under it. These are not videos and
they do not need a running kernel: when the book is built, the Julia cell runs
the TrackPad calculation once for **every** position of every slider and stores
the results — together with a small drawing routine — inside the page. Dragging a
slider therefore just redraws a result that TrackPad already computed, and the
`▶ play` button sweeps the first slider automatically.

The machinery lives in the local package `julia_env/TrackPadWidgets`, whose
`explorer` function takes the sliders, the panels, and a function of the slider
values that returns the curves to draw. Every such cell is ordinary Julia: open
the notebook, change the lattice, and the figure changes with it.

## Running the examples

All code cells are written in Julia (≥ 1.10) against the project-local
environment in `julia_env/`. To execute them interactively:

```bash
# one-time: install IJulia kernel bound to this environment
julia --project=julia_env -e 'using Pkg; Pkg.instantiate(); using IJulia; notebook()'
```

`julia_env` already tracks the TrackPad.jl development checkout at
`~/src/juliapkg/TrackPad`. If your checkout lives elsewhere, adjust it with:

```bash
julia --project=julia_env -e 'using Pkg; Pkg.develop(path="/path/to/TrackPad")'
```

To build or preview the book itself (Python side is managed with `uv`):

```bash
uv run jupyter-book start   # live-reload preview at http://localhost:3000
uv run jupyter-book build   # static HTML export into _build/
```


