# Fundamental Accelerator Physics (Jupyter Book 2)

Lecture notes for fundamental accelerator physics, built with
[Jupyter Book 2](https://mystmd.org) and simulated with
[TrackPad.jl](https://github.com/MSU-Beam-Dynamics/TrackPad.jl), the MSU Beam
Dynamics group's particle-tracking package.

Reorganized from the earlier `fundamental_AP` notes; chapters:

- **Math Preparation** — linear algebra, relativistic mechanics, ODEs & vector calculus
- **Transverse Dynamics** — multipoles, Hill's equation, transfer matrices, Twiss/FODO,
  emittance, coupling, dispersion, chromaticity
- **Longitudinal Dynamics** — RF acceleration, synchrotron motion, buckets, linac dynamics

## Repository layout

```
myst.yml                 # Jupyter Book 2 / MyST site config + table of contents
content/                 # all pages (MyST markdown)
  images/                # figures shared by chapters
julia_env/               # Julia project used by every code cell
                         #   TrackPad.jl comes from GitHub (see below)
  TrackPadWidgets/       # the one local package: self-contained interactive figures
pyproject.toml           # Python side: jupyter-book via uv
```

## Building the book

```bash
uv sync                      # install jupyter-book (first time only)
PATH="/opt/homebrew/opt/node@22/bin:$PATH" uv run jupyter-book start --execute --execute-parallel 1  # live preview → http://localhost:3000
PATH="/opt/homebrew/opt/node@22/bin:$PATH" uv run jupyter-book build --html --execute --execute-parallel 1  # static HTML → _build/html/index.html
```

Use Homebrew's Node 22 runtime for this project. Node 26 can trigger an
upstream MyST HTTP-client crash after a page has built. Install Node 22 once
if it is not already present:

```bash
brew install node@22
```

## Running the simulation cells

Code cells are Julia ≥ 1.10. The pinned environment lives in `julia_env/` and
uses CairoMakie for static, headless figures. TrackPad.jl is not registered, so
add it straight from its repository:

```bash
julia --project=julia_env -e 'using Pkg; Pkg.add(url="https://github.com/MSU-Beam-Dynamics/TrackPad.jl")'
julia --project=julia_env -e 'using Pkg; Pkg.instantiate()'
```

`Pkg.instantiate()` on its own is enough once the manifest already resolves
TrackPad. To move the book onto a newer TrackPad later:

```bash
julia --project=julia_env -e 'using Pkg; Pkg.update("TrackPad")'
```

Only if you are developing TrackPad itself does a local checkout come into it,
in which case point the environment at your working copy instead:

```bash
julia --project=julia_env -e 'using Pkg; Pkg.develop(path="/path/to/TrackPad.jl")'
```

Install the dedicated Jupyter kernel once. Every page containing Julia cells
selects this kernel, so MyST executes the cells in this pinned environment
rather than whichever Julia project happens to be current:

```bash
julia --project=julia_env -e 'using IJulia; IJulia.installkernel("Julia (TrackPad — Fundamental AP)", "--project=$(pwd())/julia_env")'
```

Confirm the installation (the identifier must be
`julia-_trackpad-_-fundamental-ap_-1.12`):

```bash
uv run jupyter kernelspec list
```

Then build or preview with execution enabled. MyST does *not* execute
`{code-cell}` blocks during an ordinary `jupyter-book start` or
`jupyter-book build`; `--execute` is required. CairoMakie renders the book's
figures headlessly through its static Cairo backend.

To force fresh outputs after changing Julia packages or a calculation:

```bash
uv run jupyter-book clean --execute
```

You can also start a notebook with this environment's kernel interactively:

```bash
julia --project=julia_env -e 'using IJulia; jupyterlab()'
```

Notes on conventions: `Beam(energy)` takes kinetic energy in eV, coordinates
are `(x, px, y, py, z, δE)` with `z = s/β₀ − ct`, and `δE = (E−E₀)/(P₀c)`.
See [`docs/src/conventions.md`](https://github.com/MSU-Beam-Dynamics/TrackPad.jl/blob/main/docs/src/conventions.md)
in the TrackPad repository for the normative reference.

## Interactive figures

Many pages carry a figure with sliders under it, built with
`julia_env/TrackPadWidgets` — the only package that does live in this
repository, `Pkg.develop`ed into `julia_env`, so `Pkg.instantiate()` above
already picks it up. Its `explorer` function
evaluates a TrackPad calculation once for every point of a small parameter
grid at *build* time and embeds the results with a small canvas renderer as
the cell's `text/html` output — the published book stays fully static and
needs no running kernel for the sliders to work. See the docstring
(`using TrackPadWidgets; ?explorer`) or any of the existing widget cells for
the pattern: `Knob` for a slider axis, `Panel` for a set of axes, `line` /
`points` for series, and `statics` for series that do not change between
frames.
