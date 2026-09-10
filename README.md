# Fundamental Accelerator Physics (Jupyter Book 2)

Accelerator Physics lecture notes for fundamental accelerator course (PHY 862 at Michigan State University), built with
[Jupyter Book 2](https://mystmd.org) and simulated with
[TrackPad.jl](https://github.com/MSU-Beam-Dynamics/TrackPad.jl), the MSU Beam
Dynamics group's particle-tracking package.

Reorganized from the earlier `fundamental_AP` notes; chapters:

- **Math Preparation** — linear algebra, relativistic mechanics, ODEs & vector calculus
- **Transverse Dynamics** — multipoles, Hill's equation, transfer matrices, C-S parameters, FODO cell, emittance, coupling, dispersion, chromaticity
- **Longitudinal Dynamics** — RF acceleration, synchrotron motion

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
uv run jupyter-book start --execute --execute-parallel 1  # live preview → http://localhost:3000
uv run jupyter-book build --html --execute --execute-parallel 1  # static HTML → _build/html/index.html
```

Use Node 22 runtime for this project. Node 26 can trigger an
upstream MyST HTTP-client crash after a page has built.  E.g. on MacOS with `brew` you may:

```bash
brew install node@22
```

## Running the simulation cells

Code cells are Julia ≥ 1.10. The pinned environment lives in `julia_env/` and
uses CairoMakie for static, headless figures. TrackPad.jl is not registered yet, so
add it straight from its repository:

```bash
julia --project=julia_env -e 'using Pkg; Pkg.add(url="https://github.com/MSU-Beam-Dynamics/TrackPad.jl")'
julia --project=julia_env -e 'using Pkg; Pkg.instantiate()'
```

