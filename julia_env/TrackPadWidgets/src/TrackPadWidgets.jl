"""
    TrackPadWidgets

Self-contained slider widgets for the *Fundamental Accelerator Physics* notes.

A widget is built by evaluating a Julia function — usually a TrackPad
calculation — once for every point of a small parameter grid, then shipping the
results plus a ~10 kB canvas renderer as the `text/html` output of the code
cell.  The published book therefore stays fully static: dragging a slider only
redraws a frame that was computed when the book was built, and no Jupyter
kernel is needed at read time.

Typical use:

```julia
using TrackPadWidgets

explorer(
    title   = "FODO cell",
    sliders = [Knob("k₁ [m⁻²]", 1.0:0.25:6.0; fmt = k -> string(round(k; digits=2)))],
    panels  = [Panel(xlabel = "s [m]", ylabel = "β [m]")],
) do k
    tw = periodic_twiss(build_ring(k), beam)
    (series   = [line(tw.s, tw.betax; label = "βₓ", color = C[1]),
                 line(tw.s, tw.betay; label = "βᵧ", color = C[2])],
     readouts = ["Qₓ" => round(tw.tunex; digits=4)])
end
```
"""
module TrackPadWidgets

export Widget, Knob, Panel, explorer, line, points, PALETTE, lattice_strip, LATTICE_COLORS

# ---------------------------------------------------------------------------
# minimal JSON writer (keeps julia_env dependency-free)
# ---------------------------------------------------------------------------

_json(io::IO, ::Nothing) = print(io, "null")
_json(io::IO, x::Bool) = print(io, x ? "true" : "false")

function _json(io::IO, x::Real)
    (isnan(x) || isinf(x)) ? print(io, "null") : print(io, x isa Integer ? x : Float64(x))
end

function _json(io::IO, s::AbstractString)
    print(io, '"')
    for c in s
        if c == '"'
            print(io, "\\\"")
        elseif c == '\\'
            print(io, "\\\\")
        elseif c == '\n'
            print(io, "\\n")
        elseif c == '\r'
            print(io, "\\r")
        elseif c == '\t'
            print(io, "\\t")
        elseif c in ('<', '>', '&')
            # keeps a literal "</script>" from ever appearing inside the payload
            print(io, "\\u", string(UInt16(c); base = 16, pad = 4))
        elseif c < ' '
            print(io, "\\u", string(UInt16(c); base = 16, pad = 4))
        else
            print(io, c)
        end
    end
    print(io, '"')
end

_json(io::IO, x::Symbol) = _json(io, String(x))

function _json(io::IO, v::Union{AbstractVector,Tuple})
    print(io, '[')
    for (i, e) in enumerate(v)
        i > 1 && print(io, ',')
        _json(io, e)
    end
    print(io, ']')
end

function _json(io::IO, d::AbstractDict)
    print(io, '{')
    first = true
    for (k, val) in d
        val === nothing && continue          # drop empty keys, the JS defaults them
        first || print(io, ',')
        first = false
        _json(io, string(k))
        print(io, ':')
        _json(io, val)
    end
    print(io, '}')
end

_json(io::IO, p::Pair) = _json(io, [p.first, p.second])

json(x) = sprint(_json, x)

# ---------------------------------------------------------------------------
# plot pieces
# ---------------------------------------------------------------------------

"""
    PALETTE

Six line colours chosen to stay legible on both the light and the dark book
theme (they are drawn on a transparent canvas over the page background).
"""
const PALETTE = ("#2563eb", "#e2711d", "#0f9d58", "#c2255c", "#7c3aed", "#0d9488")

_clean(v) = [x === nothing ? nothing :
             (isfinite(x) ? round(Float64(x); sigdigits = 5) : nothing) for x in v]

function _series(kind, x, y; panel = 1, color = PALETTE[1], label = nothing,
                 dash = false, alpha = nothing, size = nothing, width = nothing,
                 axis = :left, fill = false, fillalpha = nothing)
    axis in (:left, :right) || throw(ArgumentError("axis must be :left or :right"))
    Dict{String,Any}(
        "type" => String(kind), "panel" => panel - 1, "color" => color,
        "label" => label, "dash" => dash ? true : nothing, "alpha" => alpha,
        "size" => size, "width" => width, "axis" => axis === :right ? "right" : nothing,
        "fill" => fill ? true : nothing, "fillAlpha" => fillalpha,
        "x" => _clean(collect(x)), "y" => _clean(collect(y)),
    )
end

"""
    line(x, y; panel=1, color, label, dash=false, alpha, width, axis=:left,
         fill=false, fillalpha=0.25)

A polyline in `panel`.  A `nothing` entry breaks the line (useful for
separatrices and for orbits that leave the plot).

`fill=true` shades the band between the curve and `y = 0`, which is how an
integrand shows the integral it stands for.  Each unbroken run is shaded
separately, so a `nothing` interrupts the shading as well as the stroke.
"""
line(x, y; kw...) = _series(:line, x, y; kw...)

"""
    points(x, y; panel=1, color, label, alpha, size, axis=:left)

A marker cloud in `panel`.  Keep clouds below a few thousand points: every
frame of the grid is embedded in the page.
"""
points(x, y; kw...) = _series(:scatter, x, y; kw...)

"""
Glyph colours by element kind, keyed by the `kind` field of a
`TrackPad.LatticeGlyph`. Anything not listed falls back to `:element`.
"""
const LATTICE_COLORS = Dict{Symbol,String}(
    :quadrupole => "#c2255c",
    :bend       => "#2563eb",
    :sextupole  => "#0f9d58",
    :octupole   => "#0d9488",
    :rf_cavity  => "#e2711d",
    :solenoid   => "#7c3aed",
    :kicker     => "#7c3aed",
    :element    => "#5b6472",
)

"""
    lattice_strip(glyphs, span; panel=1, ymax=1.0, strip=0.09,
                  colors=LATTICE_COLORS, width=1.6, baseline="#9aa4b2")

Turn the `LatticeGlyph` objects returned by `TrackPad.lattice_plot_data` into
widget series: a beamline band across the top `strip` fraction of the panel's
`[0, ymax]` range. Focusing multipoles sit above the baseline and defocusing
ones below, following the sign of each glyph's `height`; bends and other
`centered` glyphs straddle it, and drifts are the baseline itself.

`span` is the length of the line, either a number or an `(s0, s1)` tuple.
The panel's `ylim` must be `(0, ymax)`, and `ymax` should exceed the data by at
least `strip` plus some headroom so the band has room of its own — the usual
choice is `ymax = datamax/(1 - strip - headroom)`.

The series carry no labels, so the beamline never enters the legend. Only the
glyph fields are read, so this stays free of any TrackPad dependency.
"""
function lattice_strip(glyphs, span; panel::Int = 1, ymax::Real = 1.0,
                       strip::Real = 0.09, colors = LATTICE_COLORS,
                       width::Real = 1.6, baseline::AbstractString = "#9aa4b2")
    s0, s1 = span isa Tuple ? span : (zero(span), span)
    base = ymax * (1 - strip/2)
    half = ymax * strip / 2
    out = Any[line([s0, s1], [base, base]; panel, color = baseline, width = 1.2, alpha = 0.7)]
    for g in glyphs
        h = Float64(g.height)
        lo, hi = g.centered ? (base - abs(h)*half/2, base + abs(h)*half/2) :
                 h >= 0     ? (base, base + h*half) : (base + h*half, base)
        col = get(colors, g.kind, get(colors, :element, "#5b6472"))
        push!(out, line([g.plot_start, g.plot_start, g.plot_end, g.plot_end, g.plot_start],
                        [lo, hi, hi, lo, lo]; panel, color = col, width))
    end
    out
end

"""
    Knob(label, values; fmt=string, init=1)

One interactive axis of the parameter grid.  `values` are the physical
parameter values handed to the `explorer` body; `fmt` turns each into the
string shown next to the slider.  `init` is the 1-based index selected when the
page loads.
"""
struct Knob
    label::String
    values::Vector{Any}
    labels::Vector{String}
    init::Int
end

function Knob(label::AbstractString, values; fmt = string, init::Int = 1)
    v = collect(values)
    Knob(String(label), Any[v...], String[string(fmt(x)) for x in v], init)
end

"""
    Panel(; xlabel, ylabel, title, xlim, ylim, xscale=:linear, yscale=:linear,
            equal=false, height=250, legend=:right)

One set of axes.  Limits default to the union of all frames, so the axes never
jump while a slider is dragged.  `equal=true` requests an equal-aspect view.

`autoscale` selects how limits are chosen for an axis left unset: `:all`
(default) takes the union over every frame, so the axes never move while a
slider is dragged; `:frame` rescales to the frame on display, which is what a
widget needs when a knob changes the size of the data by orders of magnitude.
An explicit `xlim`/`ylim` always wins over both.

`ticklabels=false` keeps the grid but drops the numeric tick labels, and gives
the freed margin back to the data — what a small inset panel wants.

`y2label`/`y2lim` add a right-hand axis.  Series sent to it with `axis=:right`
are scaled against `y2lim` (or, unset, against their own extent) instead of
`ylim`, so two quantities in different units can share one set of x values.
Right-axis series never contribute to the left axis' automatic limits.

`share` links panels that must be read against each other: every panel given
the same share name is drawn with one set of limits, wide enough to hold every
point on any of them and recomputed for the frame on display.  It applies only
where a limit is not already fixed, so combine it with `autoscale=:frame` and
leave `xlim`/`ylim` unset.

`basis` and `minwidth` control how a panel shares the row with its siblings:
panels are laid out in a wrapping flex row, so `basis="100%"` puts one on a row
of its own and a small `basis`/`minwidth` pair lets several sit side by side.

`xscale`/`yscale` may be `:linear` (default) or `:log10`.  On a logarithmic
axis pass the *physical* values to `line`/`points` and the *physical* limits to
`xlim`/`ylim` — the widget takes the logarithm itself and labels the ticks as
decades, so plots never have to show `log10(quantity)` as if it were the
quantity.  Non-positive values on a logarithmic axis are dropped, exactly as a
`nothing` entry would be.  `equal` is ignored unless both axes are linear.
"""
Base.@kwdef struct Panel
    xlabel::Union{Nothing,String} = nothing
    ylabel::Union{Nothing,String} = nothing
    title::Union{Nothing,String} = nothing
    xlim::Union{Nothing,Tuple{Float64,Float64}} = nothing
    ylim::Union{Nothing,Tuple{Float64,Float64}} = nothing
    xscale::Symbol = :linear  # :linear or :log10
    yscale::Symbol = :linear  # :linear or :log10
    autoscale::Symbol = :all  # :all (union of frames) or :frame (per frame)
    equal::Bool = false
    height::Int = 250
    legend::Symbol = :right   # :right :left :bottom :bottomleft :bottomright
    ticklabels::Bool = true   # false drops the numeric tick labels (inset panels)
    basis::Union{Nothing,String} = nothing   # CSS flex-basis, e.g. "100%" or "190px"
    minwidth::Union{Nothing,Int} = nothing   # CSS min-width in px
    share::Union{Nothing,String} = nothing   # panels sharing a name share limits
    y2label::Union{Nothing,String} = nothing        # right-hand axis label
    y2lim::Union{Nothing,Tuple{Float64,Float64}} = nothing  # right-hand axis limits
end

# ---------------------------------------------------------------------------
# the widget itself
# ---------------------------------------------------------------------------

const _JS = read(joinpath(@__DIR__, "widget.js"), String)

"""
    Widget

Wraps the generated markup.  Displaying it from a code cell writes a
`text/html` output, which the book renders in place.
"""
struct Widget
    html::String
end

Base.show(io::IO, ::MIME"text/html", w::Widget) = print(io, w.html)
Base.show(io::IO, ::MIME"text/plain", ::Widget) =
    print(io, "TrackPadWidgets.Widget (interactive; renders as HTML)")

function _limits(frames, statics, ipanel, which; log::Bool = false)
    lo, hi = Inf, -Inf
    keep(v) = v !== nothing && isfinite(v) && (!log || v > 0)
    right(s) = which == "y" && get(s, "axis", nothing) == "right"
    for s in statics
        s["panel"] == ipanel || continue
        right(s) && continue
        for v in s[which]
            keep(v) || continue
            lo = min(lo, v); hi = max(hi, v)
        end
    end
    for f in frames, s in f["series"]
        s["panel"] == ipanel || continue
        right(s) && continue
        for v in s[which]
            keep(v) || continue
            lo = min(lo, v); hi = max(hi, v)
        end
    end
    (isfinite(lo) && isfinite(hi)) || return log ? (1.0, 10.0) : (0.0, 1.0)
    if log
        l, h = log10(lo), log10(hi)
        pad = h - l < 1e-12 ? 0.5 : 0.05 * (h - l)
        return (10.0^(l - pad), 10.0^(h + pad))
    end
    if hi - lo < 1e-12 * max(1.0, abs(hi))
        pad = max(abs(hi) * 0.05, 1.0)
        return (lo - pad, hi + pad)
    end
    pad = 0.05 * (hi - lo)
    (lo - pad, hi + pad)
end

"""
    explorer(body; sliders, panels, title="", note="", play=true, interval=90)

Evaluate `body(values...)` once for every point of the grid spanned by
`sliders` and package the results as an interactive figure.

`statics` holds series that are identical in every frame (a reference curve, a
unit circle, an axis marker); they are stored once instead of once per frame,
which is usually what keeps a widget small.

`body` receives one physical value per slider and returns either a vector of
series or a named tuple `(; series, readouts)`, where `readouts` is a list of
`"name" => value` pairs shown as chips under the plot.

The number of grid points times the number of points per frame is what ends up
in the built HTML, so prefer ~20–60 slider steps and decimated curves.
"""
function explorer(body; sliders::Vector{Knob}, panels::Vector{Panel},
                  statics::Vector = Any[], title::AbstractString = "",
                  note::AbstractString = "", play::Bool = true, interval::Int = 90)
    isempty(sliders) && throw(ArgumentError("explorer needs at least one Knob"))
    dims = Tuple(length(s.values) for s in sliders)

    frames = Vector{Dict{String,Any}}(undef, prod(dims))
    for (k, I) in enumerate(CartesianIndices(dims))     # first slider varies fastest
        out = body((sliders[j].values[I[j]] for j in 1:length(sliders))...)
        ser, reads = out isa NamedTuple ?
                     (collect(out.series), collect(get(out, :readouts, Pair[]))) :
                     (collect(out), Pair[])
        frames[k] = Dict{String,Any}(
            "series" => ser,
            "readouts" => [[string(p.first), string(p.second)] for p in reads],
        )
    end

    pspec = map(enumerate(panels)) do (i, p)
        Dict{String,Any}(
            "xlabel" => p.xlabel, "ylabel" => p.ylabel, "title" => p.title,
            # A missing xlim/ylim key tells the renderer to rescale per frame.
            "xlim" => p.xlim !== nothing ? collect(p.xlim) :
                      p.autoscale === :frame ? nothing :
                      collect(_limits(frames, statics, i - 1, "x"; log = p.xscale === :log10)),
            "ylim" => p.ylim !== nothing ? collect(p.ylim) :
                      p.autoscale === :frame ? nothing :
                      collect(_limits(frames, statics, i - 1, "y"; log = p.yscale === :log10)),
            "xscale" => p.xscale === :log10 ? "log10" : nothing,
            "yscale" => p.yscale === :log10 ? "log10" : nothing,
            "equal" => p.equal ? true : nothing,
            "height" => p.height, "legend" => String(p.legend),
            "ticklabels" => p.ticklabels ? nothing : false,
            "basis" => p.basis, "minwidth" => p.minwidth, "share" => p.share,
            "y2label" => p.y2label,
            "y2lim" => p.y2lim !== nothing ? collect(p.y2lim) : nothing,
        )
    end

    spec = Dict{String,Any}(
        "title" => isempty(title) ? nothing : title,
        "note" => isempty(note) ? nothing : note,
        "play" => play ? nothing : false,
        "interval" => interval,
        "sliders" => [Dict{String,Any}("label" => s.label, "values" => s.labels,
                                       "init" => clamp(s.init, 1, length(s.values)) - 1)
                      for s in sliders],
        "panels" => pspec,
        "statics" => isempty(statics) ? nothing : collect(statics),
        "frames" => frames,
    )

    payload = json(spec)
    id = "tpw-" * string(hash(payload); base = 36)
    Widget(string("<div id=\"", id, "\"></div>\n<script>\n", _JS,
                  "\nwindow.__TPW__.mount(\"", id, "\", ", payload, ");\n</script>"))
end

end # module
