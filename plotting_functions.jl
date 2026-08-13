using Plots
using JLD2
using LaTeXStrings
using Printf

# ==============================================================================
# plotting_lib.jl
#
# All the machinery for plotting the (β, x) sweep: data access, labels,
# styling, the plot functions themselves, and saving.
#
# This file DEFINES things and runs nothing. It is included by the plot
# drivers -- run_plot.jl, run_plot_r_xi.jl and compare_runs.jl -- which are the
# files you actually run. Do not include it from anywhere else.
#
# Hardcoded because there is only one correct answer:
#   - the ROW coordinate is horizontal on a log scale, the COLUMN coordinate
#     vertical and linear
#   - the flattening convention (the row coordinate varies fastest -- see
#     flatten_grid)
#   - NaN points dropped
#   - colorbar on (GR does not infer it from marker_z)
#
# WHICH coordinates those are is per-file, not hardcoded, because there are now
# two sweeps: (β, x) from run_jc_pump_disp_asy.jl and (r, ξ) from run_r_xi.jl.
# A loaded run carries an `axes` field describing its own axes; anything
# without one -- the older .jld2 files, and the bare (β_vals, x_vals) tuples
# compare_runs.jl builds -- falls back to the original β/x behaviour. See
# axis_spec.
#
# Varies per figure: title, colorbar label, colour limits, marker size,
# colormap. Anything else is forwarded straight to Plots via kwargs.
# ==============================================================================


# ==============================================================================
# DATA ACCESS
# ==============================================================================

# ------------------------------------------------------------------------
# load_results: pull a saved run back in.
#
# Axes come from the file -- the grid that actually produced the data --
# never from config.jl, which may have been edited since the run.
# ------------------------------------------------------------------------
function load_results(fname)
    isfile(fname) || error("no such file: $fname (working directory is $(pwd()))")
    return jldopen(fname, "r") do f
        (outs        = f["outs"],
         labels      = f["labels"],
         β_vals      = f["β_vals"],
         x_vals      = f["x_vals"],
         params      = f["params"],
         status_grid = f["status_grid"],
         fname       = fname)
    end
end

# ------------------------------------------------------------------------
# load_results_rxi: the same, for a run_r_xi.jl sweep.
#
# Separate function rather than a mode flag on load_results because the two
# file formats genuinely have different keys. Reusing the names β_vals/x_vals
# to hold r/ξ would make a .jld2's field names stop describing its contents,
# which is exactly the trap this split avoids.
#
# The returned tuple carries an `axes` field, so plot_map draws it on (r, ξ)
# without any per-call configuration. ROW = r (horizontal, log), COLUMN = ξ
# (vertical, linear) -- matching run_r_xi.jl's grid layout, which in turn
# matches the (β, x) sweep's, so flatten_grid needs no change.
# ------------------------------------------------------------------------
function load_results_rxi(fname)
    isfile(fname) || error("no such file: $fname (working directory is $(pwd()))")
    return jldopen(fname, "r") do f
        r_vals = f["r_vals"]
        ξ_vals = f["ξ_vals"]
        (outs        = f["outs"],
         labels      = f["labels"],
         r_vals      = r_vals,
         ξ_vals      = ξ_vals,
         params      = f["params"],
         status_grid = f["status_grid"],
         fname       = fname,
         axes        = (xvals = r_vals,
                        yvals = ξ_vals,
                        xlab  = L"r = \kappa_1/\kappa_2",
                        ylab  = L"\xi = 4g_1g_2/(\kappa_1\kappa_2)",
                        xsc   = :log10,
                        ysc   = :identity))
    end
end

# ------------------------------------------------------------------------
# axis_spec: what to draw a grid against.
#
# A run loaded by load_results_rxi says so itself via its `axes` field. Older
# (β, x) result files predate the field, and compare_runs.jl hands plot_map a
# bare (β_vals = ..., x_vals = ...) tuple it builds inline, so the fallback
# reproduces the original hardcoded behaviour exactly. Do not make `axes`
# mandatory -- the three stored .jld2 files on disk do not have it.
# ------------------------------------------------------------------------
function axis_spec(d)
    hasproperty(d, :axes) && return d.axes
    hasproperty(d, :β_vals) && hasproperty(d, :x_vals) ||
        error("cannot determine plot axes: expected an `axes` field, or " *
              "β_vals and x_vals; got fields $(propertynames(d))")
    return (xvals = d.β_vals,
            yvals = d.x_vals,
            xlab  = L"\beta",
            ylab  = L"x",
            xsc   = :log10,
            ysc   = :identity)
end

# ------------------------------------------------------------------------
# getobs: fetch an observable grid, with a readable error on a typo.
# ------------------------------------------------------------------------
function getobs(d, key::Symbol)
    haskey(d.outs, key) ||
        error("no observable :$key -- available: " *
              join(sort(string.(collect(keys(d.outs)))), ", "))
    return d.outs[key]
end

# ------------------------------------------------------------------------
# flatten_grid: an (N_row × N_col) grid -> three flat vectors, for scatter.
#
# vec() walks column-major, so the ROW index varies fastest. `outer` on the
# row values and `inner` on the column values reproduces that ordering
# exactly; swapping them silently transposes the plot with no error, and the
# length assertions cannot catch it on a square grid -- which len = 7 on both
# axes always is.
#
# Both sweeps put their horizontal/log coordinate on the rows (β for (β, x),
# r for (r, ξ)), so this function is shared unchanged; only the names are
# generic now.
#
# NaN entries (failed grid points) are dropped, so they show as missing
# markers rather than breaking the colour range.
# ------------------------------------------------------------------------
function flatten_grid(M, row_vals, col_vals)
    N_row, N_col = size(M)
    @assert length(row_vals) == N_row "row values length $(length(row_vals)) ≠ grid rows $N_row"
    @assert length(col_vals) == N_col "column values length $(length(col_vals)) ≠ grid cols $N_col"

    rowflat = repeat(row_vals, outer = N_col)
    colflat = repeat(col_vals, inner = N_row)
    vflat   = vec(M)

    keep = findall(isfinite, vflat)
    return rowflat[keep], colflat[keep], vflat[keep]
end


# ==============================================================================
# LABELS
# ==============================================================================

fmt_num(v::Real) = @sprintf("%g", v)

# ------------------------------------------------------------------------
# param_string: whatever is FIXED for this run, as a LaTeX math fragment for
# the title. γ and γ_φ are omitted when zero -- a title listing five numbers
# of which two are always 0 is mostly noise.
#
# Which parameters are fixed depends on which sweep produced the run, so this
# branches on the shape of `params`:
#
#   (β, x) sweep   κ_1, κ_2, g_2 are held fixed and η, g_1 vary
#   (r, ξ) sweep   x, β, κ_geo are held fixed and ALL FOUR of κ_1, κ_2, g_1,
#                  g_2 vary across the grid -- so printing them would be
#                  actively wrong, not merely uninformative
#
# Detected on :κ_geo, which only the (r, ξ) params tuple has.
# ------------------------------------------------------------------------
function param_string(p)
    parts = String[]
    if hasproperty(p, :κ_geo)
        push!(parts, "x = $(fmt_num(p.x_fixed))")
        push!(parts, "\\beta = $(fmt_num(p.β_fixed))")
        push!(parts, "\\sqrt{\\kappa_1\\kappa_2} = $(fmt_num(p.κ_geo))")
    else
        push!(parts, "\\kappa_1 = $(fmt_num(p.κ_1))")
        push!(parts, "\\kappa_2 = $(fmt_num(p.κ_2))")
        push!(parts, "g_2 = $(fmt_num(p.g_2))")
    end
    p.γ_val     != 0 && push!(parts, "\\gamma = $(fmt_num(p.γ_val))")
    p.γ_phi_val != 0 && push!(parts, "\\gamma_\\phi = $(fmt_num(p.γ_phi_val))")
    return join(parts, ",\\; ")
end

const TITLE_NAME = Dict(
    :tracedist        => "\\mathrm{Trace\\; distance}",
    :concurrence_full => "\\mathrm{Concurrence\\; (full)}",
    :concurrence_adia => "\\mathrm{Concurrence\\; (adiabatic)}",
    :concurrence_diff => "\\mathrm{Concurrence\\; difference}",
    :purity_full      => "\\mathrm{Purity\\; (full)}",
    :purity_adia      => "\\mathrm{Purity\\; (adiabatic)}",
    :purity_diff      => "\\mathrm{Purity\\; difference}",
)

# ------------------------------------------------------------------------
# NOTE on :concurrence_diff -- analyze() computes C(ρ_full) - C(ρ_adia),
# two concurrences subtracted. That is NOT C(ρ_full - ρ_adia): the
# difference of two density matrices has zero trace and is not positive
# semidefinite, so it is not a state and Wootters' construction does not
# apply to it. The label reflects what is actually computed.
# ------------------------------------------------------------------------
const CBAR_LABEL = Dict(
    :tracedist        => "D(\\rho_{\\mathrm{full}}, \\rho_{\\mathrm{adia}})",
    :concurrence_full => "C(\\rho_{\\mathrm{full}})",
    :concurrence_adia => "C(\\rho_{\\mathrm{adia}})",
    :concurrence_diff => "C(\\rho_{\\mathrm{full}}) - C(\\rho_{\\mathrm{adia}})",
    :purity_full      => "\\mathrm{tr}(\\rho_{\\mathrm{full}}^2)",
    :purity_adia      => "\\mathrm{tr}(\\rho_{\\mathrm{adia}}^2)",
    :purity_diff      => "\\mathrm{tr}(\\rho_{\\mathrm{full}}^2) - \\mathrm{tr}(\\rho_{\\mathrm{adia}}^2)",
)

# ------------------------------------------------------------------------
# CLIM_SOURCE: which key's data range sets the colour limits.
#
# All three concurrence panels take their limits from concurrence_full, so
# the three are directly comparable by eye. Keys not listed fall back to
# their own range.
# ------------------------------------------------------------------------
const CLIM_SOURCE = Dict(
    :concurrence_full => :concurrence_full,
    :concurrence_adia => :concurrence_full,
    :concurrence_diff => :concurrence_full,
)


# ==============================================================================
# STYLE
# ==============================================================================

function apply_theme!()
    default(titlefontsize     = 11,
            guidefontsize     = 10,
            tickfontsize      = 9,
            legendfontsize    = 9,
            grid              = true,
            gridalpha         = 0.25,
            framestyle        = :box,
            linewidth         = 2,
            markerstrokewidth = 0.5,
            dpi               = 300)
    return nothing
end

# ------------------------------------------------------------------------
# resolve_clims: look up the source key, take its finite range.
#
# Warns if the plotted data falls outside those limits -- Plots clamps
# out-of-range values to the end colours silently, which would hide a
# negative concurrence difference rather than showing it.
# ------------------------------------------------------------------------
function resolve_clims(d, key::Symbol)
    src = get(CLIM_SOURCE, key, key)
    ref = filter(isfinite, vec(getobs(d, src)))
    isempty(ref) && error("colour-limit source :$src is entirely NaN")
    lo, hi = minimum(ref), maximum(ref)

    if src !== key
        own = filter(isfinite, vec(getobs(d, key)))
        if !isempty(own) && (minimum(own) < lo || maximum(own) > hi)
            @warn ":$key ranges [$(minimum(own)), $(maximum(own))] but is drawn on " *
                  ":$src limits [$lo, $hi]; out-of-range points are clamped to the end " *
                  "colours and will look saturated. Pass clims explicitly to override."
        end
    end

    return (lo, hi)
end

# Markers should scale with grid spacing (~1/len), or a len=21 grid overlaps
# into a solid block. The constant 50 reproduces markersize 7 at len = 7.
default_markersize(N) = clamp(round(Int, 50 / N), 2, 14)


# ==============================================================================
# THE PLOT
# ==============================================================================

# ------------------------------------------------------------------------
# plot_map: scatter of one observable over the sweep plane, marker colour
# encoding the value.
#
# Which plane that is comes from axis_spec(d): (β, x) for a run_plot.jl run
# or a bare compare_runs.jl axis tuple, (r, ξ) for a load_results_rxi one.
#
# `key` may be a Symbol naming an observable, or a raw (N_row × N_col) matrix,
# so N_1_grid and friends can be mapped with no extra code.
#
# Any extra keyword is passed straight to Plots, so per-figure tweaks need
# no changes here.
# ------------------------------------------------------------------------
function plot_map(d, key;
                  clims = nothing,
                  markersize = nothing,
                  cmap = :viridis,
                  title = nothing,
                  colorbar_title = nothing,
                  show_params = true,
                  kwargs...)

    M = key isa Symbol ? getobs(d, key) : key

    ax = axis_spec(d)
    xflat, yflat, vflat = flatten_grid(M, ax.xvals, ax.yvals)
    isempty(vflat) && error("nothing finite to plot")

    cl = if clims !== nothing
        clims
    elseif key isa Symbol
        resolve_clims(d, key)
    else
        (minimum(vflat), maximum(vflat))
    end

    ttl = if title !== nothing
        title
    elseif key isa Symbol
        name = get(TITLE_NAME, key, replace(string(key), "_" => "\\_"))
        show_params ? latexstring(name * "\\quad (" * param_string(d.params) * ")") :
                      latexstring(name)
    else
        ""
    end

    cbar = if colorbar_title !== nothing
        colorbar_title
    elseif key isa Symbol && haskey(CBAR_LABEL, key)
        latexstring(CBAR_LABEL[key])
    elseif key isa Symbol
        get(d.labels, key, string(key))
    else
        ""
    end

    ms = markersize === nothing ? default_markersize(maximum(size(M))) : markersize

    scatter(xflat, yflat;
            marker_z          = vflat,
            c                 = cmap,
            clims             = cl,
            colorbar          = true,
            markersize        = ms,
            markerstrokewidth = 0.5,
            markerstrokecolor = :black,
            xscale            = ax.xsc,
            yscale            = ax.ysc,
            xlabel            = ax.xlab,
            ylabel            = ax.ylab,
            title             = ttl,
            colorbar_title    = cbar,
            legend            = false,
            size              = (760, 520),
            kwargs...)
end

# ------------------------------------------------------------------------
# save_fig: PDF (vector) for the document.
# ------------------------------------------------------------------------
function save_fig(p, name; dir = "figures3", formats = ("pdf",))
    isdir(dir) || mkpath(dir)
    for ext in formats
        path = joinpath(dir, "$(name).$(ext)")
        savefig(p, path)
        println("  wrote $path")
    end
    return p
end