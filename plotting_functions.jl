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
# This file DEFINES things and runs nothing. It is included by plot_try.jl,
# which is the file you actually run. Do not include it from anywhere else --
# a single include point is what keeps ordering from becoming a problem.
#
# Hardcoded because there is only one correct answer:
#   - β horizontal on a log scale, x vertical linear
#   - axis labels β and x
#   - the flattening convention (β varies fastest -- see flatten_grid)
#   - NaN points dropped
#   - colorbar on (GR does not infer it from marker_z)
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
# getobs: fetch an observable grid, with a readable error on a typo.
# ------------------------------------------------------------------------
function getobs(d, key::Symbol)
    haskey(d.outs, key) ||
        error("no observable :$key -- available: " *
              join(sort(string.(collect(keys(d.outs)))), ", "))
    return d.outs[key]
end

# ------------------------------------------------------------------------
# flatten_grid: an (N_β × N_x) grid -> three flat vectors, for scatter.
#
# vec() walks column-major, so the ROW index (β) varies fastest. `outer` on
# β and `inner` on x reproduces that ordering exactly; swapping them
# silently transposes the plot with no error.
#
# NaN entries (failed grid points) are dropped, so they show as missing
# markers rather than breaking the colour range.
# ------------------------------------------------------------------------
function flatten_grid(M, β_vals, x_vals)
    N_β, N_x = size(M)
    @assert length(β_vals) == N_β "β_vals length $(length(β_vals)) ≠ grid rows $N_β"
    @assert length(x_vals) == N_x "x_vals length $(length(x_vals)) ≠ grid cols $N_x"

    βflat = repeat(β_vals, outer = N_x)
    xflat = repeat(x_vals, inner = N_β)
    vflat = vec(M)

    keep = findall(isfinite, vflat)
    return βflat[keep], xflat[keep], vflat[keep]
end


# ==============================================================================
# LABELS
# ==============================================================================

fmt_num(v::Real) = @sprintf("%g", v)

# ------------------------------------------------------------------------
# param_string: fixed physical parameters as a LaTeX math fragment for the
# title. γ and γ_φ are omitted when zero -- a title listing five numbers of
# which two are always 0 is mostly noise.
# ------------------------------------------------------------------------
function param_string(p)
    parts = String[]
    push!(parts, "\\kappa_1 = $(fmt_num(p.κ_1))")
    push!(parts, "\\kappa_2 = $(fmt_num(p.κ_2))")
    push!(parts, "g_2 = $(fmt_num(p.g_2))")
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
# plot_map: scatter of one observable over the (β, x) plane, marker colour
# encoding the value.
#
# `key` may be a Symbol naming an observable, or a raw (N_β × N_x) matrix,
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

    βflat, xflat, vflat = flatten_grid(M, d.β_vals, d.x_vals)
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

    scatter(βflat, xflat;
            marker_z          = vflat,
            c                 = cmap,
            clims             = cl,
            colorbar          = true,
            markersize        = ms,
            markerstrokewidth = 0.5,
            markerstrokecolor = :black,
            xscale            = :log10,
            xlabel            = L"\beta",
            ylabel            = L"x",
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