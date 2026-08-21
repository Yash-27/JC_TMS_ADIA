include("config_x_xi.jl")
include("plotting_functions.jl")

# ==============================================================================
# run_plot_x_xi.jl
#
# The file you run to draw a run_x_xi.jl sweep -- the (x, ξ) plane at fixed
# (β, r). All machinery lives in plotting_functions.jl; this is just the list of
# figures to make.
#
#   julia --project=. run_plot_x_xi.jl
#
# Everything is at top level, so `d` and `figs` stay available in the session
# afterwards -- e.g. figs[:tracedist], or
#   plot_map(d, :tracedist; markersize = 10, show_params = false)
#   plot_map(d, d.outs[:concurrence_full])   # or any raw (N_ξ × N_x) matrix
#
# As with the other two plot drivers, config_x_xi.jl is included ONLY to build
# the filename. Axes, parameters in the title and everything else come out of
# the .jld2. Practical consequence: editing the config to try new parameters
# breaks plotting until you re-run the sweep, because the filename it resolves
# to will not exist yet.
#
# Do not include this file and run_plot_r_xi.jl in the SAME session: both define
# `overrides(d, key)` and `flat_safe_clims(v)` at top level with identical
# signatures, so the second include silently replaces the first's methods. The
# two copies are currently identical, so nothing breaks today -- but that is a
# coincidence, not a guarantee. Run one driver per session, or lift both helpers
# into plotting_functions.jl (the clean fix for issue 2).
# ==============================================================================

const RESULTS_FILE_XXI = OUTFILE_XXI
const FIG_DIR_XXI      = "Full_vs_Adia_xxi"
const FIG_SUFFIX_XXI   = string("_xxi_b_", PARAMS_XXI.β_fixed,
                                "_r_",     PARAMS_XXI.r_fixed,
                                "_kg_",    PARAMS_XXI.κ_geo)

# Which maps to draw. Unknown keys are skipped with a warning rather than
# erroring, so a half-finished list still runs.
const MAPS_TO_PLOT_XXI = (
    :tracedist,
    :concurrence_full,
    :concurrence_adia,
    :concurrence_diff,
)

# ------------------------------------------------------------------------
# Per-key overrides.
#
# :concurrence_diff is a SIGNED quantity, and the default CLIM_SOURCE sends it
# to :concurrence_full's range, which is non-negative -- so every negative value
# clamps to the bottom viridis colour, indistinguishable from zero (known issue
# 2). The documented workaround is a symmetric range and a diverging colormap,
# applied here at the call site, as run_plot_r_xi.jl does.
#
# :concurrence_adia is NOT flat on this plane -- it is constant down each column
# but varies with x, from 0.236 at x = 0.129 to exactly 0 for x ≥ 0.485 (the
# adiabatic entanglement threshold x* = 0.48389). So flat_safe_clims will pass
# it straight through here. It is kept anyway: it costs nothing, and it is what
# keeps the panel drawable if this sweep is ever narrowed to a single x column,
# where the panel WOULD be constant and GR would drop the colorbar.
# ------------------------------------------------------------------------
function overrides(d, key)
    v = filter(isfinite, vec(get(d.outs, key, Float64[])))
    isempty(v) && return NamedTuple()

    # if key === :concurrence_diff
    #     m = maximum(abs, v)
    #     m == 0 && return NamedTuple()
    #     return (clims = (-m, m), cmap = :balance)

    # elseif key === :concurrence_adia
    #     # Its OWN range, stated explicitly. Passing clims = nothing would not
    #     # do this -- that is plot_map's default, and it routes through
    #     # resolve_clims, which CLIM_SOURCE sends straight back to
    #     # :concurrence_full. The range has to be computed here.
    #     return (clims = flat_safe_clims(v),)
    # end
    return NamedTuple()
end

# ------------------------------------------------------------------------
# flat_safe_clims: a colour range GR can actually draw, for a panel that may be
# CONSTANT.
#
# A data range ~1e-15 wide on a value of ~0.1 makes GR emit
#
#     GKS: Possible loss of precision in routine SET_WINDOW
#     GKS: Rectangle definition is invalid in routine CELLARRAY
#
# and drop the colorbar. An equality test (lo == hi) does not catch it -- the
# values differ, just in the 15th digit -- so the guard has to be RELATIVE.
#
# When the panel is flat, widen to (0, 2v): a band anchored at zero with the
# constant at its midpoint. The markers still render as one uniform colour,
# which is the honest depiction of a constant field; the widening only gives the
# colorbar readable ticks. Zero is also the natural anchor for a concurrence,
# and it keeps GR's tick labels short enough not to overprint the colorbar
# title -- the same collision si_scale works around in compare_runs.jl.
#
# A panel with genuine structure keeps its exact range, unchanged.
#
# Duplicated from run_plot_r_xi.jl rather than shared: plotting_functions.jl is
# the only file both plot drivers include, and lifting these two helpers into it
# is the clean fix for issue 2 across all three drivers at once. Until that
# happens, the copies must stay in step.
# ------------------------------------------------------------------------
function flat_safe_clims(v; rtol = 1e-9)
    lo, hi = extrema(v)
    scale  = max(abs(lo), abs(hi), 1e-12)
    hi - lo >= rtol * scale && return (lo, hi)      # real structure: leave alone

    # Flat. hi ≈ lo ≈ the constant.
    hi <= 0 && return (-1e-9, 1e-9)                 # flat at (or below) zero
    return (0.0, 2hi)
end


apply_theme!()

d = load_results_xxi(RESULTS_FILE_XXI)

println("loaded $(d.fname)")
println("  observables: ", join(sort(string.(collect(keys(d.outs)))), ", "))
println("  grid: $(length(d.ξ_vals)) ξ × $(length(d.x_vals)) x")
println("  fixed: β = $(d.params.β_fixed), r = $(d.params.r_fixed), " *
        "κ_geo = $(d.params.κ_geo)")

figs = Dict{Symbol,Any}()
for key in MAPS_TO_PLOT_XXI
    if haskey(d.outs, key)
        figs[key] = plot_map(d, key; overrides(d, key)...)
        save_fig(figs[key], string(key) * FIG_SUFFIX_XXI; dir = FIG_DIR_XXI)
    else
        @warn "skipping :$key -- not present in $(d.fname)"
    end
end

println("\ndone -- $(length(figs)) figure(s) in $FIG_DIR_XXI/")
