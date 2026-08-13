include("config_r_xi.jl")
include("plotting_functions.jl")

# ==============================================================================
# run_plot_r_xi.jl
#
# The file you run to draw a run_r_xi.jl sweep. All machinery lives in
# plotting_functions.jl; this is just the list of figures to make.
#
#   julia --project=. run_plot_r_xi.jl
#
# Everything is at top level, so `d` and `figs` stay available in the session
# afterwards -- e.g. figs[:tracedist], or
#   plot_map(d, :tracedist; markersize = 10, show_params = false)
#   plot_map(d, d.outs[:concurrence_full])   # or any raw (N_r × N_ξ) matrix
#
# As with run_plot.jl, config_r_xi.jl is included ONLY to build the filename.
# Axes, parameters in the title and everything else come out of the .jld2.
# Practical consequence: editing the config to try new parameters breaks
# plotting until you re-run the sweep, because the filename it resolves to
# will not exist yet.
# ==============================================================================

const RESULTS_FILE_RXI = OUTFILE_RXI
const FIG_DIR_RXI      = "Full_vs_Adia_rxi"
const FIG_SUFFIX_RXI   = string("_rxi_x_", PARAMS_RXI.x_fixed,
                                "_b_",     PARAMS_RXI.β_fixed,
                                "_kg_",    PARAMS_RXI.κ_geo)

# Which maps to draw. Unknown keys are skipped with a warning rather than
# erroring, so a half-finished list still runs.
const MAPS_TO_PLOT_RXI = (
    :tracedist,
    :concurrence_full,
    :concurrence_adia,
    :concurrence_diff,
)

# ------------------------------------------------------------------------
# Per-key overrides.
#
# :concurrence_diff is a SIGNED quantity, and the default CLIM_SOURCE sends it
# to :concurrence_full's range, which is non-negative -- so every negative
# value clamps to the bottom viridis colour, indistinguishable from zero
# (known issue 2). The documented workaround is a symmetric range and a
# diverging colormap, applied here at the call site.
#
# :concurrence_adia gets the same treatment for a different reason: this sweep
# predicts it is CONSTANT to machine precision, so its own data range is ~1e-15
# wide and drawing it on concurrence_full's range would render a uniform block
# that cannot distinguish "identically flat" from "merely small". Its own
# range shows the numerical noise floor, which is the informative thing.
# ------------------------------------------------------------------------
function overrides(d, key)
    v = filter(isfinite, vec(get(d.outs, key, Float64[])))
    isempty(v) && return NamedTuple()

    if key === :concurrence_diff
        m = maximum(abs, v)
        m == 0 && return NamedTuple()
        return (clims = (-m, m), cmap = :balance)

    elseif key === :concurrence_adia
        # Its OWN range, stated explicitly. Passing clims = nothing would not
        # do this -- that is plot_map's default, and it routes through
        # resolve_clims, which CLIM_SOURCE sends straight back to
        # :concurrence_full. The range has to be computed here.
        return (clims = flat_safe_clims(v),)
    end
    return NamedTuple()
end

# ------------------------------------------------------------------------
# flat_safe_clims: a colour range GR can actually draw, for a panel that is
# expected to be CONSTANT.
#
# This sweep predicts concurrence_adia is flat to machine precision, so its
# own data range is ~1e-15 wide on a value of ~0.1. GR cannot build a
# colorbar that narrow: it emits
#
#     GKS: Possible loss of precision in routine SET_WINDOW
#     GKS: Rectangle definition is invalid in routine CELLARRAY
#
# and drops the bar. An equality test (lo == hi) does not catch this -- the
# values differ, just in the 15th digit -- so the guard has to be RELATIVE.
#
# When the panel is flat, widen to (0, 2v) -- a band anchored at zero with the
# constant at its midpoint. The markers still render as one uniform colour,
# which is the honest depiction of a constant field; the widening only gives
# the colorbar readable ticks, and cannot be mistaken for structure because
# there is no colour variation to read.
#
# (0, 2v) rather than a narrow band around v for a SECOND reason, also GR's.
# The colorbar title sits at a fixed offset and the tick labels grow rightward
# into it (the collision si_scale exists to work around in compare_runs.jl).
# A 5% band around 0.1175 makes GR label ticks "0.1225" -- six characters,
# which overprints the title. Anchoring at zero widens the range enough that
# GR picks round values ("0.05", "0.10") that fit. Zero is also the natural
# anchor for a concurrence, which is non-negative by construction.
#
# A panel with genuine structure keeps its exact range, unchanged.
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

d = load_results_rxi(RESULTS_FILE_RXI)

println("loaded $(d.fname)")
println("  observables: ", join(sort(string.(collect(keys(d.outs)))), ", "))
println("  grid: $(length(d.r_vals)) r × $(length(d.ξ_vals)) ξ")
println("  fixed: x = $(d.params.x_fixed), β = $(d.params.β_fixed), " *
        "κ_geo = $(d.params.κ_geo)")

figs = Dict{Symbol,Any}()
for key in MAPS_TO_PLOT_RXI
    if haskey(d.outs, key)
        figs[key] = plot_map(d, key; overrides(d, key)...)
        save_fig(figs[key], string(key) * FIG_SUFFIX_RXI; dir = FIG_DIR_RXI)
    else
        @warn "skipping :$key -- not present in $(d.fname)"
    end
end

println("\ndone -- $(length(figs)) figure(s) in $FIG_DIR_RXI/")
