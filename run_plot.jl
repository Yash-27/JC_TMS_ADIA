include("config_jc_pump_disp_asy.jl")
include("plotting_functions.jl")

# ==============================================================================
# plot_try.jl
#
# The file you run. All machinery lives in plotting_lib.jl; this is just the
# list of figures to make.
#
# Everything here is at top level, so `d` and `figs` stay available in the
# session afterwards -- e.g. figs[:tracedist], or
#   plot_map(d, :tracedist; markersize = 10, show_params = false)
# ==============================================================================

const RESULTS_FILE = "results_k1_$(PARAMS.κ_1)_k2_$(PARAMS.κ_2)_g2_$(PARAMS.g_2).jld2"
const FIG_DIR      = "Full_vs_Adia"
const FIG_SUFFIX   = "_k1_$(PARAMS.κ_1)_k2_$(PARAMS.κ_2)_g2_$(PARAMS.g_2)"

# Which maps to draw. Add or remove keys here; unknown keys are skipped with
# a warning rather than erroring, so a half-finished list still runs.
const MAPS_TO_PLOT = (
    :tracedist,
    :concurrence_full,
    :concurrence_adia,
    :concurrence_diff,
)


apply_theme!()

d = load_results(RESULTS_FILE)

println("loaded $(d.fname)")
println("  observables: ", join(sort(string.(collect(keys(d.outs)))), ", "))
println("  grid: $(length(d.β_vals)) β × $(length(d.x_vals)) x")

figs = Dict{Symbol,Any}()
for key in MAPS_TO_PLOT
    if haskey(d.outs, key)
        figs[key] = plot_map(d, key)
        save_fig(figs[key], string(key) * FIG_SUFFIX; dir = FIG_DIR)
    else
        @warn "skipping :$key -- not present in $(d.fname)"
    end
end

println("\ndone -- $(length(figs)) figure(s) in $FIG_DIR/")