# ==============================================================================
# config.jl
#
# Everything you change between runs lives here. This file has NO dependencies
# -- it references no functions and includes no other file, so its position in
# the include order is irrelevant.
#
# Observables are selected by NAME. The names come from OBSERVABLE_REGISTRY in
# run_jc_pump_disp_asy.jl; a typo here raises an error listing the valid ones
# rather than silently computing fewer observables than you asked for.
# ==============================================================================

# --- Physical parameters and grid ---------------------------------------
# Field names match run_sweep's keyword arguments and are splatted straight in,
# so adding a parameter here requires adding it to run_sweep's signature too.
#
#   x = 4η² / (κ_1 κ_2)          -- drive strength relative to threshold
#   β = (g_1² κ_2) / (g_2² κ_1)  -- coupling asymmetry between the arms
#
# x_range is swept linearly, β_range logarithmically.
const PARAMS = (
    κ_1       = 2.0,
    κ_2       = 2.0,
    g_2       = 1.0,
    γ_val     = 0.0,
    γ_phi_val = 0.0,

    len       = 7,
    x_range   = (0.01, 0.7225),
    β_range   = (0.01, 100.0),
)

# --- Observables to compute ---------------------------------------------
# :single  observables expand automatically into _full / _adia / _diff.
# :compare observables produce one grid under their own name.
#
# Currently available: :concurrence, :purity, :tracedist
const ACTIVE_OBSERVABLES = [
    :concurrence,
    # :purity,
    :tracedist,
]

# --- Output --------------------------------------------------------------
const OUTFILE = "results.jld2"