# ==============================================================================
# config_r_xi.jl
#
# Config for the SECOND sweep: the (r, ξ) plane at fixed (x, β). Same role as
# config_jc_pump_disp_asy.jl, and equally dependency-free.
#
# The model has five rates (κ_1, κ_2, g_1, g_2, η) and one unphysical overall
# frequency scale, so exactly FOUR dimensionless numbers exist:
#
#   x = 4η² / (κ_1 κ_2)          -- drive strength relative to threshold
#   β = (g_1² κ_2) / (g_2² κ_1)  -- coupling asymmetry between the arms
#   r = κ_1 / κ_2                -- cavity-decay asymmetry
#   ξ = 4 g_1 g_2 / (κ_1 κ_2)    -- coupling strength relative to decay
#
# (x, β) is what run_jc_pump_disp_asy.jl sweeps. (r, ξ) is the exact
# complement: together the four are a complete set, so this sweep covers the
# directions the first one holds fixed, with no overlap and no gaps.
#
# ξ is the adiabatic elimination's OWN small parameter -- it is the ratio the
# elimination expands in -- so the ξ axis is the physically informative one.
#
# --- Which scale is pinned, and why it is κ_geo -------------------------
# Four dimensionless numbers do not fix five rates; one scale must be pinned.
# Pinning κ_geo = sqrt(κ_1 κ_2) rather than κ_2 gives
#
#   κ_1 = κ_geo sqrt(r)              κ_2 = κ_geo / sqrt(r)
#   η   = (κ_geo/2) sqrt(x)                     <- free of BOTH r and ξ
#   g_1 = (κ_geo/2) sqrt(ξ) (β r)^(1/4)
#   g_2 = (κ_geo/2) sqrt(ξ) (β r)^(-1/4)
#
# so κ_1 κ_2 and η are constant across the whole grid and only the ratio
# moves. Pinning κ_2 instead puts a sqrt(r) into η and breaks the r <-> 1/r
# symmetry of the solve cost -- the Liouvillian gap at r = 0.1 comes out ~3x
# smaller than at r = 10 purely because κ_1 shrank, which is a choice of units
# rather than physics. With κ_geo pinned the cost is symmetric in log r.
#
# κ_geo = 2.0 is chosen so that r = 1 reproduces run A's κ_1 = κ_2 = 2 exactly.
# ==============================================================================

# --- Physical parameters and grid ---------------------------------------
# Field names match run_sweep_rxi's keyword arguments and are splatted straight
# in, so adding a parameter here requires adding it to that signature too.
#
# r_range is swept logarithmically, ξ_range linearly.
const PARAMS_RXI = (
    x_fixed   = 0.3663,      # see note below
    β_fixed   = 0.5,         # symmetric arms: g_1/g_2 = sqrt(r), r is the only asymmetry
    κ_geo     = 2.0,         # sqrt(κ_1 κ_2), the pinned scale
    γ_val     = 0.0,
    γ_phi_val = 0.0,

    len       = 7,
    r_range   = (0.1, 10.0),  # log10; len = 7 over two decades puts r = 1 exactly at index 4
    ξ_range   = (0.1, 1.0),   # linear
)

# --- Note on x_fixed = 0.3663 -------------------------------------------
# Two constraints pick this value, and they pull in opposite directions.
#
# It is a point of the existing (β, x) sweep's x grid (linear, 0.01 -> 0.7225,
# len 7 -- so 0.3663 is the 4th). Combined with β_fixed = 1.0, which is
# likewise exactly on the existing β grid, the point (r = 1, ξ = 1) of THIS
# sweep sits on exactly the same physical parameters as run A's (β = 1,
# x = 0.3663) point, making results_k1_2.0_k2_2.0_g2_1.0.jld2 a free
# regression test. (Verified at x = 0.485: the adiabatic states agreed to the
# last bit and the full states to 1.7e-8, the iterative solver's tolerance.)
#
# The upper limit is physical, not numerical. At β = 1 the ADIABATIC model's
# concurrence goes to exactly zero for x ≳ 0.4, measured on run A's grid:
#
#     x        0.1288  0.2475  0.3663  0.4850  0.6038  0.7225
#     C_adia   0.2360  0.2073  0.1175  0.0000  0.0000  0.0000
#     C_full   0.2610  0.2866  0.2681  0.2221  0.1541  0.0636
#
# Above 0.4 the concurrence_adia panel is uniformly zero, concurrence_diff
# becomes an exact duplicate of concurrence_full, and the flatness self-test
# in run_r_xi.jl degenerates into checking that 0 stays 0. 0.3663 is the
# largest grid value that keeps all four panels independent and the self-test
# meaningful -- it pins a NONZERO constant to ~4e-15 across the whole plane.
#
# That x ≳ 0.4 region is interesting in its own right -- the elimination fails
# qualitatively there, predicting a separable state where the full model still
# has C = 0.22 -- but it is a separate run. OUTFILE_RXI keys on x, so setting
# x_fixed = 0.485 here lands beside this run rather than overwriting it.
#
# --- Note on ξ_range starting at 0.1, not 0 -----------------------------
# ξ = 0 means g_1 = g_2 = 0. H_adia and every entry of J_adia are then
# identically zero, so L_adia is the zero superoperator and its kernel is the
# ENTIRE 16-dimensional state space. steadystate.eigenvector would return an
# arbitrary member of it without warning rather than failing. The full model
# is only slightly better behaved. Keep the lower end strictly positive.
#
# If the elimination error turns out to collapse fast at small ξ, switch this
# axis to log spacing -- one line in run_r_xi.jl -- to resolve the ξ -> 0
# limit, where the elimination becomes exact.

# --- Observables to compute ---------------------------------------------
# Same registry as the (β, x) sweep -- both drivers include observables.jl.
const ACTIVE_OBSERVABLES_RXI = [
    :concurrence,
    # :purity,
    :tracedist,
]

# --- Output --------------------------------------------------------------
# Keyed on everything that defines the run, INCLUDING len and both ranges via
# their endpoints. The (β, x) sweep's OUTFILE keys on three parameters only,
# so two of its runs differing in len or range silently overwrite each other
# (known issue 7); no reason to reproduce that here.
const OUTFILE_RXI = string(
    "results_rxi",
    "_x_",  PARAMS_RXI.x_fixed,
    "_b_",  PARAMS_RXI.β_fixed,
    "_kg_", PARAMS_RXI.κ_geo,
    "_r_",  PARAMS_RXI.r_range[1], "-", PARAMS_RXI.r_range[2],
    "_xi_", PARAMS_RXI.ξ_range[1], "-", PARAMS_RXI.ξ_range[2],
    "_n_",  PARAMS_RXI.len,
    ".jld2")
