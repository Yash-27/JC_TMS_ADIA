# ==============================================================================
# config_x_xi.jl
#
# Config for the THIRD sweep: the (x, ξ) plane at fixed (β, r). Same role as
# config_jc_pump_disp_asy.jl and config_r_xi.jl, and equally dependency-free --
# it references no function and includes no file, so its position in any
# driver's include order is irrelevant.
#
# The model has five rates (κ_1, κ_2, g_1, g_2, η) and one unphysical overall
# frequency scale, so exactly FOUR dimensionless numbers exist:
#
#   x = 4η² / (κ_1 κ_2)          -- drive strength relative to threshold
#   β = (g_1² κ_2) / (g_2² κ_1)  -- coupling asymmetry between the arms
#   r = κ_1 / κ_2                -- cavity-decay asymmetry
#   ξ = 4 g_1 g_2 / (κ_1 κ_2)    -- coupling strength relative to decay
#
# --- Why this slice, given the other two sweeps already exist --------------
# Sweep 1 sweeps (β, x) at fixed κ_1, κ_2, g_2. Holding those fixed does NOT
# hold ξ fixed: with κ and g_2 pinned, ξ ∝ g_1 ∝ sqrt(β), so sweep 1 CONFOUNDS
# β with ξ and structurally cannot separate them. Sweep 2 sweeps (r, ξ) but
# pins x = 0.3663 and stops at ξ = 1.
#
# So neither existing sweep shows how the elimination error behaves in the
# (drive, coupling) plane at SYMMETRIC arms -- which is the plane the
# elimination's validity actually lives in, since ξ is its own small parameter
# and x is what its 4η² - κ_1κ_2 denominators blow up on. This sweep is that
# plane: r = β = 1, everything else varying.
#
# Two independently known landmarks sit inside this box and are on no existing
# grid:
#
#   * the FULL model's concurrence maximum, C = 0.29676 at x = 0.291,
#     ξ ≈ 2.79, β = r = 1 (files_online/07). Bracketed here by the grid points
#     x ∈ {0.2475, 0.3663} and ξ = 2.9375.
#   * the ADIABATIC entanglement threshold x* = 0.48388826 at β = 1 (root of
#     x³ - 2x² + 9x - 4, files_online/06). Crossed between columns 4 and 5, so
#     concurrence_adia is exactly zero over the right-hand third of the plane.
#
# --- What is pinned, and why it is κ_geo ---------------------------------
# Four dimensionless numbers do not fix five rates; one scale must be pinned.
# Pinning κ_geo = sqrt(κ_1 κ_2), as sweep 2 does:
#
#   κ_1 = κ_geo sqrt(r)              κ_2 = κ_geo / sqrt(r)
#   η   = (κ_geo/2) sqrt(x)
#   g_1 = (κ_geo/2) sqrt(ξ) (β r)^(1/4)
#   g_2 = (κ_geo/2) sqrt(ξ) (β r)^(-1/4)
#
# At r = β = 1 this collapses to κ_1 = κ_2 = κ_geo, η = (κ_geo/2)sqrt(x),
# g_1 = g_2 = (κ_geo/2)sqrt(ξ) -- so only η and g move across the grid, one
# per axis, which is as clean as this model gets.
#
# κ_geo = 2.0 reproduces run A's κ_1 = κ_2 = 2 exactly. That matters here for
# more than tidiness: see the note on ξ_range.
# ==============================================================================

# --- Physical parameters and grid ---------------------------------------
# Field names match run_sweep_xxi's keyword arguments and are splatted straight
# in, so adding a parameter here requires adding it to that signature too.
#
# ξ_range is swept logarithmically, x_range linearly.
const PARAMS_XXI = (
    β_fixed   = 1.0,          # symmetric couplings -- and C_adia's optimum at every x
    r_fixed   = 1.0,          # symmetric decay -- κ_1 = κ_2 = κ_geo
    κ_geo     = 2.0,          # sqrt(κ_1 κ_2), the pinned scale
    γ_val     = 0.0,
    γ_phi_val = 0.0,

    len       = 7,
    ξ_range   = (0.2, 5.0),      # log10; see the note below -- ξ = 1 is ON this grid
    x_range   = (0.01, 0.7225),  # linear; identical to sweep 1's x grid
)

# --- Note on ξ_range = (0.2, 5.0) rather than (0.1, 5.0) -----------------
# The endpoints are geometrically symmetric about 1 (0.2 = 1/5), so at len = 7
# the log-spaced grid puts ξ = 1 EXACTLY at index 4:
#
#     0.2000  0.3415  0.5833  1.0000  1.7145  2.9375  5.0000
#
# and ξ = 1 with β = r = 1 and κ_geo = 2 gives κ_1 = κ_2 = 2, g_1 = g_2 = 1,
# η = sqrt(x) -- which is precisely run A's β = 1 row. All seven x values here
# are on run A's x grid too, so THE ENTIRE ξ = 1 ROW OF THIS SWEEP IS ALREADY
# STORED in results_k1_2.0_k2_2.0_g2_1.0.jld2. That makes it a free 7-point
# regression test of the coordinate inversion, the truncation and the driver at
# once, including at the expensive dim-4096 point -- verified before this run:
#
#     x        0.0100  0.1288  0.2475  0.3663  0.4850  0.6038  0.7225
#     C_adia   0.0939  0.2360  0.2073  0.1176  0.0000  0.0000  0.0000
#     C_full   0.0942  0.2610  0.2866  0.2681  0.2221  0.1541  0.0636
#     D        0.0002  0.0152  0.0483  0.0898  0.1300  0.1605  0.1747
#
# (0.1, 5.0) at len = 7 gives 0.1, 0.192, 0.368, 0.707, 1.357, 2.605, 5.0 --
# no point in common with any stored run, and the cross-check is lost. If you
# widen this axis, keep the endpoints reciprocal (0.1, 10.0) works too, giving
# ξ = 1 at index 4 again.
#
# --- Note on ξ_range not starting at 0 ----------------------------------
# ξ = 0 means g_1 = g_2 = 0. H_adia and every entry of J_adia are then
# identically zero, so L_adia is the zero superoperator and its kernel is the
# ENTIRE 16-dimensional state space; steadystate.eigenvector would return an
# arbitrary member of it without warning rather than failing. Keep the lower
# end strictly positive.
#
# Note also that ξ = 0.2 is NOT small enough for the elimination's own linear
# regime: D ∝ ξ holds only below ξ ≈ 0.01 (measured log-log slopes fall from
# 0.97 to 0.20 across ξ = 0.003 -> 3). This grid deliberately spans the
# SATURATING region, where the interesting behaviour is; do not read a slope
# off it and extrapolate to ξ -> 0.
#
# --- Note on x_range topping out at 0.7225 ------------------------------
# x < 1 is a hard constraint -- the adiabatic model has 4η² - κ_1κ_2 in every
# denominator and diverges at x = 1. 0.7225 = 0.85² is sweep 1's ceiling and is
# kept identical here so the two grids share x values.
#
# COST WARNING: estimate_truncation puts every point of this grid at N = 14
# (dim 900) EXCEPT the x = 0.7225 column, which lands at N = 31 (dim 4096) --
# for every ξ from 0.1 to 5, since the truncation there is set by η, not g.
# That is CLAUDE.md's ~40x-per-solve truncation cliff, 7 of the 49 points, and
# it is the single cost fact governing this run. run_x_xi.jl serializes those
# points behind a real semaphore for that reason. Dropping x_range to
# (0.01, 0.60375) keeps the whole grid at dim 900 and the run at minutes.

# --- Observables to compute ---------------------------------------------
# Same registry as both other sweeps -- every driver includes observables.jl.
const ACTIVE_OBSERVABLES_XXI = [
    :concurrence,
    # :purity,
    :tracedist,
]

# --- Output --------------------------------------------------------------
# Keyed on everything that defines the run, INCLUDING len and both ranges via
# their endpoints, as OUTFILE_RXI is. Sweep 1's OUTFILE keys on three
# parameters only, so two of its runs differing in len or range silently
# overwrite each other (known issue 7); no reason to reproduce that here.
const OUTFILE_XXI = string(
    "results_xxi",
    "_b_",  PARAMS_XXI.β_fixed,
    "_r_",  PARAMS_XXI.r_fixed,
    "_kg_", PARAMS_XXI.κ_geo,
    "_xi_", PARAMS_XXI.ξ_range[1], "-", PARAMS_XXI.ξ_range[2],
    "_x_",  PARAMS_XXI.x_range[1], "-", PARAMS_XXI.x_range[2],
    "_n_",  PARAMS_XXI.len,
    ".jld2")

# Written after every completed point and deleted on success. A single dim-4096
# solve here runs for minutes, so an all-or-nothing write at the end is how 45
# minutes of solves get lost to a killed process (recorded in CLAUDE.md).
const PARTIAL_XXI = replace(OUTFILE_XXI, ".jld2" => "_partial.jld2")
