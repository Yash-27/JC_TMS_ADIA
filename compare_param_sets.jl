using LinearAlgebra
using Printf

include("truncation.jl")
include("jc_pump_disp_asy.jl")
include("observables.jl")

# ==============================================================================
# compare_param_sets.jl
#
#   julia --project=. -t auto compare_param_sets.jl
#
# Solve the full model for a list of raw parameter sets and print the pairwise
# trace distance between their rho_full.
#
# Edit SETS below. That is the whole interface -- everything else is reporting.
#
# For each set it prints the four dimensionless invariants
#
#     x = 4η²/(κ₁κ₂)      β = (g₁²κ₂)/(g₂²κ₁)
#     r = κ₁/κ₂           ξ = 4g₁g₂/(κ₁κ₂)
#
# so you can read off which pairs sit at the same point of the 4D space. Pairs
# that do should show D ≈ 0; pairs that do not should not.
# ==============================================================================

BLAS.set_num_threads(1)

# ------------------------------------------------------------------------
# THE PARAMETER SETS. Add, remove or edit rows freely.
#
# The three below are only a starting point, chosen to show both outcomes at
# once: "A" and "A x17.3" are the same physical point (all five rates scaled
# together, so all four invariants are unchanged) and must give D ≈ 0, while
# "C" is a genuinely different point and must not.
#
# η must satisfy η < sqrt(κ₁κ₂)/2, i.e. x < 1 -- above that there is no steady
# state at all and the solve will fail or return garbage. Checked below.
# ------------------------------------------------------------------------
# from_invariants: (x, β, r, ξ) -> the five rates, at a chosen overall scale.
# Same inversion run_r_xi.jl uses, with κ_geo = sqrt(κ₁κ₂) as the pinned scale.
#
# Use this rather than typing rounded decimals when a row is meant to sit at a
# SPECIFIC dimensionless point. Hardcoding 8-digit values for the row below put
# β at 1.00039 instead of 1, and that 4e-4 error alone lifted its D_adia from
# ~1e-16 to 7e-6 -- six orders above the floor, easily misread as a small real
# effect. Let the arithmetic be exact.
function from_invariants(x, β, r, ξ; κ_geo = 2.0)
    amp = (κ_geo / 2) * sqrt(ξ)
    q   = (β * r)^0.25
    return (g_1 = amp * q, g_2 = amp / q,
            κ_1 = κ_geo * sqrt(r), κ_2 = κ_geo / sqrt(r),
            η   = (κ_geo / 2) * sqrt(x))
end

const SETS = [
    (name = "A",        g_1 = 1.0,  g_2 = 1.0,  κ_1 = 2.0, κ_2 = 2.0, η = 0.6964194138592059),
    (name = "A x17.3",  g_1 = 17.3, g_2 = 17.3, κ_1 = 34.6, κ_2 = 34.6, η = 12.048055859764262),
    merge((name = "A same xβ",), from_invariants(0.485, 1.0, 3.0, 0.6)),
    (name = "B",        g_1 = 0.75, g_2 = 0.75, κ_1 = 2.5, κ_2 = 1.5, η = 0.6),
    (name = "C",        g_1 = 1.4,  g_2 = 0.6,  κ_1 = 3.0, κ_2 = 1.0, η = 0.5),
]

# "A same xβ" is the row that makes the two tables say different things. It
# carries A's (x, β) = (0.485, 1) with (r, ξ) moved from (1, 1) to (3, 0.6).
#
# Against A it must show D_adia ≈ 0 -- the adiabatic model sees only x and β --
# but D_full ≠ 0, because the full model sees all four. That single pair of
# entries is the whole 4-vs-2 claim.
#
# Note what separates it from the "A x17.3" row. The scale test is a theorem
# (L -> λL preserves ker L), so its zero was guaranteed before it ran. Nothing
# forces the full model to move HERE: that it does, and by how much, is a
# measured result.

# Fock truncation. `nothing` calls estimate_truncation for each set, which is
# what the sweeps do. Set an integer to force the same small cutoff on every
# set instead -- faster, and worth doing when you only care about whether two
# states agree rather than about their converged values.
const FORCE_N = nothing

# γ and γ_φ. Note that turning either on adds TWO more dimensionless numbers
# (γ/scale, γ_φ/scale), so the four invariants printed below stop being a
# complete description of the state.
const γ_VAL     = 0.0
const γ_PHI_VAL = 0.0


invariants(s) = (x = 4 * s.η^2 / (s.κ_1 * s.κ_2),
                 β = (s.g_1^2 * s.κ_2) / (s.g_2^2 * s.κ_1),
                 r = s.κ_1 / s.κ_2,
                 ξ = 4 * s.g_1 * s.g_2 / (s.κ_1 * s.κ_2))

# ==============================================================================
# SOLVE
# ==============================================================================
println("Solving $(length(SETS)) parameter set(s) on $(Threads.nthreads()) thread(s).\n")

@printf("%-10s %8s %8s %8s %8s %9s | %8s %9s %8s %9s | %7s %6s %8s\n",
        "set", "g_1", "g_2", "κ_1", "κ_2", "η",
        "x", "β", "r", "ξ", "N", "dim", "time")
println("-"^116)

rho_full = Vector{Matrix{ComplexF64}}(undef, length(SETS))
rho_adia = Vector{Matrix{ComplexF64}}(undef, length(SETS))
ok       = falses(length(SETS))

for (i, s) in enumerate(SETS)
    inv = invariants(s)

    if !(0 < inv.x < 1)
        @warn "set $(s.name): x = $(inv.x) is not in (0,1) -- η must be below " *
              "sqrt(κ_1 κ_2)/2 or no steady state exists. Skipping."
        continue
    end

    N_1, N_2 = if FORCE_N === nothing
        try
            estimate_truncation(s.κ_1, s.κ_2, s.η, s.g_1, s.g_2)
        catch e
            @warn "set $(s.name): truncation failed, skipping" exception = e
            continue
        end
    else
        (FORCE_N, FORCE_N)
    end

    local ρf, ρa, dt
    try
        dt = @elapsed ((ρf, ρa) = run_sim(s.g_1, s.g_2, s.κ_1, s.κ_2, s.η,
                                          γ_VAL, γ_PHI_VAL, N_1, N_2))
    catch e
        @warn "set $(s.name): solve failed, skipping" exception = e
        continue
    end

    rho_full[i] = Matrix{ComplexF64}(ρf.data)
    rho_adia[i] = Matrix{ComplexF64}(ρa.data)
    ok[i] = true

    @printf("%-10s %8.4f %8.4f %8.4f %8.4f %9.4f | %8.5f %9.5f %8.5f %9.5f | %3d,%-3d %6d %7.2fs\n",
            s.name, s.g_1, s.g_2, s.κ_1, s.κ_2, s.η,
            inv.x, inv.β, inv.r, inv.ξ,
            N_1, N_2, (N_1 + 1) * (N_2 + 1) * 4, dt)
end

live = findall(ok)
length(live) < 2 && error("need at least two successful sets to compare; got $(length(live))")

# ==============================================================================
# PAIRWISE TRACE DISTANCE
#
# calc_tracedist (observables.jl) normalizes both arguments first. That matters:
# steadystate.iterative imposes tr(rho) = 1 as an appended row of its linear
# system, so it holds only to the solver residual, and an unnormalized
# comparison would report that defect as physical disagreement.
# ==============================================================================
function print_matrix(title, states, idx, note)
    println("\n", title)
    println(note)
    @printf("%-10s", "")
    for j in idx; @printf("%12s", SETS[j].name); end
    println()
    for i in idx
        @printf("%-10s", SETS[i].name)
        for j in idx
            @printf("%12.3e", i == j ? 0.0 : calc_tracedist(states[i], states[j]))
        end
        println()
    end
end

print_matrix("D(rho_full) -- pairwise trace distance, full model", rho_full, live,
             "Rows/columns are sets. Compare against the ~1e-8 floor noted below.")

print_matrix("D(rho_adia) -- pairwise trace distance, adiabatic model", rho_adia, live,
             "The adiabatic model depends only on (x, β), so any two sets sharing\n" *
             "those two agree here even when their r and ξ differ.")

# ------------------------------------------------------------------------
# The elimination error, per set. This is the ONLY number here that compares
# the two models to each other; every entry of the two matrices above is a
# distance WITHIN one model, between two parameter sets.
#
# It is the scale against which those matrices should be read. If it is larger
# than the between-set distances, then the gap between the models dominates
# anything the parameter choice does, and the two matrices are being compared
# in a regime where the adiabatic model is not a good approximation to begin
# with. This is the :tracedist observable of the sweeps, evaluated pointwise.
# ------------------------------------------------------------------------
println("\nD(rho_full, rho_adia) -- the elimination error, one set at a time")
println("The only cross-model number here. Read the two matrices above against it.")
@printf("%-12s %12s   %s\n", "set", "D(full,adia)", "ξ  (elimination is exact as ξ -> 0)")
for i in live
    @printf("%-12s %12.3e   %.4f\n", SETS[i].name,
            calc_tracedist(rho_full[i], rho_adia[i]), invariants(SETS[i]).ξ)
end

# ==============================================================================
# HOW TO READ IT
# ==============================================================================
println("""

Reading the D(rho_full) table
-----------------------------
steadystate.iterative is not deterministic: it seeds a shadow vector from the
unseeded global RNG and stops at reltol = sqrt(eps) ~ 1.5e-8 (relative to
‖L rho_0‖). Solving the SAME parameters twice already differs by ~1e-8, so
anything at or below ~1e-7 in this table is zero. Only entries orders of
magnitude above that are physical disagreement.

Two sets agree in rho_full if and only if all FOUR of (x, β, r, ξ) match --
that is what "same point of the 4D space" means. Because the map from the five
rates to (x, β, r, ξ) plus an overall scale is a bijection, the only way two
DIFFERENT raw sets can share all four is if one is a uniform rescaling of the
other by some λ. That is a theorem (L -> λL preserves ker L), not a
measurement, so a zero here confirms the parametrisation, not the physics.

Avoid λ = 2 or 0.5 if you build such a pair: halving is exact in binary and can
give a bit-identical Liouvillian, so D = 0 exactly -- which is indistinguishable
from a test that never varied anything. Use something like 17.3.
""")
