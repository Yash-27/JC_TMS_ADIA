using Base.Threads
using Printf
using LinearAlgebra
using JLD2
using Statistics

include("config_r_xi.jl")
include("truncation.jl")
include("jc_pump_disp_asy.jl")
include("observables.jl")

# ==============================================================================
# run_r_xi.jl
#
# The SECOND sweep driver. Same shape as run_jc_pump_disp_asy.jl -- top-level
# begin block at the bottom, so `include`ing this file RUNS the sweep -- but
# over the (r, ξ) plane at fixed (x, β) instead of the (β, x) plane at fixed
# (κ_1, κ_2, g_2).
#
# Run it with:   julia --project=. -t auto run_r_xi.jl
#
# What is shared and what is not:
#
#   shared   jc_pump_disp_asy.jl (run_sim), truncation.jl
#            (estimate_truncation), observables.jl (the whole registry).
#            All three already take physical parameters point by point, so
#            none of them knows or cares which coordinates are being swept.
#
#   its own  the sweep loop below. It duplicates the two-pass scheduler from
#            run_jc_pump_disp_asy.jl rather than sharing it, because the
#            coordinate inversion and the cost model both differ, and because
#            including that file to borrow code would run its sweep.
#
# ==============================================================================
# WHAT THIS SWEEP IS EXPECTED TO SHOW -- and why that makes it a self-test
#
# Substituting the config's inversion into jc_pump_disp_asy.jl's adiabatic
# model, every r-dependence cancels ANALYTICALLY:
#
#   H_adia   ∝ η g_1 g_2 / (η² - κ_1κ_2/4).  η, g_1g_2 = ξκ_geo²/4 and
#              κ_1κ_2 = κ_geo² are all r-free, so H_adia is r-free.
#
#   J_adia[1] = Γ_sqrt_1 (g_1 ϵ_1 σ⁺_1 - g_2 σ⁻_2), with Γ_sqrt_1 ∝ κ_1√κ_2/D,
#              ϵ_1 = 2η/κ_1 and D = 4η² - κ_1κ_2 r-free. Then
#                  Γ_sqrt_1 g_1 ϵ_1 ∝ η g_1 √κ_2 ∝ (βr)^(1/4) r^(-1/4) = β^(1/4)
#                  Γ_sqrt_1 g_2     ∝ κ_1 √κ_2 g_2 ∝ r^(1/2-1/4-1/4) = r^0
#              -- both r-free. Same for J_adia[2] by the 1<->2 symmetry.
#
# and since g ∝ sqrt(ξ), every term of L_adia carries exactly one factor of ξ,
# so L_adia -> ξ L_adia and its KERNEL -- the steady state -- is unchanged.
#
# So the adiabatic steady state is exactly a function of (x, β) alone. This
# PROVES analytically, in these coordinates, what CLAUDE.md's finding recorded
# empirically from the three (β, x) runs. It holds only for γ = γ_φ = 0, which
# is what the config sets; turning either on breaks the homogeneity argument.
#
# Two consequences worth stating before the run:
#
#   1. Every _adia panel must be FLAT to machine precision across this whole
#      grid. That is a free end-to-end regression test: if concurrence_adia
#      varies by more than ~1e-13, something is miswired. The end-of-run
#      report checks it explicitly.
#
#   2. Therefore tracedist and concurrence_diff here are PURE adiabatic-
#      elimination error, with no adiabatic-side variation mixed in. That is
#      not true of the (β, x) sweep, where both models move.
# ==============================================================================

# Disable BLAS multithreading to prevent CPU oversubscription during our
# own @threads loop.
BLAS.set_num_threads(1)

# Cost is assumed to scale as dim^COST_EXPONENT; see the note on the cost
# model below for the second factor. The end-of-run report prints the measured
# scaling so this can be checked against real data.
const COST_EXPONENT = 2.0

# Points at or above BIG_DIM are memory-heavy enough to contend with each
# other; at most MAX_BIG of them run concurrently.
#
# NOTE: 2000 is inherited from the (β, x) sweep, where it was calibrated to sit
# between that grid's 900 and 2704+ dimension clusters. This grid's dimension
# distribution is NOT the same. Pass 1 prints the min/median/max before any
# solve starts -- read it on the first run and re-derive this if the clusters
# have moved.
const BIG_DIM = 2000
const MAX_BIG = 2

fmt_time(s) = s < 60  ? @sprintf("%.0fs", s) :
              s < 3600 ? @sprintf("%dm%02ds", floor(Int, s/60), round(Int, s%60)) :
                         @sprintf("%dh%02dm", floor(Int, s/3600), floor(Int, (s%3600)/60))

nanmat() = fill(ComplexF64(NaN), 4, 4)

# ------------------------------------------------------------------------
# rxi_to_physical: the coordinate inversion, derived in config_r_xi.jl.
#
# Returns a NamedTuple rather than a positional tuple ON PURPOSE.
# run_sim(g_1, g_2, κ_1, κ_2, η, ...) puts the couplings first and the decay
# rates second; estimate_truncation(κ_1, κ_2, η, g_1, g_2) does the opposite.
# Two call sites, two different orders, five same-typed Float64 arguments --
# a positional tuple here would make a silent transposition a matter of time.
# ------------------------------------------------------------------------
function rxi_to_physical(r, ξ; x, β, κ_geo)
    sr  = sqrt(r)
    amp = (κ_geo / 2) * sqrt(ξ)
    q   = (β * r)^0.25
    return (κ_1 = κ_geo * sr,
            κ_2 = κ_geo / sr,
            η   = (κ_geo / 2) * sqrt(x),
            g_1 = amp * q,
            g_2 = amp / q)
end

# ------------------------------------------------------------------------
# liouvillian_gap: the slowest linear relaxation rate, for the cost model.
#
# The (β, x) sweep's cost model uses 1/(1 - sqrt(x)), which is the gap for
# κ_1 = κ_2 with the κ/2 prefactor dropped as constant. Neither simplification
# survives here: x is FIXED, so that factor is a constant and does no ordering
# work at all, while r varies and the two decay rates are no longer equal.
#
# Setting g = 0, the drift matrix on (a_1, a_2†) is
#
#       M = -[ κ_1/2   η    ]
#            [ η       κ_2/2]
#
# with eigenvalues -( (κ_1+κ_2)/4 ± sqrt( ((κ_1-κ_2)/4)² + η² ) ), so
#
#       gap = (κ_1+κ_2)/4 - sqrt( ((κ_1-κ_2)/4)² + η² )
#           = (κ_geo/4) [ (√r + 1/√r) - sqrt( (√r - 1/√r)² + 4x ) ]
#
# Two checks: at r = 1 this collapses to (κ_geo/2)(1 - sqrt(x)), the old
# formula; and it vanishes at x = 1 for EVERY r, confirming that x alone sets
# threshold and that the r axis cannot walk into the divergence on its own.
#
# Krylov iteration count goes as 1/gap, so cost ∝ dim^COST_EXPONENT / gap.
# The gap is symmetric under r <-> 1/r, which is the whole reason κ_geo is the
# pinned scale.
# ------------------------------------------------------------------------
function liouvillian_gap(r, x, κ_geo)
    sr = sqrt(r)
    s, d = sr + 1/sr, sr - 1/sr
    return (κ_geo / 4) * (s - sqrt(d^2 + 4x))
end

function run_sweep_rxi(; x_fixed::Float64 = 0.485,
                         β_fixed::Float64 = 1.0,
                         κ_geo::Float64 = 2.0,
                         γ_val::Float64 = 0.0,
                         γ_phi_val::Float64 = 0.0,
                         len::Int = 7,
                         r_range = (0.1, 10.0),
                         ξ_range = (0.1, 1.0))

    0 < x_fixed < 1 ||
        error("x_fixed must satisfy 0 < x < 1: the adiabatic model has " *
              "4η² - κ_1κ_2 in every denominator and diverges at x = 1")
    ξ_range[1] > 0 ||
        error("ξ_range must start strictly above 0: at ξ = 0 both couplings " *
              "vanish, L_adia is identically zero and its kernel is the whole " *
              "state space (see config_r_xi.jl)")

    # r log-spaced, ξ linear -- mirroring the (β, x) sweep, where the ratio-like
    # coordinate is logarithmic and the bounded one is linear.
    r_vals = collect(exp10.(range(log10.(r_range)..., length = len)))
    ξ_vals = collect(range(ξ_range..., length = len))

    # ROW = r, COLUMN = ξ. This mirrors the (β, x) sweep's ROW = β, COLUMN = x,
    # which is what lets flatten_grid and plot_map be reused unchanged: the row
    # coordinate is the horizontal/log axis, the column coordinate the
    # vertical/linear one.
    N_r, N_ξ = length(r_vals), length(ξ_vals)
    total = N_r * N_ξ

    rho_full_grid = [nanmat() for _ in 1:N_r, _ in 1:N_ξ]
    rho_adia_grid = [nanmat() for _ in 1:N_r, _ in 1:N_ξ]

    N_1_grid    = zeros(Int, N_r, N_ξ)
    N_2_grid    = zeros(Int, N_r, N_ξ)
    dim_grid    = zeros(Int, N_r, N_ξ)
    time_grid   = fill(NaN, N_r, N_ξ)
    status_grid = fill(:ok, N_r, N_ξ)

    phys(i, j) = rxi_to_physical(r_vals[i], ξ_vals[j];
                                 x = x_fixed, β = β_fixed, κ_geo = κ_geo)

    println("Sweep on $(Threads.nthreads()) thread(s), $total grid points.")
    @printf("  fixed: x = %.4f, β = %.4f, κ_geo = %.4f  (η = %.4f throughout)\n",
            x_fixed, β_fixed, κ_geo, (κ_geo / 2) * sqrt(x_fixed))

    # ======================================================================
    # PASS 1 -- truncation only. Cheap relative to the steady-state solve,
    # and it gives the Hilbert dimension of every point before committing to
    # any of them, so total work is known rather than discovered.
    # ======================================================================
    print("Pass 1/2: estimating truncations... ")
    t_pass1 = @elapsed begin
        @threads for idx in CartesianIndices((N_r, N_ξ))
            i, j = idx.I
            try
                p = phys(i, j)
                N_1, N_2 = estimate_truncation(p.κ_1, p.κ_2, p.η, p.g_1, p.g_2)
                N_1_grid[i, j] = N_1
                N_2_grid[i, j] = N_2
                dim_grid[i, j] = (N_1 + 1) * (N_2 + 1) * 4
            catch e
                status_grid[i, j] = :truncation_failed
                @warn "truncation failed at (r=$(r_vals[i]), ξ=$(ξ_vals[j]))" exception = e
            end
        end
    end
    @printf("done in %s\n", fmt_time(t_pass1))

    # --- what we now know about the work ahead ---
    todo = [idx for idx in CartesianIndices((N_r, N_ξ)) if status_grid[idx] === :ok]
    isempty(todo) && error("every point failed truncation; nothing to solve")

    dims = [dim_grid[idx] for idx in todo]
    dmax = maximum(dims)
    mem_max = dmax^2 * 16 / 2^20          # one dense rho, MiB

    println("  Hilbert dimensions: min $(minimum(dims)), median $(round(Int, median(dims))), max $dmax")
    @printf("  largest dense rho: %.0f MiB (%d concurrent solves possible)\n",
            mem_max, Threads.nthreads())
    # Known to underestimate: steadystate.iterative holds a whole Krylov
    # subspace, not one dense rho (known issue 5). Treat as a lower bound.
    mem_max * Threads.nthreads() > 4096 &&
        @warn "peak memory could exceed 4 GiB with all threads on large points"
    dmax >= BIG_DIM || @info "no point reaches BIG_DIM = $BIG_DIM; the big/small " *
                             "queue split is inactive this run"

    # Sort expensive-first: with a dynamic queue this is longest-processing-
    # time-first, which keeps the tail from being one huge job starting last.
    gap_of(i)  = liouvillian_gap(r_vals[i], x_fixed, κ_geo)
    cost = Dict(idx => Float64(dim_grid[idx])^COST_EXPONENT /
                       max(gap_of(idx.I[1]), 1e-6)
                for idx in todo)
    order = sort(todo; by = idx -> -cost[idx])
    cost_total = sum(values(cost))

    # ======================================================================
    # PASS 2 -- the actual solves, on TWO dynamic work queues.
    #
    # A single queue lets all threads grab the largest points at once,
    # because they sort first, and steadystate.iterative's Krylov vectors are
    # dim^2 complex numbers each. Splitting into a big queue (capped at
    # MAX_BIG preferring workers) and a small one keeps every core busy
    # without several large workspaces competing.
    #
    # NOTE: MAX_BIG is a spawn-time PREFERENCE, not a concurrency limit --
    # once small_idx drains, every thread falls through to the big queue
    # (known issue 6, carried over verbatim from the (β, x) sweep).
    # ======================================================================
    println("Pass 2/2: solving...")

    big_idx   = [idx for idx in order if dim_grid[idx] >= BIG_DIM]
    small_idx = [idx for idx in order if dim_grid[idx] <  BIG_DIM]
    n_big_workers = min(MAX_BIG, Threads.nthreads(), max(length(big_idx), 1))

    @printf("  %d large point(s) (dim >= %d) on %d worker(s); %d small on the rest\n",
            length(big_idx), BIG_DIM, n_big_workers, length(small_idx))

    next_big    = Atomic{Int}(0)
    next_small  = Atomic{Int}(0)
    started     = Atomic{Int}(0)   # dispatch counter, for the start lines
    inflight    = Set{CartesianIndex{2}}()   # guarded by print_lock
    print_lock  = ReentrantLock()
    done_count  = 0            # guarded by print_lock
    done_cost   = 0.0          # guarded by print_lock
    busy_time   = 0.0          # summed per-point solve time, guarded
    t_start     = time()

    # Take from the preferred queue; fall through to the other when empty.
    function take_job(prefer_big::Bool)
        if prefer_big
            k = atomic_add!(next_big, 1) + 1
            k <= length(big_idx) && return big_idx[k]
            k = atomic_add!(next_small, 1) + 1
            k <= length(small_idx) && return small_idx[k]
        else
            k = atomic_add!(next_small, 1) + 1
            k <= length(small_idx) && return small_idx[k]
            k = atomic_add!(next_big, 1) + 1
            k <= length(big_idx) && return big_idx[k]
        end
        return nothing
    end

    @sync for w in 1:Threads.nthreads()
        prefer_big = w <= n_big_workers
        Threads.@spawn while true
            idx = take_job(prefer_big)
            idx === nothing && break
            i, j = idx.I

            p = phys(i, j)
            N_1, N_2 = N_1_grid[i, j], N_2_grid[i, j]

            # Announce on dispatch, not just on completion: the expensive
            # points run first and can take minutes.
            n_started = atomic_add!(started, 1) + 1
            lock(print_lock) do
                push!(inflight, idx)
                @printf("  start  %2d/%-2d  r=%7.3f  ξ=%.4f  N₁=%2d N₂=%2d  dim=%5d  (%d running)\n",
                        n_started, length(order), r_vals[i], ξ_vals[j], N_1, N_2,
                        dim_grid[i, j], length(inflight))
            end

            dt = @elapsed try
                rho_full, rho_adia = run_sim(p.g_1, p.g_2, p.κ_1, p.κ_2, p.η,
                                             γ_val, γ_phi_val, N_1, N_2)
                @assert size(rho_full.data) == (4, 4) "rho_full is not a 2-qubit state"
                @assert size(rho_adia.data) == (4, 4) "rho_adia is not a 2-qubit state"
                rho_full_grid[i, j] = Matrix{ComplexF64}(rho_full.data)
                rho_adia_grid[i, j] = Matrix{ComplexF64}(rho_adia.data)
            catch e
                status_grid[i, j] = :solve_failed
                lock(print_lock) do
                    @warn "solve failed at (r=$(r_vals[i]), ξ=$(ξ_vals[j]), dim=$(dim_grid[i,j]))" exception = e
                end
            end
            time_grid[i, j] = dt

            lock(print_lock) do
                delete!(inflight, idx)
                done_count += 1
                done_cost  += cost[idx]
                busy_time  += dt

                # ETA: fit one scale factor from work already done, apply it
                # to the remaining cost, divide by thread count. Unreliable
                # until a few points of each size have finished.
                rem_cost = cost_total - done_cost
                eta = rem_cost * (busy_time / done_cost) / Threads.nthreads()
                pct = 100 * done_cost / cost_total

                @printf("  done   %2d/%-2d  r=%7.3f  ξ=%.4f  dim=%5d  %8.2fs  |  %5.1f%%  ETA %s  (%d running)\n",
                        done_count, length(order), r_vals[i], ξ_vals[j],
                        dim_grid[i, j], dt, pct, fmt_time(eta), length(inflight))
            end
        end
    end

    wall = time() - t_start

    # ======================================================================
    # REPORT
    #
    # Wrapped in try/catch on purpose. Everything below is cosmetic, and it
    # sits before the return -- an exception here would discard a completed
    # sweep because a summary line had a bug. The report may fail, the
    # results may not.
    # ======================================================================
    try
    n_ok = count(==(:ok), status_grid)
    println("\nSweep complete. $n_ok/$total points solved.")
    @printf("  pass 1 (truncation): %s\n", fmt_time(t_pass1))
    @printf("  pass 2 (solves):     %s wall, %s summed across threads\n",
            fmt_time(wall), fmt_time(busy_time))
    @printf("  parallel efficiency: %.0f%% of %d threads\n",
            100 * busy_time / (wall * Threads.nthreads()), Threads.nthreads())

    ts = filter(isfinite, vec(time_grid))
    if !isempty(ts)
        @printf("  per-point: min %.2fs, median %.2fs, max %.2fs\n",
                minimum(ts), median(ts), maximum(ts))

        pts = [(Float64(dim_grid[idx]), time_grid[idx])
               for idx in CartesianIndices(time_grid)
               if isfinite(time_grid[idx]) && time_grid[idx] > 0.05 && dim_grid[idx] > 0]
        if length(pts) ≥ 4
            lx = log.(first.(pts)); ly = log.(last.(pts))
            slope = sum((lx .- mean(lx)) .* (ly .- mean(ly))) / sum((lx .- mean(lx)).^2)
            @printf("  measured scaling: time ~ dim^%.2f (COST_EXPONENT is %.2f)\n",
                    slope, COST_EXPONENT)
        end

        # vec() is essential: collect(CartesianIndices(...)) of an (N_r × N_ξ)
        # grid is itself a MATRIX, and sort on a multidimensional array demands
        # a `dims` keyword.
        worst = sort(vec(collect(CartesianIndices(time_grid)));
                     by = idx -> isfinite(time_grid[idx]) ? -time_grid[idx] : 0.0)[1:min(3, total)]
        println("  slowest points:")
        for idx in worst
            i, j = idx.I
            isfinite(time_grid[idx]) || continue
            @printf("    r=%7.3f  ξ=%.4f  dim=%5d  %6.2fs\n",
                    r_vals[i], ξ_vals[j], dim_grid[i, j], time_grid[idx])
        end
    end

    n_bad = total - n_ok
    n_bad > 0 && @warn "$n_bad point(s) failed -- see status_grid"
    catch e
        @warn "report failed (results are unaffected and returned below)" exception = e
    end

    return (full = rho_full_grid,
            adia = rho_adia_grid,
            r_vals = r_vals,
            ξ_vals = ξ_vals,
            N_1_grid = N_1_grid,
            N_2_grid = N_2_grid,
            dim_grid = dim_grid,
            time_grid = time_grid,
            status_grid = status_grid,
            params = (; x_fixed, β_fixed, κ_geo, γ_val, γ_phi_val, len, r_range, ξ_range))
end

# ------------------------------------------------------------------------
# check_adia_flat: the self-test described in the header.
#
# L_adia depends on (x, β) alone, both of which are FIXED here, so every
# adiabatic observable must be constant across the entire grid. Anything
# above a loose numerical floor means the sweep is miswired -- most likely
# the coordinate inversion, since that is the only new code on the path.
#
# Reported, not thrown: a violation is worth seeing next to the data that
# shows it, and the states are already saved by the time this runs.
# ------------------------------------------------------------------------
function check_adia_flat(outs; tol = 1e-10)
    adia_keys = sort([k for k in keys(outs) if endswith(String(k), "_adia")])
    isempty(adia_keys) && return nothing

    println("\nAdiabatic flatness check (L_adia depends on (x, β) alone, both fixed):")
    for k in adia_keys
        v = filter(isfinite, vec(outs[k]))
        isempty(v) && (@printf("  %-22s  all NaN\n", String(k)); continue)
        spread = maximum(v) - minimum(v)
        ok = spread < tol
        @printf("  %-22s  spread = %.3e   %s\n", String(k), spread,
                ok ? "flat (as predicted)" : "NOT FLAT")
        ok || @warn "$k varies by $spread across the (r, ξ) grid, but the " *
                    "adiabatic Liouvillian is analytically independent of both. " *
                    "Suspect rxi_to_physical, or a nonzero γ / γ_φ (which breaks " *
                    "the homogeneity argument)."
    end
    return nothing
end


# ==============================================================================
# MAIN
# ==============================================================================
begin
    obs_list = select_observables(ACTIVE_OBSERVABLES_RXI)
    labels   = obs_labels(obs_list)

    res  = run_sweep_rxi(; PARAMS_RXI...)
    outs = analyze(res, obs_list)

    summarize(outs)
    check_adia_flat(outs)

    # `params` is the parameter set that actually produced these states --
    # read plot axes and titles from here, not from config_r_xi.jl, which may
    # have been edited since.
    full_data   = res.full
    adia_data   = res.adia
    r_vals      = res.r_vals
    ξ_vals      = res.ξ_vals
    params      = res.params
    N_1_grid    = res.N_1_grid
    N_2_grid    = res.N_2_grid
    dim_grid    = res.dim_grid
    time_grid   = res.time_grid
    status_grid = res.status_grid

    @save OUTFILE_RXI outs labels full_data adia_data r_vals ξ_vals params N_1_grid N_2_grid dim_grid time_grid status_grid
    println("\nSaved to $OUTFILE_RXI")
end
