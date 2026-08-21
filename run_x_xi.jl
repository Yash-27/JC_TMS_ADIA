using Base.Threads
using Printf
using LinearAlgebra
using JLD2
using Statistics

include("config_x_xi.jl")
include("truncation.jl")
include("jc_pump_disp_asy.jl")
include("observables.jl")

# ==============================================================================
# run_x_xi.jl
#
# The THIRD sweep driver. Same shape as the other two -- a top-level begin
# block at the bottom, so `include`ing this file RUNS the sweep -- over the
# (x, ξ) plane at fixed (β, r) = (1, 1).
#
# Run it with:   julia --project=. -t auto run_x_xi.jl
#
# What is shared and what is not:
#
#   shared   jc_pump_disp_asy.jl (run_sim / run_adia), truncation.jl
#            (estimate_truncation), observables.jl (the whole registry).
#            All three take physical parameters point by point, so none of them
#            knows or cares which coordinates are being swept.
#
#   its own  the sweep loop below. Like run_r_xi.jl it duplicates the two-pass
#            scheduler rather than sharing it -- the coordinate inversion and
#            the cost model differ, and including either other driver to borrow
#            code would run its sweep. That is known issue 8, and this file
#            makes it THREE copies. It also means the semaphore fix below
#            exists in this copy only.
#
# ==============================================================================
# WHAT THIS SWEEP IS EXPECTED TO SHOW -- and the two self-tests that come free
#
# 1. ρ_adia IS A FUNCTION OF x ALONE HERE, so every _adia panel must be
#    constant DOWN EACH COLUMN (across ξ) while varying across columns.
#
#    Substituting the inversion into the adiabatic model, g ∝ sqrt(ξ) makes
#    every term of L_adia carry exactly one factor of ξ, so L_adia -> ξ L_adia
#    and its KERNEL -- the steady state -- is unchanged. β and r are fixed, so
#    x is the only coordinate ρ_adia can see. check_adia_x_only asserts this at
#    the end of the run and additionally re-derives each column value from a
#    direct run_adia call, which costs microseconds.
#
#    This is the same statement run_r_xi.jl's check_adia_flat makes, in the
#    complementary direction: there the whole grid is one constant, here each
#    column is.
#
# 2. THE ξ = 1 ROW IS ALREADY STORED. With β = r = 1 and κ_geo = 2, ξ = 1 gives
#    κ_1 = κ_2 = 2, g_1 = g_2 = 1, η = sqrt(x) -- run A's β = 1 row, on run A's
#    own x grid. ξ_range = (0.2, 5.0) is chosen so ξ = 1 lands exactly on grid
#    index 4 (see config_x_xi.jl). So row 4 must reproduce
#    results_k1_2.0_k2_2.0_g2_1.0.jld2's β = 1 row: the adiabatic values to the
#    last bit, the full ones to ~1.5e-8 (steadystate.iterative's unseeded-RNG
#    floor). That is a regression test of the inversion, the truncation and this
#    driver at once -- including at the expensive dim-4096 point.
#
# Neither test is imposed anywhere in the code; both are consequences.
#
# Expected structure, from independently known results:
#
#   * concurrence_adia is exactly 0 for x ≥ 0.485 (three of seven columns) --
#     the adiabatic entanglement threshold x* = 0.48388826 at β = 1. There
#     concurrence_diff becomes an exact duplicate of concurrence_full. That is
#     the analytic boundary, not a bug; tracedist stays fully informative.
#   * concurrence_full should peak near x ≈ 0.29, ξ ≈ 2.9 at C ≈ 0.297 -- the
#     full model's global optimum, which no existing run covers.
# ==============================================================================

# Disable BLAS multithreading to prevent CPU oversubscription during our
# own @threads loop.
BLAS.set_num_threads(1)

# Cost is assumed to scale as dim^COST_EXPONENT; see the note on the cost
# model below for the second factor. The end-of-run report prints the measured
# scaling so this can be checked against real data.
const COST_EXPONENT = 2.0

# ------------------------------------------------------------------------
# BIG_DIM / MAX_BIG -- and why this driver enforces the cap for real.
#
# Measured over this exact grid (Pass 1 prints the histogram; read it and
# check): 42 of the 49 points sit at N = 14, dim = 900, and the seven points of
# the x = 0.7225 column sit at N = 31, dim = 4096. The truncation there is set
# by η, not by g, so it is dim 4096 for EVERY ξ from 0.1 to 5. BIG_DIM = 2000
# therefore separates the two clusters cleanly, exactly as it does in sweep 1.
#
# What does NOT carry over is sweep 1 and 2's treatment of MAX_BIG as a
# spawn-time preference. There, once the small queue drains every thread falls
# through to the big queue (known issue 6). Here that is not a tail
# inefficiency, it is an out-of-memory risk: one dense ρ at dim 4096 is 268 MiB
# and steadystate.iterative holds a whole Krylov subspace of them (20-30
# vectors, known issue 5), so a single big solve can sit in the several-GiB
# range. Eight threads doing that at once does not fit in 16 GiB.
#
# So the big solves go through a real counting semaphore. The queue preference
# is kept -- it starts the expensive points early, which is what keeps the tail
# short -- and the semaphore is what actually bounds concurrency.
# ------------------------------------------------------------------------
const BIG_DIM = 2000
const MAX_BIG = 1

# ------------------------------------------------------------------------
# KRYLOV_VECTORS: how many dim²-sized complex vectors one steadystate.iterative
# solve is holding at its peak, for the memory estimate in Pass 1.
#
# Issue 5 records that the old estimate -- ONE dense ρ -- is wrong by "at least
# 6-8x, probably much more", because restarted GMRES keeps its whole subspace.
# This constant makes the estimate honest instead of a documented lower bound,
# and it is CALIBRATED, not guessed: at dim 4096 one dense ρ is 256 MiB, and a
# single solve at (x = 0.7225, ξ = 1.357) measured 5.04 GiB peak RSS
# (/usr/bin/time -l, -t 1). 5.04 GiB / 256 MiB = 20.2.
#
# So 20 reproduces the one point where it has been measured. It is one
# measurement at one dimension -- treat it as an order-of-magnitude figure, and
# re-measure if the solver or its options change.
# ------------------------------------------------------------------------
const KRYLOV_VECTORS = 20

fmt_time(s) = s < 60  ? @sprintf("%.0fs", s) :
              s < 3600 ? @sprintf("%dm%02ds", floor(Int, s/60), round(Int, s%60)) :
                         @sprintf("%dh%02dm", floor(Int, s/3600), floor(Int, (s%3600)/60))

nanmat() = fill(ComplexF64(NaN), 4, 4)

# ------------------------------------------------------------------------
# xxi_to_physical: the coordinate inversion, derived in config_x_xi.jl.
#
# Identical algebra to run_r_xi.jl's rxi_to_physical with the roles of the free
# and the fixed coordinates swapped: there (r, ξ) vary and x is a keyword, here
# (ξ, x) vary and r is. Kept as a separate function rather than reusing that
# one so this driver is self-contained -- including run_r_xi.jl would run its
# sweep.
#
# Returns a NamedTuple rather than a positional tuple ON PURPOSE.
# run_sim(g_1, g_2, κ_1, κ_2, η, ...) puts the couplings first and the decay
# rates second; estimate_truncation(κ_1, κ_2, η, g_1, g_2) does the opposite.
# Two call sites, two different orders, five same-typed Float64 arguments --
# a positional tuple here would make a silent transposition a matter of time.
#
# At r = β = 1 this reduces to κ_1 = κ_2 = κ_geo, g_1 = g_2 = (κ_geo/2)sqrt(ξ),
# η = (κ_geo/2)sqrt(x): one rate per axis and nothing else moving. The general
# r and β are carried anyway so the file works off that slice unchanged.
# ------------------------------------------------------------------------
function xxi_to_physical(ξ, x; β, r, κ_geo)
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
# and at r = 1, which is where this sweep runs, that is (κ_geo/2)(1 - sqrt(x)).
# Krylov iteration count goes as 1/gap, so cost ∝ dim^COST_EXPONENT / gap.
#
# Unlike sweep 2 -- where x is fixed and this factor is a constant that does no
# ordering work at all -- x is the swept coordinate here, so the gap factor
# spans 1/(1 - sqrt(0.01)) = 1.11 to 1/(1 - sqrt(0.7225)) = 6.67 across the
# grid and genuinely orders the queue. It multiplies the 20.7x from dim² on the
# last column, which is why that column dominates the run.
#
# Copied from run_r_xi.jl rather than shared, for the reason in the header.
# ------------------------------------------------------------------------
function liouvillian_gap(r, x, κ_geo)
    sr = sqrt(r)
    s, d = sr + 1/sr, sr - 1/sr
    return (κ_geo / 4) * (s - sqrt(d^2 + 4x))
end

# ------------------------------------------------------------------------
# write_checkpoint: dump everything known so far to PARTIAL_XXI.
#
# Called under the print lock after every completed point. CLAUDE.md records 45
# minutes of solves lost when a long run wrapped in a shell `timeout` was
# killed before its single end-of-run write; a dim-4096 solve here runs for
# minutes, so an all-or-nothing write is the same trap. Writing 49 4x4 matrices
# plus a few integer grids costs milliseconds.
#
# Never throws: a failed checkpoint must not take a running sweep with it.
# The states are stored raw (no observables) -- analyze() is instant, so
# recovering from a partial file means loading it and calling analyze.
# ------------------------------------------------------------------------
function write_checkpoint(path, full, adia, ξ_vals, x_vals, params,
                          dim_grid, time_grid, status_grid)
    try
        jldsave(path;
                full_data = full, adia_data = adia,
                ξ_vals, x_vals, params,
                dim_grid, time_grid, status_grid)
    catch e
        @warn "checkpoint write to $path failed (the sweep is unaffected)" exception = e
    end
    return nothing
end

function run_sweep_xxi(; β_fixed::Float64 = 1.0,
                         r_fixed::Float64 = 1.0,
                         κ_geo::Float64 = 2.0,
                         γ_val::Float64 = 0.0,
                         γ_phi_val::Float64 = 0.0,
                         len::Int = 7,
                         ξ_range = (0.2, 5.0),
                         x_range = (0.01, 0.7225))

    0 < x_range[1] && x_range[2] < 1 ||
        error("x_range must lie strictly inside (0, 1): the adiabatic model " *
              "has 4η² - κ_1κ_2 in every denominator and diverges at x = 1")
    ξ_range[1] > 0 ||
        error("ξ_range must start strictly above 0: at ξ = 0 both couplings " *
              "vanish, L_adia is identically zero and its kernel is the whole " *
              "state space (see config_x_xi.jl)")

    # ξ log-spaced, x linear -- ξ is a ratio spanning more than a decade, x is
    # bounded above by 1 and shared with sweep 1's linear grid.
    ξ_vals = collect(exp10.(range(log10.(ξ_range)..., length = len)))
    x_vals = collect(range(x_range..., length = len))

    # ROW = ξ, COLUMN = x. This mirrors the other two sweeps' ROW = β / ROW = r,
    # i.e. the row coordinate is the horizontal/log axis and the column
    # coordinate the vertical/linear one, which is what lets flatten_grid and
    # plot_map be reused unchanged.
    N_ξ, N_x = length(ξ_vals), length(x_vals)
    total = N_ξ * N_x

    rho_full_grid = [nanmat() for _ in 1:N_ξ, _ in 1:N_x]
    rho_adia_grid = [nanmat() for _ in 1:N_ξ, _ in 1:N_x]

    N_1_grid    = zeros(Int, N_ξ, N_x)
    N_2_grid    = zeros(Int, N_ξ, N_x)
    dim_grid    = zeros(Int, N_ξ, N_x)
    time_grid   = fill(NaN, N_ξ, N_x)
    status_grid = fill(:ok, N_ξ, N_x)

    params = (; β_fixed, r_fixed, κ_geo, γ_val, γ_phi_val, len, ξ_range, x_range)

    phys(i, j) = xxi_to_physical(ξ_vals[i], x_vals[j];
                                β = β_fixed, r = r_fixed, κ_geo = κ_geo)

    println("Sweep on $(Threads.nthreads()) thread(s), $total grid points.")
    @printf("  fixed: β = %.4f, r = %.4f, κ_geo = %.4f  (κ_1 = %.4f, κ_2 = %.4f throughout)\n",
            β_fixed, r_fixed, κ_geo, κ_geo * sqrt(r_fixed), κ_geo / sqrt(r_fixed))
    println("  ξ (rows, log10): ", join((@sprintf("%.4f", v) for v in ξ_vals), "  "))
    println("  x (cols, linear): ", join((@sprintf("%.4f", v) for v in x_vals), "  "))
    flush(stdout)

    # ======================================================================
    # PASS 1 -- truncation only. Cheap relative to the steady-state solve,
    # and it gives the Hilbert dimension of every point before committing to
    # any of them, so total work is known rather than discovered.
    # ======================================================================
    print("Pass 1/2: estimating truncations... ")
    flush(stdout)
    t_pass1 = @elapsed begin
        @threads for idx in CartesianIndices((N_ξ, N_x))
            i, j = idx.I
            try
                p = phys(i, j)
                N_1, N_2 = estimate_truncation(p.κ_1, p.κ_2, p.η, p.g_1, p.g_2)
                N_1_grid[i, j] = N_1
                N_2_grid[i, j] = N_2
                dim_grid[i, j] = (N_1 + 1) * (N_2 + 1) * 4
            catch e
                status_grid[i, j] = :truncation_failed
                @warn "truncation failed at (ξ=$(ξ_vals[i]), x=$(x_vals[j]))" exception = e
            end
        end
    end
    @printf("done in %s\n", fmt_time(t_pass1))

    # --- what we now know about the work ahead ---
    todo = [idx for idx in CartesianIndices((N_ξ, N_x)) if status_grid[idx] === :ok]
    isempty(todo) && error("every point failed truncation; nothing to solve")

    dims = [dim_grid[idx] for idx in todo]
    dmax = maximum(dims)
    mem_rho  = dmax^2 * 16 / 2^20            # one dense rho, MiB
    mem_peak = mem_rho * KRYLOV_VECTORS      # what a solve actually costs

    println("  Hilbert dimensions: min $(minimum(dims)), median $(round(Int, median(dims))), max $dmax")
    # The per-dimension histogram, not just the summary: BIG_DIM has to sit
    # between the clusters, and that is only checkable by seeing them.
    for d in sort(unique(dims))
        @printf("    dim %5d : %2d point(s)%s\n", d, count(==(d), dims),
                d >= BIG_DIM ? "   <- big queue" : "")
    end
    # One literal format string: @printf takes the format at macro-expansion
    # time and rejects a `*`-concatenated one with "First argument to @printf
    # after io must be a format string".
    @printf("  largest solve: one dense rho %.0f MiB, estimated peak %.1f GiB; %d concurrent max => %.1f GiB\n",
            mem_rho, mem_peak / 1024, MAX_BIG, mem_peak * MAX_BIG / 1024)
    mem_peak * MAX_BIG > 8192 &&
        @warn "estimated peak for the big queue alone is " *
              @sprintf("%.1f GiB", mem_peak * MAX_BIG / 1024) *
              "; lower MAX_BIG or cap x_range below the truncation cliff"
    dmax >= BIG_DIM || @info "no point reaches BIG_DIM = $BIG_DIM; the big/small " *
                             "queue split and the semaphore are both inactive this run"
    flush(stdout)

    # Sort expensive-first: with a dynamic queue this is longest-processing-
    # time-first, which keeps the tail from being one huge job starting last.
    gap_of(j)  = liouvillian_gap(r_fixed, x_vals[j], κ_geo)
    cost = Dict(idx => Float64(dim_grid[idx])^COST_EXPONENT /
                       max(gap_of(idx.I[2]), 1e-6)
                for idx in todo)
    order = sort(todo; by = idx -> -cost[idx])
    cost_total = sum(values(cost))

    # ======================================================================
    # PASS 2 -- the actual solves, on TWO dynamic work queues plus a
    # semaphore.
    #
    # The queues are sweep 1 and 2's: big (dim >= BIG_DIM, preferred by
    # MAX_BIG workers at spawn time) and small, both sorted expensive-first.
    # @threads splits into contiguous chunks and the expensive points are
    # contiguous (the whole last x column), so a plain @threads would hand one
    # thread all of them.
    #
    # The semaphore is what this driver adds. Unlike the other two, the cap on
    # concurrent big solves is ENFORCED, not merely preferred -- see the note
    # on BIG_DIM above for why that matters at dim 4096 on 16 GiB.
    # ======================================================================
    println("Pass 2/2: solving...")

    big_idx   = [idx for idx in order if dim_grid[idx] >= BIG_DIM]
    small_idx = [idx for idx in order if dim_grid[idx] <  BIG_DIM]
    n_big_workers = min(MAX_BIG, Threads.nthreads(), max(length(big_idx), 1))

    @printf("  %d large point(s) (dim >= %d), at most %d at a time; %d small on the rest\n",
            length(big_idx), BIG_DIM, MAX_BIG, length(small_idx))
    flush(stdout)

    big_sem     = Base.Semaphore(MAX_BIG)
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
    # The fall-through is safe in BOTH directions here: a thread that ends up
    # with a big job blocks on the semaphore rather than allocating alongside
    # the others, so the tail serializes instead of thrashing.
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
            is_big = dim_grid[i, j] >= BIG_DIM

            # Announce on dispatch, not just on completion: the expensive
            # points run first and can take many minutes.
            n_started = atomic_add!(started, 1) + 1
            lock(print_lock) do
                @printf("  queued %2d/%-2d  ξ=%.4f  x=%.4f  N₁=%2d N₂=%2d  dim=%5d%s\n",
                        n_started, length(order), ξ_vals[i], x_vals[j], N_1, N_2,
                        dim_grid[i, j], is_big ? "   (waiting for a big slot)" : "")
                flush(stdout)
            end

            # Only the big solves are gated. try/finally so a thrown solve
            # cannot leak a permit and deadlock the tail.
            is_big && Base.acquire(big_sem)
            try
                lock(print_lock) do
                    push!(inflight, idx)
                    @printf("  start  %2d/%-2d  ξ=%.4f  x=%.4f  dim=%5d  (%d running)\n",
                            n_started, length(order), ξ_vals[i], x_vals[j],
                            dim_grid[i, j], length(inflight))
                    flush(stdout)
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
                        @warn "solve failed at (ξ=$(ξ_vals[i]), x=$(x_vals[j]), dim=$(dim_grid[i,j]))" exception = e
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
                    # until a few points of each size have finished, and it does
                    # not know about the semaphore -- so it under-estimates the
                    # big-column tail, which cannot use every thread.
                    rem_cost = cost_total - done_cost
                    eta = rem_cost * (busy_time / done_cost) / Threads.nthreads()
                    pct = 100 * done_cost / cost_total

                    @printf("  done   %2d/%-2d  ξ=%.4f  x=%.4f  dim=%5d  %8.2fs  |  %5.1f%%  ETA %s  (%d running)\n",
                            done_count, length(order), ξ_vals[i], x_vals[j],
                            dim_grid[i, j], dt, pct, fmt_time(eta), length(inflight))
                    flush(stdout)

                    write_checkpoint(PARTIAL_XXI, rho_full_grid, rho_adia_grid,
                                     ξ_vals, x_vals, params,
                                     dim_grid, time_grid, status_grid)
                end
            finally
                is_big && Base.release(big_sem)
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

        # Per-dimension cluster timings. With two clusters 4.5x apart in dim
        # this is the number to quote when budgeting the next run.
        for d in sort(unique(dim_grid[idx] for idx in CartesianIndices(time_grid)
                            if isfinite(time_grid[idx])))
            tv = [time_grid[idx] for idx in CartesianIndices(time_grid)
                  if isfinite(time_grid[idx]) && dim_grid[idx] == d]
            @printf("    dim %5d : %2d point(s), median %8.2fs, total %s\n",
                    d, length(tv), median(tv), fmt_time(sum(tv)))
        end

        pts = [(Float64(dim_grid[idx]), time_grid[idx])
               for idx in CartesianIndices(time_grid)
               if isfinite(time_grid[idx]) && time_grid[idx] > 0.05 && dim_grid[idx] > 0]
        if length(pts) ≥ 4
            lx = log.(first.(pts)); ly = log.(last.(pts))
            slope = sum((lx .- mean(lx)) .* (ly .- mean(ly))) / sum((lx .- mean(lx)).^2)
            @printf("  measured scaling: time ~ dim^%.2f (COST_EXPONENT is %.2f)\n",
                    slope, COST_EXPONENT)
        end

        # vec() is essential: collect(CartesianIndices(...)) of an (N_ξ × N_x)
        # grid is itself a MATRIX, and sort on a multidimensional array demands
        # a `dims` keyword.
        worst = sort(vec(collect(CartesianIndices(time_grid)));
                     by = idx -> isfinite(time_grid[idx]) ? -time_grid[idx] : 0.0)[1:min(3, total)]
        println("  slowest points:")
        for idx in worst
            i, j = idx.I
            isfinite(time_grid[idx]) || continue
            @printf("    ξ=%.4f  x=%.4f  dim=%5d  %6.2fs\n",
                    ξ_vals[i], x_vals[j], dim_grid[i, j], time_grid[idx])
        end
    end

    n_bad = total - n_ok
    n_bad > 0 && @warn "$n_bad point(s) failed -- see status_grid"
    catch e
        @warn "report failed (results are unaffected and returned below)" exception = e
    end

    return (full = rho_full_grid,
            adia = rho_adia_grid,
            ξ_vals = ξ_vals,
            x_vals = x_vals,
            N_1_grid = N_1_grid,
            N_2_grid = N_2_grid,
            dim_grid = dim_grid,
            time_grid = time_grid,
            status_grid = status_grid,
            params = params)
end

# ------------------------------------------------------------------------
# check_adia_x_only: the first self-test described in the header.
#
# L_adia depends on (x, β) alone and β is fixed, so every adiabatic observable
# must be constant DOWN EACH COLUMN of this grid (columns are x, rows are ξ)
# and must vary ACROSS columns -- a panel that is flat in both directions would
# mean x is not reaching the model at all, which a column-only check cannot
# see.
#
# Then, for each _adia key whose observable is in OBSERVABLE_REGISTRY, the
# column value is re-derived from a direct run_adia call at that x. That is a
# stronger statement than internal consistency: it checks the sweep's
# coordinate inversion against the model itself, and it costs microseconds
# because run_adia is a 4x4 dense eigensolve with no cavities.
#
# Reported, not thrown: a violation is worth seeing next to the data that shows
# it, and the states are already saved by the time this runs.
# ------------------------------------------------------------------------
function check_adia_x_only(outs, x_vals, params; tol = 1e-10)
    adia_keys = sort([k for k in keys(outs) if endswith(String(k), "_adia")])
    isempty(adia_keys) && return nothing

    by_name = Dict(o.name => o for o in OBSERVABLE_REGISTRY)

    println("\nAdiabatic self-test (ρ_adia depends on x alone; β, r fixed):")
    for k in adia_keys
        M = outs[k]

        # --- constant down each column (across ξ) ---
        # worst_j starts at 1, not 0: if every column is bit-identical the
        # spread never exceeds 0.0, worst_j is never assigned, and indexing
        # x_vals with 0 (or printing NaN for it) would report a nonexistent
        # column for what is the ideal outcome.
        worst, worst_j = 0.0, 1
        col_vals = fill(NaN, size(M, 2))
        for j in axes(M, 2)
            v = filter(isfinite, M[:, j])
            isempty(v) && continue
            col_vals[j] = first(v)
            s = maximum(v) - minimum(v)
            s > worst && ((worst, worst_j) = (s, j))
        end
        finite_cols = filter(isfinite, col_vals)
        if isempty(finite_cols)
            @printf("  %-22s  all NaN\n", String(k))
            continue
        end

        ok = worst < tol
        @printf("  %-22s  ξ-spread = %.3e (worst column x = %.4f)   %s\n",
                String(k), worst, x_vals[worst_j],
                ok ? "constant in ξ (as predicted)" : "NOT CONSTANT")
        ok || @warn "$k varies by $worst along ξ at fixed x, but L_adia -> ξ L_adia " *
                    "leaves the steady state unchanged. Suspect xxi_to_physical, " *
                    "or a nonzero γ / γ_φ (which breaks the homogeneity argument)."

        # --- and NOT flat across columns, or x is not reaching the model ---
        x_spread = maximum(finite_cols) - minimum(finite_cols)
        @printf("  %-22s  x-spread = %.3e   %s\n", "", x_spread,
                x_spread > tol ? "varies with x (as it must)" :
                                 "SUSPICIOUS -- flat in x too")

        # --- re-derive each column from run_adia directly ---
        name = Symbol(String(k)[1:end-5])          # strip "_adia"
        o = get(by_name, name, nothing)
        (o === nothing || o.kind !== :single) && continue

        worst_direct = 0.0
        for j in axes(M, 2)
            isfinite(col_vals[j]) || continue
            # ξ = 1 is a deliberate arbitrary probe, NOT a grid point (the grid
            # only contains ξ = 1 when the range endpoints are reciprocal). That
            # is what makes this an independent re-derivation rather than a
            # restatement: the ξ-spread check above already established the
            # column is constant in ξ, so if the model is wired correctly ANY ξ
            # must reproduce it, including one the sweep never solved at.
            p = xxi_to_physical(1.0, x_vals[j];
                                β = params.β_fixed, r = params.r_fixed,
                                κ_geo = params.κ_geo)
            ρ = run_adia(p.g_1, p.g_2, p.κ_1, p.κ_2, p.η,
                         params.γ_val, params.γ_phi_val)
            worst_direct = max(worst_direct, abs(real(o.f(ρ.data)) - col_vals[j]))
        end
        @printf("  %-22s  vs direct run_adia: worst |Δ| = %.3e   %s\n", "",
                worst_direct, worst_direct < tol ? "agrees" : "DISAGREES")
        worst_direct < tol ||
            @warn "$k disagrees with a direct run_adia call by $worst_direct -- " *
                  "the sweep and the model are not seeing the same parameters"
    end
    return nothing
end

# ------------------------------------------------------------------------
# check_against_run_A: the second self-test described in the header.
#
# With β = r = 1 and κ_geo = 2, ξ = 1 means κ_1 = κ_2 = 2, g_1 = g_2 = 1,
# η = sqrt(x) -- exactly run A's β = 1 row. If ξ = 1 is on this grid and run
# A's .jld2 is on disk, every observable along that row must reproduce it:
# the adiabatic values to the last bit (dense eigensolver, no RNG) and the full
# ones to ~1.5e-8 (steadystate.iterative seeds a shadow vector from the
# unseeded global RNG and stops at reltol = sqrt(eps)).
#
# Skipped silently when ξ = 1 is not on the grid or the file is gone -- it is a
# free bonus, not a prerequisite.
# ------------------------------------------------------------------------
function check_against_run_A(outs, ξ_vals, x_vals, params;
                             fname = "results_k1_2.0_k2_2.0_g2_1.0.jld2",
                             tol_full = 1e-6, tol_adia = 1e-12)
    (params.β_fixed ≈ 1.0 && params.r_fixed ≈ 1.0 && params.κ_geo ≈ 2.0) || return nothing
    i1 = findfirst(v -> isapprox(v, 1.0; rtol = 1e-9), ξ_vals)
    i1 === nothing && return nothing
    isfile(fname) || return nothing

    ref = jldopen(fname, "r") do f
        (outs = f["outs"], β_vals = f["β_vals"], x_vals = f["x_vals"])
    end
    jb = findfirst(v -> isapprox(v, 1.0; rtol = 1e-9), ref.β_vals)
    jb === nothing && return nothing

    # Only compare x values both grids actually have.
    cols = [(j, findfirst(v -> isapprox(v, x_vals[j]; rtol = 1e-9), ref.x_vals))
            for j in eachindex(x_vals)]
    cols = [(j, jr) for (j, jr) in cols if jr !== nothing]
    isempty(cols) && return nothing

    # NOTE the two grids are indexed differently and it is not symmetric:
    # THIS sweep is [ξ, x] (row = ξ), run A is [β, x] (row = β). So the
    # reference entry is ref[β_index, x_index] = [jb, jr], NOT [jr, jb].
    # Getting it backwards reads a different β at a different x and reports
    # mismatches of order the observable itself -- which is exactly what it did
    # on the first run of this check.

    println("\nξ = 1 row vs run A's β = 1 row ($(length(cols))/$(length(x_vals)) x values in common):")
    println("  (identical physical parameters: κ_1 = κ_2 = 2, g_1 = g_2 = 1, η = √x)")
    for k in sort(collect(keys(outs)))
        haskey(ref.outs, k) || continue
        worst = 0.0
        for (j, jr) in cols
            a, b = outs[k][i1, j], ref.outs[k][jb, jr]
            (isfinite(a) && isfinite(b)) || continue
            worst = max(worst, abs(a - b))
        end
        tol = endswith(String(k), "_adia") ? tol_adia : tol_full
        @printf("  %-22s  worst |Δ| = %.3e   %s\n", String(k), worst,
                worst <= tol ? "reproduces run A" : "MISMATCH (tol $(tol))")
        worst <= tol ||
            @warn "$k does not reproduce run A along ξ = 1 (worst |Δ| = $worst). " *
                  "These are the same physical parameters, so this is a real " *
                  "disagreement -- suspect xxi_to_physical or the truncation."
    end
    return nothing
end


# ==============================================================================
# MAIN
# ==============================================================================
begin
    obs_list = select_observables(ACTIVE_OBSERVABLES_XXI)
    labels   = obs_labels(obs_list)

    res  = run_sweep_xxi(; PARAMS_XXI...)
    outs = analyze(res, obs_list)

    summarize(outs)
    check_adia_x_only(outs, res.x_vals, res.params)
    check_against_run_A(outs, res.ξ_vals, res.x_vals, res.params)

    # `params` is the parameter set that actually produced these states --
    # read plot axes and titles from here, not from config_x_xi.jl, which may
    # have been edited since.
    full_data   = res.full
    adia_data   = res.adia
    ξ_vals      = res.ξ_vals
    x_vals      = res.x_vals
    params      = res.params
    N_1_grid    = res.N_1_grid
    N_2_grid    = res.N_2_grid
    dim_grid    = res.dim_grid
    time_grid   = res.time_grid
    status_grid = res.status_grid

    @save OUTFILE_XXI outs labels full_data adia_data ξ_vals x_vals params N_1_grid N_2_grid dim_grid time_grid status_grid
    println("\nSaved to $OUTFILE_XXI")

    # The checkpoint has served its purpose; leaving it behind invites loading
    # a stale partial instead of the real thing.
    if isfile(PARTIAL_XXI)
        rm(PARTIAL_XXI)
        println("Removed checkpoint $PARTIAL_XXI")
    end
end
