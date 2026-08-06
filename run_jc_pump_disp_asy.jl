using Base.Threads
using Printf
using LinearAlgebra
using JLD2
using Dates
using Statistics

include("config_jc_pump_disp_asy.jl")
include("truncation.jl")
include("jc_pump_disp_asy.jl")

# ==============================================================================
# run_jc_pump_disp_asy.jl
#
# Heavy lifting only. Nothing tunable lives here -- physical parameters, grid
# size, observable selection and output filename all come from config.jl.
#
# Two sections, separated by the ===== banner below:
#
#   1. run_sweep(...)  -- generates a dimensionless (β, x) grid, converts each
#                         point to physical parameters, estimates the Fock
#                         truncation, solves both models, and stores the two
#                         reduced 2-qubit states as raw 4x4 matrices.
#
#   2. observables     -- pure functions of those stored matrices. Nothing in
#                         this section reaches back into section 1 except
#                         through the returned `res` object, so it can be cut
#                         into its own file later with no edits.
# ==============================================================================

# Disable BLAS multithreading to prevent CPU oversubscription during our
# own @threads loop.
BLAS.set_num_threads(1)


# ==============================================================================
# run_sweep with timing / ETA
#
# Drop-in replacement for run_sweep in run_jc_pump_disp_asy.jl. Same
# signature, same return value plus two extra fields (time_grid, dim_grid).
#
# Three changes beyond the timing itself, all needed to make an ETA mean
# anything:
#
# 1. TWO PASSES. estimate_truncation is cheap relative to the steady-state
#    solve, so it runs for every point FIRST. That gives the Hilbert-space
#    dimension of every solve before any of them start -- so the total work
#    is known up front rather than discovered as it goes.
#
# 2. COST-WEIGHTED ETA. Points differ by orders of magnitude in cost, so
#    "fraction of points done" is a poor progress measure. Cost is modelled
#    as dim^COST_EXPONENT with a single scale factor fitted from the points
#    already finished, and the ETA is over remaining COST, not count.
#
# 3. DYNAMIC WORK QUEUE, LONGEST FIRST. @threads splits the index range into
#    contiguous chunks, and the expensive points are contiguous (the whole
#    last x column), so one thread would get all of them. An atomic work
#    queue with points sorted expensive-first balances the load instead.
# ==============================================================================

# Cost is assumed to scale as dim^COST_EXPONENT. dim^2 is the memory
# scaling of a dense density matrix and a reasonable first guess for the
# iterative solve; the true exponent depends on Krylov iteration counts,
# which I have not measured. The fitted scale factor absorbs the constant,
# and the end-of-run report prints the measured scaling so this can be
# tuned against real data.
const COST_EXPONENT = 2.0

# Points at or above BIG_DIM are memory-heavy enough to contend with each
# other; at most MAX_BIG of them run concurrently. dim 2000 sits between
# the 900 cluster and the 2704+ cluster of the current grid.
const BIG_DIM = 2000
const MAX_BIG = 2

fmt_time(s) = s < 60  ? @sprintf("%.0fs", s) :
              s < 3600 ? @sprintf("%dm%02ds", floor(Int, s/60), round(Int, s%60)) :
                         @sprintf("%dh%02dm", floor(Int, s/3600), floor(Int, (s%3600)/60))

nanmat() = fill(ComplexF64(NaN), 4, 4)

function run_sweep(; κ_1::Float64 = 2.0,
                     κ_2::Float64 = 2.0,
                     g_2::Float64 = 1.0,
                     γ_val::Float64 = 0.0,
                     γ_phi_val::Float64 = 0.0,
                     len::Int = 7,
                     x_range = (0.01, 0.7225),
                     β_range = (0.01, 100.0))

    x_vals = collect(range(x_range..., length = len))
    β_vals = collect(exp10.(range(log10.(β_range)..., length = len)))

    N_β, N_x = length(β_vals), length(x_vals)
    total = N_β * N_x

    rho_full_grid = [nanmat() for _ in 1:N_β, _ in 1:N_x]
    rho_adia_grid = [nanmat() for _ in 1:N_β, _ in 1:N_x]

    N_1_grid    = zeros(Int, N_β, N_x)
    N_2_grid    = zeros(Int, N_β, N_x)
    dim_grid    = zeros(Int, N_β, N_x)
    time_grid   = fill(NaN, N_β, N_x)
    status_grid = fill(:ok, N_β, N_x)

    η_of(x)  = sqrt(x * κ_1 * κ_2 / 4)
    g1_of(β) = g_2 * sqrt(β * κ_1 / κ_2)

    println("Sweep on $(Threads.nthreads()) thread(s), $total grid points.")

    # ======================================================================
    # PASS 1 -- truncation only. Cheap, and it tells us the size of every
    # solve before committing to any of them.
    # ======================================================================
    print("Pass 1/2: estimating truncations... ")
    t_pass1 = @elapsed begin
        @threads for idx in CartesianIndices((N_β, N_x))
            i, j = idx.I
            try
                N_1, N_2 = estimate_truncation(κ_1, κ_2, η_of(x_vals[j]), g1_of(β_vals[i]), g_2)
                N_1_grid[i, j] = N_1
                N_2_grid[i, j] = N_2
                dim_grid[i, j] = (N_1 + 1) * (N_2 + 1) * 4
            catch e
                status_grid[i, j] = :truncation_failed
                @warn "truncation failed at (β=$(β_vals[i]), x=$(x_vals[j]))" exception = e
            end
        end
    end
    @printf("done in %s\n", fmt_time(t_pass1))

    # --- what we now know about the work ahead ---
    todo = [idx for idx in CartesianIndices((N_β, N_x)) if status_grid[idx] === :ok]
    isempty(todo) && error("every point failed truncation; nothing to solve")

    dims = [dim_grid[idx] for idx in todo]
    dmax = maximum(dims)
    mem_max = dmax^2 * 16 / 2^20          # dense rho, MiB

    println("  Hilbert dimensions: min $(minimum(dims)), median $(round(Int, median(dims))), max $dmax")
    @printf("  largest dense rho: %.0f MiB (%d concurrent solves possible)\n",
            mem_max, Threads.nthreads())
    mem_max * Threads.nthreads() > 4096 &&
        @warn "peak memory could exceed 4 GiB with all threads on large points"

    # Sort expensive-first. With a dynamic queue this is the longest-
    # processing-time-first heuristic, which keeps the tail from being one
    # huge job starting last.
    # Cost model. Two factors, both needed:
    #
    #   dim^COST_EXPONENT   -- work per Krylov iteration
    #   1/(1 - sqrt(x))     -- number of iterations
    #
    # The second is the important one and is derived, not fitted. Setting
    # g = 0, the drift matrix on (a_1, a_2†) is
    #
    #       M = -[ κ_1/2   η    ]
    #            [ η       κ_2/2]
    #
    # whose eigenvalues for κ_1 = κ_2 = κ are -(κ/2 ± η). The slowest
    # relaxation rate -- the Liouvillian gap -- is therefore κ/2 - η, and
    # with η = sqrt(x κ_1 κ_2)/2 = sqrt(x) κ/2 that is
    #
    #       gap = (κ/2)(1 - sqrt(x))
    #
    # Krylov iteration count scales as 1/gap, so cost picks up a factor
    # 1/(1 - sqrt(x)): critical slowing down as the OPO approaches
    # threshold. This dominates dimension in practice -- measured times at
    # fixed dim=900 rise by 2.4x from x=0.01 to x=0.37, matching the
    # predicted 2.28x.
    #
    # The κ/2 prefactor is a constant across the sweep and is absorbed into
    # the fitted scale factor, so only the (1 - sqrt(x)) part appears here.
    gap_factor(x) = 1.0 / max(1.0 - sqrt(x), 1e-3)   # guard x -> 1
    cost = Dict(idx => Float64(dim_grid[idx])^COST_EXPONENT * gap_factor(x_vals[idx.I[2]])
                for idx in todo)
    order = sort(todo; by = idx -> -cost[idx])
    cost_total = sum(values(cost))

    # ======================================================================
    # PASS 2 -- the actual solves, on TWO dynamic work queues.
    #
    # A single queue lets all threads grab the largest points at once,
    # because they sort first. steadystate.iterative holds Krylov vectors
    # of dim^2 complex numbers each -- 268 MiB apiece at dim 4096 -- so
    # four concurrent large solves thrash memory. Measured: identical work
    # cost 2h06m of summed CPU time on one run and 3h33m on another, purely
    # from contention. Total CPU time is scheduling-independent, so that
    # difference is stalling, not load imbalance.
    #
    # Splitting into a big queue (capped at MAX_BIG concurrent workers) and
    # a small queue keeps every core busy without four large workspaces
    # competing. Workers fall through to the other queue when theirs is
    # empty, so nothing idles at the end.
    # ======================================================================
    println("Pass 2/2: solving...")

    big_idx   = [idx for idx in order if dim_grid[idx] >= BIG_DIM]
    small_idx = [idx for idx in order if dim_grid[idx] <  BIG_DIM]
    n_big_workers = min(MAX_BIG, Threads.nthreads(), max(length(big_idx), 1))

    @printf("  %d large point(s) (dim >= %d) on %d worker(s); %d small on the rest\n",
            length(big_idx), BIG_DIM, n_big_workers, length(small_idx))

    next_big    = Atomic{Int}(0)
    next_small  = Atomic{Int}(0)
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

            η   = η_of(x_vals[j])
            g_1 = g1_of(β_vals[i])
            N_1, N_2 = N_1_grid[i, j], N_2_grid[i, j]

            # Announce on dispatch, not just on completion. The expensive
            # points run first and can take minutes, so without this the
            # display sits unchanged long enough to look hung.
            lock(print_lock) do
                push!(inflight, idx)
                @printf("  start  ..../%-2d  β=%8.3f  x=%.4f  N₁=%2d N₂=%2d  dim=%5d  (%d running)\n",
                        length(order), β_vals[i], x_vals[j], N_1, N_2,
                        dim_grid[i, j], length(inflight))
            end

            dt = @elapsed try
                rho_full, rho_adia = run_sim(g_1, g_2, κ_1, κ_2, η,
                                             γ_val, γ_phi_val, N_1, N_2)
                @assert size(rho_full.data) == (4, 4) "rho_full is not a 2-qubit state"
                @assert size(rho_adia.data) == (4, 4) "rho_adia is not a 2-qubit state"
                rho_full_grid[i, j] = Matrix{ComplexF64}(rho_full.data)
                rho_adia_grid[i, j] = Matrix{ComplexF64}(rho_adia.data)
            catch e
                status_grid[i, j] = :solve_failed
                lock(print_lock) do
                    print("\r\x1B[K")
                    @warn "solve failed at (β=$(β_vals[i]), x=$(x_vals[j]), dim=$(dim_grid[i,j]))" exception = e
                end
            end
            time_grid[i, j] = dt

            lock(print_lock) do
                delete!(inflight, idx)
                done_count += 1
                done_cost  += cost[idx]
                busy_time  += dt

                # ETA: fit one scale factor c from work already done
                # (c = busy_time / done_cost), apply it to the remaining
                # cost, then divide by thread count for wall-clock.
                # Unreliable until a few points of each size have finished.
                rem_cost = cost_total - done_cost
                eta = rem_cost * (busy_time / done_cost) / Threads.nthreads()
                pct = 100 * done_cost / cost_total

                # No separate progress bar: it was drawn with \r and no
                # newline, so the next start/done line overwrote it every
                # time. The numbers live on the done line instead, which
                # also leaves a readable scrollback of how long each point
                # took -- more useful here than a bar, given the spread.
                @printf("  done   %2d/%-2d  β=%8.3f  x=%.4f  dim=%5d  %8.2fs  |  %5.1f%%  ETA %s  (%d running)\n",
                        done_count, length(order), β_vals[i], x_vals[j],
                        dim_grid[i, j], dt, pct, fmt_time(eta), length(inflight))
            end
        end
    end

    wall = time() - t_start

    # ======================================================================
    # REPORT
    #
    # Wrapped in try/catch on purpose. Everything below is cosmetic, and an
    # exception here would propagate out of run_sweep BEFORE the return --
    # discarding a completed sweep because a summary line had a bug. That
    # is exactly what happened once. Never again: the report may fail, the
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

        # Measured cost scaling, so COST_EXPONENT can be checked against
        # reality rather than left as a guess. Slope of log(time) vs
        # log(dim) over points that took long enough to time meaningfully.
        pts = [(Float64(dim_grid[idx]), time_grid[idx])
               for idx in CartesianIndices(time_grid)
               if isfinite(time_grid[idx]) && time_grid[idx] > 0.05 && dim_grid[idx] > 0]
        if length(pts) ≥ 4
            lx = log.(first.(pts)); ly = log.(last.(pts))
            slope = sum((lx .- mean(lx)) .* (ly .- mean(ly))) / sum((lx .- mean(lx)).^2)
            @printf("  measured scaling: time ~ dim^%.2f (COST_EXPONENT is %.2f)\n",
                    slope, COST_EXPONENT)
        end

        # Where the time actually went.
        # vec() is essential: collect(CartesianIndices(...)) of an (N_β × N_x)
        # grid is itself an (N_β × N_x) MATRIX, and sort on a multidimensional
        # array demands a `dims` keyword.
        worst = sort(vec(collect(CartesianIndices(time_grid)));
                     by = idx -> isfinite(time_grid[idx]) ? -time_grid[idx] : 0.0)[1:min(3, total)]
        println("  slowest points:")
        for idx in worst
            i, j = idx.I
            isfinite(time_grid[idx]) || continue
            @printf("    β=%8.3f  x=%.4f  dim=%5d  %6.2fs\n",
                    β_vals[i], x_vals[j], dim_grid[i, j], time_grid[idx])
        end
    end

    n_bad = total - n_ok
    n_bad > 0 && @warn "$n_bad point(s) failed -- see status_grid"
    catch e
        @warn "report failed (results are unaffected and returned below)" exception = e
    end

    return (full = rho_full_grid,
            adia = rho_adia_grid,
            β_vals = β_vals,
            x_vals = x_vals,
            N_1_grid = N_1_grid,
            N_2_grid = N_2_grid,
            dim_grid = dim_grid,
            time_grid = time_grid,
            status_grid = status_grid,
            params = (; κ_1, κ_2, g_2, γ_val, γ_phi_val, len, x_range, β_range))
end


# ==============================================================================
# SECTION 2 -- OBSERVABLES
#
# Everything below is a pure function of the stored 4x4 matrices. To add an
# observable: write the function, add one entry to OBSERVABLE_REGISTRY, and
# name it in config.jl. Nothing else changes.
# ==============================================================================

# σ_y ⊗ σ_y, worked out explicitly. With σ_y = [0 -i; i 0] the off-diagonal
# blocks are -i σ_y = [0 -1; 1 0] and +i σ_y = [0 1; -1 0], so the whole
# thing is real -- no need to carry complex arithmetic or a sparse operator.
const YY = Float64[0  0  0 -1
                   0  0  1  0
                   0  1  0  0
                  -1  0  0  0]

# ------------------------------------------------------------------------
# calc_concurrence: Wootters' construction.
#
#   ρ̃ = (σ_y ⊗ σ_y) ρ* (σ_y ⊗ σ_y),  R = ρ ρ̃,
#   λ_i = sqrt of the eigenvalues of R, sorted descending,
#   C = max(0, λ_1 - λ_2 - λ_3 - λ_4).
#
# R is not Hermitian, so eigvals returns complex numbers; analytically they
# are real and non-negative (R is similar to √ρ ρ̃ √ρ, which is PSD), and
# abs() absorbs the numerical imaginary part without ever throwing.
#
# The leading normalization guards against a trace defect from the iterative
# steady-state solver leaking into C.
# ------------------------------------------------------------------------
function calc_concurrence(ρ_data)
    ρ = ρ_data / tr(ρ_data)
    ρ_conj  = conj(ρ)
    temp    = YY * ρ_conj
    ρ_tilde = temp * YY
    R       = ρ * ρ_tilde
    evals   = eigvals(R)
    λ       = sort!(sqrt.(abs.(evals)), rev = true)
    return max(0.0, λ[1] - λ[2] - λ[3] - λ[4])
end

# ------------------------------------------------------------------------
# calc_purity: tr(ρ²), normalization-independent.
# ------------------------------------------------------------------------
calc_purity(ρ_data) = real(tr(ρ_data * ρ_data)) / real(tr(ρ_data))^2

# ------------------------------------------------------------------------
# calc_tracedist: D = ½‖ρ - σ‖₁ = ½ Σ|λ_i| over eigenvalues of the
# difference. Both states are normalized first, so a trace defect in either
# solver output cannot masquerade as physical disagreement. Hermitian()
# routes to the symmetric eigensolver, returning real eigenvalues directly.
# ------------------------------------------------------------------------
function calc_tracedist(ρ_full, ρ_adia)
    Δ = ρ_full / tr(ρ_full) - ρ_adia / tr(ρ_adia)
    return 0.5 * sum(abs, eigvals(Hermitian(Δ)))
end

# ------------------------------------------------------------------------
# Obs: an observable is data, not a code path.
#
#   kind = :single   -> f(ρ),        auto-expanded to _full / _adia / _diff
#   kind = :compare  -> f(ρf, ρa),   stored under `name` alone
# ------------------------------------------------------------------------
struct Obs
    name::Symbol
    kind::Symbol
    f::Function
    label::String
end

# The registry lives here, next to the function definitions it references.
# config.jl selects from it by name only, so config has no dependencies.
const OBSERVABLE_REGISTRY = [
    Obs(:concurrence, :single,  calc_concurrence, "Concurrence"),
    Obs(:purity,      :single,  calc_purity,      "Purity  tr(rho^2)"),
    Obs(:tracedist,   :compare, calc_tracedist,   "Trace distance"),
]

# ------------------------------------------------------------------------
# select_observables: resolve names from config.jl into registry entries.
# Errors loudly on an unknown name rather than silently computing fewer
# observables than asked for.
# ------------------------------------------------------------------------
function select_observables(names, registry = OBSERVABLE_REGISTRY)
    by_name = Dict(o.name => o for o in registry)
    selected = Obs[]
    for n in names
        haskey(by_name, n) ||
            error("unknown observable :$n -- available: " *
                  join(string.(getfield.(registry, :name)), ", "))
        push!(selected, by_name[n])
    end
    isempty(selected) && @warn "no observables selected"
    return selected
end

# Keys each observable contributes to the output Dict.
function obs_keys(o::Obs)
    if o.kind === :single
        return [Symbol(o.name, :_full), Symbol(o.name, :_adia), Symbol(o.name, :_diff)]
    elseif o.kind === :compare
        return [o.name]
    else
        error("unknown observable kind: $(o.kind)")
    end
end

# Plot labels, generated from the same source as the data so the two cannot
# drift apart.
function obs_labels(obs_list)
    d = Dict{Symbol,String}()
    for o in obs_list
        if o.kind === :single
            d[Symbol(o.name, :_full)] = o.label * " (full)"
            d[Symbol(o.name, :_adia)] = o.label * " (adiabatic)"
            d[Symbol(o.name, :_diff)] = "delta " * o.label
        else
            d[o.name] = o.label
        end
    end
    return d
end

# ------------------------------------------------------------------------
# analyze: loops grid points × observables and fills one Float64 matrix per
# generated key. Single-threaded on purpose -- 49 points of 4x4 linear
# algebra is microseconds, and threading would only reintroduce the
# whole-loop-teardown failure mode.
#
# Grids are prefilled with NaN, so any point that throws simply stays NaN
# and shows up as a gap on the heatmap instead of killing the run.
# ------------------------------------------------------------------------
function analyze(res, obs_list)
    N_β, N_x = size(res.full)

    all_keys = Symbol[]
    for o in obs_list
        append!(all_keys, obs_keys(o))
    end
    outs = Dict(k => fill(NaN, N_β, N_x) for k in all_keys)

    n_failed = 0

    for idx in CartesianIndices(res.full)
        ρf = res.full[idx]
        ρa = res.adia[idx]

        for o in obs_list
            try
                if o.kind === :single
                    vf = real(o.f(ρf))
                    va = real(o.f(ρa))
                    outs[Symbol(o.name, :_full)][idx] = vf
                    outs[Symbol(o.name, :_adia)][idx] = va
                    outs[Symbol(o.name, :_diff)][idx] = vf - va
                else
                    outs[o.name][idx] = real(o.f(ρf, ρa))
                end
            catch e
                n_failed += 1
                @warn "observable :$(o.name) failed at index $(Tuple(idx))" exception = e
            end
        end
    end

    n_failed > 0 && @warn "$n_failed observable evaluation(s) failed; those entries are NaN"

    return outs
end

# ------------------------------------------------------------------------
# summarize: quick console readout so a run tells you something before you
# get around to plotting.
# ------------------------------------------------------------------------
function summarize(outs)
    println("\nObservable summary (NaN entries excluded):")
    for k in sort(collect(keys(outs)))
        v = filter(isfinite, vec(outs[k]))
        if isempty(v)
            @printf("  %-22s  all NaN\n", String(k))
        else
            @printf("  %-22s  min = %+.6e   max = %+.6e\n", String(k), minimum(v), maximum(v))
        end
    end
end


# ==============================================================================
# MAIN
# ==============================================================================
begin
    obs_list = select_observables(ACTIVE_OBSERVABLES)
    labels   = obs_labels(obs_list)

    res  = run_sweep(; PARAMS...)
    outs = analyze(res, obs_list)

    summarize(outs)

    # Saving the states alongside the observables costs ~25 kB and means
    # Section 2 can later be lifted into its own file and re-run without
    # repeating the sweep. `params` is the parameter set that actually
    # produced these states -- read plot axes from here, not from config.jl,
    # which may have been edited since.
    full_data   = res.full
    adia_data   = res.adia
    β_vals      = res.β_vals
    x_vals      = res.x_vals
    params      = res.params
    N_1_grid    = res.N_1_grid
    N_2_grid    = res.N_2_grid
    dim_grid    = res.dim_grid
    time_grid   = res.time_grid
    status_grid = res.status_grid

    @save OUTFILE outs labels full_data adia_data β_vals x_vals params N_1_grid N_2_grid dim_grid time_grid status_grid
    println("\nSaved to $OUTFILE")
end