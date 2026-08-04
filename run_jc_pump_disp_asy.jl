using Base.Threads
using Printf
using LinearAlgebra
using JLD2

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
# SECTION 1 -- SWEEP
# ==============================================================================

nanmat() = fill(ComplexF64(NaN), 4, 4)

# ------------------------------------------------------------------------
# run_sweep: every parameter is a keyword argument, supplied by splatting
# PARAMS from config.jl. Grids and storage are function-locals, which also
# removes the boxed-global access the @threads closure used to incur.
#
# Returns a NamedTuple carrying the states together with the parameters that
# produced them, so the two can never desync.
# ------------------------------------------------------------------------
function run_sweep(; κ_1::Float64,
                     κ_2::Float64,
                     g_2::Float64,
                     γ_val::Float64,
                     γ_phi_val::Float64,
                     len::Int,
                     x_range,
                     β_range)

    # --- Dimensionless grid ---
    x_vals = collect(range(x_range..., length = len))
    β_vals = collect(exp10.(range(log10.(β_range)..., length = len)))

    N_β = length(β_vals)
    N_x = length(x_vals)

    # --- Result storage ---
    # Raw 4x4 ComplexF64 matrices, not Operator objects: concretely typed,
    # and they survive serialization independently of the QuantumOptics
    # version. Prefilled with NaN so a failed grid point leaves a readable
    # marker rather than an #undef that throws on access downstream.
    rho_full_grid = [nanmat() for _ in 1:N_β, _ in 1:N_x]
    rho_adia_grid = [nanmat() for _ in 1:N_β, _ in 1:N_x]

    N_1_grid    = zeros(Int, N_β, N_x)
    N_2_grid    = zeros(Int, N_β, N_x)
    status_grid = fill(:ok, N_β, N_x)

    total_iters = N_β * N_x
    update_interval = max(1, total_iters ÷ 100)
    counter = Atomic{Int}(0)
    print_lock = ReentrantLock()
    bar_width = 30

    println("Starting multithreaded sweep on $(Threads.nthreads()) threads...")
    println("Sweeping β ∈ [$(β_vals[1]), $(β_vals[end])] and x ∈ [$(x_vals[1]), $(x_vals[end])]")
    println("Total grid points: $total_iters\n")

    grid = CartesianIndices((N_β, N_x))

    @threads for idx in grid
        i, j = idx.I
        β = β_vals[i]
        x = x_vals[j]

        # Convert dimensionless -> physical parameters
        η   = sqrt(x * κ_1 * κ_2 / 4)
        g_1 = g_2 * sqrt(β * κ_1 / κ_2)

        N_1 = 0
        N_2 = 0

        # A throw inside @threads aborts the whole loop and takes every
        # completed point with it. Contain it per-point instead.
        try
            N_1, N_2 = estimate_truncation(κ_1, κ_2, η, g_1, g_2)
            N_1_grid[i, j] = N_1
            N_2_grid[i, j] = N_2

            rho_full, rho_adia = run_sim(g_1, g_2, κ_1, κ_2, η,
                                         γ_val, γ_phi_val, N_1, N_2)

            # Both must land in the same 2-qubit basis for the compare
            # observables to mean anything.
            @assert size(rho_full.data) == (4, 4) "rho_full is not a 2-qubit state"
            @assert size(rho_adia.data) == (4, 4) "rho_adia is not a 2-qubit state"

            rho_full_grid[i, j] = Matrix{ComplexF64}(rho_full.data)
            rho_adia_grid[i, j] = Matrix{ComplexF64}(rho_adia.data)
        catch e
            status_grid[i, j] = :failed
            lock(print_lock) do
                print("\r\x1B[K")
                @warn "grid point (β=$β, x=$x) failed" exception = e
            end
        end

        # Progress bar
        current = atomic_add!(counter, 1) + 1
        if current % update_interval == 0 || current == total_iters
            lock(print_lock) do
                pct = (current / total_iters) * 100
                filled = round(Int, pct / 100 * bar_width)
                bar = "[" * repeat("=", filled) * ">" * repeat(" ", bar_width - filled) * "]"
                print("\r\x1B[K")
                @printf("%s %5.1f%% | β=%.3f  x=%.4f  g₁=%.3f  η=%.4f  N₁=%d  N₂=%d",
                        bar, pct, β, x, g_1, η, N_1, N_2)
            end
        end
    end

    n_failed = count(==(:failed), status_grid)
    println("\n\nSweep complete. ", total_iters - n_failed, "/", total_iters, " points solved.")
    n_failed > 0 && @warn "$n_failed grid point(s) failed -- see status_grid"

    return (full = rho_full_grid,
            adia = rho_adia_grid,
            β_vals = β_vals,
            x_vals = x_vals,
            N_1_grid = N_1_grid,
            N_2_grid = N_2_grid,
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
if abspath(PROGRAM_FILE) == @__FILE__

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
    status_grid = res.status_grid

    @save OUTFILE outs labels full_data adia_data β_vals x_vals params N_1_grid N_2_grid status_grid
    println("\nSaved to $OUTFILE")
end