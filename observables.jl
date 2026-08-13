# ==============================================================================
# observables.jl
#
# Pure functions of the stored 4x4 reduced states, plus the registry that turns
# a name in a config file into a computed grid.
#
# Lifted verbatim out of Section 2 of run_jc_pump_disp_asy.jl, which had always
# been written to allow it: nothing here reaches back into the sweep except
# through the returned `res` object. It moved because there are now TWO sweep
# drivers -- run_jc_pump_disp_asy.jl over (β, x) and run_r_xi.jl over (r, ξ) --
# and neither can include the other to borrow these definitions, since both run
# a full sweep on include.
#
# This file DEFINES things and runs nothing. It requires LinearAlgebra (tr,
# eigvals, Hermitian) and Printf (@printf) to be in scope; both drivers `using`
# them before including it.
#
# To add an observable: write the function, add one entry to
# OBSERVABLE_REGISTRY, and name it in the config file of whichever sweep you
# want it in. Nothing else changes -- and it becomes available to both sweeps.
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
# The config files select from it by name only, so they have no dependencies.
const OBSERVABLE_REGISTRY = [
    Obs(:concurrence, :single,  calc_concurrence, "Concurrence"),
    Obs(:purity,      :single,  calc_purity,      "Purity  tr(rho^2)"),
    Obs(:tracedist,   :compare, calc_tracedist,   "Trace distance"),
]

# ------------------------------------------------------------------------
# select_observables: resolve names from a config file into registry entries.
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
#
# Axis-agnostic: it reads only res.full / res.adia and their shape, which is
# why both the (β, x) and the (r, ξ) sweep can share it unchanged.
# ------------------------------------------------------------------------
function analyze(res, obs_list)
    N_row, N_col = size(res.full)

    all_keys = Symbol[]
    for o in obs_list
        append!(all_keys, obs_keys(o))
    end
    outs = Dict(k => fill(NaN, N_row, N_col) for k in all_keys)

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
