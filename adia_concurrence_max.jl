using LinearAlgebra
using Printf

include("jc_pump_disp_asy.jl")
include("observables.jl")

# ==============================================================================
# adia_concurrence_max.jl
#
#   julia --project=. adia_concurrence_max.jl
#   julia --project=. adia_concurrence_max.jl x_lo x_hi β_lo β_hi len
#
# Where in the (x, β) plane is the ADIABATIC model's concurrence largest?
#
# Answer it yourself, from this repo's own equations: every number below comes
# from run_adia in jc_pump_disp_asy.jl -- the same H_adia and J_adia the sweeps
# solve -- so if the code says something different from the claim, the code
# wins. There is no closed form in the loop and no .jld2 is read (except for one
# optional cross-check at the end).
#
# NO SWEEP, NO SOLVE OF THE FULL MODEL, NO THREADS NEEDED. The adiabatic
# Liouvillian is 4x4; a point costs microseconds, so the scan below is
# thousands of points in a second or two. That is what the run_adia extraction
# bought -- before it, one point meant a dim-900 Krylov solve as well.
#
# What it prints, in order:
#   1. a coarse (β, x) table of C_adia, with the largest cell marked
#   2. a golden-section refinement of that cell -> (x*, β*, C*)
#   3. four independent checks that the answer is real and not an artifact
#
# x < 1 is required (the adiabatic denominators carry 4η² - κ₁κ₂ and blow up at
# x = 1), and C_adia is zero above x* = 0.4839 anyway, so the default window
# stops well short.
# ==============================================================================

BLAS.set_num_threads(1)

# ------------------------------------------------------------------------
# Search window. Command-line args override, so you can zoom in without
# editing the file:
#
#   julia --project=. adia_concurrence_max.jl 0.10 0.20 0.8 1.25 21
#
# Defaults span the whole entangled region: x up to the x* = 0.4839 threshold
# and β over two decades, log-spaced (β and 1/β are physically equivalent
# arms, so a log axis is the symmetric one).
# ------------------------------------------------------------------------
const X_LO, X_HI = length(ARGS) >= 2 ? (parse(Float64, ARGS[1]), parse(Float64, ARGS[2])) : (0.01, 0.48)
const Β_LO, Β_HI = length(ARGS) >= 4 ? (parse(Float64, ARGS[3]), parse(Float64, ARGS[4])) : (0.1, 10.0)
const LEN        = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 25

# ------------------------------------------------------------------------
# The raw parameters the dimensionless point is realized with.
#
# ρ_adia depends only on (x, β) -- proven in CLAUDE.md, and check 2 below
# re-measures it -- so these three numbers cannot affect any answer here.
# They are run A's values, which makes the optional check 4 possible.
# ------------------------------------------------------------------------
const Κ_1, Κ_2, G_2 = 2.0, 2.0, 1.0
const Γ, Γ_ΦI = 0.0, 0.0     # the (x, β)-only argument needs both zero

# ------------------------------------------------------------------------
# C_of: the single function under test. (x, β) -> concurrence of ρ_adia.
#
# Inversion is sweep 1's, verbatim (run_jc_pump_disp_asy.jl:107-108):
#     η_of(x)  = sqrt(x κ₁ κ₂ / 4)      g1_of(β) = g₂ sqrt(β κ₁ / κ₂)
# then run_adia solves, and calc_concurrence from observables.jl scores it --
# the same Wootters routine the sweeps report.
# ------------------------------------------------------------------------
η_of(x)  = sqrt(x * Κ_1 * Κ_2 / 4)
g1_of(β) = G_2 * sqrt(β * Κ_1 / Κ_2)

function C_of(x, β; κ_1 = Κ_1, κ_2 = Κ_2, g_2 = G_2)
    η   = sqrt(x * κ_1 * κ_2 / 4)
    g_1 = g_2 * sqrt(β * κ_1 / κ_2)
    ρ   = run_adia(g_1, g_2, κ_1, κ_2, η, Γ, Γ_ΦI)
    return calc_concurrence(ρ.data)
end

# ------------------------------------------------------------------------
# gss_max: golden-section maximization of f on [a, b].
#
# Derivative-free and needs no packages, which is the point -- it makes the
# refinement below auditable by eye. It assumes ONE interior maximum on the
# bracket; that is true here (the table in step 1 shows it), but if you widen
# the window to include the C = 0 dead zone, read the table first rather than
# trusting the refinement, since a flat zero region has no unique argmax.
# ------------------------------------------------------------------------
function gss_max(f, a, b; tol = 1e-12, maxit = 300)
    φ = (sqrt(5.0) - 1) / 2
    c, d = b - φ * (b - a), a + φ * (b - a)
    fc, fd = f(c), f(d)
    for _ in 1:maxit
        if fc > fd
            b, d, fd = d, c, fc
            c = b - φ * (b - a); fc = f(c)
        else
            a, c, fc = c, d, fd
            d = a + φ * (b - a); fd = f(d)
        end
        (b - a) < tol && break
    end
    m = (a + b) / 2
    return (m, f(m))
end

# best_x: at fixed β, the optimal x and the C there. Used both on its own and
# as the inner problem of the 2D refinement.
best_x(β; lo = X_LO, hi = X_HI) = gss_max(x -> C_of(x, β), lo, hi)

# ==============================================================================
# 1. COARSE SCAN -- the picture. Rows are β (log-spaced), columns are x
#    (linear), matching the grid convention everything else in the repo uses.
# ==============================================================================
x_vals = collect(range(X_LO, X_HI; length = LEN))
β_vals = exp10.(range(log10(Β_LO), log10(Β_HI); length = LEN))

C = [C_of(x, β) for β in β_vals, x in x_vals]
imax = argmax(C)

println("="^78)
println("ADIABATIC CONCURRENCE over (x, β) -- from run_adia, this repo's own model")
println("="^78)
@printf("realized at κ₁ = %.4g, κ₂ = %.4g, g₂ = %.4g, γ = γ_φ = %.4g\n", Κ_1, Κ_2, G_2, Γ)
@printf("window: x ∈ [%.4g, %.4g] linear, β ∈ [%.4g, %.4g] log, %d × %d = %d points\n",
        X_LO, X_HI, Β_LO, Β_HI, LEN, LEN, LEN^2)

# Table, thinned to at most ~9 rows and ~9 columns so it stays on screen. The
# scan itself is at full resolution; only the printing is subsampled.
row_step = max(1, cld(LEN, 9))
col_step = max(1, cld(LEN, 9))
rows = unique(vcat(1:row_step:LEN, imax[1], LEN))  |> sort
cols = unique(vcat(1:col_step:LEN, imax[2], LEN))  |> sort

println("\nC_adia   (rows = β, cols = x; * marks the grid maximum)")
print("   β \\ x ")
for j in cols; @printf("%9.4f", x_vals[j]); end
println()
for i in rows
    @printf("%8.4g ", β_vals[i])
    for j in cols
        mark = (i == imax[1] && j == imax[2]) ? "*" : " "
        @printf("%8.4f%s", C[i, j], mark)
    end
    println()
end

@printf("\ngrid maximum: C_adia = %.10f  at  x = %.6f, β = %.6f\n",
        C[imax], x_vals[imax[2]], β_vals[imax[1]])
@printf("zero cells (C_adia = 0 exactly): %d of %d\n", count(iszero, C), length(C))

# ==============================================================================
# 2. REFINE -- the grid maximum is only as good as the grid spacing, so
#    golden-section it. Outer search over β, inner search over x, each inner
#    call a full 1D optimization: this finds a genuine 2D maximum, not the best
#    x along one guessed β.
# ==============================================================================
println("\n" * "="^78)
println("REFINEMENT (golden section, nested: best x at each β, then best β)")
println("="^78)

β_star, C_star = gss_max(β -> best_x(β)[2], Β_LO, Β_HI)
x_star, _      = best_x(β_star)

@printf("β* = %.9f\n", β_star)
@printf("x* = %.9f\n", x_star)
@printf("C* = %.10f\n", C_star)
@printf("\nη*/η_threshold = sqrt(x*) = %.6f   (so the optimum is %.1f%% of the way to threshold)\n",
        sqrt(x_star), 100 * sqrt(x_star))
@printf("realized as: η = %.6f, g₁ = %.6f, g₂ = %.6f, κ₁ = κ₂ = %.6f\n",
        η_of(x_star), g1_of(β_star), G_2, Κ_1)

# ==============================================================================
# 3. CHECKS. A bare argmax is not evidence -- these four are what make it one.
# ==============================================================================
println("\n" * "="^78)
println("CHECKS")
println("="^78)

# (a) Is it actually a peak? Step off it in both coordinates and both
#     directions. All four neighbours must come out LOWER. This is the check
#     that a converged-looking optimizer on a monotone function would fail.
println("\n(a) local peak test -- all four neighbours must be lower")
for (lbl, xx, bb) in (("x* - 1%", 0.99x_star, β_star), ("x* + 1%", 1.01x_star, β_star),
                      ("β* - 1%", x_star, 0.99β_star), ("β* + 1%", x_star, 1.01β_star))
    c = C_of(xx, bb)
    @printf("    %-9s  C = %.10f   ΔC = %+.3e  %s\n",
            lbl, c, c - C_star, c < C_star ? "lower  ok" : "HIGHER -- not a maximum!")
end

# (b) Does the answer depend on the raw parameters? It must not: with
#     γ = γ_φ = 0 the adiabatic Liouvillian is homogeneous in the couplings and
#     ρ_adia is a function of (x, β) alone (CLAUDE.md). Three unrelated (κ₁, κ₂,
#     g₂) at the SAME (x*, β*) must give the same C to ~1e-13. If they do not,
#     the search was over a different surface for each choice and step 2 is
#     meaningless.
println("\n(b) (x, β)-only invariance -- three raw realizations of the same point")
ref = C_of(x_star, β_star)
for (κ1, κ2, g2) in ((2.0, 2.0, 1.0), (2.5, 1.5, 0.75), (7.3, 0.41, 3.1))
    c = C_of(x_star, β_star; κ_1 = κ1, κ_2 = κ2, g_2 = g2)
    @printf("    κ₁ = %5.2f  κ₂ = %5.2f  g₂ = %5.2f   C = %.13f   Δ = %+.2e\n",
            κ1, κ2, g2, c, c - ref)
end

# (c) The arm-swap symmetry C(x, β) = C(x, 1/β). It is exact, it was not
#     imposed anywhere above, and it forces any maximum over β to sit at
#     β = 1 or to come in a mirror pair -- so β* = 1 is the symmetric answer
#     rather than a coincidence.
println("\n(c) arm-swap symmetry  C(x, β) = C(x, 1/β)  at x = x*")
for β in (0.5, 0.8, 1.25, 2.0, 3.0)
    c1, c2 = C_of(x_star, β), C_of(x_star, 1 / β)
    @printf("    β = %5.3f: C = %.10f    1/β = %6.3f: C = %.10f   Δ = %+.2e\n",
            β, c1, 1 / β, c2, c1 - c2)
end

# (d) Optional: agreement with the recorded sweep. Run A's grid has a column at
#     x = 0.12875 and a row at β = 1, so its stored concurrence_adia there must
#     match C_of at the same point. This is the one place a .jld2 is read, and
#     it is what ties this script to the numbers already in the repo. Skipped
#     silently if the file has been deleted or moved.
const RUN_A = "results_k1_2.0_k2_2.0_g2_1.0.jld2"
println("\n(d) against the recorded sweep, $RUN_A")
if isfile(RUN_A)
    using JLD2
    d = JLD2.load(RUN_A)
    Cg, βg, xg = d["outs"][:concurrence_adia], d["β_vals"], d["x_vals"]
    worst = maximum(abs(Cg[i, j] - C_of(xg[j], βg[i]))
                    for i in eachindex(βg), j in eachindex(xg) if isfinite(Cg[i, j]))
    ig = argmax(Cg)
    @printf("    stored grid max      C = %.13f  at x = %.5f, β = %.5f\n",
            Cg[ig], xg[ig[2]], βg[ig[1]])
    @printf("    recomputed here      C = %.13f  at the same grid point\n",
            C_of(xg[ig[2]], βg[ig[1]]))
    @printf("    worst |stored - recomputed| over all %d grid points: %.2e\n",
            length(Cg), worst)
    @printf("    and C* = %.10f from step 2 %s the stored grid max\n",
            C_star, C_star >= Cg[ig] ? "≥" : "< (!!)")
    println("    (the stored grid has no column at x*, so it should land just below C*)")
else
    println("    file not found -- skipped")
end

println("\n" * "="^78)
@printf("VERDICT: C_adia is maximal at x = %.6f, β = %.6f, where C = %.10f\n",
        x_star, β_star, C_star)
println("="^78)
