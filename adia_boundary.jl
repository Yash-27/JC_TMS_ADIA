using LinearAlgebra
using Printf

include("jc_pump_disp_asy.jl")
include("observables.jl")
include("plotting_functions.jl")

# ==============================================================================
# adia_boundary.jl
#
#   julia --project=. adia_boundary.jl
#   julia --project=. adia_boundary.jl x_lo x_hi N_x N_u
#
# THE ENTANGLEMENT BOUNDARY OF THE ADIABATIC MODEL, drawn and checked.
#
# file_concurrence/concurrence_adia.md proves (its §3.5) that C_adia > 0 exactly
# on
#           u ∈ [2, u_c(x))     and     x < x_c
#
# where u = β + 1/β ≥ 2 and u_c(x) is the closed form of its §4.1. This script
# puts that claim next to this repo's own numbers: it maps C_adia over the
# (x, u) plane with run_adia -- the same H_adia and J_adia every sweep solves --
# and overlays the analytic curve. If the coloured region does not stop exactly
# on the curve, the .md is wrong, or this repo is, and the disagreement is the
# result.
#
# NO SWEEP, NO FULL MODEL, NO THREADS. The adiabatic Liouvillian is 4x4 and one
# point costs ~120 µs measured, so the default 220 x 180 = 39,600-point grid is
# about 5 seconds. That is what the run_adia extraction bought.
#
# Three things it produces:
#   1. the figure, Adia_boundary/adia_boundary.pdf
#   2. a numerically bisected boundary, found WITHOUT reference to the closed
#      form, and the worst disagreement between the two
#   3. six checks (below) that make the agreement evidence rather than a
#      redrawing of the same formula twice
#
# WHY THE AXES ARE WHAT THEY ARE. u_c ~ 1/(4x²) as x -> 0, so u_c spans 2 to
# ~10⁴ over the interesting range and the y axis has to be log. x cannot start
# at 0 -- u_c diverges there and every denominator below carries 8x²q1 -- so the
# default x_lo is 0.02, where u_c = 577.75. On log-log the small-x asymptote is
# a straight line of slope -2, which is check (d).
# ==============================================================================

BLAS.set_num_threads(1)

# ------------------------------------------------------------------------
# Window. CLI args override so you can zoom without editing the file:
#
#   julia --project=. adia_boundary.jl 0.10 0.55 300 240
#
# x_hi defaults to x_c + 0.1, which is the point: three of the columns sit
# ABOVE the threshold and must come out uniformly zero.
# ------------------------------------------------------------------------
const X_LO = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 0.02
const N_X  = length(ARGS) >= 3 ? parse(Int,     ARGS[3]) : 220
const N_U  = length(ARGS) >= 4 ? parse(Int,     ARGS[4]) : 180

# ------------------------------------------------------------------------
# Raw parameters the dimensionless point is realized with.
#
# ρ_adia depends only on (x, β) -- proven in CLAUDE.md under the (r, ξ) sweep,
# and re-measured as check (f) below -- so these cannot affect any answer here.
# Run A's values, which is what lets check (f) quote run A's stored number.
# ------------------------------------------------------------------------
const Κ_1, Κ_2, G_2 = 2.0, 2.0, 1.0
const Γ, Γ_ΦI       = 0.0, 0.0        # the (x, β)-only argument needs both zero

# ------------------------------------------------------------------------
# C_of: (x, β) -> concurrence of ρ_adia. The single function under test.
#
# Inversion is sweep 1's, verbatim (run_jc_pump_disp_asy.jl:107-108). These
# five lines are duplicated from adia_concurrence_max.jl rather than included,
# because that file has no main() guard -- including it would run its whole
# scan. Same reason observables.jl exists. If they ever disagree, this one is
# wrong: adia_concurrence_max.jl is the older and better-checked copy.
# ------------------------------------------------------------------------
function C_of(x, β; κ_1 = Κ_1, κ_2 = Κ_2, g_2 = G_2)
    η   = sqrt(x * κ_1 * κ_2 / 4)
    g_1 = g_2 * sqrt(β * κ_1 / κ_2)
    ρ   = run_adia(g_1, g_2, κ_1, κ_2, η, Γ, Γ_ΦI)
    return calc_concurrence(ρ.data)
end

# u = β + 1/β with β ≥ 1 inverts to this. u < 2 is unphysical (no real β), and
# u = 2 is β = 1 exactly. The clamp absorbs roundoff at the u = 2 endpoint,
# where u^2 - 4 can come out as a small negative.
β_of_u(u) = (u + sqrt(max(u^2 - 4, 0.0))) / 2

C_of_u(x, u) = C_of(x, β_of_u(u))

# ==============================================================================
# THE CLOSED FORMS UNDER TEST -- concurrence_adia.md §3.3, §4.1, §10.7
#
# Transcribed from the .md and nothing else. No numerical constant from that
# file is hardcoded below: x_c and x* are re-derived from their polynomials, so
# a typo in either root shows up as a failed check rather than as agreement.
# ==============================================================================

# Auxiliary polynomials (§1.1) and the four bridging identities (§4.1).
q1_of(x) = (x - 1)^4 + 4x
p1_of(x) = (x - 1)^4 + 4x * (x + 1)
p2_of(x) = (x^2 - 1)^2 - 4x * (x + 1)

# u_c(x), §4.1: the parabola Φ's positive root, i.e. the sudden-death point.
#   u_c = [ B̃ + sqrt(B̃² + 16 x² q1 C̃0) ] / (8 x² q1)
# Returns NaN where the radicand goes negative (off the physical domain).
function uc_closed(x)
    q1, p1, p2 = q1_of(x), p1_of(x), p2_of(x)
    lam = q1 + 4x            # λ,  a + c = x λ
    kap = q1 - 4x            # κ,  a - c = x κ   ( = (1-x)^4 )
    sig = p1 + p2 - 2q1      # σ,  b     = x σ   ( = 4x(x²-2x-1) )
    tau = p1 + p2 - q1 + 4x  # τ,  g²    = x κ τ ( = (1-x²)²   )

    Btil  = kap * tau - x * sig * lam
    C0til = 2 * kap * tau - x * (kap^2 + sig^2)

    disc = Btil^2 + 16 * x^2 * q1 * C0til
    disc < 0 && return NaN
    return (Btil + sqrt(disc)) / (8 * x^2 * q1)
end

# C_adia(β = 1, x), §10.7 -- the closed form for the whole u = 2 edge.
C_at_u2(x) = (2 * sqrt(x) * (1 - x) - x * (1 + x)) / (2 * (x^2 + 1))

# ------------------------------------------------------------------------
# u*(x), §10.3 -- where G_u changes sign, i.e. the internal structure of §10's
# PROOF rather than of the physics. L := S'D - 2SD' is linear in u with positive
# slope, so it has exactly one zero, at the rational function below. NO SQUARE
# ROOT, unlike u_c.
#
# Why it is worth drawing: it splits the domain into the three regions the proof
# actually uses. For x > x_dagger, u*(x) < 2 -- off the physical domain -- so
# G_u ≥ 0 throughout and dC/du < 0 follows immediately from F_u < 0 (negative
# minus non-negative). Only the sliver x < x_dagger, u ∈ [2, u*(x)) has both
# derivatives negative and needed the degree-174 discriminant. That sliver is
# region 3, and on this figure it is visibly tiny.
# ------------------------------------------------------------------------
V2_10_3(x) = x^7 - 7x^6 + 19x^5 - 21x^4 + 7x^3 + 23x^2 + 5x + 5
V3_10_3(x) = x^9 - 13x^8 + 68x^7 - 160x^6 + 214x^5 - 154x^4 - 28x^3 - 56x^2 + x - 1
u_star(x)  = -V3_10_3(x) / (2x * V2_10_3(x))

# §10.4: L(2,x) changes sign at the unique root of V in (0,1). Equivalently the
# x at which u*(x) crosses 2 -- which is check (g) below, and a free test of the
# V2/V3 transcription against the independent V.
V_10_4(x) = x^7 - 7x^6 + 25x^5 - 27x^4 + 51x^3 + 3x^2 + 19x - 1

# A sign-change bisector. Used three times below, always on a function with one
# root in the bracket, so it needs nothing cleverer.
function bisect(f, a, b; tol = 1e-15, maxit = 200)
    fa = f(a)
    for _ in 1:maxit
        m = (a + b) / 2
        fm = f(m)
        if (fa < 0) == (fm < 0)
            a, fa = m, fm
        else
            b = m
        end
        (b - a) <= tol * max(1.0, abs(b)) && break
    end
    return (a + b) / 2
end

# x_c, §3.3: the unique real root of x³ - 2x² + 9x - 4, which lies in (0,1).
const X_C = bisect(x -> x^3 - 2x^2 + 9x - 4, 0.0, 1.0)

# x*, §10.7: (t*)² where t* is the root of R in (0,1). This is the argmax of
# C_adia along the u = 2 edge -- NOT the threshold. The .md and (since this
# session) CLAUDE.md both use x* for the argmax and x_c for the threshold.
const T_STAR = bisect(t -> t^6 + t^5 - 3t^4 - 2t^3 - 3t^2 - t + 1, 0.0, 1.0)
const X_STAR = T_STAR^2

const X_HI = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : X_C + 0.1

# x_dagger, §10.4: the unique real root of V in (0,1). Below it, region 3 exists.
const X_DAGGER = bisect(V_10_4, 0.0, 1.0)

# ==============================================================================
# 1. THE GRID -- C_adia over (x, u), from run_adia
#
# ROWS are u and COLUMNS are x, so the grid is [N_u, N_x]. That is the repo's
# convention read the usual way (see CLAUDE.md "Conventions"): the row
# coordinate is the one drawn on a log axis. Here the log axis is vertical
# rather than horizontal, so flatten_grid's output is handed to scatter with the
# two vectors swapped -- see the plot section, where that swap is the one place
# this script departs from plot_map.
# ==============================================================================
x_vals = collect(range(X_LO, X_HI; length = N_X))

# The u axis has to contain the whole curve, so its top comes from the tallest
# point of the boundary, which is at x_lo. 1.15 leaves headroom above it.
const U_TOP = 1.15 * uc_closed(X_LO)
u_vals = exp10.(range(log10(2.0), log10(U_TOP); length = N_U))

println("="^78)
println("ADIABATIC ENTANGLEMENT BOUNDARY -- run_adia vs concurrence_adia.md §4.1")
println("="^78)
@printf("realized at κ₁ = %.4g, κ₂ = %.4g, g₂ = %.4g, γ = γ_φ = %.4g\n", Κ_1, Κ_2, G_2, Γ)
@printf("x ∈ [%.4g, %.4g] linear × %d    u ∈ [2, %.4g] log × %d    = %d points\n",
        X_LO, X_HI, N_X, U_TOP, N_U, N_X * N_U)
@printf("x_c = %.12f   (root of x³-2x²+9x-4, §3.3)\n", X_C)
@printf("x*  = %.12f   (t*² with t* the root of R, §10.7 -- the ARGMAX, not the threshold)\n", X_STAR)

t_grid = @elapsed C = [C_of_u(x, u) for u in u_vals, x in x_vals]
@printf("\ngrid solved in %.2f s  (%.0f µs/point)\n", t_grid, 1e6 * t_grid / length(C))

n_zero = count(iszero, C)
@printf("C_adia = 0 exactly at %d of %d points (%.1f%%)\n",
        n_zero, length(C), 100 * n_zero / length(C))
@printf("C_adia max over the grid: %.10f\n", maximum(C))

# ==============================================================================
# 2. THE NUMERICAL BOUNDARY -- bisected, with the closed form kept out of it
#
# For each x column: is C > 0 at u = 2? If not, the whole column is dead and the
# boundary does not exist there. If so, bracket the crossing by DOUBLING u until
# C hits exactly zero, then geometrically bisect. Nothing here reads uc_closed,
# which is the point -- otherwise the comparison in check (c) would be a formula
# against itself.
#
# C vanishes LINEARLY in (u_c - u) (N = vg - √S is linear near its root), and
# calc_concurrence clamps at exactly 0 via max(0.0, ...), so the crossing is a
# clean step and bisection resolves it to the last few digits.
# ==============================================================================
function u_boundary_numeric(x)
    C_of_u(x, 2.0) > 0 || return NaN          # x ≥ x_c: no entangled region

    lo, hi = 2.0, 4.0
    while C_of_u(x, hi) > 0                   # bracket by doubling
        lo = hi
        hi *= 2
        hi > 1e9 && return NaN
    end
    for _ in 1:200                            # geometric bisection: u spans decades
        mid = sqrt(lo * hi)
        if C_of_u(x, mid) > 0
            lo = mid
        else
            hi = mid
        end
        (hi - lo) <= 1e-13 * hi && break
    end
    return (lo + hi) / 2
end

t_bis = @elapsed u_num = [u_boundary_numeric(x) for x in x_vals]
u_ana = [uc_closed(x) for x in x_vals]

@printf("boundary bisected for %d columns in %.2f s\n", count(isfinite, u_num), t_bis)

# ==============================================================================
# 3. CHECKS. Six of them, and they are the reason to run this rather than to
#    look at the picture.
# ==============================================================================
println("\n" * "="^78)
println("CHECKS")
println("="^78)

# (a) The closed form against the three values the .md itself records in §4.4,
#     which came from a completely different route (root-finding on N).
#
#     Tolerance is RELATIVE and loose (1e-4) because §4.4's table is printed to
#     5-6 significant figures, which is all it can be tested to. Note its first
#     row: it quotes 577.76 for the closed form, but this reproduces 577.75148,
#     which rounds to 577.75. Both of this script's independent routes (§4.1
#     here and the bisection in check (b)) agree on 577.751476, and §4.4's own
#     second column -- the older root-finding result -- says "~577.75" too. So
#     the 577.76 is a last-digit slip in that one table cell of the .md, not a
#     disagreement about u_c. Nothing depends on it; it is a display value.
println("\n(a) uc_closed against concurrence_adia.md §4.4's stored values")
for (x, want) in ((0.02, 577.76), (0.10, 18.4304), (0.40, 2.16302))
    got = uc_closed(x)
    rel = abs(got - want) / want
    @printf("    x = %.2f   u_c = %10.5f   .md says %10.5f   rel = %.2e %s\n",
            x, got, want, rel, rel < 1e-4 ? "ok" : "MISMATCH")
end
println("    (the x = 0.02 row: .md prints 577.76, this gives 577.75148 -> 577.75.")
println("     Last-digit slip in that cell; §4.4's own root-finding column agrees with us.)")

# (b) THE LOAD-BEARING CHECK. The numerically bisected boundary against the
#     closed form, column by column. These two share no code.
println("\n(b) bisected boundary vs closed form, over all live columns")
live = findall(i -> isfinite(u_num[i]) && isfinite(u_ana[i]), eachindex(u_num))
if isempty(live)
    println("    no live columns -- widen the window")
else
    rel = [abs(u_num[i] - u_ana[i]) / u_ana[i] for i in live]
    worst = argmax(rel)
    @printf("    live columns: %d   (x from %.4f to %.4f)\n",
            length(live), x_vals[live[1]], x_vals[live[end]])
    @printf("    worst relative error: %.3e  at x = %.5f  (u_num = %.9f, u_ana = %.9f)\n",
            rel[worst], x_vals[live[worst]], u_num[live[worst]], u_ana[live[worst]])
    @printf("    median relative error: %.3e\n", sort(rel)[cld(length(rel), 2)])
    println(maximum(rel) < 1e-6 ? "    -> the two agree; §4.1 reproduces run_adia" :
                                  "    -> DISAGREEMENT: investigate before trusting either")
end

# (c) The curve must meet the physical floor u = 2 exactly at x_c. Solve
#     u_c(x) = 2 numerically and compare against the cubic's root -- two
#     independent routes to the same threshold.
println("\n(c) u_c(x) = 2 must happen exactly at x_c")
x_from_curve = bisect(x -> uc_closed(x) - 2.0, X_LO, 0.99)
@printf("    u_c(x_c)              = %.12f   (want exactly 2)\n", uc_closed(X_C))
@printf("    root of u_c(x) - 2    = %.12f\n", x_from_curve)
@printf("    x_c from the cubic    = %.12f   Δ = %+.2e\n", X_C, x_from_curve - X_C)

# (d) The small-x asymptote of §4.5, u_c = 1/(4x²) - 1/x + 4 + O(x). On log-log
#     this is the straight line of slope -2 that the left end of the figure
#     should lie on.
println("\n(d) small-x asymptote  u_c ≈ 1/(4x²) - 1/x + 4   (§4.5)")
for x in (0.005, 0.01, 0.02, 0.05)
    a = 1 / (4x^2) - 1 / x + 4
    @printf("    x = %.3f   u_c = %11.4f   asymptote = %11.4f   ratio = %.6f\n",
            x, uc_closed(x), a, uc_closed(x) / a)
end

# (e) Above x_c the model must give zero at EVERY u, not merely at u = 2. This
#     is the "no revival" half of §3.5 -- the half that needed x_v < x_c -- and
#     the grid columns beyond x_c are a direct test of it.
println("\n(e) no revival: every column with x > x_c must be identically zero")
dead_cols = findall(>(X_C), x_vals)
if isempty(dead_cols)
    println("    none in window (raise x_hi above x_c to exercise this)")
else
    worst_dead = maximum(maximum(@view C[:, j]) for j in dead_cols)
    @printf("    %d columns above x_c, %d points each\n", length(dead_cols), N_U)
    @printf("    largest C_adia anywhere in them: %.3e %s\n",
            worst_dead, worst_dead == 0 ? "(exactly zero)" : "(NOT zero -- revival?)")
end

# (f) Tie back to numbers already in the repo, so this script is anchored to
#     something outside itself. The u = 2 edge is §10.7's closed form, and its
#     argmax is x*; run A's stored grid maximum sits on the same edge.
println("\n(f) the u = 2 edge against §10.7 and against run A's stored value")
@printf("    C_of(x*, β=1)      = %.13f   §10.7 C(2,x*) = %.13f   Δ = %+.2e\n",
        C_of(X_STAR, 1.0), C_at_u2(X_STAR), C_of(X_STAR, 1.0) - C_at_u2(X_STAR))
@printf("    C_of(0.12875, β=1) = %.13f   run A stored  = %.13f   Δ = %+.2e\n",
        C_of(0.12875, 1.0), 0.2360436749507588, C_of(0.12875, 1.0) - 0.2360436749507588)
@printf("    C_at_u2(x_c)       = %.3e   (the edge must vanish at the threshold)\n",
        C_at_u2(X_C))

# (g) The proof's own geometry, drawn on the figure below. u*(x) comes from V2
#     and V3 (§10.3); x_dagger comes from V (§10.4). Nothing ties those three
#     polynomials together in the .md, so u*(x_dagger) = 2 is a genuine
#     cross-check on all three transcriptions at once.
println("\n(g) region structure: u*(x) must cross u = 2 exactly at x_dagger")
@printf("    x_dagger (root of V, §10.4) = %.15f\n", X_DAGGER)
@printf("    .md quotes                  = %.15f\n", 0.051842973189252369)
@printf("    u*(x_dagger)                = %.12f   (want exactly 2)\n", u_star(X_DAGGER))
x_from_ustar = bisect(x -> u_star(x) - 2.0, 1e-4, 0.4)
@printf("    root of u*(x) - 2           = %.15f   Δ vs x_dagger = %+.2e\n",
        x_from_ustar, x_from_ustar - X_DAGGER)
@printf("    u* at x = 0.10 (> x_dagger) = %.4f  -> below 2, so region 3 is empty there\n",
        u_star(0.10))

# ==============================================================================
# 4. THE FIGURE
#
# Departs from plot_map in exactly one way, and it is deliberate: plot_map puts
# the ROW coordinate on the horizontal log axis, and here the log coordinate (u)
# belongs on the VERTICAL axis with x horizontal. So flatten_grid is reused for
# what it is good at -- pairing every cell with its coordinates and dropping
# non-finite ones -- and its two outputs are handed to scatter swapped.
#
# Dead cells are drawn, in grey, rather than dropped. An absence of markers
# would be ambiguous between "computed, came out zero" and "never sampled", and
# the whole point of the figure is that the zero region was measured.
# ==============================================================================
apply_theme!()

uflat, xflat, cflat = flatten_grid(C, u_vals, x_vals)
alive = cflat .> 0

p = scatter(xflat[.!alive], uflat[.!alive];
            marker = (:square, 2.4, :grey85), markerstrokewidth = 0,
            label  = "C = 0 (measured)",
            yscale = :log10,
            xlims  = (X_LO, X_HI),
            ylims  = (2.0, U_TOP),
            yticks = ([2, 3, 5, 10, 20, 50, 100, 200, 500],
                      ["2", "3", "5", "10", "20", "50", "100", "200", "500"]),
            xlabel = L"x \;=\; 4\eta^2/(\kappa_1\kappa_2)",
            ylabel = L"u \;=\; \beta + 1/\beta",
            # Plain string, not a latexstring: a "§" inside math mode renders
            # with a stray overline over the rest of the line under GR.
            title  = "Adiabatic entanglement region: run_adia vs the closed form\n" *
                     "colour = C_adia,  red curve = u_c(x)",
            size   = (900, 640),
            legend = :topright,
            colorbar = true)

scatter!(p, xflat[alive], uflat[alive];
         marker = (:square, 2.4), markerstrokewidth = 0,
         marker_z = cflat[alive], c = :viridis,
         clims = (0.0, maximum(cflat)),
         colorbar_title = L"C_{\mathrm{adia}}",
         label = "C > 0 (measured)")

# The analytic curve, sampled far finer than the grid so it reads as a curve
# rather than as a second scatter. Only the physical part, u_c ≥ 2.
x_curve = range(X_LO, min(X_HI, X_C); length = 2000)
u_curve = uc_closed.(x_curve)
keep    = findall(u -> isfinite(u) && u >= 2.0, u_curve)
plot!(p, x_curve[keep], u_curve[keep];
      lw = 2.5, lc = :red, label = "u_c(x) closed form (§4.1)")

# The bisected boundary, thinned so the markers stay readable on top of it.
step = max(1, cld(length(live), 45))
plot_idx = live[1:step:end]
scatter!(p, x_vals[plot_idx], u_num[plot_idx];
         marker = (:circle, 4.5, :white), markerstrokewidth = 1.4,
         markerstrokecolor = :red,
         label = "bisected from run_adia")

vline!(p, [X_STAR]; ls = :dash, lc = :black, lw = 1.6,
       label = "x* = $(round(X_STAR, digits = 5))  (argmax on the u = 2 edge)")
vline!(p, [X_C]; ls = :dashdot, lc = :darkorange, lw = 2.0,
       label = "x_c = $(round(X_C, digits = 5))  (threshold)")

# x* is a property of the u = 2 edge, not of the whole column, so mark the point
# it actually names: the global maximum of C_adia over this entire plane.
scatter!(p, [X_STAR], [2.0];
         marker = (:star5, 11, :red), markerstrokewidth = 1.0,
         markerstrokecolor = :black,
         label = "max C_adia = $(round(C_at_u2(X_STAR), digits = 6))")

save_fig(p, "adia_boundary"; dir = "Adia_boundary")

println("\n" * "="^78)
if !isempty(live)
    @printf("VERDICT: closed form and run_adia agree on the boundary to %.1e (worst, relative)\n",
            maximum(abs(u_num[i] - u_ana[i]) / u_ana[i] for i in live))
end
println("="^78)
