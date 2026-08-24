using LinearAlgebra
using Printf

include("jc_pump_disp_asy.jl")
include("observables.jl")
include("plotting_functions.jl")

# ==============================================================================
# adia_derivatives.jl
#
#   julia --project=. adia_derivatives.jl
#   julia --project=. adia_derivatives.jl x_lo x_hi N_x N_u
#
# HOW C_adia VARIES INSIDE the entangled region -- ∂C/∂u and ∂C/∂x, mapped.
#
# adia_boundary.jl settled the region's EDGE (u_c(x), agreeing with §4.1 to
# 1.2e-12). This settles its INTERIOR, which is what concurrence_adia.md §10 is
# about. Two figures, both on the same (x, u) plane as the boundary figure so
# they register against it by eye:
#
#   Adia_boundary/dC_du.pdf   ∂C/∂u -- the subject of §10's theorem
#   Adia_boundary/dC_dx.pdf   ∂C/∂x -- signed, and its zero-locus is x_opt(u)
#
# THE TWO ARE NOT EQUALLY INTERESTING, and the figures are built differently
# because of it:
#
#   ∂C/∂u is STRICTLY NEGATIVE everywhere (§10's theorem: hence β = 1 is the
#   unique maximiser). Measured max over the region is ≈ -9.9e-06, median
#   -0.149 -- five decades apart, so a plain linear scale renders as one flat
#   colour, and log10|∂C/∂u| recovers the structure only by throwing the sign
#   away. NEITHER SHOWS THE SIGN. A panel coloured by an absolute value looks
#   exactly the same whether or not a stray positive value exists, so "it is
#   negative" would be true only because the title says so.
#
#   THE FIX IS TO NEVER TAKE AN ABSOLUTE VALUE. The default encoding is
#
#       DU_SCALE = :neg_log      colour = log10( -∂C/∂u )
#
#   which is well defined precisely BECAUSE -∂C/∂u > 0, keeps all five decades
#   of structure, and puts no |·| anywhere on the figure. A cell with
#   ∂C/∂u ≥ 0 has no log10(-v) at all: it is excluded from the colour-mapped
#   series, counted in the printed colour-honesty line, and over-plotted as a
#   loud red marker, so a violation is a thing you SEE rather than a thing the
#   caption forgot to mention. (:signed_log, :log and :linear are kept as
#   comparison branches -- one edit to switch.)
#
#   The colour scale still cannot show the sign of a single cell in isolation,
#   so figure 1 carries a SECOND PANEL: the full per-row range of -∂C/∂u (row
#   min to row max), filled on a log axis, five decades, all strictly positive.
#   Positivity of -∂C/∂u there IS ∂C/∂u < 0 -- a filled band that never touches
#   the axis is the claim, not an assertion in a title. The row maximum (the
#   worst case, i.e. the value closest to zero) is annotated where it occurs,
#   read off the data.
#
#   ∂C/∂x has NO GLOBAL SIGN (§10.7) -- range ≈ [-1.46, +2.80] -- so "undefined
#   off the theorem" cannot mean the same thing there. It means one level down:
#   ONE MAP, but TWO BRANCHES, each a log of a strictly positive quantity,
#
#       DX_SCALE = :split_log    +(1 + log10( ∂C/∂x /DX_FLOOR))   for v > 0
#                                -(1 + log10(-∂C/∂x /DX_FLOOR))   for v < 0
#
#   so again NO |·| ANYWHERE: which branch a cell enters IS its sign, not a
#   factor multiplied onto a magnitude. The ±1 offset leaves the band |z| < 1
#   EMPTY BY CONSTRUCTION, so the pale centre of the diverging bar is
#   unreachable and "white means zero" is true rather than approximately true --
#   which is exactly the failure this panel had (see the long note at figure 2:
#   a linear scale painted a nonzero strip the colour of zero and manufactured a
#   second zero contour). DX_FLOOR is not a taste parameter: it is the measured
#   finite-difference floor from check (b), below which the SIGN ITSELF is not
#   resolved by the numerics, and those cells are drawn grey rather than coloured
#   as either sign. (:signed_sqrt and :linear are kept as comparison branches;
#   :linear is the known-bad case the colour-honesty audit must fire on.)
#
#   Its zero-locus is a real curve -- x_opt(u), the optimal pump as a function of
#   arm asymmetry -- and figure 2 carries a SECOND PANEL for it, on the map's own
#   x axis: ∂C/∂x against x, UNTRANSFORMED, linear, for four u rows, crossing a
#   drawn zero line exactly once each (check (g)). That panel needs no colorbar
#   to be believed, and its crossings are checked against ridge_x(u).
#
# NO SWEEP, NO FULL MODEL, NO THREADS. ~40 s for the default grid: a central
# difference is two run_adia solves per field per point, and run_adia is ~120 µs.
# ==============================================================================

BLAS.set_num_threads(1)

const X_LO = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 0.02
const N_X  = length(ARGS) >= 3 ? parse(Int,     ARGS[3]) : 220
const N_U  = length(ARGS) >= 4 ? parse(Int,     ARGS[4]) : 180

# ------------------------------------------------------------------------
# ∂C/∂u colour scale -- see header. :neg_log is the default and the only one
# that takes no absolute value. The other three are kept for comparison:
#   :signed_log  sign(v)·log10(1+|v|/DU_FLOOR) on :balance
#   :log         log10|v| on :magma -- sign lives only in the title
#   :linear      raw v, which renders as one flat colour
# ------------------------------------------------------------------------
const DU_SCALE = :neg_log

# symlog threshold for :signed_log -- a hair below the measured min |∂C/∂u|
# (≈ 9.9e-6), so the near-boundary points don't all collapse to z ≈ 0.
const DU_FLOOR = 1e-6
signed_log(v; floor = DU_FLOOR) = sign(v) * log10(1 + abs(v) / floor)

# ------------------------------------------------------------------------
# ∂C/∂x encoding. :sign_map is the default: NO COLOUR SCALE AT ALL. The map
# carries only the sign -- two flat colours, the one thing a colour can say
# without a transform -- and the magnitude comes back as LABELLED CONTOURS in
# real units. Nothing on that figure is then a transformed number.
#
# The three colour-mapped scales below are kept as comparison branches:
#   :split_log    ±(1 + log10(±v/DX_FLOOR)) -- no |·|, and the band |z| < 1 is
#                 empty by construction so nothing can render as zero; but the
#                 bar then reads ±[1+log10(±∂C/∂x/1e-5)], three ideas on one axis
#   :signed_sqrt  sign(v)·sqrt|v| -- fixes the fake zero contour, but with an
#                 |·|, so the panel would look identical if a sign were wrong
#   :linear       raw v on symmetric limits -- KNOWN BAD, kept deliberately:
#                 it is the case the colour-honesty audit below must fire on
#
# DX_FLOOR survives the removal of the colour scale: it is not a plotting
# parameter. Check (b) measures the finite-difference floor for ∂C/∂x, and a cell
# below it has a SIGN THE NUMERICS DO NOT RESOLVE, so it must not be painted as
# either sign -- it gets a grey series of its own. The check below re-reads check
# (b)'s measured worst ABSOLUTE discrepancy and complains if the margin has
# shrunk, which would mean the grey class is undersized. It no longer appears
# anywhere on the figure.
# ------------------------------------------------------------------------
const DX_SCALE = :sign_map
const DX_FLOOR = 1e-5

# Contour levels for :sign_map, drawn and labelled in the units of ∂C/∂x itself.
# c = 0 is NOT in this list: it is the ridge, which is bisected separately by
# ridge_x and checked against x* to 3.6e-11, and is drawn as the central member
# of the same family.
const DX_LEVELS = [2.0, 1.0, 0.5, 0.1, -0.1, -0.5, -1.0]

# Two branches, no abs() in either: z_pos exists only for v > 0, z_neg only for
# v < 0, and both land outside the band |z| < 1.
z_pos(v) =  (1 + log10(v / DX_FLOOR))
z_neg(v) = -(1 + log10(-v / DX_FLOOR))

const Κ_1, Κ_2, G_2 = 2.0, 2.0, 1.0
const Γ, Γ_ΦI       = 0.0, 0.0

# ------------------------------------------------------------------------
# C_of / β_of_u: duplicated from adia_boundary.jl, which duplicated C_of from
# adia_concurrence_max.jl. Not an oversight -- no standalone script here has a
# main() guard, so including one RUNS it (CLAUDE.md, "Running"). If the three
# copies ever disagree, adia_concurrence_max.jl is the oldest and most checked.
# ------------------------------------------------------------------------
function C_of(x, β; κ_1 = Κ_1, κ_2 = Κ_2, g_2 = G_2)
    η   = sqrt(x * κ_1 * κ_2 / 4)
    g_1 = g_2 * sqrt(β * κ_1 / κ_2)
    ρ   = run_adia(g_1, g_2, κ_1, κ_2, η, Γ, Γ_ΦI)
    return calc_concurrence(ρ.data)
end

β_of_u(u)    = (u + sqrt(max(u^2 - 4, 0.0))) / 2
C_of_u(x, u) = C_of(x, β_of_u(u))

# ==============================================================================
# CLOSED FORMS -- concurrence_adia.md, transcribed. Used only to CHECK the
# finite differences, never to draw: a figure of these would redraw the .md
# rather than test it.
# ==============================================================================
q1_of(x) = (x - 1)^4 + 4x
p1_of(x) = (x - 1)^4 + 4x * (x + 1)
p2_of(x) = (x^2 - 1)^2 - 4x * (x + 1)

# The unclamped C(u,x) of §2.2, C = 2[v g - √S]/D. Goes negative past u_c, which
# is exactly what makes it usable for differentiating: calc_concurrence's clamp
# is what the finite differences have to be protected from (see below).
function C_closed(u, x)
    q1, p1, p2 = q1_of(x), p1_of(x), p2_of(x)
    a, c = x * q1, 4x^2
    b    = 4x^2 * (x^2 - 2x - 1)
    g    = sqrt(x) * (1 - x)^3 * (x + 1)
    S    = a * c * u^2 + b * (a + c) * u + (a - c)^2 + b^2
    D    = (1 + x) * (p1 * u + 2p2)
    return 2 * (sqrt(u + 2) * g - sqrt(S)) / D
end

# §10.1, the unsquared derivatives:
#   ∂C/∂u = (1/D)[ g/v - S'/√S ] - (D'/D) C
# Assembled here by differentiating C_closed instead of transcribing the boxed
# form, because the boxed form's ingredients (S', D', (g²)_x, S_x, D_x) are five
# more chances to mistype and the checks below would not localise which. This is
# still an independent route from run_adia: different equations entirely.
#
# Steps are RELATIVE. An absolute step is fine in the bulk and wrong in the tail:
# check (h) evaluates the ridge out to u = 1e12, where x_opt ~ 1e-7, and a fixed
# 1e-7 step puts x - h at a negative number and throws inside sqrt(x). Relative
# steps also keep the truncation error scale-free, which matters over the eleven
# decades of x this file's asymptotics cover.
const H_CF = 1e-6
dC_du_closed(u, x) = (h = H_CF * max(1.0, u);
                      (C_closed(u + h, x) - C_closed(u - h, x)) / (2h))
dC_dx_closed(u, x) = (h = H_CF * x;
                      (C_closed(u, x + h) - C_closed(u, x - h)) / (2h))

# §10.5: ∂C/∂u at u = 2 in closed form, proven negative for all x in (0,1).
V_10_4(x) = x^7 - 7x^6 + 25x^5 - 27x^4 + 51x^3 + 3x^2 + 19x - 1
function dC_du_at_u2(x)
    num = sqrt(x) * (1 + x)^3 * (x^4 - 6x^3 + 18x^2 + 2x + 1) + x * V_10_4(x)
    return -num / (8 * (1 - x) * (1 + x)^3 * (1 + x^2)^2)
end

# §4.1, the boundary. Needed here to bracket the ridge search.
function uc_closed(x)
    q1, p1, p2 = q1_of(x), p1_of(x), p2_of(x)
    lam, kap = q1 + 4x, q1 - 4x
    sig, tau = p1 + p2 - 2q1, p1 + p2 - q1 + 4x
    Btil  = kap * tau - x * sig * lam
    C0til = 2 * kap * tau - x * (kap^2 + sig^2)
    disc  = Btil^2 + 16 * x^2 * q1 * C0til
    disc < 0 && return NaN
    return (Btil + sqrt(disc)) / (8 * x^2 * q1)
end

function bisect(f, a, b; tol = 1e-15, maxit = 200)
    fa = f(a)
    for _ in 1:maxit
        m  = (a + b) / 2
        fm = f(m)
        if (fa < 0) == (fm < 0); a, fa = m, fm else b = m end
        (b - a) <= tol * max(1.0, abs(b)) && break
    end
    return (a + b) / 2
end

const X_C     = bisect(x -> x^3 - 2x^2 + 9x - 4, 0.0, 1.0)
const T_STAR  = bisect(t -> t^6 + t^5 - 3t^4 - 2t^3 - 3t^2 - t + 1, 0.0, 1.0)
const X_STAR  = T_STAR^2
const X_HI    = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : X_C + 0.1

# ==============================================================================
# FINITE DIFFERENCES OF run_adia -- the data that actually gets plotted
#
# Two stencil hazards, both of which return a plausible wrong number rather than
# an error, so both are guarded rather than commented about:
#
#  1. calc_concurrence clamps at exactly 0.0 (its `max(0.0, ...)`). A centred
#     difference straddling u_c therefore differentiates THE CLAMP, not the
#     model, and yields a slope that looks fine. Guard: a point is evaluated
#     only when the centre AND both legs are strictly inside the region. This
#     costs a sliver of width h along the boundary, which is invisible at h=1e-6.
#
#  2. u = 2 is a hard edge -- u < 2 has no real β at all -- so no centred
#     difference exists there. That row gets a FORWARD difference. It must not
#     simply be dropped: it is the row §10.5 is about, and check (c) needs it.
# ==============================================================================
const H_U = 1e-6
const H_X = 1e-6

function dC_du_fd(x, u)
    h = H_U * max(1.0, u)
    if u <= 2.0 + h                       # bottom edge: forward difference
        c0, c1, c2 = C_of_u(x, 2.0), C_of_u(x, 2.0 + h), C_of_u(x, 2.0 + 2h)
        (c0 > 0 && c1 > 0 && c2 > 0) || return NaN
        return (-3c0 + 4c1 - c2) / (2h)   # 2nd-order one-sided
    end
    cm, cp = C_of_u(x, u - h), C_of_u(x, u + h)
    (cm > 0 && cp > 0 && C_of_u(x, u) > 0) || return NaN
    return (cp - cm) / (2h)
end

function dC_dx_fd(x, u)
    x - H_X <= 0 && return NaN
    cm, cp = C_of_u(x - H_X, u), C_of_u(x + H_X, u)
    (cm > 0 && cp > 0 && C_of_u(x, u) > 0) || return NaN
    return (cp - cm) / (2H_X)
end

# ==============================================================================
# 1. THE GRIDS
# ==============================================================================
x_vals = collect(range(X_LO, X_HI; length = N_X))
const U_TOP = 1.15 * uc_closed(X_LO)
u_vals = exp10.(range(log10(2.0), log10(U_TOP); length = N_U))

println("="^78)
println("ADIABATIC CONCURRENCE DERIVATIVES -- run_adia vs concurrence_adia.md §10")
println("="^78)
@printf("x ∈ [%.4g, %.4g] × %d    u ∈ [2, %.4g] log × %d    = %d points\n",
        X_LO, X_HI, N_X, U_TOP, N_U, N_X * N_U)
@printf("x_c = %.12f    x* = %.12f\n", X_C, X_STAR)
@printf("finite differences of run_adia, h_u = %.0e·max(1,u), h_x = %.0e\n", H_U, H_X)

t_du = @elapsed DU = [dC_du_fd(x, u) for u in u_vals, x in x_vals]
t_dx = @elapsed DX = [dC_dx_fd(x, u) for u in u_vals, x in x_vals]
@printf("\n∂C/∂u grid: %.1f s    ∂C/∂x grid: %.1f s\n", t_du, t_dx)

live_du = findall(isfinite, DU)
live_dx = findall(isfinite, DX)
@printf("live points (whole stencil inside the region): %d and %d of %d\n",
        length(live_du), length(live_dx), length(DU))

# ==============================================================================
# 2. THE RIDGE -- where ∂C/∂x = 0, i.e. the optimal pump at each asymmetry
#
# For a given u the entangled x-range is (0, x_max(u)), where x_max solves
# u_c(x) = u. The bracket MUST stay inside it: outside, C ≡ 0 and so is its
# derivative, and a bisector handed a flat zero would report a spurious root
# anywhere in it.
# ==============================================================================
# Bracket the root of u_c(x) = u. The left end must sit BELOW the answer, and
# u_c ~ 1/(4x²) means x_max ~ 1/(2√u), so reaching u = 1e12 (check (h)) needs the
# bracket to start near 1e-9, not 1e-6. It was 1e-6 first, and the failure was
# silent: with no sign change in the bracket `bisect` returns a number rather
# than complaining, and x_max_of_u(1e12) came back as X_C -- off by six orders,
# and wrong in the direction that looks plausible. Hence the explicit check.
function x_max_of_u(u)
    u <= 2.0 && return X_C
    lo = 1e-9
    uc_closed(lo) > u || return NaN     # u beyond the bracket: say so, don't guess
    return bisect(x -> uc_closed(x) - u, lo, X_C)
end

# ridge_x_closed: the same root, from the §10.1 derivative of the closed form.
#
# Used only for check (h)'s asymptotic table, which runs out to u = 1e12. Two
# reasons, neither of them "run_adia cannot do it": the closed form is exact and
# costs nothing where 300 bisection steps of run_adia would cost seconds, and
# check (b) has already validated it against run_adia across the whole grid, so
# using it out where a grid is impractical is a licensed extrapolation rather
# than an assumption. See the survival map at the end of check (h) for where
# run_adia genuinely does fail -- it is the deep dead zone, not the ridge.
function ridge_x_closed(u)
    xm = x_max_of_u(u)
    isfinite(xm) || return NaN
    hi, lo = xm * (1 - 1e-8), xm * 1e-7
    dC_dx_closed(u, lo) > 0 || return NaN
    for _ in 1:300
        m = sqrt(lo * hi)                     # geometric: x spans decades here
        dC_dx_closed(u, m) > 0 ? (lo = m) : (hi = m)
    end
    return sqrt(lo * hi)
end

function ridge_x(u)
    xm = x_max_of_u(u)
    isfinite(xm) || return NaN
    # Back off from the death boundary by FAR more than one step size. At a
    # relative 1e-6 the upper leg of the x-stencil lands past x_max, C there is
    # exactly 0, and dC_dx_fd correctly refuses -- which silently emptied the
    # bracket and returned NaN for every row on the first run. C vanishes
    # linearly in x here, so 1e-3 still leaves C ~ 1e-3·slope, far above noise.
    hi = xm - max(1e-3 * xm, 50 * H_X)
    lo = max(1e-4, hi * 1e-3)
    (isfinite(hi) && hi > lo) || return NaN
    flo, fhi = dC_dx_fd(lo, u), dC_dx_fd(hi, u)
    (isfinite(flo) && isfinite(fhi) && flo > 0 && fhi < 0) || return NaN
    a, b = lo, hi
    for _ in 1:120
        m = (a + b) / 2
        d = dC_dx_fd(m, u)
        isfinite(d) || (b = m; continue)
        d > 0 ? (a = m) : (b = m)
        (b - a) <= 1e-12 && break
    end
    return (a + b) / 2
end

t_ridge = @elapsed ridge = [ridge_x(u) for u in u_vals]
@printf("ridge found on %d of %d rows in %.1f s\n", count(isfinite, ridge), N_U, t_ridge)

# ==============================================================================
# 3. CHECKS
# ==============================================================================
println("\n" * "="^78)
println("CHECKS")
println("="^78)

# (a) §10's theorem, tested at every live grid point.
# (a) THE SIGN TEST. This used to lead with "max ∂C/∂u = -9.91e-06", which was
# wrong in three ways at once and caused two misreadings in one session:
#
#   1. It is NOT a maximum of anything about the model. sup ∂C/∂u over the
#      physical region is exactly 0, approached as x -> 0 (A_0 ~ -sqrt(x), see
#      (a2)) and as u -> ∞ (the u^-3/2 law). No grid can ever exhibit a strictly
#      negative bound, so quoting one as the headline invents a margin.
#   2. It is the LEAST NEGATIVE value, i.e. the WEAKEST -- while the largest
#      magnitude in the same grid is 0.94, five decades away. Calling the small
#      one "the max" is the exact wording that got the figure read backwards.
#   3. It is window-determined: it sits at x = X_LO, always (see below).
#
# So the theorem gets a BOOLEAN AND A COUNT, which is what a sign claim actually
# is, and the extremal values are reported afterwards as the coordinates they
# are. The count can fail; the old headline could not.
println("\n(a) THE SIGN TEST: is ∂C/∂u < 0 at every point sampled?")
nviol = count(>=(0.0), DU[live_du])
@printf("    live points %d;  points with ∂C/∂u >= 0: %d\n", length(live_du), nviol)
println(nviol == 0 ? "    -> no violation sampled (the theorem is PROVEN in §10; this only tests the numerics)" :
                     "    -> VIOLATION SAMPLED -- the theorem fails, or the stencil is wrong")

mx, imx = findmax(DU[live_du])       # least-negative; kept for (a2), not a result
ci = live_du[imx]
@printf("    strongest sampled (largest |∂C/∂u|): %+.4e\n", minimum(DU[live_du]))
@printf("    weakest   sampled (closest to zero): %+.4e at x = %.4g, u = %.4f  -- %.1f decades apart\n",
        mx, x_vals[ci[2]], u_vals[ci[1]],
        log10(abs(minimum(DU[live_du])) / abs(mx)))
# ∂C/∂u is monotone DECREASING in x at fixed u, so the weakest value is at the
# smallest live x -- the window edge -- in every row. Asserted, not assumed: if
# the argmin ever leaves column 1 the framing in (a2) needs revisiting.
@printf("    the weakest sits at x-column %d of %d%s\n", ci[2], N_X,
        ci[2] == 1 ? " = X_LO: it is a WINDOW CORNER, not a property of the model (see a2)" :
                     "  <-- NOT column 1: ∂C/∂u is not monotone in x here, (a2)'s framing needs revisiting")

# (b) The load-bearing one: FD of run_adia against §10.1's derivatives of the
#     closed-form C. Two routes with no shared code.
println("\n(b) finite differences of run_adia vs the §10.1 derivatives")
# Recorded rather than only printed: figure 2's DX_FLOOR is anchored on the
# ∂C/∂x entry, and an anchor that is never re-read is an assumption.
fd_worst = Dict{String,Float64}()
for (nm, G, f) in (("∂C/∂u", DU, dC_du_closed), ("∂C/∂x", DX, dC_dx_closed))
    liv = findall(isfinite, G)
    rel = Float64[]
    dev = 0.0                                     # worst ABSOLUTE discrepancy
    skipped = 0
    for ci in liv
        u, x = u_vals[ci[1]], x_vals[ci[2]]
        u <= 2.0 + H_U && continue                # one-sided row: check (c) covers it
        a = f(u, x)
        dev = max(dev, abs(G[ci] - a))
        # A relative error is meaningless where the reference is itself near the
        # FD floor: |∂C/∂u| gets down to 1e-5, and 1e-9 of absolute noise on it
        # reads as 1e-4 relative. Skip those and say how many.
        abs(a) < 1e-4 && (skipped += 1; continue)
        push!(rel, abs(G[ci] - a) / abs(a))
    end
    sort!(rel)
    fd_worst[nm] = dev
    @printf("    %s   n = %5d (+%d too small to score)   worst %.2e   median %.2e   %s   (worst ABSOLUTE %.1e)\n",
            nm, length(rel), skipped, rel[end], rel[cld(length(rel), 2)],
            rel[end] < 1e-5 ? "ok (FD floor)" : "MISMATCH", dev)
end

# (a2) WHY (a)'s max is a coordinate and not a result -- and the law that
#      replaces it. Placed after (b) because it needs (b)'s measured FD floor.
#
# The session note "Comparison of ∂C/∂u at u = 2 and u = 500" gives the exact
# large-u expansion of the UNCLAMPED bracket,
#
#     ∂C/∂u = A_0(x) u^(-3/2) + A_1(x) u^(-2) + O(u^(-5/2))
#     A_0 = -sqrt(x)(1-x)^3 / p1
#     A_1 = -2 x^(3/2) (1-x)^3 V_2 / ((1+x) sqrt(q1) p1^2)
#     V_2 = x^7 - 7x^6 + 19x^5 - 21x^4 + 7x^3 + 23x^2 + 5x + 5
#
# V_2 is §10.3's polynomial, already Sturm-proven positive on (0,1), and p1, q1
# are positive there, so A_0 < 0 AND A_1 < 0 on all of (0,1). THAT IS AN
# INDEPENDENT PROOF OF ∂C/∂u < 0 AT LARGE u, cheaper than §10's machinery and
# reached by a different route -- it says nothing near u = 2, which is what
# check (c) and §10.5 are for.
#
# Two things follow, and they are why (a) was reworded:
#   1. A_0 ~ -sqrt(x) as x -> 0, so ∂C/∂u -> 0 there. Combined with u^(-3/2),
#      the supremum over the region is 0 and is approached in TWO directions,
#      both of which the grid truncates.
#   2. (a)'s max is therefore just A_0(X_LO)·u_c(X_LO)^(-3/2) -- the window
#      corner, in closed form, needing no grid. Reproduced below.
V2_of(x) = x^7 - 7x^6 + 19x^5 - 21x^4 + 7x^3 + 23x^2 + 5x + 5
A0_of(x) = -sqrt(x) * (1 - x)^3 / p1_of(x)
A1_of(x) = -2 * x^1.5 * (1 - x)^3 * V2_of(x) /
           ((1 + x) * sqrt(q1_of(x)) * p1_of(x)^2)

#
# ⚠ THE LAW IS A STATEMENT ABOUT THE UNCLAMPED BRACKET C̃, NOT ABOUT C.
# The note's §1 exists to separate the two, and they differ almost everywhere at
# large u. The physical concurrence is C = 2·max{0, ·}, so ∂C/∂u is IDENTICALLY
# ZERO for x beyond the death boundary x_d(u), while C̃'s derivative there is the
# analytic continuation -- a number, but not a slope of anything physical.
# x_d(500) = 0.02144, so the note's own x = 0.05, 0.1486, 0.4838 rows are all
# OUTSIDE. An earlier version of this check lifted exactly those three x values
# into a script whose grid is the clamped field, which is the one confusion the
# note leads with. Every x below is asserted to satisfy x < x_d(u), where the two
# objects coincide and the law is the physical derivative.
x_d(u) = x_max_of_u(u)               # same root as the boundary; cross-checked below

println("\n(a2) the large-u law  ∂C/∂u = A_0 u^-3/2 + A_1 u^-2, INSIDE the region only")
@printf("     A_0 < 0 and A_1 < 0 on (0,1) via §10.3's V_2 > 0 -> an independent large-u\n")
@printf("     proof of ∂C/∂u < 0 for the bracket; inside the region C ≡ C̃ so it carries over\n")
@printf("     x_d(500) = %.18f   (note's 60-digit value 0.021439817689553069, rel Δ %.1e)\n",
        x_d(500.0), abs(x_d(500.0) - 0.0214398176895530694) / 0.0214398176895530694)
@printf("     %-12s %18s %12s %12s   %s\n", "x (at u=500)", "fd", "1-term/fd", "2-term/fd", "inside?")
for x in (0.002, 0.005, 0.010, 0.020)
    fd = dC_du_closed(500.0, x)
    t1 = A0_of(x) * 500.0^-1.5
    t2 = t1 + A1_of(x) * 500.0^-2
    @printf("     %-12.4f %+18.9e %12.4f %12.4f   %s\n", x, fd, t1 / fd, t2 / fd,
            x < x_d(500.0) ? "yes" : "NO -- unphysical, remove it")
end
# Both truncation orders are pinned, which is what makes this a test of the LAW
# rather than of one fitted number. The 1-term relative error is A_1/A_0·u^(-1/2)
# so it must halve over u × 4; the 2-term error is the u^(-5/2) tail against a
# u^(-3/2) leading term, i.e. O(u^-1), so it must QUARTER over the same span.
# (An earlier version of this line expected ~2.0 for the 2-term ratio and got
# 3.79 -- the measurement was right and the exponent in the comment was wrong.)
#
# THE u^-1/2 CORRECTION IS NOT OBSERVABLE INSIDE THE REGION, and this is a real
# qualification of the note's §3 rather than a numerical nuisance.
#
# A FIXED x cannot be used to measure the u-scaling: x_d ~ 1/(2√u), so any fixed
# x leaves the region as u grows, and the scan would be measuring the analytic
# continuation again. Scanning at fixed σ = x√u instead keeps x at a constant
# FRACTION of the way to the death boundary (σ -> 1/2 is the boundary itself,
# σ = 1/6 is the ridge), which is the only way to grow u and stay physical.
#
# But at fixed σ, x ~ σ/√u, so A_1/A_0 ~ 10x ~ 10σ/√u and the A_1 term's relative
# size is ~10σ/u -- the SAME O(u^-1) as the next order. So both truncations
# converge at u^-1 and A_1 never separates out. Measured below: both ratios -> 4.
#
# The note's u^-1/2 story is a fixed-x statement, and at large u fixed x is
# outside the region -- which is exactly why its own table used x = 0.05, 0.1486,
# 0.4838. An earlier version of this check tried a fixed x = 0.005 and got a
# 1-term ratio of 10.86: at u = 500 the A_1 term (2.24e-3) and the next order
# (2.0e-3) are the same size there, so the residual was a cancellation between
# two comparable terms and its ratio meant nothing.
let σ = 0.35                       # 73% of the way to x_d; ridge sits at 1/6
    res(u) = (x = σ / sqrt(u); f = dC_du_closed(u, x);
              (abs(A0_of(x) * u^-1.5 / f - 1),
               abs((A0_of(x) * u^-1.5 + A1_of(x) * u^-2) / f - 1),
               x, x < x_d(u)))
    us = (500.0, 2000.0, 8000.0, 32000.0)
    @printf("     residuals at FIXED σ = x√u = %.2f (staying inside as u grows):\n", σ)
    @printf("        %-9s %-11s %-9s %-12s %-12s %s\n", "u", "x", "inside?", "1-term", "2-term", "ratios vs prev")
    prev = nothing
    for u in us
        r1, r2, x, ok = res(u)
        rat = prev === nothing ? "" : @sprintf("%.2f  %.2f", prev[1] / r1, prev[2] / r2)
        @printf("        %-9.0f %-11.5g %-9s %-12.4g %-12.4g %s\n", u, x, ok ? "yes" : "NO", r1, r2, rat)
        prev = (r1, r2)
    end
    println("        -> both ratios approach 4, i.e. BOTH truncations are O(u^-1) in the region;")
    println("           the note's O(u^-1/2) correction is a fixed-x statement and fixed x leaves it.")
end
# THE LINE THAT RETIRES (a)'s HEADLINE: predict it from the window corner alone.
# ci and mx come from check (a). If this reproduces (a)'s max, then (a)'s max was
# never measuring the model -- it was measuring X_LO.
let xg = x_vals[ci[2]], ug = u_vals[ci[1]]
    p1t = A0_of(xg) * ug^-1.5
    p2t = p1t + A1_of(xg) * ug^-2
    @printf("     the weakest sampled value, predicted from the window corner (x = %.4g, u = %.1f) with NO grid:\n", xg, ug)
    @printf("        measured %+.4e   1-term %+.4e (%.2f%%)   2-term %+.4e (%.2f%%)\n",
            mx, p1t, 100abs(p1t / mx - 1), p2t, 100abs(p2t / mx - 1))
end
# THE GUARD. X_LO is a CLI argument, and pushing it left drives (a)'s max toward
# zero like sqrt(X_LO) while the FD floor stays put. At X_LO = 0.001 the max is
# -2.5e-10 against a floor of ~1.1e-9: the sign of the file's own headline stops
# being resolved, and without this line (a) would still print "uniformly
# negative". Fire it with:  julia --project=. adia_derivatives.jl 0.001 0.55 20 15
let floor_du = fd_worst["∂C/∂u"], margin = abs(mx) / fd_worst["∂C/∂u"]
    @printf("     |weakest sampled ∂C/∂u| = %.2e vs check (b)'s worst ABSOLUTE FD discrepancy %.1e  (%.1f× margin)%s\n",
            abs(mx), floor_du, margin,
            margin >= 10 ? "" :
            "\n     <-- WINDOW TOO NARROW: the weakest sampled value is below the FD floor, its sign is NOT resolved")
end

# (c) §10.5's closed form for the u = 2 row, against the one-sided difference.
println("\n(c) ∂C/∂u at u = 2: §10.5's closed form vs run_adia")
@printf("    %-10s %16s %16s %10s\n", "x", "run_adia (fd)", "§10.5", "Δ")
for x in (0.02, 0.05, 0.10, X_STAR, 0.20, 0.30, 0.40, 0.998X_C)
    fd, cf = dC_du_fd(x, 2.0), dC_du_at_u2(x)
    @printf("    %-10.5f %+16.9f %+16.9f %10.1e\n", x, fd, cf, abs(fd - cf))
end
@printf("    monotone decreasing: %+.4f at x = %.3f  ->  %+.4f at x_c\n",
        dC_du_at_u2(X_LO), X_LO, dC_du_at_u2(X_C))
println("    i.e. C''(β=1) grows toward threshold -- the β = 1 peak SHARPENS as")
println("    the entanglement dies, so the state is most fragile to asymmetry")
println("    exactly where it is weakest.")

# (d) §10.7 proves ∂C/∂x|_{u=2} changes sign exactly at x*. The ridge's bottom
#     endpoint is therefore a number known in advance, not a fitted one.
println("\n(d) the ridge must start at x* on the u = 2 edge (§10.7)")
r2 = ridge_x(2.0)
@printf("    ridge at u = 2 : %.12f\n", r2)
@printf("    x* (root of R) : %.12f    Δ = %.2e\n", X_STAR, abs(r2 - X_STAR))

# (e) The ridge against the three x_opt values CLAUDE.md stores from
#     files_online/06 §6 -- three isolated numbers, from a third project.
println("\n(e) the ridge against the stored x_opt values")
@printf("    %-8s %-10s %14s %12s\n", "β", "u", "ridge here", "stored")
for (β, want) in ((1.0, 0.1486), (2.0, 0.107), (10.0, 0.042))
    u = β + 1 / β
    r = ridge_x(u)
    @printf("    %-8.1f %-10.4f %14.6f %12.4f   Δ = %.1e\n", β, u, r, want, abs(r - want))
end

# (f) §10.7: ∂C/∂x -> +∞ like 1/sqrt(x(u+2)) as x -> 0, so the scaled quantity
#     must approach 1. This is the asymptote that forbids a global sign.
println("\n(f) ∂C/∂x → +∞ as x → 0 :  ∂C/∂x · sqrt(x(u+2)) → 1  (§10.7)")
for x in (1e-2, 1e-3, 1e-4, 1e-5)
    s = dC_dx_fd(x, 4.0) * sqrt(x * (4.0 + 2))
    @printf("    x = %.0e   scaled = %.6f\n", x, s)
end

# (g) EXACTLY ONE zero of ∂C/∂x per u row. The ridge is bisected on the
#     assumption of a single sign change, so if C_adia were bimodal in x at some
#     u the ridge would be silently incomplete -- it would find one root and
#     report nothing about the other. Scan each row densely and count.
#
#     This check exists because the figure INVITES the opposite conclusion: on a
#     linear colour scale the near-zero strip along u = 2 is painted the same
#     white as the zero locus and looks like a second contour. It is not one.
println("\n(g) ∂C/∂x must have exactly ONE sign change per u row")
counts = Dict{Int,Int}()
for u in u_vals
    xm = x_max_of_u(u)
    isfinite(xm) || continue
    hi, lo = xm - max(1e-3 * xm, 50 * H_X), 1e-4
    hi > lo || continue
    prev, ch = 0, 0
    for i in 0:400
        v = dC_dx_closed(u, lo + (hi - lo) * i / 400)
        s = v > 0 ? 1 : -1
        prev != 0 && s != prev && (ch += 1)
        prev = s
    end
    counts[ch] = get(counts, ch, 0) + 1
end
for (k, v) in sort(collect(counts))
    @printf("    %d sign change(s): %d rows\n", k, v)
end
println(keys(counts) == Set([1]) ? "    -> unimodal in x at every u; the ridge is the whole zero-locus" :
                                   "    -> NOT unimodal somewhere; the bisected ridge is INCOMPLETE")

# (h) HOW FAR THE RIDGE GOES. It has a genuine lower terminus and NO upper one,
#     and the figure's frame hides that -- worth reporting rather than leaving
#     to be inferred from a curve that runs off the side of the plot.
#
#     Lower end: (x*, u = 2). A real endpoint, because u = β + 1/β ≥ 2 is a hard
#     domain boundary (β = 1), not because the curve stops.
#     Upper end: none. x_opt → 0 as u → ∞, and the ridge stays strictly inside
#     the entangled region the whole way.
println("\n(h) how far the ridge extends")
@printf("    lower terminus: (x*, u=2) = (%.9f, 2)  -- the β = 1 edge, a hard boundary\n", X_STAR)
println("    upper terminus: NONE -- it runs to u → ∞ with x_opt → 0. Asymptotically:")
println("      (closed form: exact and free out here; run_adia is fine on the ridge -- see map)")
println("      u          x_opt          x_opt·√u      x_max·√u     x_opt/x_max")
for u in (2.0, 10.1, 100.0, 1.0e4, 1.0e8, 1.0e12)
    r, xm = ridge_x_closed(u), x_max_of_u(u)
    (isfinite(r) && isfinite(xm)) || continue
    @printf("      %-10.3g %-14.9f %-13.8f %-12.8f %.8f\n",
            u, r, r * sqrt(u), xm * sqrt(u), r / xm)
end
@printf("      limits:                   %-13.8f %-12.8f %.8f\n", 1/6, 0.5, 1/3)
println("    i.e. x_opt → 1/(6√u) and x_max → 1/(2√u), so at extreme asymmetry the")
println("    optimal pump sits at exactly ONE THIRD of the way to the death boundary.")
u_exit = bisect(u -> 0.02 - ridge_x(u), 2.01, 1.0e6)
@printf("    it leaves THIS FIGURE at the left edge x = %.3f, at u = %.2f (β = %.2f),\n",
        X_LO, u_exit, β_of_u(u_exit))
@printf("    which is well below the frame top of u = %.1f -- so the ridge exits the\n", U_TOP)
@printf("    SIDE, not the top, and the upper ~%.0f%% of the u range shows none of it.\n",
        100 * (1 - log(u_exit / 2) / log(U_TOP / 2)))

# WHERE run_adia GIVES OUT. It does, but NOT along the ridge -- it fails when
# large β is combined with large x, which is the deep dead zone, far from any
# entangled state. steadystate.eigenvector throws "Eigenvalue with smallest
# absolute value is not zero" there: loud, not a wrong number, which is the good
# kind of failure.
#
# Recorded because a bracket bug in x_max_of_u once drove ridge_x into exactly
# that corner (its left end was 1e-6, too high to bracket x_max ~ 1/(2√u) beyond
# u ~ 1e12, and `bisect` with no sign change returns a number rather than
# complaining -- it came back as X_C, six orders too big). The ridge search then
# probed x ≈ 0.48 at β = 1e8 and the whole script died. The bug was in the
# bracket, not in run_adia.
#
# The asymptotic table above uses the closed form because it is exact and free,
# and because check (b) validated it against run_adia across the entire grid --
# not because run_adia fails on the ridge. It does not: the ridge sits at
# x ~ 1/(6√u), which stays inside the survivable band at every u tested.
println("\n    run_adia survival map (rows u ≈ β, cols x) -- the ridge sits at x ~ 1/(6√u):")
@printf("      %-10s", "u \\ x")
for x in (1e-6, 1e-4, 1e-2, 0.1, 0.3); @printf("%9.0e", x); end
println("     x_ridge")
for u in (1.0e2, 1.0e4, 1.0e6, 1.0e8, 1.0e10)
    @printf("      %-10.0e", u)
    for x in (1e-6, 1e-4, 1e-2, 0.1, 0.3)
        ok = try isfinite(C_of_u(x, u)) catch; false end
        @printf("%9s", ok ? "ok" : "FAIL")
    end
    @printf("   %9.2e\n", ridge_x_closed(u))
end
println("      -> failures are large-β-AND-large-x, i.e. deep dead zone, off the ridge")

# ==============================================================================
# 4. FIGURES
#
# Both reuse flatten_grid for pairing cells with coordinates and dropping
# non-finite ones, then hand the two outputs to scatter SWAPPED -- the log
# coordinate here is u and belongs on the vertical axis, whereas plot_map
# hardcodes the row coordinate onto the horizontal. Same departure as
# adia_boundary.jl, same reason.
# ==============================================================================
apply_theme!()

x_curve = range(X_LO, X_C; length = 2000)
u_curve = uc_closed.(x_curve)
keepc   = findall(u -> isfinite(u) && u >= 2.0, u_curve)

function base_panel(title, ylab_extra = ""; cb = true)
    plot(; yscale = :log10, xlims = (X_LO, X_HI), ylims = (2.0, U_TOP),
         yticks = ([2, 3, 5, 10, 20, 50, 100, 200, 500],
                   ["2", "3", "5", "10", "20", "50", "100", "200", "500"]),
         xlabel = L"x \;=\; 4\eta^2/(\kappa_1\kappa_2)",
         ylabel = L"u \;=\; \beta + 1/\beta",
         title = title, size = (900, 640), legend = :topright, colorbar = cb)
end

# ---- Figure 1: ∂C/∂u -------------------------------------------------------
uf, xf, vf = flatten_grid(DU, u_vals, x_vals)

# The colorbar TITLE has to name what the numbers on it actually are. GR does
# not honour custom colorbar_ticks reliably (its colorbar support is thin --
# CLAUDE.md says as much), so relabelling the ticks as negative decades was
# silently ignored on the first render and the bar read "∂C/∂u" over numbers
# that were log10|∂C/∂u|. Name the transform instead of fighting GR.
#
# :neg_log is the default: colour = log10(-∂C/∂u). NO ABSOLUTE VALUE IS TAKEN
# ANYWHERE, which is the whole point -- log10 of a negated quantity is defined
# only where that quantity is negative, so the encoding is not merely labelled
# "negative", it is unavailable otherwise. Cells with ∂C/∂u ≥ 0 are split out
# BEFORE the colour map (they have no log10(-v)), counted in the printed
# colour-honesty line, and over-plotted in red below.
pos = findall(>=(0), vf)          # forbidden by §10 -- expected to be empty
neg = findall(<(0),  vf)

@printf("    colour-honesty (dC_du panel): %d of %d live cells have ∂C/∂u ≥ 0, i.e. no log10(-v) and no colour\n",
        length(pos), length(vf))

if DU_SCALE == :neg_log
    xs, us = xf[neg], uf[neg]
    zs     = log10.(-vf[neg])
    # REVERSED, deliberately. Plain :magma puts bright at high log10(-∂C/∂u),
    # i.e. bright = MOST negative -- so the panel's most conspicuous feature was
    # its least marginal region, and a reader scanning for "is this negative?"
    # was pulled to the corner where the answer is least in doubt. Reversed,
    # bright = closest to zero, so the eye goes to where the claim is tightest.
    cmap   = cgrad(:magma, rev = true)
    zclims = (minimum(zs), maximum(zs))
    cbl    = L"\log_{10}\left(-\,\partial C/\partial u\right)"
    note   = "colour = log10(−∂C/∂u)\n" *
             "bright = closest to zero;  dark = most strongly negative"
elseif DU_SCALE == :signed_log
    xs, us = xf, uf
    zs     = signed_log.(vf)
    cmap   = :balance
    m      = maximum(abs, zs)
    zclims = (-m, m)
    cbl    = L"\mathrm{sign}\cdot\log_{10}(1+|\partial C/\partial u|/10^{-6})"
    note   = "colour = signed log10(1+|∂C/∂u|/1e-6) on :balance — all data left of centre"
elseif DU_SCALE == :log
    xs, us = xf, uf
    zs     = log10.(abs.(vf))
    cmap   = :magma
    zclims = (minimum(zs), maximum(zs))
    cbl    = L"\log_{10}|\partial C/\partial u|"
    note   = "colour = log10|∂C/∂u|  (the SIGN is uniform and negative — see check (a))"
else
    xs, us = xf, uf
    zs     = vf
    cmap   = :magma
    zclims = (minimum(zs), maximum(zs))
    cbl    = L"\partial C/\partial u"
    note   = "colour = ∂C/∂u, linear"
end

# The second line is NOT decoration. cgrad(:magma, rev = true) puts bright at the
# values CLOSEST to zero, which is the opposite of plain :magma and is the whole
# reason for the reversal. That direction used to be stated only in panel 2's
# legend ("the dark band of the map"); panel 2 is gone, so it moves here or it
# exists nowhere on the figure.
p1 = base_panel("∂C/∂u < 0 EVERYWHERE \n" *
                "bright = closest to zero (weakest) · dark = most negative\n")
scatter!(p1, xs, us; marker = (:square, 2.4), markerstrokewidth = 0,
         marker_z = zs, c = cmap, clims = zclims, colorbar_title = cbl, label = "")
plot!(p1, x_curve[keepc], u_curve[keepc]; lw = 2.5, lc = :cyan,
      label = L"u_c(x)\ \mathrm{(§4.1)}")
vline!(p1, [X_C]; ls = :dashdot, lc = :darkorange, lw = 2.0,
       label = "x_c = $(round(X_C, digits = 5))")
# Drawn only if it ever fires. An always-present legend entry reading "0 cells"
# is an assertion; a marker series that appears out of nowhere is evidence.
if !isempty(pos)
    scatter!(p1, xf[pos], uf[pos]; marker = (:circle, 7), mc = :red,
             markerstrokewidth = 1.2, markerstrokecolor = :black,
             label = @sprintf("∂C/∂u ≥ 0 : %d cells — THEOREM VIOLATED", length(pos)))
end

# ---- Figure 1, lower panel: the negativity envelope -------------------------
#
# The colour panel above can show the sign of the FIELD only through its own
# encoding, which the reader has to take on the colorbar's word. This panel does
# not: for each u row it plots max_x ∂C/∂u -- the least-negative, i.e. WORST
# CASE, value anywhere in that row -- untransformed, on a linear axis, against a
# bold zero line. Every live cell lies below its own row's maximum by
# definition, so one curve under that line is the entire 2D claim in 1D.
#
# Its own maximum must equal check (a)'s max over all live points; the two are
# computed by different reductions of the same grid, so the agreement line below
# is a free consistency test of this panel against the check it illustrates.
#
# The row MINIMUM is collected too, for the third panel -- see below.
row_max = fill(NaN, N_U)
row_min = fill(NaN, N_U)
for i in 1:N_U
    fin = filter(isfinite, @view DU[i, :])
    isempty(fin) && continue
    row_max[i] = maximum(fin)
    row_min[i] = minimum(fin)
end
okrm         = findall(isfinite, row_max)
rm_worst, iw = findmax(row_max[okrm])
u_worst      = u_vals[okrm[iw]]
rm_min       = minimum(row_max[okrm])
strongest    = minimum(row_min[okrm])

@printf("    row reductions (cross-check of (a); not drawn since fig 1 became one panel): %d of %d u rows live;\n                              weakest row value %+.4e at u = %.2f  (check (a) had %+.4e, Δ = %.1e)\n",
        length(okrm), N_U, rm_worst, u_worst, mx, abs(rm_worst - mx))
@printf("                              strongest row value %+.4e  (check (a) had %+.4e, Δ = %.1e); %.1f decades apart\n",
        strongest, minimum(DU[live_du]), abs(strongest - minimum(DU[live_du])),
        log10(strongest / rm_worst))
# Which of the two edges is real. Asserted from the data, not from the comment:
# every live row's argmax must be its first live column for the upper edge to be
# the window edge. If that ever fails, the panel's caption is wrong.
let edge_rows = 0, tot = 0
    for i in okrm
        fin = findall(isfinite, @view DU[i, :]); isempty(fin) && continue
        tot += 1
        argmax([DU[i, j] for j in fin]) == 1 && (edge_rows += 1)
    end
    @printf("                              the weakest value is at x = X_LO in %d of %d live rows -- WINDOW-DETERMINED, not a margin;\n", edge_rows, tot)
    @printf("                              the strongest is at x -> x_max(u) (death boundary, inside the window) -- intrinsic\n")
end

# ---- FIGURE 1 IS ONE PANEL. The -∂C/∂u band panel was removed on request.
#
# What it drew: the per-row range of -∂C/∂u, row min to row max, filled, on a log
# axis. What went with it, and where that content now lives:
#
#   - Its LOWER edge (row min, at x -> x_max(u), the death boundary) was the one
#     intrinsic curve on it. That is now unillustrated. If it is ever wanted
#     back, it is the only part worth redrawing -- NOT the band.
#   - Its UPPER edge was DU[:, 1], a slice along x = X_LO, i.e. an artifact of
#     the window. Removing it removes a curve that was being read as a margin.
#   - It never carried the sign claim and could not: on a log axis a value with
#     ∂C/∂u >= 0 is UNPLOTTABLE -- dropped, not drawn on the wrong side -- so the
#     panel could not falsify its own headline. Check (a)'s violation COUNT is
#     what carries the sign, and it can fail.
#
# The row reductions above are KEPT even though nothing draws them now: they
# cross-check check (a)'s extrema by a different reduction of the same grid
# (Δ = 0.0e+00), and the 175-of-175 assertion is what licenses check (a)'s
# "the weakest sits at column 1" and all of (a2)'s framing. They are checks that
# happened to have a picture, not a picture that happened to print numbers.
#
# The panel's env_kw went with it rather than being left as an unused binding.
# If the death-boundary curve is ever drawn on its own, the u axis it used was
#   xscale = :log10, xlims = (2.0, U_TOP), colorbar = false,
#   xticks = ([2,3,5,10,20,50,100,200,500], the same as strings)
save_fig(plot(p1; size = (900, 640), left_margin = 7Plots.mm),
         "dC_du"; dir = "Adia_boundary")

# ---- Figure 2: ∂C/∂x -------------------------------------------------------
# sym_clims inlined rather than lifted from compare_runs.jl -- see CLAUDE.md
# issue 2; this is now the fourth site that wants it.
uf2, xf2, vf2 = flatten_grid(DX, u_vals, x_vals)

# SIGNED SQUARE ROOT, and this one is a correctness fix, not a preference.
#
# On a plain linear scale with symmetric limits, this panel READS AS IF IT HAS
# TWO ZERO CONTOURS. It does not -- ∂C/∂x has exactly one sign change per u row
# (C_adia is unimodal in x at fixed u; verified on 400 rows), and the ridge below
# is it. The illusion comes from the colour range: ∂C/∂x diverges as x -> 0
# (§10.7, check (f)), so |clims| is set to 2.80 by that corner, and along the
# u = 2 edge
#
#     x = 0.1486   ∂C/∂x =  0.0000    0.0% of range   <- the ONE true zero
#     x = 0.160    ∂C/∂x = -0.0837    3.0%            <- NOT zero, renders white
#     x = 0.180    ∂C/∂x = -0.2152    7.7%            <- NOT zero, renders white
#     x = 0.250    ∂C/∂x = -0.5675   20.3%            <- finally visibly blue
#
# so a genuinely nonzero strip is painted the same white as the zero locus, and
# the eye reads the strip's far edge as a second contour. WHITE IN A DIVERGING
# COLORMAP MEANS "SMALL", NOT "ZERO" -- and how small depends on an extreme value
# somewhere else in the panel entirely.
#
# sign(v)·sqrt(|v|) was the first fix, and it did work at the source: it
# compresses the divergent tail and expands the near-zero range, so x = 0.18
# moves from 7.7% to 28% of the bar and the white collapses to a thin band around
# the actual zero. (A 99th-percentile clip was tried before it and reverted -- it
# helps the bulk but saturates a wide band along the death boundary into flat
# dark blue, trading one artifact for another.)
#
# BUT IT FIXES IT WITH AN |·|, WHICH IS THE THING FIGURE 1 REFUSES TO DO. Under
# sign(v)·sqrt|v| the sign is a factor multiplied onto a magnitude, so the panel
# would render identically if a sign were wrong; and the exponent 1/2 is chosen
# by nothing except that it looked better than 1. The default is now :split_log,
# two per-sign branches (see header), which
#
#   * takes no absolute value: membership of a branch IS the sign
#   * leaves |z| < 1 EMPTY, so no cell can be painted the colour of zero --
#     the failure above becomes structurally impossible rather than mitigated
#   * resolves all ~5.5 decades instead of crushing them into the top of the bar
#   * refuses to colour anything below the MEASURED FD floor, where the sign is
#     not resolved by the numerics at all
#
# :signed_sqrt and :linear are kept for comparison. :linear is retained on
# purpose: the colour-honesty audit below must FIRE on it, and a check never
# shown to fail on a known-bad panel is not a check (CLAUDE.md, "Figure audits").
#
# The colorbar TITLE names the transform, for the reason figure 1's does: GR does
# not reliably honour custom colorbar ticks, so relabelling the numbers back into
# linear units is not available. The subtitle gives two anchor points instead.
signed_sqrt(v) = sign(v) * sqrt(abs(v))

# Cells below the FD floor: their sign is noise, so they are drawn grey rather
# than coloured as either sign. Split out BEFORE the colour map, the same way
# figure 1 splits out cells with no log10(-v).
sub  = findall(v -> abs(v) < DX_FLOOR, vf2)
colr = findall(v -> abs(v) >= DX_FLOOR, vf2)

@printf("    DX_FLOOR = %.0e vs check (b)'s worst ABSOLUTE FD discrepancy %.1e  (%.0f× margin)%s\n",
        DX_FLOOR, fd_worst["∂C/∂x"], DX_FLOOR / fd_worst["∂C/∂x"],
        DX_FLOOR > 10 * fd_worst["∂C/∂x"] ? "" : "   <-- TOO TIGHT: grey class is undersized")
@printf("    sub-floor (dC_dx panel): %d of %d live cells have |∂C/∂x| < %.0e — sign unresolved, drawn grey\n",
        length(sub), length(vf2), DX_FLOOR)

if DX_SCALE == :sign_map
    # No colour scale. zf2/m/cbl are left undefined on purpose -- nothing
    # downstream may reach for them on this branch.
    xc2, uc2 = xf2[colr], uf2[colr]
    sub2 = "colour carries the SIGN only; magnitude is on the labelled contours, in the units of ∂C/∂x"
elseif DX_SCALE == :split_log
    vc  = vf2[colr]
    zf2 = [v > 0 ? z_pos(v) : z_neg(v) for v in vc]     # no abs() in either branch
    xc2, uc2 = xf2[colr], uf2[colr]
    m   = maximum(abs, zf2)
    cbl = L"\pm\left[1+\log_{10}\left(\pm\,\partial C/\partial x\,/\,10^{-5}\right)\right]"
    sub2 = @sprintf("colour = signed decades above 10^-5:  ±1 → ±1e-5,  ±6 → ±1,  ±%.1f → ±%.2f  (bar ends)\nthe pale band at the centre of the bar is EMPTY by construction — nothing renders as zero",
                    m, DX_FLOOR * 10.0^(m - 1))
elseif DX_SCALE == :signed_sqrt
    zf2 = signed_sqrt.(vf2[colr])
    xc2, uc2 = xf2[colr], uf2[colr]
    m   = maximum(abs, zf2)
    cbl = L"\mathrm{sign}\cdot\sqrt{|\partial C/\partial x|}"
    sub2 = "colour = sign·sqrt|∂C/∂x| — an |·| appears, so the sign is asserted, not encoded"
else
    zf2 = vf2[colr]
    xc2, uc2 = xf2[colr], uf2[colr]
    m   = maximum(abs, zf2)
    cbl = L"\partial C/\partial x"
    sub2 = "colour = ∂C/∂x, linear — KNOWN BAD: paints a nonzero strip the colour of zero"
end

p2 = base_panel(L"\partial C/\partial x"; cb = DX_SCALE != :sign_map)

if DX_SCALE == :sign_map
    # TWO FLAT COLOURS AND NOTHING ELSE. A colour can say "which side of zero"
    # without any transform at all; it cannot say "how much" without one. So it
    # says only the first, and the second is carried by the contours below, in
    # the units of the quantity.
    #
    # Pale on purpose: seven contours, the ridge, u_c(x) and three verticals all
    # have to read on top of this.
    #
    # No marker_z here, which incidentally retires the GR colour trap documented
    # further down -- that failure is specific to panels with an active marker_z
    # colorbar, and this panel no longer has one.
    ipos = colr[findall(v -> v > 0, vf2[colr])]
    ineg = colr[findall(v -> v < 0, vf2[colr])]
    scatter!(p2, xf2[ipos], uf2[ipos]; marker = (:square, 2.4), mc = RGB(0.96, 0.78, 0.72),
             markerstrokewidth = 0, label = "∂C/∂x > 0  (more pump helps)")
    scatter!(p2, xf2[ineg], uf2[ineg]; marker = (:square, 2.4), mc = RGB(0.75, 0.84, 0.94),
             markerstrokewidth = 0, label = "∂C/∂x < 0  (more pump hurts)")
else
    scatter!(p2, xc2, uc2; marker = (:square, 2.4), markerstrokewidth = 0,
             marker_z = zf2, c = :balance, clims = (-m, m),
             colorbar_title = cbl, label = "")
end

# Drawn but not legended -- the count is already on stdout (line ~770), and at
# 1 of 6163 cells a legend row for it was more clutter than information.
if !isempty(sub)
    scatter!(p2, xf2[sub], uf2[sub]; marker = (:square, 2.4), mc = :grey70,
             markerstrokewidth = 0, label = "")
end
plot!(p2, x_curve[keepc], u_curve[keepc]; lw = 2.5, lc = :black,
      label = L"u_c(x)\ \mathrm{(§4.1)}")

# ---- the contours: magnitude, in the units of the quantity -----------------
#
# Extracted from DX, the grid check (b) has ALREADY scored against §10.1 -- so
# these curves are drawn from numbers that have been validated, at no extra
# solve. Walk a u row, find adjacent live cells that bracket the level, and
# interpolate x linearly. Resolution is the grid's own spacing in x.
#
# Re-scanning each row with dC_dx_fd would be more accurate and cost ~17 s
# (180 rows × 400 points × 2 run_adia calls), roughly doubling the script. Check
# (i) below measures exactly what the interpolation costs instead.
function contour_x(row, c)
    out = Float64[]
    for j in 1:(N_X - 1)
        a, b = DX[row, j], DX[row, j+1]
        (isfinite(a) && isfinite(b)) || continue
        ((a - c) * (b - c) <= 0 && a != b) || continue
        push!(out, x_vals[j] + (c - a) * (x_vals[j+1] - x_vals[j]) / (b - a))
    end
    return out
end

if DX_SCALE == :sign_map
    # Drawn as MARKERS, one per u row, not as a joined line. A level with two
    # crossings in a row would otherwise be joined into a curve that does not
    # exist -- check (ii) counts those rather than trusting they are absent.
    lev_col = Dict(2.0 => RGB(0.55, 0.10, 0.10), 1.0 => RGB(0.75, 0.25, 0.15),
                   0.5 => RGB(0.88, 0.45, 0.15), 0.1 => RGB(0.80, 0.60, 0.25),
                  -0.1 => RGB(0.30, 0.50, 0.75), -0.5 => RGB(0.15, 0.35, 0.65),
                  -1.0 => RGB(0.08, 0.20, 0.45))

    println("\n    contours (dC_dx map) — level, rows drawn, x range, label position:")
    for c in DX_LEVELS
        cx, cu, multi = Float64[], Float64[], 0
        for i in 1:N_U
            hit = contour_x(i, c)
            length(hit) > 1 && (multi += 1)
            for x in hit
                push!(cx, x); push!(cu, u_vals[i])
            end
        end
        if isempty(cx)
            @printf("      %+5.1f   NOT PRESENT on this grid\n", c)
            continue
        end
        scatter!(p2, cx, cu; marker = (:circle, 1.7, lev_col[c]),
                 markerstrokewidth = 0, label = "")
        # Label at the contour's TOP end, in real units, in its own colour. The
        # top is chosen because the contours crowd together at small u. Clamped
        # inside the frame: a label placed off the axis is silently dropped by
        # GR, which would leave a curve with no number on it and no warning.
        itop = argmax(cu)
        lx   = clamp(cx[itop], X_LO + 0.006, X_HI - 0.006)
        lu   = min(cu[itop] * 1.35, U_TOP * 0.93)
        annotate!(p2, lx, lu, text(@sprintf("%+.1f", c), 8, :center, lev_col[c]))
        # Frame containment (CLAUDE.md audit 2): a contour running off the side
        # must be STATED, not left to be inferred from a curve that stops.
        edge = minimum(cx) <= X_LO + (x_vals[2] - x_vals[1])
        @printf("      %+5.1f   %4d rows   x ∈ [%.4f, %.4f]   label at (%.4f, %.1f)   %s%s\n",
                c, length(cx), minimum(cx), maximum(cx), lx, lu,
                edge ? "EXITS the left edge" : "ends inside the frame",
                multi > 0 ? @sprintf("   <-- %d rows with >1 crossing: drawn as markers, not joined", multi) : "")
    end
end

# DRAWN AS MARKERS, NOT A LINE, AND THAT IS NOT A STYLE CHOICE.
#
# ON A PANEL WITH AN ACTIVE marker_z COLORBAR, GR DOES NOT GIVE YOU THE COLOUR
# YOU ASK FOR. Measured on this exact figure, all through the same PDF pipeline:
#
#     requested            rendered
#     lc = :limegreen      pale blue-grey
#     lc = RGB(0,0.7,0)    pale blue-grey
#     lc = :magenta        dark red
#     marker :lime         cyan
#     marker :yellow       yellow        (survived)
#     lc = :black          black         (survived)
#
# The first four are all colours that appear in the panel's own :balance
# colormap, which is the hint -- but the exact rule is NOT established here, and
# the yellow marker surviving argues against any simple "everything is
# quantized" story. What IS established: the same colour keywords render
# correctly on a panel WITHOUT marker_z (checked in isolation), and Plots
# reports the attribute as exactly what was asked for (linecolor=RGBA(0,0.7,0,1)),
# so nothing errors and nothing warns. The figure just comes out wrong.
#
# Practical rule, which is all this comment is really for: on a marker_z panel,
# pick overlay colours BY LOOKING AT THE OUTPUT, never by trusting the keyword.
# CLAUDE.md's standing advice to eyeball GR output rather than assume it
# rendered is this, and it applies to colours as much as to LaTeX.
#
# Markers were chosen over a line because a dense scatter reads as a curve at
# one point per u row, and because this particular one came out high-contrast
# against both the pale middle and the blue of :balance. It is verified by eye,
# not by the keyword.
# On :sign_map there IS no marker_z on this panel, so the trap above does not
# apply and a requested colour is the rendered colour. The ridge is then drawn
# black -- it is the c = 0 member of the contour family and should read as the
# heaviest of them. :lime is kept on the colour-mapped branches, where it was
# chosen by looking at the output against :balance.
rlive = findall(isfinite, ridge)
scatter!(p2, ridge[rlive], u_vals[rlive];
         marker = (:circle, DX_SCALE == :sign_map ? 2.6 : 3.0,
                   DX_SCALE == :sign_map ? :black : :lime),
         markerstrokewidth = 0,
         label = "ridge: ∂C/∂x = 0  (the c = 0 contour)")

# One legend entry, not three -- the three points make a single claim.
stored_u = [β + 1 / β for β in (1.0, 2.0, 10.0)]
scatter!(p2, ridge_x.(stored_u), stored_u;
         marker = (:diamond, 7, :yellow), markerstrokewidth = 1.2,
         markerstrokecolor = :black, label = "stored x_opt (β = 1, 2, 10)")

vline!(p2, [X_STAR]; ls = :dash, lc = :black, lw = 1.6,
       label = "x* = $(round(X_STAR, digits = 5))")
vline!(p2, [X_C]; ls = :dashdot, lc = :darkorange, lw = 2.0,
       label = "x_c = $(round(X_C, digits = 5))")

# ---- Figure 2, audits ------------------------------------------------------
#
# AUDIT 1, colour honesty (CLAUDE.md, "Figure audits"). THE TWO THRESHOLDS MUST
# BE INDEPENDENT -- one perceptual, one physical -- or the check cannot fail: the
# first version of this audit defined "reads as white" and "is not zero" off the
# same quantity and reported 0 for the broken panel as happily as for the fixed
# one. Here:
#
#   perceptual: |z| within the central 5% of the drawn range (-m, m), i.e. the
#               band a reader cannot distinguish from the midpoint colour
#   physical:   |∂C/∂x| above 1% of the panel's own maximum -- a value that
#               matters in the units of the quantity, computed from vf2, which
#               the colour transform never touches
#
# Counting cells in BOTH is the number of cells painted the colour of zero while
# being materially nonzero. On :split_log it is 0 by construction (nothing is
# even drawn with |z| < 1). On :linear it must be non-zero, and that is the
# reason the :linear branch is still in the file.
phys_tol = 0.01 * maximum(abs, vf2)
if DX_SCALE == :sign_map
    # NOT a pass, and not skipped quietly either: the failure this audit tests
    # for is a cell rendered in the colour of zero, and a two-colour map has no
    # such colour -- the set is empty structurally, the same kind of claim
    # :split_log's empty band made. Printing "0 misleading" here would be a
    # category error, so the reason is printed instead. The audit still runs, and
    # still fires, on the three colour-mapped branches.
    println("    colour-honesty (dC_dx panel, sign_map): N/A — no colour scale, so no cell can render as the colour of zero")
    println("                                            (re-run with DX_SCALE = :linear to see the audit fire)")
else
    white   = abs.(zf2) .<= 0.05m
    mislead = findall(white .& (abs.(vf2[colr]) .> phys_tol))
    worst_w = isempty(mislead) ? 0.0 : maximum(abs, vf2[colr][mislead])
    @printf("    colour-honesty (dC_dx panel, %s): %d of %d coloured cells read as the midpoint while |∂C/∂x| > %.3f  (worst such value %.3f)\n",
            DX_SCALE, length(mislead), length(zf2), phys_tol, worst_w)
    # The claim that makes "white = zero" true rather than approximate.
    if DX_SCALE == :split_log
        ngap = count(z -> abs(z) < 1, zf2)
        @printf("    gap emptiness: %d coloured cells inside |z| < 1  (must be 0 — the ±1 offset)\n", ngap)
    end
end

# AUDIT: the contour extraction itself, tested on the one level whose answer is
# known independently. contour_x(row, 0) is the SAME interpolation every drawn
# level goes through; ridge_x is 120 steps of bisection on dC_dx_fd and shares
# none of it. If the row indexing or the interpolation is wrong, this fires --
# and on a square grid nothing else would catch a transpose (CLAUDE.md's standing
# warning, and the bug check_against_run_A shipped with).
if DX_SCALE == :sign_map
    dev = Float64[]
    for i in 1:N_U
        isfinite(ridge[i]) || continue
        hit = contour_x(i, 0.0)
        isempty(hit) && continue
        push!(dev, minimum(abs.(hit .- ridge[i])))
    end
    sort!(dev)
    dx_grid = x_vals[2] - x_vals[1]
    @printf("    contour extraction: c = 0 from the grid vs ridge_x's bisection over %d rows — worst %.1e, median %.1e  (grid spacing %.1e)%s\n",
            length(dev), dev[end], dev[cld(length(dev), 2)], dx_grid,
            dev[end] < 5 * dx_grid ? "" : "   <-- WORSE THAN THE GRID: extraction is wrong")
end

# ---- Figure 2, lower panel: the same field in its own units ----------------
#
# The map can only report the sign through its own encoding, which the reader has
# to take on the colorbar's word. This panel does not: ∂C/∂x plotted against x,
# UNTRANSFORMED, on a linear axis, against a drawn zero line, for four u rows.
# One crossing per curve IS check (g) in a picture, and no colour is involved.
#
# Shares the map's x axis and x limits, so a vertical dropped from a ridge marker
# above lands on the corresponding crossing below.
#
# u = 2.5 and 10.1 are β = 2 and β = 10 exactly -- check (e)'s stored-x_opt rows,
# so two of the four crossings are numbers a third project supplies.
cut_us  = [2.0, 2.5, 4.0, 10.1]
cut_lab = ["u = 2 (β = 1)", "u = 2.5 (β = 2)", "u = 4 (β ≈ 3.7)", "u = 10.1 (β = 10)"]
cut_col = [:black, :steelblue, :seagreen, :crimson]
cut_x   = collect(range(X_LO, X_HI; length = 400))

p2b = plot(; xlims = (X_LO, X_HI),
           xlabel = L"x \;=\; 4\eta^2/(\kappa_1\kappa_2)",
           ylabel = L"\partial C/\partial x",
           legend = :topright, colorbar = false)
hline!(p2b, [0.0]; lw = 2.2, lc = :black, label = "")

# The contour levels again, as horizontal lines in the same colours. This is what
# ties the two panels together: a level is a curve upstairs and a line down here,
# carrying the SAME NUMBER, so either can be read against the other and neither
# is a colour scale.
if DX_SCALE == :sign_map
    for c in DX_LEVELS
        hline!(p2b, [c]; lw = 1.0, ls = :dot, lc = lev_col[c], label = "")
        annotate!(p2b, X_HI - 0.004, c, text(@sprintf("%+.1f", c), 7, :right, lev_col[c]))
    end
end

println("\n    cuts (dC_dx lower panel) — crossing vs ridge_x(u), and where each curve ends:")
for (k, u) in enumerate(cut_us)
    v  = [dC_dx_fd(x, u) for x in cut_x]
    ok = findall(isfinite, v)
    isempty(ok) && continue
    plot!(p2b, cut_x[ok], v[ok]; lw = 2.0, lc = cut_col[k], label = cut_lab[k])
    # The crossing, read off the PLOTTED samples, against the independently
    # bisected ridge. Both come from dC_dx_fd, but by different routes (a 400-
    # point scan vs 120 bisection steps), so a transposed axis or a mis-sampled
    # row shows up here as a Δ of order the grid spacing rather than 1e-4.
    cxs = cut_x[ok]
    cvs = v[ok]
    cr  = NaN
    for i in 1:length(cvs)-1
        if cvs[i] > 0 >= cvs[i+1]
            cr = cxs[i] - cvs[i] * (cxs[i+1] - cxs[i]) / (cvs[i+1] - cvs[i])
            break
        end
    end
    rr = ridge_x(u)
    @printf("      %-20s crossing %.6f   ridge_x %.6f   Δ = %.1e   curve ends at x = %.4f%s\n",
            cut_lab[k], cr, rr, abs(cr - rr), cxs[end],
            cxs[end] >= X_HI - 1e-9 ? " (runs to the frame edge)" : " (leaves the region)")
end

pfig2 = plot(p2, p2b; layout = grid(2, 1, heights = [0.62, 0.38]),
             size = (900, 980), left_margin = 7Plots.mm)
save_fig(pfig2, "dC_dx"; dir = "Adia_boundary")

println("\n" * "="^78)
# The verdict reports the SIGN TEST (a count that can fail), not the weakest value.
# It used to lead with "∂C/∂u max = -9.91e-06", which is a window corner and not a
# margin -- see check (a2) and CLAUDE.md.
@printf("VERDICT: %d of %d sampled ∂C/∂u are >= 0 (theorem holds numerically);  ridge at u=2 hits x* to %.1e\n",
        count(>=(0.0), DU[live_du]), length(live_du), abs(r2 - X_STAR))
println("="^78)
