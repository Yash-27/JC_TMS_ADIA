using JLD2
using LinearAlgebra
using LaTeXStrings
using Printf               # @printf in the summary block
using Plots.PlotMeasures   # mm, for right_margin

include("plotting_functions.jl")

# ==============================================================================
# compare_runs.jl
#
# Trace distance between two different runs, comparing each model against
# ITSELF across parameter sets:
#
#     D_full(β, x) = ½ ‖ ρ_full^A(β, x)  −  ρ_full^B(β, x) ‖₁
#     D_adia(β, x) = ½ ‖ ρ_adia^A(β, x)  −  ρ_adia^B(β, x) ‖₁
#
# This is NOT the :tracedist stored inside each .jld2. That one is
# full-vs-adiabatic within a single parameter set. These are full-vs-full
# and adia-vs-adia ACROSS two parameter sets, so they measure how much each
# model's two-qubit steady state moves when the cavity decay rates change,
# at fixed dimensionless coordinates (β, x).
#
# Two further panels difference one scalar observable the same way, once per
# model:
#
#     ΔC_full(β, x) = C(ρ_full^A(β, x))  −  C(ρ_full^B(β, x))
#     ΔC_adia(β, x) = C(ρ_adia^A(β, x))  −  C(ρ_adia^B(β, x))
#
# SIGNED, deliberately. Trace distance is a norm and cannot be negative;
# a concurrence difference can, and the sign is the content -- it says
# which parameter set entangles the two qubits more at that (β, x). Taking
# |ΔC| would throw exactly that away and leave a map that only repeats what
# the D panels already show. So ΔC gets a symmetric colour range and a
# diverging colormap, which is what CLAUDE.md issue 2 prescribes for a
# signed quantity and what the in-file :concurrence_diff panel gets wrong.
#
# ΔC_adia against ΔC_full is the sharper version of the D_adia/D_full
# comparison below: it asks whether the elimination reproduces not just some
# change in the state but the change in the one quantity being studied. The
# two can disagree -- D_adia can be tiny while ΔC_adia is not, if the
# adiabatic states move only along directions concurrence is blind to.
#
# Reading the two panels together: where D_adia tracks D_full, the
# adiabatic model responds to the κ asymmetry the same way the full model
# does. Where they part company, the elimination is not carrying the κ
# dependence faithfully -- which the in-file :tracedist cannot show, since
# it never varies the parameters.
#
# WHAT THE AXES DO AND DO NOT HOLD FIXED
#   β and x are dimensionless; the raw parameters behind them are not shared
#   between the runs. run_sweep converts back via
#       η  = sqrt(x κ_1 κ_2 / 4)        g_1 = g_2 sqrt(β κ_1 / κ_2)
#   so at a matched (β, x) the two runs sit at different η and different
#   g_1. The figures answer "how much does the state depend on the κ
#   asymmetry at fixed dimensionless coordinates", not "... at fixed drive".
#
# All four maps are immune to the unverified relative-sign issue in the
# adiabatic model (CLAUDE.md issue 1). For the full panels, neither operand
# touches J_adia at all. For the adiabatic trace distance, if the issue is
# real the correction is the same fixed local σ_z rotation on qubit 1 in
# both runs, and trace distance is invariant under a unitary applied to both
# arguments. ΔC_adia is safer still: concurrence is a local-unitary
# invariant, so each operand is individually unchanged by the rotation and
# the argument does not even need the two runs to share it. The issue would
# still shift the in-file :tracedist; it cannot shift any of these.
#
# Run with:   julia --project=. compare_runs.jl
#
# Everything is at top level, so D_grids and figs stay live in the session:
#   figs[:adia]
#   plot_map((β_vals = A.β_vals, x_vals = A.x_vals), D_grids[:adia])
# ==============================================================================


# ==============================================================================
# WHICH TWO RUNS, AND WHICH MODELS
# ==============================================================================

# (κ_1, κ_2, g_2) for each run -- the three numbers OUTFILE keys on.
const RUN_A = (κ_1 = 2.0, κ_2 = 2.0, g_2 = 1.0)
const RUN_B = (κ_1 = 2.5, κ_2 = 1.5, g_2 = 1.0)

const FIG_DIR = "Full_vs_Full"

# One row per figure. `key` selects which stored state grid is compared,
# the rest is labelling and per-figure styling:
#
#   kw        extra keywords forwarded straight to plot_map -> scatter.
#             Safe for anything plot_map does not already set explicitly;
#             right_margin and xticks are not in that list. Do NOT put
#             xscale, xlabel, ylabel, c, size, legend, colorbar,
#             markerstrokewidth or markerstrokecolor here -- plot_map sets
#             those itself and whether the splat overrides or errors on the
#             duplicate is unverified.
#
#   decade_x  label the β axis at powers of ten as 10^n rather than
#             letting GR fall back to its 1E-02 form.
#
# Drop a row to skip that figure.
const MODELS = (
    (key = :full, word = "\\mathrm{full}",      tag = "full_full",
     kw = NamedTuple(), decade_x = false),

    (key = :adia, word = "\\mathrm{adiabatic}", tag = "adia_adia",
     kw = (right_margin = 8mm,), decade_x = true),
)

resultsfile(r) = "results_k1_$(r.κ_1)_k2_$(r.κ_2)_g2_$(r.g_2).jld2"


# ==============================================================================
# DATA ACCESS
# ==============================================================================

# ------------------------------------------------------------------------
# load_states: pull both stored state grids and their axes back in.
#
# load_results in plotting_functions.jl deliberately does not carry
# full_data or adia_data -- it only needs the observable grids -- so this
# reads the file directly rather than extending it. status_grid comes along
# because a failed solve must not be compared against anything.
#
# The two grids are keyed :full / :adia so MODELS can select one by symbol.
#
# `outs` is the observable grids the sweep already computed. The concurrence
# panel reads :concurrence_full straight out of it rather than recomputing
# Wootters' construction here -- unlike tracedist below, which has to be
# duplicated because it is never stored per-point. analyze() prefills every
# observable grid with NaN and only overwrites on success, so a NaN in there
# already means "no value", consistent with the rest of the project.
#
# Axes come from the file, never from config.jl, which may have been edited
# since either run.
# ------------------------------------------------------------------------
function load_states(fname)
    isfile(fname) || error("no such file: $fname (working directory is $(pwd()))")
    return jldopen(fname, "r") do f
        (states      = (full = f["full_data"], adia = f["adia_data"]),
         outs        = f["outs"],
         β_vals      = f["β_vals"],
         x_vals      = f["x_vals"],
         params      = f["params"],
         status_grid = f["status_grid"],
         fname       = fname)
    end
end

# ------------------------------------------------------------------------
# assert_same_grid: point-by-point comparison is only meaningful if the two
# runs were swept over the identical (β, x) mesh.
#
# Checked elementwise, not by length or endpoints: two runs with the same
# len and the same range endpoints but a different spacing rule would pass
# a length check and then silently compare mismatched physical points.
# ------------------------------------------------------------------------
function assert_same_grid(A, B)
    size(A.states.full) == size(B.states.full) ||
        error("grid size mismatch: $(size(A.states.full)) in $(A.fname) vs " *
              "$(size(B.states.full)) in $(B.fname)")
    A.β_vals == B.β_vals ||
        error("β_vals differ between $(A.fname) and $(B.fname) -- points are not comparable")
    A.x_vals == B.x_vals ||
        error("x_vals differ between $(A.fname) and $(B.fname) -- points are not comparable")
    return nothing
end


# ==============================================================================
# THE OBSERVABLE
# ==============================================================================

# ------------------------------------------------------------------------
# tracedist: D = ½‖ρ − σ‖₁ = ½ Σ|λ_i| over eigenvalues of the difference.
#
# Same construction as calc_tracedist in run_jc_pump_disp_asy.jl, repeated
# here on purpose: that file has a top-level begin block that runs the
# whole sweep on include, so it cannot be included just to borrow three
# lines.
#
# Both states are normalized first, so a trace defect left by the iterative
# solver in either run cannot masquerade as physical disagreement.
# Hermitian() routes to the symmetric eigensolver, returning real
# eigenvalues directly.
# ------------------------------------------------------------------------
function tracedist(ρ, σ)
    Δ = ρ / tr(ρ) - σ / tr(σ)
    return 0.5 * sum(abs, eigvals(Hermitian(Δ)))
end

# ------------------------------------------------------------------------
# cross_tracedist: the (N_β × N_x) grid of D between the two runs, for one
# model -- :full compares ρ_full against ρ_full, :adia compares ρ_adia
# against ρ_adia. It never crosses the two, which is what the in-file
# :tracedist already does.
#
# NaN means "no data", as everywhere else in this project: a point where
# either run's solve did not return :ok becomes a gap in the heatmap rather
# than a number computed from a NaN-filled matrix. (analyze() in the sweep
# skips this check and lets eigvals throw once per failed point -- not
# repeated here.)
# ------------------------------------------------------------------------
function cross_tracedist(A, B, key::Symbol)
    SA, SB = A.states[key], B.states[key]
    D = fill(NaN, size(SA))
    for j in axes(D, 2), i in axes(D, 1)
        (A.status_grid[i, j] === :ok && B.status_grid[i, j] === :ok) || continue
        D[i, j] = tracedist(SA[i, j], SB[i, j])
    end
    return D
end

# ------------------------------------------------------------------------
# cross_concurrence: the (N_β × N_x) grid of C_A − C_B for one model --
# :concurrence_full compares the full model against itself across the two
# runs, :concurrence_adia the adiabatic model against itself. Like
# cross_tracedist, it never crosses the two models; that is what the in-file
# :concurrence_diff already does.
#
# Signed, and A − B in that order to match the A:(...) B:(...) convention in
# the title -- positive means run A is the more entangled of the two at that
# point.
#
# Note what this is NOT: it is a difference of two concurrences, not the
# concurrence of a difference. Same distinction CLAUDE.md draws for the
# in-file :concurrence_diff, and for the same reason -- ρ^A − ρ^B has zero
# trace and is not PSD, so Wootters' construction does not apply to it.
#
# Gated on status_grid like cross_tracedist, AND on both operands being
# finite: analyze() does not consult status_grid (CLAUDE.md issue 7), so an
# observable can be NaN at a point whose solve was perfectly :ok. Either
# way the point drops out rather than propagating a NaN into the colour
# range.
#
# Missing key is an error, not a NaN grid: it means one of the two runs was
# swept without :concurrence in ACTIVE_OBSERVABLES, and a blank figure is a
# worse answer than a message saying so.
# ------------------------------------------------------------------------
function cross_concurrence(A, B, key::Symbol = :concurrence_full)
    for r in (A, B)
        haskey(r.outs, key) || error(
            "$(r.fname) has no :$key -- it was swept without :concurrence in " *
            "ACTIVE_OBSERVABLES. Stored observables: $(sort(collect(keys(r.outs))))")
    end
    CA, CB = A.outs[key], B.outs[key]
    D = fill(NaN, size(CA))
    for j in axes(D, 2), i in axes(D, 1)
        (A.status_grid[i, j] === :ok && B.status_grid[i, j] === :ok) || continue
        (isfinite(CA[i, j]) && isfinite(CB[i, j]))                   || continue
        D[i, j] = CA[i, j] - CB[i, j]
    end
    return D
end


# ==============================================================================
# LABELS
# ==============================================================================

# A parameter goes inside the per-run parentheses only if it distinguishes
# A from B. Anything the two runs share is stated once, outside them --
# repeating it under both labels is the param_string failure mode: it reads
# as a difference and isn't one. (γ and γ_φ stay out entirely; at 0 they say
# nothing about either run.)
κ_string(p) = "\\kappa_1{=}$(fmt_num(p.κ_1)),\\; \\kappa_2{=}$(fmt_num(p.κ_2))"
g_string(p) = "g_2{=}$(fmt_num(p.g_2))"

# Titles and colorbar labels are built here rather than registered in
# TITLE_NAME / CBAR_LABEL because plot_map only consults those tables for a
# Symbol key, and a Symbol key would also drag in show_params -- which
# prints a single parameter block. There are two parameter sets in these
# figures and one of them would be a lie.
#
# g_2 is equal in both runs in every comparison done so far, so it lands in
# the shared trailing block. The branch is not dead code: OUTFILE keys on
# g_2, so two runs differing only in it are a legal pair to compare, and
# then it is the whole point of the figure and has to sit inside A:(...)
# and B:(...).
function fig_title(m, A, B)
    shared_g = A.params.g_2 == B.params.g_2
    per(p)   = κ_string(p) * (shared_g ? "" : ",\\; " * g_string(p))
    return latexstring(
        "\\mathrm{Trace\\; distance\\; between\\; }" * m.word * "\\mathrm{\\; models}" *
        "\\quad A:(" * per(A.params) * ")" *
        "\\quad B:(" * per(B.params) * ")" *
        (shared_g ? "\\quad " * g_string(A.params) : ""))
end

fig_cbar(m, factor = "") =
    latexstring("D(\\rho^{A}_{" * m.word * "}, \\rho^{B}_{" * m.word * "})" *
                (isempty(factor) ? "" : "\\; " * factor))

# The concurrence panels reuse the same shared/per-run split, so all four
# figures stay readable side by side. They take the same `m` row as the
# trace-distance panels -- see the concurrence loop for why MODELS drives
# both families.
function conc_title(m, A, B)
    shared_g = A.params.g_2 == B.params.g_2
    per(p)   = κ_string(p) * (shared_g ? "" : ",\\; " * g_string(p))
    return latexstring(
        "\\mathrm{Concurrence\\; difference,\\; }" * m.word * "\\mathrm{\\; model}" *
        "\\quad A:(" * per(A.params) * ")" *
        "\\quad B:(" * per(B.params) * ")" *
        (shared_g ? "\\quad " * g_string(A.params) : ""))
end

conc_cbar(m, factor = "") =
    latexstring("(C(\\rho^{A}_{" * m.word * "}) - C(\\rho^{B}_{" * m.word * "}))" *
                (isempty(factor) ? "" : "\\; " * factor))

# ------------------------------------------------------------------------
# sym_clims: a symmetric range (-m, m) about zero for a signed quantity.
#
# The reason is CLAUDE.md issue 2. On a range taken from the data, zero
# lands wherever the data happens to put it, and a diverging colormap's
# neutral midpoint then marks some arbitrary nonzero value -- so a panel
# with, say, all-positive-but-one values reads as half red half blue and
# the sign is unrecoverable by eye. Forcing the range symmetric pins the
# midpoint colour to ΔC = 0, which is the one value that means "the two
# runs agree here".
#
# Returns nothing if there is nothing finite, so the caller can fall back to
# plot_map's own default rather than passing clims = (0, 0).
# ------------------------------------------------------------------------
function sym_clims(D)
    finite = filter(isfinite, vec(D))
    isempty(finite) && return nothing
    m = maximum(abs, finite)
    (m == 0 || !isfinite(m)) && return nothing
    return (-m, m)
end

# ------------------------------------------------------------------------
# si_scale: pull a common power of ten out of the data, so the colorbar
# ticks read 0, 0.5, 1 ... instead of 0, 5E-15, 1E-14 ...
#
# This fixes a collision, not just an aesthetic. GR draws colorbar_title at
# a FIXED offset to the right of the bar, and the tick labels grow rightward
# from that same bar, so a wide label like "1.25E-14" runs straight through
# the rotated title. right_margin cannot separate them -- it translates bar,
# ticks and title together, leaving their spacing untouched (checked at 8mm
# and 25mm: identical overlap). Narrower labels are the only lever.
#
# The exponent then reappears in the title, set in LaTeX like the rest of
# the figure rather than in GR's E-notation.
#
# It goes AFTER the quantity, as "D(...) x 10^-14" -- read as "multiply the
# tick by this to recover D", the usual convention on a scaled axis. As a
# prefix, "10^-14 D(...)" would literally name the plotted quantity as
# 10^-14 times D, which is the reciprocal of what is drawn.
#
# Returns (scaled_grid, latex_suffix). An exponent of 0 is a no-op with an
# empty suffix, so data of order 1 passes through untouched -- as does
# all-NaN or all-zero data.
# ------------------------------------------------------------------------
function si_scale(D)
    finite = filter(isfinite, vec(D))
    isempty(finite) && return D, ""
    m = maximum(abs, finite)
    (m == 0 || !isfinite(m)) && return D, ""
    e = floor(Int, log10(m))
    e == 0 && return D, ""
    return D ./ 10.0^e, "\\times 10^{$e}"
end

# ------------------------------------------------------------------------
# decade_ticks: β ticks at whole powers of ten, labelled 10^n.
#
# GR's own log-axis labels come out as 1E-02, which reads badly next to the
# LaTeX in the rest of the figure. Ticks sit at the decades, not at the
# swept β values -- those are log-spaced but land off the decades (0.0464,
# 0.215, ...) and would make a cluttered axis.
#
# The range is taken from the data, so this still does the right thing if
# β_range changes.
# ------------------------------------------------------------------------
function decade_ticks(vals)
    lo, hi = extrema(vals)
    exps   = floor(Int, log10(lo)):ceil(Int, log10(hi))
    return ([10.0^e for e in exps], [latexstring("10^{$e}") for e in exps])
end

run_pair_tag() = "_A_k1_$(RUN_A.κ_1)_k2_$(RUN_A.κ_2)" *
                 "_B_k1_$(RUN_B.κ_1)_k2_$(RUN_B.κ_2)" *
                 "_g2_$(RUN_A.g_2)"

fig_name(m)      = "tracedist_" * m.tag * run_pair_tag()
conc_fig_name(m) = "concurrence_diff_$(m.key)" * run_pair_tag()

# :full -> :concurrence_full, :adia -> :concurrence_adia. The sweep names
# the _full/_adia pair by appending the model to the observable (the :single
# expansion in run_jc_pump_disp_asy.jl), so MODELS' own key is exactly the
# suffix needed here.
conc_key(m) = Symbol(:concurrence_, m.key)


# ==============================================================================
# LINE CUTS
# ==============================================================================
#
# A heatmap answers "where is it large"; a line cut answers "how does it grow".
# Holding β fixed and reading D against x is the direction that matters here --
# x is the drive relative to threshold, and the adiabatic model's denominators
# all carry 4η² − κ_1κ_2, so the interesting behaviour is the approach to x = 1.
#
# CUT_β is a REQUESTED value, snapped to the nearest swept β. It is not asserted
# to be on the grid: β_range or len can change, and a cut that silently moves to
# the nearest available point is more useful than one that errors. The chosen
# values are printed so the snap is never invisible.

const CUT_β   = 1.0      # requested β for the first curve (10^0)
const CUT_N   = 2        # that point plus the next (CUT_N - 1) in increasing β
const CUT_KEY = :full    # which D grid to cut -- :full or :adia

# ------------------------------------------------------------------------
# nearest_β: index of the swept β closest to a requested value.
#
# Matched in LOG space, because β_range is swept logarithmically. The swept
# values are 0.01, 0.0464, 0.215, 1.0, 4.64, 21.5, 100 -- on a linear metric
# a request for β = 1 sits far closer to 0.215 and 4.64 than those points are
# to each other in any meaningful sense, and requests anywhere below 1 would
# all collapse toward the same few small values. Log distance is the metric
# the grid was built with, so it is the one to search in.
# ------------------------------------------------------------------------
nearest_β(β_vals, β) = argmin(abs.(log10.(β_vals) .- log10(β)))

# ------------------------------------------------------------------------
# model_by_key: look up a MODELS row by its key rather than by position, so
# reordering MODELS or dropping a row cannot silently cut the wrong grid.
# ------------------------------------------------------------------------
function model_by_key(k::Symbol)
    i = findfirst(m -> m.key === k, MODELS)
    if i === nothing
        have = join(string.(getfield.(MODELS, :key)), ", ")
        error("no MODELS row with key :$k -- have $have")
    end
    return MODELS[i]
end

# ------------------------------------------------------------------------
# plot_linecut: D against x, one curve per β row.
#
# Not built on plot_map -- that draws a scatter over the (β, x) plane with a
# colour axis and forces legend = false, none of which applies to a 1D cut.
# The two share save_fig and apply_theme! and nothing else.
#
# NaN still means "no data": non-finite points are dropped per curve rather
# than breaking the line, consistent with the heatmaps. A row that is
# entirely NaN is skipped with a warning instead of contributing an empty
# legend entry.
#
# Markers are on because the grid is coarse (len = 7). Without them the
# reader cannot tell which points are data and which are interpolation --
# and with a NaN dropped mid-curve, a bare line would hide the gap entirely.
# ------------------------------------------------------------------------
function plot_linecut(x_vals, D, rows, β_vals; kwargs...)
    p     = plot()
    drawn = 0
    for i in rows
        row  = D[i, :]
        keep = isfinite.(row)
        if !any(keep)
            @warn "line cut at β = $(β_vals[i]) is entirely NaN -- skipped"
            continue
        end
        plot!(p, x_vals[keep], row[keep];
              label      = latexstring("\\beta = " * fmt_num(β_vals[i])),
              marker     = :circle,
              markersize = 4)
        drawn += 1
    end
    drawn == 0 && error("nothing finite to plot in any requested row")

    plot!(p;
          xlabel = L"x",
          legend = :topleft,
          size   = (700, 480),
          kwargs...)
    return p
end


# ==============================================================================
# RUN
# ==============================================================================

apply_theme!()

A = load_states(resultsfile(RUN_A))
B = load_states(resultsfile(RUN_B))
assert_same_grid(A, B)

println("A: $(A.fname)")
println("B: $(B.fname)")
println("  grid: $(length(A.β_vals)) β × $(length(A.x_vals)) x")

# plot_map only reads β_vals and x_vals when handed a raw matrix, so a
# two-field NamedTuple is a complete `d` for these figures. Nothing in
# plotting_functions.jl needed changing.
d = (β_vals = A.β_vals, x_vals = A.x_vals)

D_grids = Dict{Symbol,Matrix{Float64}}()
figs    = Dict{Symbol,Any}()

for m in MODELS
    D = cross_tracedist(A, B, m.key)
    D_grids[m.key] = D            # stored unscaled -- si_scale is presentation only

    Dplot, factor = si_scale(D)

    kw = m.decade_x ? (; m.kw..., xticks = decade_ticks(A.β_vals)) : m.kw

    figs[m.key] = plot_map(d, Dplot;
                           title          = fig_title(m, A, B),
                           colorbar_title = fig_cbar(m, factor),
                           kw...)
    save_fig(figs[m.key], fig_name(m); dir = FIG_DIR)
end

# --- Concurrence difference, signed, one panel per model ----------------------
#
# Driven by MODELS rather than a table of its own: "which models exist" is
# one fact, and a second tuple listing :full and :adia again would be a
# second place to edit when a model is added or dropped. Only `key` and
# `word` are taken from the row -- `kw`, `decade_x` and `tag` are styling
# tuned for the trace-distance panels and do not carry over, because the two
# families have different needs. There, D_adia sits at 1e-14 and D_full at
# 1e-1, so only the adiabatic panel had label-width trouble; here both
# panels are the same order of magnitude, so both get the same treatment.
#
# clims comes from the SCALED grid, not the raw one: si_scale divides the
# data but plot_map applies clims to what it is handed, so limits taken
# before scaling would be off by that power of ten and clamp everything to
# one end colour.
#
# :balance is diverging with a neutral midpoint. Unlike the trace-distance
# panels, this is not viridis -- a sequential ramp on signed data is
# precisely the failure CLAUDE.md issue 2 describes.
#
# Each panel takes its OWN symmetric range. The alternative -- one range
# shared across both -- was rejected because the point of the pair is
# whether the two models agree in sign and structure, and a shared range set
# by whichever panel happens to be larger would flatten the smaller one to
# uniform white. The printed maxima below are what to compare magnitudes by.
for m in MODELS
    key = conc_key(m)
    ΔC  = cross_concurrence(A, B, key)
    D_grids[key] = ΔC         # stored unscaled, like the trace-distance grids

    ΔCplot, factor = si_scale(ΔC)
    cl = sym_clims(ΔCplot)

    figs[key] = plot_map(d, ΔCplot;
                         title          = conc_title(m, A, B),
                         colorbar_title = conc_cbar(m, factor),
                         cmap           = :balance,
                         xticks         = decade_ticks(A.β_vals),
                         right_margin   = 8mm,
                         (cl === nothing ? (;) : (; clims = cl))...)
    save_fig(figs[key], conc_fig_name(m); dir = FIG_DIR)
end

# --- Line cut: D vs x at two adjacent β -------------------------------------
#
# Reuses fig_title / fig_cbar so the cut carries the same A:(...) B:(...)
# header as the heatmap it comes from, and fig_cbar doubles as the y label --
# it already names exactly the plotted quantity.
#
# The row range is clamped rather than wrapped: asking for CUT_N points from
# the last β should give the one that exists, not silently fold around to
# β = 0.01 and draw two curves four decades apart under adjacent-looking
# labels.
let m = model_by_key(CUT_KEY)
    D    = D_grids[m.key]
    i0   = nearest_β(A.β_vals, CUT_β)
    rows = i0:min(i0 + CUT_N - 1, length(A.β_vals))

    @printf("\nline cut of :%s at β = %s (requested %g)\n",
            m.key, join([fmt_num(A.β_vals[i]) for i in rows], ", "), CUT_β)

    figs[:linecut] = plot_linecut(A.x_vals, D, rows, A.β_vals;
                                  title  = fig_title(m, A, B),
                                  ylabel = fig_cbar(m))
    save_fig(figs[:linecut], "tracedist_linecut_" * m.tag * run_pair_tag();
             dir = FIG_DIR)
end

# Same try/catch discipline as the report block in run_sweep: a bug in a
# summary line must not discard work that already succeeded. The figures
# are on disk by this point, but the live figs and D_grids are worth
# keeping too.
#
# The two panels are NOT drawn on a shared colour scale -- each takes its
# own range, since plot_map falls back to the data range for a raw matrix.
# The printed maxima below are what to check before reading the two
# figures against each other by eye.
#
# The extremum is reported by largest MAGNITUDE and printed with its sign.
# For the two trace-distance grids that is identical to the maximum, since
# a norm is non-negative; for the signed ΔC grid it is the only sensible
# choice, as the strongest disagreement between the runs may well be
# negative and reporting max() there would name some near-zero point
# instead.
function report(label, D, axes_from)
    finite = filter(isfinite, vec(D))
    nbad   = length(D) - length(finite)
    println("\n$label over $(length(finite)) point(s)" *
            (nbad > 0 ? " ($nbad skipped)" : ""))
    isempty(finite) && return
    iext = argmax(map(v -> isfinite(v) ? abs(v) : -Inf, D))
    # %g, not %f: these values can come out at machine precision, and
    # %.4f renders every one of them as a flat 0.0000.
    @printf("  min       %.4g\n", minimum(finite))
    @printf("  max       %.4g\n", maximum(finite))
    @printf("  mean      %.4g\n", sum(finite) / length(finite))
    @printf("  largest   %.4g  at β = %g, x = %g\n",
            D[iext], axes_from.β_vals[iext[1]], axes_from.x_vals[iext[2]])
end

try
    for m in MODELS
        report("tracedist :$(m.key)", D_grids[m.key], A)
    end
    for m in MODELS
        report("$(conc_key(m)) (signed, A - B)", D_grids[conc_key(m)], A)
    end
catch err
    @warn "summary failed (figures are already saved)" err
end
