# CLAUDE.md

Julia code for a cavity-QED reservoir-engineering study: two qubits, each dispersively
coupled to its own cavity, with the two cavities linked by a two-mode-squeezing (TMS)
drive. The question the code answers is how well the **adiabatically-eliminated 2-qubit
model** reproduces the **full 4-body model** across a 2D parameter sweep.

## Running

Two sweeps, each with its own config and plot driver, sharing everything below them.

```bash
# sweep 1 -- the (β, x) plane at fixed κ_1, κ_2, g_2
julia --project=. -t auto run_jc_pump_disp_asy.jl
julia --project=. run_plot.jl
julia --project=. compare_runs.jl        # two (β, x) runs against each other

# sweep 2 -- the (r, ξ) plane at fixed x, β
julia --project=. -t auto run_r_xi.jl
julia --project=. run_plot_r_xi.jl

# one-off parameter-set tools, no sweep
julia --project=. -t auto compare_param_sets.jl
julia --project=. scale_invariance_test.jl [a b c]

# where the adiabatic concurrence peaks -- adiabatic model only, seconds, no threads
julia --project=. adia_concurrence_max.jl [x_lo x_hi β_lo β_hi len]
```

Both `run_jc_pump_disp_asy.jl` and `run_r_xi.jl` have a top-level `begin` block at the
bottom that executes the whole sweep on include. There is no `main()` guard —
`include`ing either one *runs* it. That is why the observables live in `observables.jl`
rather than inside a runner: neither runner can include the other to borrow code.

Threads matter a lot: the sweep is thread-parallel and `-t auto` (or `-t N`) is the
difference between minutes and hours. `BLAS.set_num_threads(1)` is set deliberately at
include time to stop BLAS from oversubscribing inside the `@threads` loop — don't remove it.

## Pipeline

```
        SHARED CORE -- knows nothing about which sweep is running
        ├─ truncation.jl        estimate_truncation(κ_1,κ_2,η,g_1,g_2) -> (N_1,N_2)
        ├─ jc_pump_disp_asy.jl  run_sim(g_1,g_2,κ_1,κ_2,η,γ,γ_φ,N_1,N_2)
        │                          -> (ρ_full_reduced, ρ_adia)
        │                       run_adia(g_1,g_2,κ_1,κ_2,η,γ,γ_φ) -> ρ_adia
        │                          the adiabatic half alone, no truncation arg
        └─ observables.jl       Obs / OBSERVABLE_REGISTRY / select_observables /
                                analyze / summarize
                 │
    ┌────────────┴────────────────────────────┐
    │                                         │
config_jc_pump_disp_asy.jl              config_r_xi.jl
run_jc_pump_disp_asy.jl                 run_r_xi.jl
   sweeps (β, x) at fixed κ, g_2           sweeps (r, ξ) at fixed x, β
   -> results_k1_*_k2_*_g2_*.jld2          -> results_rxi_*.jld2
    │                                         │
    ├─ run_plot.jl                            └─ run_plot_r_xi.jl
    └─ compare_runs.jl  (two .jld2)
                 │
        plotting_functions.jl   plot_map / save_fig / axis_spec,
                                included by all three plot drivers

    STANDALONE -- no sweep, no .jld2
        compare_param_sets.jl   pairwise D between named parameter sets
        scale_invariance_test.jl  randomized fibre + perturbation control
        adia_concurrence_max.jl   argmax of C_adia over (x, β), via run_adia
```

**`run_adia` is the only definition of the adiabatic model.** It was extracted out of
`run_sim`, which now calls it, so anything wanting `ρ_adia` alone skips the dim-900 Krylov
solve sitting next to it and costs microseconds instead of seconds. The extraction is
behaviour-preserving — verified: `run_sim`'s `ρ_adia` and a direct `run_adia` call at the
same parameters agree to `0.0` exactly, and `adia_concurrence_max.jl` check (d) reproduces
all 49 stored `concurrence_adia` values of run A with worst error `0.00e+00`.

Include order is fixed by each runner's include block. Both configs have no dependencies
by design and can sit anywhere in that order. `observables.jl` needs `LinearAlgebra` and
`Printf` in scope; every driver `using`s them before including it.

## Dimensionless coordinates

With `γ = γ_φ = 0` and all detunings zero the model has **five rates** (`κ_1`, `κ_2`,
`g_1`, `g_2`, `η`) and one unphysical overall frequency scale, so exactly **four**
dimensionless numbers exist:

| symbol | definition | meaning | swept by |
|---|---|---|---|
| `x` | `4η² / (κ_1 κ_2)` | drive relative to the linear OPO threshold | sweep 1, linear |
| `β` | `(g_1² κ_2) / (g_2² κ_1)` | Purcell-rate asymmetry between the arms | sweep 1, log10 |
| `r` | `κ_1 / κ_2` | cavity-decay asymmetry | sweep 2, log10 |
| `ξ` | `4 g_1 g_2 / (κ_1 κ_2)` | coupling strength relative to decay | sweep 2, linear |

The two sweeps are exact complements: together the four coordinates are a complete set,
so sweep 2 covers precisely the directions sweep 1 holds fixed, with no overlap and no
gaps. **`ξ` is the adiabatic elimination's own small parameter** — the elimination is
the `ξ → 0` limit — which is why it is the physically informative axis.

`(g_1, g_2, κ_1, κ_2, η) → (x, β, r, ξ, scale)` is a **bijection**. Two raw parameter
sets therefore agree on all four coordinates *if and only if* one is a uniform rescaling
of the other. This is used in two places below and is worth remembering: there is no
second, independent way to land on the same point.

**Sweep 1** converts back via `η_of(x) = sqrt(x κ_1 κ_2 / 4)` and
`g1_of(β) = g_2 sqrt(β κ_1 / κ_2)`; `κ_1`, `κ_2`, `g_2` stay fixed for a run. Note this
means **`ξ` is not fixed along sweep 1** — with `κ` and `g_2` held, `ξ ∝ g_1 ∝ sqrt(β)`,
and for run A (`κ_1 = κ_2 = 2`, `g_2 = 1`) it is exactly `ξ = sqrt(β)`. So sweep 1
confounds β with ξ and structurally cannot separate them. Sweep 2 at fixed β is what
deconfounds them.

**`x < 1` is a hard constraint.** The adiabatic model has `4η² − κ_1κ_2` in every
denominator, so it diverges at `x = 1`. Sweep 1's `x_range` tops out at `0.7225 = 0.85²`.
Anything approaching 1 is also where the Liouvillian gap closes and the iterative solve
gets slow (see the cost model).

## The (r, ξ) sweep — `run_r_xi.jl`

Sweeps `(r, ξ)` at fixed `(x, β)`, with `κ_geo = sqrt(κ_1 κ_2)` as the pinned scale:

```
κ_1 = κ_geo sqrt(r)              κ_2 = κ_geo / sqrt(r)
η   = (κ_geo/2) sqrt(x)                  <- free of BOTH r and ξ
g_1 = (κ_geo/2) sqrt(ξ) (β r)^(1/4)
g_2 = (κ_geo/2) sqrt(ξ) (β r)^(-1/4)
```

**Pin `κ_geo`, not `κ_2`.** Pinning `κ_2` puts a `sqrt(r)` into η and breaks the
`r ↔ 1/r` symmetry of the solve cost — the gap at `r = 0.1` comes out ~3× smaller than
at `r = 10` purely because `κ_1` shrank, which is units, not physics. `κ_geo = 2.0` is
chosen so `r = 1` reproduces run A's `κ_1 = κ_2 = 2`.

**Grid layout mirrors sweep 1 exactly**: ROW = `r` (horizontal, log10), COLUMN = `ξ`
(vertical, linear), matching ROW = β / COLUMN = x. That is what lets `flatten_grid` and
`plot_map` be reused unchanged.

**Cost model differs.** Sweep 1 uses `dim² / (1 − sqrt(x))`, which collapses to a
constant when x is fixed, and whose derivation assumed `κ_1 = κ_2`. `run_r_xi.jl` uses
the general two-mode gap

```
gap = (κ_geo/4) [ (√r + 1/√r) − sqrt( (√r − 1/√r)² + 4x ) ]
```

which reduces to `(κ_geo/2)(1 − sqrt(x))` at `r = 1`, vanishes at `x = 1` for every `r`
(confirming x alone sets threshold), and is symmetric under `r ↔ 1/r`.

`BIG_DIM = 2000` is inherited from sweep 1 and calibrated against *its* dimension
clusters. Read Pass 1's histogram before trusting it here — the recorded run came out
uniformly `dim = 900`, so the big/small split never engaged.

### The adiabatic model is analytically independent of r and ξ

Substituting the inversion into `jc_pump_disp_asy.jl:95,105-106`, every `r` cancels:

- `H_adia ∝ η g_1g_2 / (η² − κ_1κ_2/4)` — η, `g_1g_2 = ξκ_geo²/4` and `κ_1κ_2 = κ_geo²`
  are all r-free.
- `Γ_sqrt_1 · g_1 ϵ_1 ∝ η g_1 √κ_2 ∝ (βr)^(1/4) r^(-1/4) = β^(1/4)` — r cancels.
- `Γ_sqrt_1 · g_2 ∝ κ_1 √κ_2 g_2 ∝ r^(1/2−1/4−1/4) = r⁰` — r cancels. Same for `J_adia[2]`.

And since `g ∝ sqrt(ξ)`, every term of `L_adia` carries exactly one factor of ξ, so
`L_adia → ξ L_adia` and its **kernel is unchanged**.

So `ρ_adia` is exactly a function of `(x, β)` alone. This **proves analytically**, in
these coordinates, the finding recorded further down that was previously only measured
across the three (β, x) runs. It holds only for `γ = γ_φ = 0`; turning either on breaks
the homogeneity argument, since those rates do not scale with g.

`check_adia_flat` in `run_r_xi.jl` asserts this at the end of every run — every `_adia`
observable must be constant to ~1e-13 across the whole grid. It is a free end-to-end
regression test of the coordinate inversion.

### Recorded run: `results_rxi_x_0.3663_b_0.5_kg_2.0_r_0.1-10.0_xi_0.1-1.0_n_7.jld2`

`x = 0.3663`, `β = 0.5`, `κ_geo = 2`, `len = 7`, `r ∈ [0.1,10]` log, `ξ ∈ [0.1,1]`
linear. All 49 points `:ok`, all `dim = 900`, 19.1 min of solve time.

**Two of the four panels are degenerate in this run, and it is the choice of β that did
it.** `concurrence_adia` is **exactly 0.000000e+00 at every point** — spread `0.0e+00`,
not merely small — so `concurrence_diff` is bit-identical to `concurrence_full`. The
adiabatic model predicts no entanglement anywhere on this plane.

That is not a bug, it is the analytic entanglement boundary. The adiabatic concurrence
depends on β only through `b = β + 1/β ≥ 2`, and is nonzero only for `b < b_max(x)`; the
`b = 2` (β = 1) threshold is `x* = 0.48388826`, the root of `x³ − 2x² + 9x − 4`. At
β = 0.5 (`b = 2.5`) the boundary has already been crossed by `x = 0.3663`. Measured on
run A's grid at `x = 0.3663`:

| β | 0.215 | 0.4 | **0.5** | 0.6 | 0.8 | **1.0** | 1.5 | 2.0 |
|---|---|---|---|---|---|---|---|---|
| `C_adia` | 0 | 0 | **0** | 0.006 | 0.092 | **0.118** | 0.041 | 0 |

The nonzero window at that x is roughly `β ∈ (0.55, 1.9)`. **To get four independent
panels, use `β = 1.0` at `x = 0.3663` (`C_adia = 0.1175`), or keep `β = 0.5` and drop x
to ≈ 0.2475 (`C_adia = 0.1035`).** Both pairs sit on run A's grid, which preserves the
free cross-check described under `compare_param_sets.jl`.

**`tracedist` is fully informative even so**, and is the real content of this run:

```
tracedist   (rows = r, cols = ξ)
 r \ ξ    0.100    0.250    0.400    0.550    0.700    0.850    1.000
 0.100   0.0723   0.1305   0.1689   0.1973   0.2193   0.2370   0.2516
 0.464   0.0419   0.0770   0.0993   0.1155   0.1285   0.1393   0.1488
 1.000   0.0460   0.0820   0.1027   0.1162   0.1257   0.1327   0.1382
 4.642   0.0757   0.1162   0.1348   0.1456   0.1527   0.1578   0.1618
10.000   0.0921   0.1284   0.1434   0.1529   0.1602   0.1666   0.1725
```

Three things to read off it:

- **The error grows monotonically with ξ at every r**, as the elimination requires. At
  `r = 1` it runs 0.046 → 0.138 over the decade.
- **It is sub-linear in ξ, not linear.** Log-log slopes at `r = 1` are 0.63 (ξ 0.1→0.25),
  0.48, then 0.39 — decreasing, i.e. saturating. The `D ∝ ξ` prediction is a small-ξ
  statement and **this grid is outside that regime**; `ξ = 0.1` is already too large.
  Reaching the linear regime needs `ξ ≲ 0.01` on a log axis.
- **It is not symmetric under `r ↔ 1/r`** (0.2516 at `r = 0.1` vs 0.1725 at `r = 10`,
  at ξ = 1). That is correct and is a consequence of `β ≠ 1`: the exact arm-swap symmetry
  of the full Liouvillian acts as `β → 1/β` *and* `r → 1/r` together, so `r ↔ 1/r`
  alone is a symmetry only on the `β = 1` slice. **At `β = 1` every panel of this sweep
  would be mirror-symmetric about `r = 1`**, which is a second free self-test but also
  makes half the r axis redundant. `β = 0.5` buys non-redundant r coverage — the one
  advantage it has over `β = 1`.

## Adding an observable

Three edits, nothing else:

1. Write `calc_foo(ρ)` (or `calc_foo(ρ_full, ρ_adia)`) in `observables.jl`.
2. Add one `Obs(:foo, :single|:compare, calc_foo, "Label")` row to `OBSERVABLE_REGISTRY`.
3. Name `:foo` in `ACTIVE_OBSERVABLES` (sweep 1) or `ACTIVE_OBSERVABLES_RXI` (sweep 2).

Because `observables.jl` is shared, step 1 makes the observable available to **both**
sweeps at once.

`:single` auto-expands into `foo_full`, `foo_adia`, `foo_diff` (= full − adia).
`:compare` produces one grid under `:foo`. `select_observables` errors loudly on an
unknown name — it is called *before* the sweep so a typo fails in seconds, not hours.

Observables are pure functions of the stored 4×4 matrices. That independence is what
allowed Section 2 of `run_jc_pump_disp_asy.jl` to be lifted verbatim into
`observables.jl` — which had to happen once there were two runners, since including
either one runs a full sweep. `analyze` reads only `res.full` / `res.adia` and their
shape, so it is axis-agnostic and both sweeps share it unchanged.

## Plotting

`run_plot.jl` is the single-run entry point (`compare_runs.jl` is the two-run one; both
`include` `plotting_functions.jl`, and nothing else should). It runs at top level so `d`
and `figs` stay live in the session:

```julia
figs[:tracedist]
plot_map(d, :tracedist; markersize = 10, show_params = false)
plot_map(d, d.outs[:concurrence_full])   # or any raw (N_β × N_x) matrix
```

`run_plot.jl` includes the config **only to build the filename** —
`results_k1_$(κ_1)_k2_$(κ_2)_g2_$(g_2).jld2`. Everything else (axes, parameters in the
title) comes out of the `.jld2`. Practical consequence: **editing the config to try new
parameters breaks plotting until you re-run the sweep**, because `run_plot.jl` will go
looking for a results file that doesn't exist yet.

Three **(β, x)** results files are currently on disk, all `len = 7`, same `x_range` and
`β_range`, all points `:ok`, all carrying the same four observables. The letters are the
names used throughout this file and in `compare_runs.jl`. (A fourth file,
`results_rxi_*.jld2`, is from the *other* sweep — different key names, different loader,
not comparable to these and not readable by `compare_runs.jl`.)

| | file | κ_1 | κ_2 | g_2 |
|---|---|---|---|---|
| **A** | `results_k1_2.0_k2_2.0_g2_1.0.jld2` | 2.0 | 2.0 | 1.0 |
| **B** | `results_k1_2.5_k2_1.5_g2_1.0.jld2` | 2.5 | 1.5 | 1.0 |
| **C** | `results_k1_2.5_k2_1.5_g2_0.75.jld2` | 2.5 | 1.5 | 0.75 |

B and C differ **only in `g_2`**; A differs from both in κ. That is what makes the three
useful together — see the finding below, which needs a g_2-only pair to state.

The config sits at C's parameters, so `run_plot.jl` currently resolves to C.
`compare_runs.jl` uses all three.

Hardcoded on purpose in `plot_map`, because there is only one right answer:

- the ROW coordinate horizontal on log10, the COLUMN coordinate vertical linear
- the flattening convention — `vec` is column-major, so the **row** index varies
  fastest; `outer` on the row values and `inner` on the column values reproduces that.
  Swapping them silently transposes the figure with no error, and the length assertions
  won't catch it on a square grid, which `len = 7` on both axes always is.
- NaN points dropped rather than plotted
- colorbar forced on (GR does not infer it from `marker_z`)

**Which** coordinates those are is per-file, not hardcoded, since there are now two
sweeps. `axis_spec(d)` returns `d.axes` when the loaded run carries one — `load_results_rxi`
attaches `(xvals, yvals, xlab, ylab, xsc, ysc)` for the (r, ξ) plane — and otherwise
falls back to the original β/x behaviour. The fallback is load-bearing: the three stored
(β, x) `.jld2` files predate the field, and `compare_runs.jl` hands `plot_map` a bare
`(β_vals = …, x_vals = …)` tuple it builds inline. Both still work untouched.

`param_string` branches on the same distinction, detecting `:κ_geo`. For a (β, x) run it
prints `κ_1, κ_2, g_2`; for an (r, ξ) run those all vary per grid point, so printing them
would be actively wrong, and it prints the fixed `x, β, sqrt(κ_1κ_2)` instead.

`run_plot_r_xi.jl` carries two per-key overrides that `run_plot.jl` does not.
`:concurrence_diff` gets `sym_clims` + `:balance` (issue 2's prescription, applied at the
call site). `:concurrence_adia` gets `flat_safe_clims`, which exists because this sweep
predicts that panel is *constant*: a data range ~1e-15 wide makes GR emit
`Rectangle definition is invalid in routine CELLARRAY` and drop the colorbar entirely.
The guard is **relative** — an `lo == hi` test does not catch a spread in the 15th digit —
and widens a flat panel to `(0, 2v)`. Anchoring at zero rather than a narrow band around
`v` also keeps GR's tick labels short enough not to overprint the colorbar title, the
same collision `si_scale` works around in `compare_runs.jl`.

Per-figure: `title`, `colorbar_title`, `clims`, `markersize`, `cmap`, `show_params`.

Three lookup tables drive the labelling, all keyed by observable symbol: `TITLE_NAME`,
`CBAR_LABEL`, `CLIM_SOURCE`. A new observable needs entries in the first two or it falls
back to a raw symbol name rendered in math mode. `CLIM_SOURCE` maps a key to whichever
key's data range sets its colour limits — all three concurrence panels point at
`concurrence_full` so they're comparable by eye. **See issue 2 for why that is wrong for
the `_diff` panel.**

Note `:concurrence_diff` is `C(ρ_full) − C(ρ_adia)`, two concurrences subtracted — *not*
`C(ρ_full − ρ_adia)`. The difference of two density matrices has zero trace and isn't
PSD, so Wootters' construction doesn't apply to it. The labels say what is actually
computed; keep it that way.

Figures go to `Full_vs_Adia/` (what `run_plot.jl` passes), `Full_vs_Full/`
(`compare_runs.jl`) and `Full_vs_Adia_rxi/` (`run_plot_r_xi.jl`), PDF only. `save_fig`'s
own default is `figures3/`, and `figures/`, `figures2/` also exist from earlier runs —
six destinations, none tracked by git.

## Comparing runs — `compare_runs.jl`

The second entry point. It compares each model against *itself* across parameter sets,
which the single-run figures structurally cannot do.

`PAIRS` at the top lists the comparisons, each a `(label, a, b)` row naming two of the
`RUN_A` / `RUN_B` / `RUN_C` triples. All three unordered pairs are currently listed, so
one invocation emits a complete set: any two of A, B, C can be read against each other.

**A and B are the two SLOTS of a pair, not two fixed runs.** C sits in the `a` slot of
both its pairs, and `RUN_A` sits in the `b` slot of `:C_vs_A`. The slots are not
interchangeable — every difference is `a − b` and every title reads `A:(…) B:(…)`, so
swapping a row flips the sign of both its ΔC panels.

**Five figures per pair, 15 in total.** Four heatmaps on the (β, x) mesh:

| grid key | quantity | signed? |
|---|---|---|
| `:full` | `½‖ρ_full^A − ρ_full^B‖₁` | no — a norm |
| `:adia` | `½‖ρ_adia^A − ρ_adia^B‖₁` | no — a norm |
| `:concurrence_full` | `C(ρ_full^A) − C(ρ_full^B)` | **yes** |
| `:concurrence_adia` | `C(ρ_adia^A) − C(ρ_adia^B)` | **yes** |

plus one **line cut** of `D` against x at fixed β, which answers "how does it grow"
where a heatmap answers "where is it large". x is the direction that matters, since the
adiabatic denominators all carry `4η² − κ_1κ_2` and the interesting behaviour is the
approach to `x = 1`. Controlled by four constants, none of them per-pair on purpose —
the cut is a way of reading one `D` grid, not a property of the comparison, and two pairs
cut at different β could not be read against each other:

- `CUT_KEY` — which `D` grid to cut (`:full`)
- `CUT_β` — a *requested* β, snapped to the nearest swept value in **log** space (the
  grid is log-spaced, so a linear metric would collapse every request below 1 onto the
  same few points). The snapped values are printed, so the snap is never invisible.
- `CUT_N` / `CUT_N_LOW` — how many grid points to take above and below the snapped β.
  Clamped at both ends rather than wrapped, so asking past either end gives the points
  that exist instead of folding four decades around.

This is **not** the in-file `:tracedist` / `:concurrence_diff`, which are full-vs-adia
*within* one run. These never cross the two models; they cross the two *parameter sets*.
The distinction is load-bearing enough that the figure titles say `full models`, plural —
singular "trace distance, full" reads as the in-file observable, which is a different
quantity.

`MODELS` declares `:full` and `:adia` once and drives both figure families. The
trace-distance loop uses the whole row; the concurrence loop takes only `key` and `word`
(`conc_key(m) = Symbol(:concurrence_, m.key)` — which works because the sweep's `:single`
expansion builds those names the same way). Drop a row to skip both its panels, in every
pair.

`run_pair(pr, D_grids, figs, pair_axes)` does one pair and **mutates** the three
containers it is handed, which is what keeps them live at top level after the script
finishes. `D_grids` and `figs` are both keyed `(pair_label, grid_key)`, where `grid_key`
is `:full` / `:adia` for the trace-distance grids and `:concurrence_full` /
`:concurrence_adia` for the signed ones; `figs` additionally holds `(label, :linecut)`,
which cannot clash because no `MODELS` key is `:linecut`.

```julia
figs[(:C_vs_A, :adia)]
figs[(:C_vs_A, :linecut)]
plot_map(pair_axes[:C_vs_A], D_grids[(:C_vs_A, :adia)])
```

`pair_axes` is keyed by label alone and holds the `d` each pair was drawn on — needed
because `assert_same_grid` only promises the two runs *of one pair* share a mesh. Do not
rename it `axes`: `Base.axes` is called inside `cross_tracedist` and `cross_concurrence`,
and a top-level `axes = Dict(…)` in `Main` makes those a name-resolution error.

Things that are load-bearing:

- **The two axes hold (β, x) fixed, not the raw parameters.** At matched (β, x) the two
  runs sit at different η whenever they differ in κ, and at a different `g_1` whenever
  they differ in **either** κ or `g_2`. So the figures never answer "at fixed drive".
  *Which* parameter a figure is about is per pair, and it is not always κ:
  - `:A_vs_B` is κ-only — it moves η and `g_1` together, and is about the κ asymmetry
  - `:C_vs_B` is `g_2`-only — η is untouched and only `g_1` moves, so it is about the
    coupling scale at fixed drive *and* fixed asymmetry
  - `:C_vs_A` differs in both and cannot separate the two causes; read it as one
    combined displacement, not as evidence about either parameter
- **`pair_tag` must stay unique per pair.** It builds the filename suffix, and `g_2` moves
  inside the per-run block only when the two runs disagree about it. An earlier version
  took `g_2` from the `a` run alone and appended it once — correct only while every pair
  shared it, and it gave `:C_vs_A` and `:C_vs_B` *the same tag*, silently overwriting five
  PDFs. The shared-`g_2` form is unchanged, so `:A_vs_B` still regenerates the files it
  first wrote rather than landing beside them.
- **Signed differences get `sym_clims` + `:balance`**, not viridis on a data range. This
  is issue 2's prescription, applied here rather than inherited — on an asymmetric range
  a diverging colormap's neutral midpoint lands on some arbitrary nonzero value and the
  sign becomes unreadable.
- **`clims` is computed from the `si_scale`d grid, not the raw one.** `plot_map` applies
  `clims` to what it is handed; limits taken before scaling are off by that power of ten
  and clamp everything to one end colour.
- **Every panel takes its own colour range**, across models *and* across pairs.
  Deliberate — see the finding below for why a shared range would hide the main result.
  No panel's colours are comparable to any other's; compare magnitudes via the printed
  summary, never by eye across panels.
- **Titles name the differing parameters once, with values per slot** — `(κ₁,κ₂): A (2,2)
  B (2.5,1.5)  both: g₂ = 1`, not the names repeated under each slot. Shared parameters
  are labelled `both:`; without that label they sit one `\quad` past the B values, the
  same gap that separates A from B, and read as belonging to B. Brackets go on
  unconditionally, including a pair differing in one parameter — a lone unbracketed title
  doesn't look tidier, just unlike its neighbours.
- **`si_scale` exists to fix a collision, not for looks.** GR pins `colorbar_title` at a
  fixed offset from the bar and the tick labels grow rightward into it, so `1.25E-14`
  overprints the title. `right_margin` cannot fix this — it translates bar, ticks and
  title together (checked at 8mm and 25mm, identical overlap). Narrower labels are the
  only lever; the exponent moves into the LaTeX title.
- **`assert_same_grid` compares `β_vals`/`x_vals` elementwise**, not by length or
  endpoints — two runs with the same `len` and endpoints but a different spacing rule
  would pass a length check and then silently compare mismatched physical points.
- **`cross_concurrence` gates on `status_grid` *and* on both operands being finite.** The
  second check is not redundant: `analyze` does not consult `status_grid` (issue 7), so
  an observable can be NaN at a point whose solve was `:ok`.
- All four maps are **immune to issue 1**. The full panels never touch `J_adia`. For
  `D_adia` the correction would be the same fixed local σ_z on qubit 1 in both runs —
  fixed because it is a convention of how `J_adia` is written, not a function of the
  parameters, so it is shared whichever parameter the pair disagrees about — and trace
  distance is invariant under a unitary applied to both arguments. `ΔC_adia` is safer
  still — concurrence is a local-unitary invariant, so each operand is individually
  unchanged and the argument doesn't need the two runs to share the rotation.

### Finding: the adiabatic steady state is a function of (β, x) alone

**Now proven analytically** — see "The adiabatic model is analytically independent of r
and ξ" above, which derives it by direct substitution in the (r, ξ) coordinates and
covers κ and `g_2` in one argument. Everything below is the earlier empirical route, kept
because it is the measurement the proof was checked against, and because the numbers are
quoted elsewhere in this file.

Measured across all three pairs, so it covers **every** way these runs differ — κ
alone, `g_2` alone, and both together:

| pair | differs in | `D_full` max | `D_adia` max | `ΔC_full` extremum |
|---|---|---|---|---|
| `:A_vs_B` | κ | 0.1608 | 2.0e-14 | −0.0601 @ β=0.215 |
| `:C_vs_A` | κ and g_2 | 0.0624 | 2.1e-14 | −0.0567 @ β=1 |
| `:C_vs_B` | g_2 | 0.1017 | 2.0e-14 | −0.0542 @ β=1 |

Every adiabatic column is at machine precision — structureless, signs scattered — and the
three stored maxima agree to 14 digits:

```
concurrence_adia   A 0.2360436749507588   B 0.2360436749507590   C 0.2360436749507577
concurrence_full   A 0.2866               B 0.2598               C 0.2614
```

So at fixed (β, x) the adiabatic steady state does not move at all. The elimination has
absorbed the raw parameters completely into the two dimensionless coordinates. The full
model has **not** — it shifts by up to 0.060 in concurrence, ~20% of peak.

**The `g_2` half of this is provable, not just measured**, and worth knowing because it
means no future run can disturb it. `H_adia ∝ g_1 g_2` and each `J_adia` entry is linear
in g (`jc_pump_disp_asy.jl:95,105-106`), so the dissipator is also quadratic and the
*entire* adiabatic Liouvillian is homogeneous of degree 2 in the couplings — given
`γ = γ_φ = 0`, there is no other term. At fixed β, `g_1 = g_2 sqrt(β κ_1/κ_2)` scales with
`g_2`, so both couplings scale together, `L → λ²L`, and the kernel — the steady state — is
exactly unchanged. **This argument dies if `γ` or `γ_φ` is ever turned on**, since those
rates do not scale with g; the κ half was only ever empirical.

Two things to read off the table:

- **`:A_vs_B` is the largest displacement**, so κ is the stronger lever on the full model,
  not `g_2`. (An earlier version of this section cited 0.060 as A-vs-B's shift; that is
  the *concurrence* figure, and comparing it against another pair's *trace distance* is
  what made `g_2` look dominant. Compare like with like.)
- **`ΔC_full` is a weak proxy for the displacement.** C and B have nearly equal peak
  concurrence (0.2614 vs 0.2598) while `D_full` reaches 0.10 — the states diverge
  substantially along directions concurrence cannot see. That is the `D`/`ΔC` disagreement
  the file's header anticipates, actually occurring.

All three `D_full` maxima sit at the same point, β = 4.64, x = 0.7225 — the top of the x
range, just right of β = 1. That is the approach to threshold, where the elimination's
`4η² − κ_1κ_2` denominators come closest to blowing up, so it is where the two models
should disagree most, and they do.

This is also why the panels are not on a shared colour scale — on a range set by
`ΔC_full`, the entire adiabatic panel would be uniform white, and *identically zero*
would be indistinguishable from *merely small*.

**Practical consequence: one point settles anything about the adiabatic model.** There is
no parameter regime to hunt for, which is what makes issue 1's check cheap.

## Where the adiabatic concurrence peaks — `adia_concurrence_max.jl`

```
C_adia is maximal at  x = 0.148550,  β = 1,  where  C = 0.2376153969
```

Measured with this repo's own `run_adia` + `calc_concurrence`, no closed form in the loop.
`η/η_th = sqrt(x) = 0.3854`, so the optimum sits at 38.5% of threshold — far below it, and
far below the `x* = 0.4839` point where adiabatic entanglement dies entirely.

The script is a coarse (β, x) table, then a **nested** golden section (best x at each β,
then best β) so what it finds is a 2D maximum and not the best x along a guessed β. CLI
args re-aim the window without editing the file: `adia_concurrence_max.jl 0.10 0.20 0.8
1.25 9` zooms on the peak. Runtime is seconds — that is `run_adia`'s doing, and it is why
this is a scan and not a sweep.

Four checks make the argmax evidence rather than an assertion, and they are the reason to
run this rather than read the number here:

- **(a) local peak test.** ±1% in each coordinate, all four neighbours lower. A converged
  optimizer sitting on a monotone function fails this; a maximum does not.
- **(b) (x, β)-only invariance.** Three unrelated `(κ_1, κ_2, g_2)` at the same `(x*, β*)`
  give `C` identical to `3e-16`. Without this the search would be over a different surface
  per realization and the answer would be meaningless. It is also an independent
  re-measurement of the finding above, on a third route.
- **(c) arm-swap symmetry** `C(x, β) = C(x, 1/β)` to `1e-15`, imposed nowhere in the code.
  This is what makes `β* = 1` structural rather than coincidental: any β-maximum must sit
  at the fixed point or come in a mirror pair.
- **(d) against the recorded sweep.** All 49 stored `concurrence_adia` values of run A
  reproduced with worst error `0.00e+00`, and run A's grid maximum
  (`0.2360436749508` at β = 1, x = 0.12875) lands just below `C*` — as it must, since the
  grid has no column at `x = 0.1486`. Skipped silently if the `.jld2` is gone.

**β = 1 is optimal at every x**, not just at the peak: `C_adia` depends on β only through
`b = β + 1/β ≥ 2` and decreases monotonically in `b` throughout the entangled region.
Asymmetry also *lowers* the optimal pump — `x_opt` = 0.1486 at β = 1, 0.107 at β = 2,
0.042 at β = 10. Both statements are `files_online/06_concurrence.md` §6; the closed form
there agrees with `run_adia` to 13 digits at both `x*` and run A's `x = 0.12875` grid
point, which is a stronger cross-check on `H_adia`/`J_adia` than anything in this repo
alone, since the two were derived independently.

For contrast, the **full** model peaks elsewhere and higher: `C = 0.29675` at `x = 0.291`,
`β = r = 1`, `ξ = 2.75` (`files_online/07`), i.e. 25% above the adiabatic maximum at
roughly twice the pump. Note `ξ_range = (0.1, 1.0)` in `config_r_xi.jl` excludes it.

## Parameter counting — `compare_param_sets.jl`, `scale_invariance_test.jl`

Two standalone scripts. Neither sweeps, neither writes a `.jld2`, neither is on the path
of anything else. They exist to establish that `ρ_full` is a function of exactly the four
dimensionless coordinates.

### `compare_param_sets.jl`

Edit the `SETS` list at the top — named `(g_1, g_2, κ_1, κ_2, η)` rows — and it prints
each set's four invariants, then three things: a pairwise `D(ρ_full)` matrix, a pairwise
`D(ρ_adia)` matrix, and a per-set `D(ρ_full, ρ_adia)` column.

Rows meant to sit at a *specific* dimensionless point must go through
`from_invariants(x, β, r, ξ)`, not hand-typed decimals. Typing 8-digit values for the
"same (x,β)" row put β at 1.00039 instead of 1, and that 4e-4 error alone lifted its
`D_adia` from ~1e-16 to 7e-6 — six orders above the floor, easily misread as a small real
effect.

Recorded output (`x = 0.485` slice, run A's parameters as the reference):

| pair | differs in | `D_adia` | `D_full` |
|---|---|---|---|
| A vs A×17.3 | overall scale only | 1.5e-16 | 2–5e-08 *(floor)* |
| **A vs "A same xβ"** | **r, ξ only** | **3.8e-16** | **5.9e-02** |
| A vs B | all four | 6.3e-02 | 6.4e-02 |

The middle row is the result: two parameter sets with identical `(x, β)` but `(r, ξ)`
moved from `(1,1)` to `(3,0.6)`. The adiabatic model cannot distinguish them at all;
the full model separates them by 6%. And **5.9e-02 is 92% of the 6.4e-02 you get from
moving all four**, so `r` and `ξ` are first-order levers, not corrections — you cannot
drop them and call the residue a perturbation.

**The two zero-rows are not the same kind of claim.** A vs A×17.3 is a *theorem*:
`L → λL` preserves `ker L`, so it was guaranteed before it ran, and confirms the
parametrisation code rather than the physics. Nothing forced the full model to move in
the middle row — that one is *measured*.

Three traps, all of which bit during development:

- **`D_full` has a noise floor of ~1e-8, `D_adia` does not.**
  `steadystate.iterative` seeds a shadow vector with `rand()` from the *unseeded global
  RNG* (`steadystate_iterative.jl:38`) and stops at `reltol = sqrt(eps) ≈ 1.5e-8`, so
  solving identical parameters twice already differs at that level. Across four runs the
  A vs A×17.3 entry came out 2.1, 3.7, 3.8, 4.7 ×10⁻⁸ while every other entry was stable
  to four significant figures. `ρ_adia` goes through a dense eigensolver with no RNG, so
  1.529e-16 is reproducible to the last digit. **Anything ≤ ~1e-7 in the `D_full` table
  is zero.** The two models are not "differently invariant"; they are measured with
  different instruments.
- **Never use λ = 2 or 0.5.** Halving is exact in binary, so a dyadic ratio can give a
  bit-identical Liouvillian and `D = 0.00e+00` — indistinguishable from a test that
  varied nothing at all. Use 17.3.
- **A null test needs a positive control.** `D ≈ 0` alone cannot distinguish "the state
  is scale-invariant" from "this code ignores its input". `scale_invariance_test.jl`
  section (b) supplies it: +1% in one invariant at fixed scale gives
  `D = 1.9e-3, 3.2e-4, 9.9e-5, 1.7e-4` for `x, β, r, ξ` — five orders above the null.

The elimination error, `D(ρ_full, ρ_adia)` per set, is the scale the other two matrices
should be read against. At the recorded points it is **0.10–0.15**, i.e. *larger* than
any of the between-set distances. So those tables compare two models that are already far
apart — `ξ ≈ 1` is nowhere near the `ξ → 0` limit the elimination is derived in.

### `scale_invariance_test.jl`

The randomized version: draws several base parameter sets, rescales each by non-dyadic λ,
reports `D_fibre` against `D_self`. Sections `a` (fibre), `b` (positive control),
`c` (that `estimate_truncation` is scale-invariant too, so both members of a pair are
solved in the same Hilbert space).

Two scope notes. It compares the **reduced** 2-qubit state, because `run_sim` traces out
the cavities before returning and this file deliberately does not modify production code
to do otherwise — the joint-state comparison would be strictly stronger, since two states
can agree after tracing while the joint states differ. Consequently its section (b)
numbers are **not** comparable to the joint-state values in
`files_online/parameter_counting_test.py`. Section (c) is also weaker than it looks: every
draw returned `(14, 14)`, which is `N_FLOOR + K_FLOOR_REGIME`, so a function returning a
constant would pass it too.

## `files_online/` — parallel analytical project

A separate, self-contained project on the *same* physics: same Hamiltonian, same TMS
drive. Where this repo is a numerical sweep, that one is closed-form — exact steady
state, exact concurrence, parameter counting — with a strict `[NUM]`/`[HAND]`/`[OPEN]`
evidence legend. Seven notes plus six NumPy scripts. **It has its own `claude.md`; that
file governs that project, not this one.** Its notes reference `notes/` and `scripts/`
subdirectories but everything is flat in the folder.

**⚠ `ξ` means the reciprocal there.** This is the highest-risk thing about the two
projects coexisting:

| quantity | this repo | `files_online` |
|---|---|---|
| `4g_1g_2/(κ_1κ_2)` | **ξ** | **μ**, also written **𝒜** |
| `κ_1κ_2/(4g_1g_2)` | `1/ξ` | **ξ** |
| `κ_1/κ_2` | `r` | `α`, also written `r` |

So adiabatic elimination is `ξ → 0` in this repo's convention and `ξ ≫ 1` in theirs.
`config_r_xi.jl`, `run_r_xi.jl`, the axis labels and the `.jld2` filenames all use the
first sense. Reconcile before either becomes a figure.

Results there that bear directly on this repo:

- **The `x* = 0.48388826` entanglement threshold** (root of `x³ − 2x² + 9x − 4`) for the
  adiabatic model at β = 1, and the `b = β + 1/β < b_max(x)` boundary at general β. This
  is the closed form behind the `C_adia = 0` dead zone documented under the (r, ξ) sweep —
  the empirical scan there just rediscovers it.
- **The full model's concurrence maximum**: `C = 0.29675` at `x = 0.291`, `β = r = 1`,
  `𝒜 = 2.75` (i.e. **this repo's ξ = 2.75**), 25% above the adiabatic maximum of 0.2376.
  Note the current `ξ_range = (0.1, 1.0)` **excludes it**.
- **The exact arm-swap symmetry** `C(x,β,r,𝒜) = C(x,1/β,1/r,𝒜)`, which is what makes
  the β = 1 slice mirror-symmetric in `r ↔ 1/r`.
- **`D ∝ μ` with α entering at first order**, i.e.
  `ρ_full = ρ_adia(x,β) + ξ·f(x,β,r) + O(ξ²)`. The (r, ξ) sweep is a direct 2D map of
  that correction term, where they have only slopes at two α values — but see the
  recorded run: at `ξ ≥ 0.1` the measured scaling is already sub-linear, so the map is
  outside the regime where that expansion is the whole story.
- Their basis convention is **qubit 1 = second tensor factor**; `jc_pump_disp_asy.jl`
  embeds qubit 1 at position 3 and qubit 2 at position 4, the *opposite* ordering. Both
  agree that concurrence, purity and trace distance are invariant under the exchange, so
  nothing currently computed here is affected — but any per-qubit population would be.

## Conventions worth not breaking

- **Grids are `[N_β, N_x]`** — β is the *row* index, x the *column*. Every report,
  observable matrix and plot assumes this.
- **NaN means "no data"**, everywhere. `nanmat()` prefills the state grids; observable
  grids prefill with `NaN`. A failed point becomes a gap in the heatmap instead of
  killing the run. Don't replace with zeros.
- **`res.params` is the truth for plot axes**, not `config.jl` — the config may have been
  edited after the run that produced a given `.jld2`.
- **The report block in `run_sweep` is wrapped in `try/catch` on purpose.** It sits before
  the `return`, and a bug in a summary line once discarded a completed sweep. Keep any new
  reporting inside that block.
- `.gitignore` is deny-by-default (`*` then `!*.jl`). Results, figures and `Manifest.toml`
  are not tracked. A new non-`.jl` file needs an explicit `!` line. **Git does not descend
  into an ignored directory**, so the `!*.jl` exception never reaches subdirectories —
  `files_online/` is entirely untracked and is *not* on GitHub, notes and `.py` scripts
  included. Only 15 files are tracked: the top-level `.jl`, `CLAUDE.md`, `Project.toml`,
  `.gitignore`.
- **`origin` is `https://github.com/Yash-27/JC_TMS_ADIA.git`** — this working directory
  *is* that repo, so there is nothing to clone into a subfolder. **GitHub's tree is 14
  files**: `.gitignore`, `CLAUDE.md`, `Project.toml`, 11 `.jl`. No `.py`, no
  `files_online/`, no `.jld2` — and `origin/main` is an ancestor of local HEAD, so
  everything there is already here. **The `.py` files and the second `claude.md` are in
  `files_online/`, which is local-only**; if you expect to find them on GitHub, they are
  not there yet.

## Session start — `.claude/session-status.sh`

A `SessionStart` hook in `.claude/settings.json` runs it before anything else, and its
stdout becomes Claude's opening context. Three sections:

1. **GIT** — `git fetch`, then ahead/behind vs upstream (falling back to `origin/main`),
   a specific flag when the working copy's `CLAUDE.md` has drifted from the pushed one,
   and the uncommitted-file list capped at 20.
2. **CONTENT** — `shasum` inventory of `files_online/`, diffed against
   `.claude/.content-manifest` from the previous session, reported as ADDED / REMOVED /
   MODIFIED. This section exists because **git cannot see that directory at all** — it is
   ignored, and git does not descend into ignored directories.
3. **READ LIST** — names the changed files and instructs Claude to open them, plus the
   reminder that `files_online/claude.md` governs that project and inverts `ξ`.

**It reports changes rather than pasting content, deliberately.** `files_online/` is
~138 KB (~35k tokens) and `CLAUDE.md` is another 57 KB the harness already loads every
session; dumping both unconditionally would spend a third of a context window re-reading
unchanged text. Steady-state cost is ~220 tokens. After an edit, Claude is pointed
straight at what moved and reads it in full.

**`fetch` only moves remote-tracking refs** — never merges, rebases, or touches the
working tree, so the hook cannot lose an edit; pulling is left to you. It never exits
non-zero, so no network / no remote / no `files_online/` still starts a session. The only
file it writes is its own manifest. Verified against all four paths: first run, unchanged,
modified, added+removed. `.claude/` is itself ignored by the deny-by-default rule, so none
of this is tracked or pushed.

## Scheduler in `run_sweep`

Two passes, and the reasons are non-obvious:

1. **Pass 1 runs `estimate_truncation` for every point first.** It's cheap relative to the
   steady-state solve, and it gives the Hilbert dimension of every solve up front so total
   work is known rather than discovered.
2. **Pass 2 solves on two atomic work queues**, big (`dim ≥ BIG_DIM = 2000`, capped at
   `MAX_BIG = 2` preferring workers) and small, both sorted expensive-first. `@threads`
   splits into contiguous chunks and the expensive points are contiguous (the whole last
   x column), so plain `@threads` hands one thread all of them.

Cost model is `dim^COST_EXPONENT × 1/(1 − sqrt(x))`. The second factor is *derived*, not
fitted: the Liouvillian gap is `(κ/2)(1 − sqrt(x))` and Krylov iteration count goes as
1/gap. It dominates dimension in practice. The end-of-run report prints the measured
`time ~ dim^slope` so `COST_EXPONENT` can be checked against reality.

## Known issues

Ordered roughly by how much they'd change a published number — but **the numbers are
stable identifiers, not a live ranking**. Both this file and `compare_runs.jl` refer to
"issue 1", "issue 2", "issue 7" in prose, so a demoted or fixed item keeps its number and
gets a note rather than being renumbered. Issue 4 is the current example: it was filed as
a correctness bug, is now known not to be one, and stays at 4.

### 1. Possible relative-sign mismatch between the full and adiabatic models — UNVERIFIED

Eliminating the cavities from `H_full` in `jc_pump_disp_asy.jl` gives steady-state cavity
amplitudes

```
a_1 ∝ (κ_2/2) g_1 σ⁻_1  +  η g_2 σ⁺_2
a_2 ∝  η g_1 σ⁺_1       + (κ_1/2) g_2 σ⁻_2
```

i.e. a **plus** between the two terms. `J_adia` in `run_sim` uses
`(g_1 ϵ_1 σ⁺_1 − g_2 σ⁻_2)` — a **minus**. The flip is consistent across both `J_adia`
entries *and* `H_adia` (which is odd in η), so the adiabatic model is internally
self-consistent; it just appears to correspond to `η → −η` relative to `H_full`.

`η → −η` in the full model is a local σ_z rotation on qubit 1, so **if this is real**:

- `concurrence_full`, `concurrence_adia`, `purity` are **unaffected** (local-unitary invariants)
- `tracedist` is **wrong** — it would be comparing `ρ_full` against `U ρ_adia U†`

This was derived by hand and **not confirmed numerically** — the check was interrupted
before it ran. Before trusting or "fixing" anything here, run the comparison: solve one
point both ways and see whether `tracedist` drops when the relative sign in both `J_adia`
entries is flipped to `+`. If it does not, the derivation above is wrong and this section
should be deleted.

The check is cheaper than it looks, and **one point settles it for the whole sweep** —
see the (β, x)-only finding above. Since the adiabatic steady state does not depend on
κ_1, κ_2 or `g_2` separately, there is no need to scan for a parameter regime where the
discrepancy shows up; any single (β, x) in any of the three runs is representative.

Note also that **no** `compare_runs.jl` figure can detect this, by construction — that is
the point of the immunity argument in that section, and it covers the line cut too, which
is just a cut through the immune `D_full` grid. It has to be tested against the in-file
`:tracedist`.

### 2. `concurrence_diff` is drawn on a colour scale that cannot show it

`CLIM_SOURCE` sends `:concurrence_diff` to `:concurrence_full`'s range. A concurrence is
non-negative, so those limits are `[0, C_max]`. The difference is *signed* and typically
much smaller in magnitude. Result:

- every negative value clamps to the bottom viridis colour, indistinguishable from zero
- all the real structure is compressed into a few percent of the ramp

The comment above `resolve_clims` says the warning exists precisely so "a negative
concurrence difference" isn't hidden — and then the configured default hides it. The
warning fires on essentially every run, which trains you to ignore it.

A signed quantity wants a symmetric range `(-m, m)` with `m = maximum(abs, finite)` and a
diverging colormap (`:balance`, `:RdBu`), not viridis on a non-negative source. Same
applies to `:purity_diff`, which currently falls back to its own range — better, but still
viridis on a signed quantity. There is no per-key colormap table to mirror `CLIM_SOURCE`;
adding one is the clean fix.

`compare_runs.jl` already does it the right way for its own signed panels (`sym_clims` +
`cmap = :balance`, both verified to render correctly through GR). That's the pattern to
lift into `plotting_functions.jl` when this gets fixed properly.

Until then, override at the call site:

```julia
m = maximum(abs, filter(isfinite, vec(d.outs[:concurrence_diff])))
plot_map(d, :concurrence_diff; clims = (-m, m), cmap = :balance)
```

### 3. `truncation.jl` constants are marked as placeholders

```julia
const K_FLOOR_REGIME   = 7    # "replace with your actual printed values"
const K_PHYSICS_REGIME = 17
const K_CLASSIFY_REF   = 22.0
```

The source comment says these should be replaced with calibrated values from the k-sweep
diagnostic. Unclear whether that was ever done, and **the diagnostic script is not in the
repo**, so there is currently no way to re-derive them. Every truncation in every sweep
depends on them.

**They also set the dimension clusters the scheduler is tuned against**, which is not
obvious from either file:

| regime | k | `N ≥` | `dim = (N_1+1)(N_2+1)·4` |
|---|---|---|---|
| floor-dominated | `K_FLOOR_REGIME = 7` | `N_FLOOR + 7 = 14` | 15·15·4 = **900** |
| physics-dominated | `K_PHYSICS_REGIME = 17` | ≈ 25 | 26·26·4 = **2704** |

Those are exactly the "900 cluster" and "2704+ cluster" cited in
`run_jc_pump_disp_asy.jl` to justify `BIG_DIM = 2000` sitting between them. So
**recalibrating these constants invalidates `BIG_DIM`**, and nothing in either file links
the two. Re-derive `BIG_DIM` from the new dimension histogram if these ever change.

### 4. `k_logic`'s stated justification is false — but the conclusion holds by accident

The comment claims that because both branches increase with k, whichever wins at
`K_CLASSIFY_REF = 22` also wins at any smaller k. **That reasoning is a non-sequitur** —
monotonicity of each branch says nothing about their *ordering*, which depends on the
slopes. Write

```
f(k) = physics − floor = (mu − N_FLOOR) + k(sqrt(var) − 1)
```

so the ordering flips iff `f` changes sign between k = 22 and the k actually used.

**The dangerous direction cannot occur.** Under-truncation needs `sqrt(var) < 1` (so `f`
decreases) together with `f(7) > 0` — classified floor-dominated at k=22, so it gets
k=7, but the physics branch actually wins there and the point receives 7σ of margin
where 17σ was intended. That requires

```
f(7) > 0  ⟹  mu + 7 sqrt(var) > 14  ⟹  mu > 7    (since sqrt(var) < 1)
```

But `photon_number_moments` returns `Var(n) = n(n+1) + |s|²`, so `mu > 7 ⟹ var > 56 ⟹
sqrt(var) > 7.5`, contradicting `sqrt(var) < 1`. The regime is unreachable.

The opposite flip (`sqrt(var) > 1`, classified physics at k=22 but floor-dominated at
k=17) *can* happen, and there the `max()` inside `N_from_k` picks the larger floor term —
over-truncation. Wasteful, never wrong.

**So this is a documentation defect, not a correctness bug.** It was previously recorded
here as something that "can under-truncate"; it cannot. What makes it safe is the
variance relation `Var ≥ n(n+1)` — an invariant that lives in a *different function*
(`photon_number_moments`) and is nowhere marked as load-bearing. Replace the Gaussian
closure with anything permitting small variance at large mean and `k_logic` starts
under-truncating silently, with no test and no warning to catch it.

If you touch either function: the coupling is the whole safety argument. Keep
`Var ≥ n(n+1)`, or replace `k_logic`'s classify-then-evaluate scheme with something that
doesn't depend on it.

### 5. Memory warning underestimates by at least 6–8×, probably much more

```julia
mem_max = dmax^2 * 16 / 2^20   # one dense ρ
```

`steadystate.iterative` holds a Krylov subspace — several vectors of `dim²` complex
numbers, not one. The 4 GiB warning threshold fires far later than it should.

The "6–8×" is likely optimistic: restarted GMRES keeps its whole subspace, typically
20–30 vectors. At dim 4096 one vector is 268 MiB (the figure the scheduler comment
quotes), so a single large solve could plausibly sit in the several-GiB range on its own —
before `MAX_BIG` fails to stop several running at once (issue 6). The two compound.

### 6. `MAX_BIG` is not actually enforced

Confirmed in `take_job`. The `prefer_big = false` branch:

```julia
k = atomic_add!(next_small, 1) + 1
k <= length(small_idx) && return small_idx[k]
k = atomic_add!(next_big, 1) + 1          # <-- falls through
k <= length(big_idx) && return big_idx[k]
```

`n_big_workers` controls **preference at spawn time, not concurrency**. Once `small_idx`
drains, every thread takes big jobs — exactly the memory contention the split exists to
prevent. It only bites in the tail, but the measured 2h06m → 3h33m CPU-time blowup cited
in the source comment is this failure mode.

The `prefer_big = true` fall-through to the small queue is fine and should stay — that
one just keeps big workers busy after the big queue empties. **The fix is a counting
semaphore on in-flight big solves, not deleting the fall-through**, which would idle
threads in the tail and trade one problem for another.

### 8. `run_r_xi.jl` duplicates sweep 1's scheduler, and the two can drift

`run_sweep_rxi` is a near-copy of `run_sweep`: same two-pass structure, same atomic
big/small work queues, same `take_job`, same ETA fit, same `try/catch` report block —
roughly 200 lines. Only three things genuinely differ: the coordinate inversion, the cost
model's gap factor, and the printf labels.

It was written that way deliberately, to keep zero risk to a working sweep, and the
observables were shared instead (`observables.jl`) because those *could* be lifted
without touching the scheduler. But the consequence is real: **a fix to the scheduler in
one file does not reach the other.** Issue 6 (`MAX_BIG` not enforced) is now present in
both copies, and any future fix has to be applied twice.

The clean resolution is a `sweep_core.jl` taking `point_params(i,j) -> (κ_1,κ_2,g_1,g_2,η)`
and `cost_factor(i,j)` as closures, with both drivers calling it. That is a real
refactor of a file that currently works, so it should be done deliberately rather than
opportunistically.

### 7. Smaller things

- `run_sim`'s docstring says `rho_ss_full` is returned **untraced**. It isn't — the last
  line is `ptrace(rho_ss_full, (1,2))` and the caller asserts 4×4. Stale doc.
- `analyze` doesn't consult `res.status_grid`, so failed points feed NaN matrices into
  `eigvals`, which throws and produces one warning per point per observable. Harmless,
  noisy, and it hides real observable failures in the noise.
- `OUTFILE` keys on `κ_1`, `κ_2`, `g_2` only. Two runs differing in `len`, `x_range` or
  `β_range` silently overwrite each other.
- `time_grid[i,j]` is written even when the solve failed, so failed-and-fast points can
  appear in the "slowest points" table.
- `Γ_sqrt_1`/`Γ_sqrt_2` are purely imaginary; the `2im` is a global phase on a jump
  operator and drops out of the dissipator. Only the relative sign inside the parenthesis
  matters (see issue 1).
- With `γ = γ_phi = 0` the last four `J_full` entries are exact zeros but still built as
  full-dimension sparse operators. Minor waste.
- `steadystate.iterative` assumes a *unique* steady state. With `γ = γ_phi = 0` that is an
  unchecked assumption; a degenerate Liouvillian kernel would silently return an arbitrary
  member of it. This applies to the **full** solve only — the adiabatic model goes through
  `steadystate.eigenvector` (`jc_pump_disp_asy.jl:120`), which has the same exposure by a
  different route: it would pick one kernel eigenvector, also without warning.
- `Dates` is imported in `run_jc_pump_disp_asy.jl` and unused.
- **`config_r_xi.jl`'s comments are stale after `β_fixed` was changed to 0.5.** The
  inline comment still reads "symmetric arms: `g_1/g_2 = sqrt(r)`, r is the only
  asymmetry", which describes β = 1, not 0.5. The longer note below `PARAMS_RXI` still
  claims "combined with `β_fixed = 1.0`, which is likewise exactly on the existing β
  grid" — β = 0.5 is **not** on that grid (`0.01, 0.0464, 0.215, 1, 4.64, 21.5, 100`),
  so the run-A cross-check that note promises does not currently exist. Fix the comments
  or move β back to 1.0; do not trust the note as written.
- `print("\r\x1B[K")` in the solve-failure catch is a leftover from a removed progress bar.
- `load_results` pulls `status_grid` out of the `.jld2` and **nothing on the single-run
  path ever reads it**. It's the one thing that would let a figure distinguish a failed
  solve from a failed observable — right now both render as an absent marker.
  `compare_runs.jl` does read it (via its own `load_states`, not `load_results`), so the
  pattern to copy exists.
- Two headers name the wrong file: `plotting_functions.jl` calls itself `plotting_lib.jl`;
  `run_plot.jl` calls itself `plot_try.jl` and says the machinery lives in
  `plotting_lib.jl`. (`plot_try.jl` was a byte-identical copy of `plotting_functions.jl`
  plus a stale driver, and has been deleted.) The "included by `plot_try.jl`" line in
  `plotting_functions.jl` has since been corrected to name the three real drivers; the
  self-naming line has not.
- `plot_map`'s docstring promises "any extra keyword is passed straight to Plots", but
  `xscale`, `xlabel`, `ylabel`, `c`, `size`, `legend`, `colorbar`, `markerstrokewidth`,
  `markerstrokecolor` are set explicitly in the `scatter` call *and* reachable through
  `kwargs...`. Whether the splat overrides or errors on the duplicate is **still
  unverified** — check before relying on it. Keywords *not* in that list do pass through
  cleanly: `compare_runs.jl` sends `right_margin` and `xticks` that way on every run.
  Note `cmap` and `clims` are named parameters of `plot_map`, so they are set through the
  signature and never reach the splat.
- `plot_map(d, M)` with a raw matrix gives `ttl = ""` and `cbar = ""`, so the advertised
  "map `N_1_grid` with no extra code" path produces an untitled, unlabelled PDF.
- `markerstrokewidth = 0.5` is set twice, in `apply_theme!` and again in `scatter`.
- PDF output goes through GR with `colorbar_title` as a `latexstring`. GR's colorbar-title
  support is thin and LaTeX there is a known weak spot — eyeball one output rather than
  assuming it rendered.

## Physics notes for anyone editing the equations

- `truncation.jl` solves 32 real unknowns (4 real + 14 complex second-order moments) by
  homotopy continuation: η is ramped 0 → η_target in 60 steps, each warm-starting from the
  last. At η = 0 the exact solution is the trivial vacuum/ground state, which is the seed.
- Its moment equations appear to use the opposite η sign convention from
  `jc_pump_disp_asy.jl` (compare the `+2η·Re(s12)` in `rn1` against the full Hamiltonian).
  This is **harmless** — the truncation only consumes `n_1`, `n_2`, `|s11|`, `|s22|`, all
  invariant under `η → −η`.
- `photon_number_moments` is a Gaussian/Wick closure at 4th order:
  `Var(n) = n(n+1) + |s|²`. It assumes `⟨a⟩ = 0`, which holds because the moment equations
  carry no first-moment variables. It is an approximation, consistent with but distinct
  from the 2nd-order cumulant truncation used for the moments themselves.
- Concurrence is Wootters': `R = ρ ρ̃` is not Hermitian, so `eigvals` returns complex
  numbers; they are analytically real and non-negative and `abs()` absorbs the numerical
  imaginary part. The leading `ρ / tr(ρ)` guards against a trace defect from the iterative
  solver leaking into C.
- `YY` (σ_y ⊗ σ_y) is real, so it's stored as a `Float64` matrix — no complex arithmetic
  needed there.
