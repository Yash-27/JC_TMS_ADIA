# CLAUDE.md

Julia code for a cavity-QED reservoir-engineering study: two qubits, each dispersively
coupled to its own cavity, with the two cavities linked by a two-mode-squeezing (TMS)
drive. The question the code answers is how well the **adiabatically-eliminated 2-qubit
model** reproduces the **full 4-body model** across a 2D parameter sweep.

## Running

Three sweeps, each with its own config and plot driver, sharing everything below them.

```bash
# sweep 1 -- the (β, x) plane at fixed κ_1, κ_2, g_2
julia --project=. -t auto run_jc_pump_disp_asy.jl
julia --project=. run_plot.jl
julia --project=. compare_runs.jl        # two (β, x) runs against each other

# sweep 2 -- the (r, ξ) plane at fixed x, β
julia --project=. -t auto run_r_xi.jl
julia --project=. run_plot_r_xi.jl

# sweep 3 -- the (x, ξ) plane at fixed β = r = 1   (~1h45m: see the cost note)
julia --project=. -t auto run_x_xi.jl
julia --project=. run_plot_x_xi.jl

# one-off parameter-set tools, no sweep
julia --project=. -t auto compare_param_sets.jl


# where the adiabatic concurrence peaks -- adiabatic model only, seconds, no threads
julia --project=. adia_concurrence_max.jl [x_lo x_hi β_lo β_hi len]

# the entanglement boundary u_c(x), measured and drawn against the closed form
julia --project=. adia_boundary.jl [x_lo x_hi N_x N_u]     # ~10 s, no threads

# how C_adia varies INSIDE that region -- ∂C/∂u and ∂C/∂x, mapped
julia --project=. adia_derivatives.jl [x_lo x_hi N_x N_u]  # ~40 s, no threads
```

All three of `run_jc_pump_disp_asy.jl`, `run_r_xi.jl` and `run_x_xi.jl` have a top-level
`begin` block at the bottom that executes the whole sweep on include. There is no
`main()` guard — `include`ing any one of them *runs* it. That is why the observables live
in `observables.jl` rather than inside a runner: no runner can include another to borrow
code.

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
    ┌────────────┴───────────────┬─────────────────────────┐
    │                            │                         │
config_jc_pump_disp_asy.jl  config_r_xi.jl          config_x_xi.jl
run_jc_pump_disp_asy.jl     run_r_xi.jl             run_x_xi.jl
  sweeps (β, x)               sweeps (r, ξ)           sweeps (x, ξ)
  at fixed κ, g_2             at fixed x, β           at fixed β = r = 1
  -> results_k1_*_k2_*_g2_*   -> results_rxi_*        -> results_xxi_*
    │                            │                         │
    ├─ run_plot.jl               └─ run_plot_r_xi.jl        └─ run_plot_x_xi.jl
    └─ compare_runs.jl  (two .jld2)
                 │
        plotting_functions.jl   plot_map / save_fig / axis_spec / flatten_grid,
                                included by all six plot drivers

    STANDALONE -- no sweep, no .jld2
        compare_param_sets.jl   pairwise D between named parameter sets
        adia_concurrence_max.jl   argmax of C_adia over (x, β), via run_adia
        adia_boundary.jl        u_c(x) measured vs closed form; own figure
        adia_derivatives.jl     ∂C/∂u and ∂C/∂x maps + the x_opt(u) ridge

    PROSE -- no code, own sections below
        files_online/           closed-form parallel project; INVERTS ξ
        file_concurrence/       closed-form C_adia: x_c, u_c(x), C(2,x); ξ-free
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
| `x` | `4η² / (κ_1 κ_2)` | drive relative to the linear OPO threshold | sweep 1, linear; **sweep 3, linear** |
| `β` | `(g_1² κ_2) / (g_2² κ_1)` | Purcell-rate asymmetry between the arms | sweep 1, log10 |
| `r` | `κ_1 / κ_2` | cavity-decay asymmetry | sweep 2, log10 |
| `ξ` | `4 g_1 g_2 / (κ_1 κ_2)` | coupling strength relative to decay | sweep 2, linear; **sweep 3, log10** |

Sweeps 1 and 2 are exact complements: together the four coordinates are a complete set,
so sweep 2 covers precisely the directions sweep 1 holds fixed, with no overlap and no
gaps. **`ξ` is the adiabatic elimination's own small parameter** — the elimination is
the `ξ → 0` limit — which is why it is the physically informative axis.

**Sweep 3 deliberately overlaps both**, and that is its point: it is the `(x, ξ)` plane at
`β = r = 1`, i.e. the two coordinates that actually control the elimination, with both
asymmetries switched off. Sweep 1 cannot deliver it (ξ is confounded with β there — see
below) and sweep 2 pins x. The three planes intersect in lines, not regions, and those
lines are used as regression tests: sweep 3's `ξ = 1` row *is* run A's `β = 1` row.

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
`b = 2` (β = 1) threshold is `x_c = 0.48388826`, the root of `x³ − 2x² + 9x − 4`. At
β = 0.5 (`b = 2.5`) the boundary has already been crossed by `x = 0.3663`. **Both halves
are proven** in `file_concurrence/concurrence_adia.md` §3.3–3.5 — `b` is that file's `u`,
and its `u_c(x)` is `b_max(x)` in closed form, so the boundary can be drawn rather than
scanned for. Measured on run A's grid at `x = 0.3663`:

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

## The (x, ξ) sweep — `run_x_xi.jl`

Sweeps `(x, ξ)` at fixed `β = r = 1`, same `κ_geo = sqrt(κ_1 κ_2) = 2` pinning as sweep 2.
At `r = β = 1` the inversion collapses to one rate per axis and nothing else moving:

```
κ_1 = κ_2 = κ_geo = 2                   <- constant across the whole grid
η   = (κ_geo/2) sqrt(x)                 <- the x axis, and only it
g_1 = g_2 = (κ_geo/2) sqrt(ξ)           <- the ξ axis, and only it
```

**Why this plane, given the other two exist.** Sweep 1 confounds β with ξ (`ξ ∝ sqrt(β)`
at fixed κ and `g_2`) and structurally cannot separate them; sweep 2 pins `x = 0.3663`
and stops at `ξ = 1`. Neither shows the error in the *(drive, coupling)* plane at
symmetric arms — which is where the elimination's validity actually lives, since `ξ` is
its own small parameter and `x` is what its `4η² − κ_1κ_2` denominators blow up on. Two
independently known landmarks sit inside this box and on no existing grid: the **full
model's optimum** (`C = 0.29676` at `x = 0.291, ξ ≈ 2.79, β = r = 1`) and the **adiabatic
entanglement threshold** `x_c = 0.48389` (proven in `file_concurrence` §3.5).

**Grid layout: ROW = ξ (horizontal, log10), COLUMN = x (vertical, linear)** — grids are
`[N_ξ, N_x]`. This keeps the repo-wide "row coordinate is the horizontal/log axis"
convention, so `flatten_grid` and `plot_map` are reused unchanged.

`ξ_range = (0.2, 5.0)` rather than `(0.1, 5.0)`, and the reason is not aesthetic:
**reciprocal endpoints put `ξ = 1` exactly on the grid** at len 7 (`0.2, 0.3415, 0.5833,
1.0, 1.7145, 2.9375, 5.0`). `(0.1, 10.0)` would work too; `(0.1, 5.0)` would not.

### Two self-tests, both free, both consequences rather than assertions

**1. `ρ_adia` is a function of x alone here**, so every `_adia` panel must be constant
*down each column* (across ξ) and must vary *across* columns. Same argument as sweep 2's
`check_adia_flat` in the complementary direction: `g ∝ sqrt(ξ)` makes every term of
`L_adia` carry one factor of ξ, so `L_adia → ξ L_adia` and the kernel is unchanged; β and
r are fixed, leaving x. `check_adia_x_only` asserts both halves — the second half matters,
since a panel flat in *both* directions would mean x is not reaching the model at all, and
a column-only check cannot see that. It then re-derives each column value from a direct
`run_adia` call, which costs microseconds and checks the sweep's inversion against the
model rather than against itself.

**2. The `ξ = 1` row is already on disk.** With `β = r = 1` and `κ_geo = 2`, `ξ = 1` gives
`κ_1 = κ_2 = 2`, `g_1 = g_2 = 1`, `η = sqrt(x)` — run A's `β = 1` row, and all seven x
values are on run A's own x grid. So `check_against_run_A` compares row 4 against
`results_k1_2.0_k2_2.0_g2_1.0.jld2` directly: adiabatic to the last bit, full to ~1.5e-8
(`steadystate.iterative`'s unseeded-RNG floor). **This is the load-bearing test** — it
validates the coordinate inversion, the truncation and the driver in one shot, including at
the expensive dim-4096 point. Verified on a 3×3 pilot grid: `concurrence_adia` `0.0e+00`,
`concurrence_full` `3.3e-08`, `tracedist` `1.7e-08`.

**Mind the two grids' index orders when reading that check.** This sweep is `[ξ, x]`, run A
is `[β, x]`, so the reference entry is `ref[β_index, x_index]`. The first version of the
check had it transposed and reported mismatches of order the observable itself — which is
what the check is for, but the trap is worth naming: both grids are square at len 7, so no
assertion catches it.

### Cost — measured, and the reason this driver enforces `MAX_BIG`

**No `results_xxi_*.jld2` exists yet** — the driver and its self-tests are verified (a 3x3
pilot grid, plus the checks below) but the full 7x7 run has not been completed and stored.

Pass 1 on the full grid, measured: **42 points at dim 900, 7 at dim 4096**, the seven being
the whole `x = 0.7225` column. The truncation there is set by η, not g, so it is dim 4096 for
**every ξ from 0.1 to 5** — raising ξ to 5 costs nothing in Hilbert dimension. This is the
truncation cliff documented further down, and it is the single cost fact governing the run.

Measured single-threaded, `-t 1`, one point each:

| point | dim | time | peak RSS |
|---|---|---|---|
| `x = 0.6038, ξ = 5` | 900 | 16.4 s | — |
| `x = 0.7225, ξ = 1.36` | 4096 | **849 s** | **5.04 GiB** |

**5 GiB for one solve is why this driver has a real counting semaphore** (`Base.Semaphore`,
`MAX_BIG = 1`) instead of sweep 1 and 2's spawn-time preference. Two concurrent big solves
would be 10 GiB before the small queue's contribution, on a 16 GiB machine — known issue 6
stops being a tail inefficiency and becomes an OOM. The queue preference is kept (it starts
the expensive points early, which is what keeps the tail short); the semaphore is what
bounds concurrency. **The fix lives in this driver only** — see issue 8.

Budget: 42 × ~16 s ≈ 11 min of CPU, a couple of minutes wall on 8 threads, plus 7 × 14 min
serialized ≈ **1 h 45 m total, essentially all of it the last x column**. Dropping
`x_range` to `(0.01, 0.60375)` keeps the whole grid at dim 900 and the run at minutes; that
is the trade if you want a fast turnaround.

Two operational things this driver adds and the other two lack: `flush(stdout)` after every
progress line, and a `PARTIAL_XXI` checkpoint of both state grids written after **every**
completed point (deleted on success). CLAUDE.md records 45 minutes of solves lost to
buffered stdout under a shell `timeout`; with single points in the 14-minute range, an
all-or-nothing write at the end is the same trap. Recovering from a checkpoint means loading
it and calling `analyze` — the states are stored raw, no observables.

### Expected structure

- `concurrence_adia` is **exactly 0 for x ≥ 0.485** (three of seven columns) — the
  adiabatic entanglement threshold at β = 1. There `concurrence_diff` becomes an exact
  duplicate of `concurrence_full`, as in sweep 2's recorded run. Analytic boundary, not a
  bug; `tracedist` stays fully informative.
- `concurrence_full` should peak near `x ≈ 0.29, ξ ≈ 2.9` at `C ≈ 0.297`.
- `tracedist` should grow monotonically in ξ at every x, **sub-linearly** — `D ∝ ξ` holds
  only below `ξ ≈ 0.01` and this grid is deliberately in the saturating region. Do not read
  a slope off it and extrapolate to `ξ → 0`.

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

`run_plot.jl` is the single-run entry point and `compare_runs.jl` the two-run one. Five
files now `include` `plotting_functions.jl` — those two, `run_plot_r_xi.jl`,
`run_plot_x_xi.jl`, and `adia_boundary.jl`, which is not a sweep driver at all and uses
only `apply_theme!`, `flatten_grid` and `save_fig`. Nothing else should.
`run_plot.jl` runs at top level so `d` and `figs` stay live in the session:

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

**Which** coordinates those are is per-file, not hardcoded, since there are now three
sweeps. `axis_spec(d)` returns `d.axes` when the loaded run carries one — `load_results_rxi`
attaches `(xvals, yvals, xlab, ylab, xsc, ysc)` for the (r, ξ) plane, `load_results_xxi`
for the (x, ξ) plane — and otherwise falls back to the original β/x behaviour. The fallback
is load-bearing: the three stored (β, x) `.jld2` files predate the field, and
`compare_runs.jl` hands `plot_map` a bare `(β_vals = …, x_vals = …)` tuple it builds inline.
Both still work untouched.

`param_string` branches on the same distinction. It **used to** detect `:κ_geo` alone,
which was sufficient while (r, ξ) was the only dimensionless-coordinate sweep; the (x, ξ)
params tuple has `κ_geo` too, so that test now sends it down the (r, ξ) branch and reaches
for a nonexistent `p.x_fixed`. The branch is therefore keyed on **which coordinate is
fixed**: `:x_fixed` → the (r, ξ) form, `:r_fixed` → the (x, ξ) form. Adding a fourth sweep
means adding a branch here, and the failure mode is a `no field` error rather than a wrong
label — which is the right way round.

- (β, x) run → prints `κ_1, κ_2, g_2`
- (r, ξ) run → prints `x, β, sqrt(κ_1κ_2)`; κ and g all vary per grid point, so printing
  them would be actively wrong
- (x, ξ) run → prints `β, r, sqrt(κ_1κ_2)`; κ_1 and κ_2 happen to be constant too (r is
  fixed), but η and g move, one per axis

`run_plot_r_xi.jl` carries two per-key overrides that `run_plot.jl` does not.
`:concurrence_diff` gets `sym_clims` + `:balance` (issue 2's prescription, applied at the
call site). `:concurrence_adia` gets `flat_safe_clims`, which exists because this sweep
predicts that panel is *constant*: a data range ~1e-15 wide makes GR emit
`Rectangle definition is invalid in routine CELLARRAY` and drop the colorbar entirely.
The guard is **relative** — an `lo == hi` test does not catch a spread in the 15th digit —
and widens a flat panel to `(0, 2v)`. Anchoring at zero rather than a narrow band around
`v` also keeps GR's tick labels short enough not to overprint the colorbar title, the
same collision `si_scale` works around in `compare_runs.jl`.

`run_plot_x_xi.jl` carries **byte-identical copies of both helpers**, for the same reason
the schedulers are duplicated: `plotting_functions.jl` is the only file both plot drivers
include, and the drivers cannot include each other. On that plane `:concurrence_adia` is
*not* flat (it varies with x), so `flat_safe_clims` passes it through unchanged and only
earns its keep if the sweep is ever narrowed to one x column. **Lifting `overrides` and
`flat_safe_clims` into `plotting_functions.jl` is the clean fix for issue 2** and would
serve all three drivers at once; until then the two copies must stay in step.

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
(`compare_runs.jl`), `Full_vs_Adia_rxi/` (`run_plot_r_xi.jl`), `Full_vs_Adia_xxi/`
(`run_plot_x_xi.jl`) and `Adia_boundary/` (`adia_boundary.jl`), PDF only. `save_fig`'s own
default is `figures3/`, and `figures/`, `figures2/` also exist from earlier runs — eight
destinations, none tracked by git.
Note the `_xxi` filename suffix carries only `β, r, κ_geo`, so two (x, ξ) runs differing
in `len` or either range overwrite each other's PDFs (the `.jld2` names do not collide).

### Figure audits — a figure is not verified because someone looked at it

**This applies to every figure in this repo, not to any one script.** The numerics here are
checked hard: every sweep and scan carries a battery of pass/fail lines, and those batteries
have caught real errors, including one in `files_online`'s own tables. **None of that
machinery was ever pointed at the figures**, whose entire QA was "eyeball it" — and a
figure was consequently wrong for several sessions in a way that survived repeated looking.

The failure was not the artifact. It was writing into this file that a colour scale
"**hides nothing**". Looking at a panel can support "I did not notice anything"; it cannot
support a claim about all 6163 of its cells. **Do not state a negative about a figure
without a test that could fail.**

Three audits, all cheap, all pure functions of data already computed:

**1. Colour honesty.** For any colour-mapped panel, count the cells that *render as*
background/midpoint but whose *value* is materially non-zero. **The two thresholds must be
independent** — one perceptual (fraction of the colorbar), one physical (a value that
matters in the units of the quantity). Measured on the `∂C/∂x` panel:

| colour variable | cells reading as white | of those, MISLEADING | worst value painted white |
|---|---|---|---|
| `v` (linear) | 534 | **47** | 0.168 |
| `sign(v)·sqrt(\|v\|)` | 30 | **0** | 0.000 |
| `:split_log` | 0 — the band is empty as a *set* | **0** | 0.000 |
| **`:sign_map`** (the shipped default) | **N/A — there is no colour scale** | — | — |

**This audit runs inside `adia_derivatives.jl` on every invocation** and the `:linear` branch
is kept precisely so it can be re-fired — it still reports 3 misleading cells on a coarse
20×15 grid, worst value 0.138 painted as zero.

**Read the last two rows carefully, they are different kinds of claim.** `:split_log` earns
its zero structurally: no cell is *drawn* within `|z| < 1`, so "reads as white" is empty as a
set rather than empty as a count. `:sign_map` does not have a colour scale at all, so the
failure mode does not exist for it — the script prints that reason rather than a `0`, because
printing `0 misleading cells` for a panel with no midpoint colour would be a category error
dressed as a pass. **Neither row licenses the phrase "the figure hides nothing"**; they
license exactly what they say.

The linear panel appeared to have **two** zero contours where the data has one. **White in
a diverging colormap means "small", and how small is set by an extreme value elsewhere in
the panel** — here the `x → 0` divergence pushed `|clims|` to 2.80, so `−0.215` rendered
identically to `0`. The fix for a divergent or long-tailed field is a *transform*, not a
clip: a percentile clip trades the fake-zero artifact for a saturated band somewhere else.

**2. Frame containment.** Report, per drawn curve, what fraction of it lies inside the axis
limits. A curve leaving the frame must be *stated*, never left to be inferred from a line
that runs off the side. On the `∂C/∂x` panel the ridge is inside for **101 of 180 rows** —
it exits the left edge at `u = 52.6`, so the upper 44% of the plotted range shows none of
it, which is invisible unless reported. Same check on the boundary figure puts `u*(x)`
inside for only 6% of columns, which is correct and expected but worth knowing.

**3. Label–encoding agreement.** If the colour variable is transformed, the colorbar
**title** must name the transform. This is not stylistic: **GR does not reliably honour
custom `colorbar_ticks`**, so relabelling the numbers back into original units is not
available, and a bar titled `∂C/∂u` over `log10|∂C/∂u|` values is simply false. Shipped
that way once. Related and equally load-bearing: **on a panel with an active `marker_z`
colorbar, GR may not give the colour that was requested** — `:limegreen` and
`RGB(0,0.7,0)` both rendered pale blue-grey, `:magenta` rendered dark red, a `:lime` marker
rendered cyan, while `:yellow` markers and `lc = :black` survived. Plots reports the
attribute as exactly what was asked for, so nothing errors and nothing warns. Pick overlay
colours **by looking at the output**, and treat the keyword as a request.

**The trap that makes all three feel done when they are not:** a check whose two thresholds
are the same number **cannot fail**. The first version of audit 1 defined "reads as white"
and "is not zero" off the same quantity, reported `0 misleading` for both the broken and
the fixed panel, and looked like a clean pass. **Before trusting a new check, confirm it
fires on the known-bad case.** Every audit above is stated with the number it produced on
the version that was wrong, for exactly that reason.

Corollary for this file: **"eyeball it" is not a verification step, it is the record that
none was performed.** Where a figure property is genuinely only inspected, say so, rather
than logging it among the checks.

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
far below the `x_c = 0.4839` point where adiabatic entanglement dies entirely.

**Both numbers now have exact closed forms**, from `file_concurrence/concurrence_adia.md`
§10.7, and the script's measured values match them to the last printed digit:

- `x* = (t*)²  = 0.1485501998979953`, where `t*` is the root of
  `R(t) = t⁶ + t⁵ − 3t⁴ − 2t³ − 3t² − t + 1` on (0,1)
- `C* = C(2,x*) = 0.23761539689930783`, from
  `C_adia(β=1,x) = [2√x(1−x) − x(1+x)] / [2(x²+1)]`

So the scan is now a *check* on the algebra rather than the only source of the number.
Keep running it: it is seconds, and it exercises `run_adia` end to end.

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
  **The stored value is also `C(2, 0.12875) = 0.23604367495075856` exactly**, so this
  check now runs against algebra as well as against itself.

**β = 1 is optimal at every x — PROVEN**, not just at the peak and no longer only sampled.
`file_concurrence/concurrence_adia.md` §10 establishes `dC/du < 0` on the entire physical
window for every `x ∈ (0, x_c)`, where `u = b = β + 1/β`, so β = 1 is the *unique* global
maximiser. The older route — `C_adia` depends on β only through `b ≥ 2` and decreases
monotonically in `b` — is `files_online/06_concurrence.md` §6 and remains a second
independent derivation. Asymmetry also *lowers* the optimal pump: `x_opt` = 0.1486 at
β = 1, 0.107 at β = 2, 0.042 at β = 10 (`files_online/06` §6; §10.7's closed form is the
`u = 2` slice only, so it does not cover the β ≠ 1 entries). `files_online`'s closed form
agrees with `run_adia` to 13 digits at both `x*` and run A's `x = 0.12875` grid point —
a stronger cross-check on `H_adia`/`J_adia` than anything in this repo alone, since the
two were derived independently, and `file_concurrence` is now a third such route.

For contrast, the **full** model peaks elsewhere and higher: `C = 0.29675` at `x = 0.291`,
`β = r = 1`, `ξ = 2.75` (`files_online/07`), i.e. 25% above the adiabatic maximum at
roughly twice the pump. `config_r_xi.jl`'s `ξ_range = (0.1, 1.0)` excludes it, but
**sweep 3 brackets it**: `ξ = 2.9375` and `x ∈ {0.2475, 0.3663}` are grid points there.

## The full model's concurrence maximum — measured, script not retained

```
C_full is maximal at  x = 0.2910,  β = 1,  r = 1,  ξ ≈ 2.79,  C = 0.296758
realized as g₁ = g₂ = 1.66938, κ₁ = κ₂ = 2, η = 0.53944
```

**Independent confirmation of `files_online/07`**, which has `C = 0.29675` at `x = 0.291`,
`β = r = 1`, `𝒜 = 2.75` — reached there analytically, reached here by solving the full
cavity Liouvillian on a grid. Two codebases, agreeing to 5 digits. Compare the adiabatic
optimum below: **the full model gets 25% more entanglement at roughly twice the pump**,
and `config_r_xi.jl`'s `ξ_range = (0.1, 1.0)` excludes `ξ ≈ 2.8` entirely, so sweep 2 has
never run near this point. **Sweep 3 does** — its `(x, ξ)` grid at `β = r = 1` reaches
`ξ = 5` and brackets the optimum, which is a third independent route to the same number.

**How many digits are real.** Two independent runs gave `C = 0.29675767` and `0.29675774`
(Δ = 7e-08), the fitted `x` differing by 1e-06 and `ξ` by 6e-06 — that is
`steadystate.iterative`'s unseeded-RNG floor at `reltol ≈ 1.5e-8`, the same one documented
under `compare_param_sets.jl`. `β = r = 1` are exact (arm-swap symmetry, not a fit); `x` is
good to ~3 figures; **`ξ` is soft** — ±5% in ξ costs only 2.7e-06 in C, so anything from
ξ ≈ 2.6 to 2.9 is indistinguishable. That is why `files_online`'s `𝒜 = 2.75` and the
numerical `2.787` are the same answer.

Structure around the optimum, all measured:

| coordinate | −5% | +5% |
|---|---|---|
| ξ | −2.6e-06 | −5.7e-05 |
| x | −2.2e-04 | −2.0e-04 |
| r | −6.4e-04 | −5.8e-04 |
| β | −8.5e-04 | −7.7e-04 |

β and r are sharply peaked at 1 and mirror-symmetric — at `x = 0.29, ξ = 2.75`,
`C_full` runs 0.0886 / 0.1956 / **0.2968** / 0.1956 / 0.0886 over β = 0.25 … 4, and
0.1095 / 0.2108 / **0.2968** / 0.2108 / 0.1095 over r = 0.25 … 4. The arm-swap symmetry
`C(x,β,r,ξ) = C(x,1/β,1/r,ξ)` held to `1.1e-08`, the solver floor. **Truncation was
verified**: re-solving the optimum at `N = 20` instead of 14 moved `C` by `−3e-09`.

### The adiabatic model at that point

```
C_adia = 0.17943675   vs   C_full = 0.29675767     ratio 1.654
D(ρ_full, ρ_adia) = 0.0779     purity 0.750 vs 0.633
```

**The elimination loses 39.5% of the entanglement that is actually there**, at the very
point one would design around. Three things make this more than a single number:

- **`ξ` and `r` are invisible to `ρ_adia`** — measured at the same `(x, β)` with
  `(r, ξ)` = (1, 0.01), (1, 25), (7.3, 2.79), (0.05, 0.4): identical to **13 digits**.
  Including at ξ = 25, two orders past where the elimination is derived. ξ *is* its
  domain of validity, and the model cannot see ξ, so it cannot tell it is being misused.
- **The point is not near a maximum of the adiabatic surface**: `∂C_adia/∂x = −0.717`
  while `∂C_adia/∂β = 0` exactly (β = 1 is stationary by symmetry). `C_adia` there is also
  24.5% below its own peak — **the two models want different pumps**, x = 0.291 vs 0.149.
- **The error grows steeply in x while the full model sits flat**: across ±10% in x,
  `C_full − C_adia` runs 0.0970 → 0.1173 → 0.1387 while `C_full` moves only
  0.2959 → 0.2968 → 0.2960.

**Method, if this is ever rebuilt.** Cost forces the shape: one full solve is ~12 s at
dim 900, so a sequential 1D optimizer chain is hopeless (hundreds of *serial* solves,
single-threaded). What worked: a (β, r) probe to establish β = r = 1 is optimal, a coarse
(x, ξ) grid on that slice under `@threads`, a zoom grid, then a **2D quadratic fit to
points already paid for** to get a sub-grid vertex for the price of one confirming solve —
then a 4D ±5% peak test and a truncation check. ~11 min on 4 threads.

## Validity of the elimination — framing and cost, script not retained

**The obvious question is degenerate and must be refused.** `D(ρ_full, ρ_adia) → 0` as
`ξ → 0` at *every* `(x, β, r)`, because that limit **is** the elimination. So `argmin D`
is the whole `ξ = 0` hyperplane — and useless in practice, since `ξ → 0` sends the
effective qubit-qubit rates to zero with it and the steady state is reached arbitrarily
slowly. The meaningful object is the **level set** `{D ≤ ε}`, described by its boundary
`ξ_max(ε; x, β, r)`: the largest coupling whose error stays within ε.

Four things established before the mapping run was stopped:

- **`D` is monotone and smooth in ξ, and saturates.** Log-log slopes at `(x,β,r)=(0.3,1,1)`
  run 0.97 → 0.93 → 0.83 → 0.68 → 0.44 → 0.20 across ξ = 0.003 → 3. `D ∝ ξ` therefore
  holds only below **ξ ≈ 0.01**; extrapolating that slope upward badly overestimates the
  error, so any boundary must be *bracketed by measurements*, never extrapolated.
- **The relative concurrence error is ~10× D**, because `|ΔC| ≲ 2D` but `C_full ≈ 0.18 < 1`.
  Measured at `(0.3, 1, 1)`: `D = 2.5e-3` while `|ΔC|/C = 2.4e-2` at ξ = 0.01. **The
  observable you would quote needs ξ roughly 5–10× smaller than the trace-distance
  criterion suggests.**
- **`|ΔC|/C_full` is non-monotone in ξ.** At `(0.3, 1, 4)` it rises to 0.105 at ξ = 0.15
  and falls back to 0.047 by ξ = 0.5, because `C_full` crosses the ξ-independent `C_adia`
  coming down. So a validity boundary must be the **first upcrossing**, not the last point
  that happens to satisfy the bound — the latter would call that cell valid to ξ = 0.5
  while hiding a 10% error inside.
- **Above `x_c = 0.4839` the concurrence criterion is not merely violated, it is
  undefined.** `C_adia = 0` exactly there, so at small ξ both models give ~0 (a 0/0
  relative error) and at any larger ξ the error is 100% by construction. That is a
  *qualitative* failure of the elimination, and no tolerance is ever met at any coupling.

Sample boundaries actually measured, at `x = 0.3, β = 1`:
`ξ_max(D ≤ 0.01, 0.02, 0.05) = 0.047, 0.111, 0.483` at `r = 1`, and
`0.035, 0.087, 0.362` at `r = 4` — asymmetry tightens the requirement.

### The truncation cliff — a repo-wide cost fact worth keeping

`estimate_truncation` holds everything up to `x = 0.62` in the `N = 14 / dim = 900`
cluster, then jumps:

| x | 0.50 | 0.55 | 0.60 | **0.62** | **0.65** | 0.68 | 0.70 |
|---|---|---|---|---|---|---|---|
| N | 14 | 14 | 14 | **14** | **24** | 27 | 29 |
| dim | 900 | 900 | 900 | **900** | **2500** | 3136 | 3600 |

Crossing it costs roughly **40× per solve** — 16× the state-vector length on a Liouvillian
gap 2.8× smaller. This is the 900-vs-2704 clustering of issue 3, and it is the biggest
cost discontinuity in the repo. Measured consequence: a 22-cell grid reaching `x = 0.70`
ran **17 minutes without completing a single cell**, while the same grid capped at
`x = 0.62` returned cells in 360–510 s. **Any new sweep touching x ≳ 0.63 must budget for
this explicitly rather than discover it.**

A second operational lesson from the same run: Julia buffers stdout when redirected, so a
long job wrapped in a shell `timeout` can be killed having produced **nothing** — 45
minutes of solves lost that way. Long runs need per-item `flush(stdout)` and incremental
persistence, not an all-or-nothing write at the end.

## Parameter counting — `compare_param_sets.jl`

A standalone script. It does not sweep, does not write a `.jld2`, and is not on the path
of anything else. It exists to establish that `ρ_full` is a function of exactly the four
dimensionless coordinates. (A second script, a randomized version of the same test, was
deliberately deleted — its results are not used anywhere.)

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
  is scale-invariant" from "this code ignores its input". There is no longer a script that
  supplies one, so if you add a row expected to give `D ≈ 0`, add a second row that
  perturbs one invariant by ~1% and must give `D` well above the `~1e-8` floor. Without
  that pairing a passing null test proves nothing.

The elimination error, `D(ρ_full, ρ_adia)` per set, is the scale the other two matrices
should be read against. At the recorded points it is **0.10–0.15**, i.e. *larger* than
any of the between-set distances. So those tables compare two models that are already far
apart — `ξ ≈ 1` is nowhere near the `ξ → 0` limit the elimination is derived in.

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

- **The `x_c = 0.48388826` entanglement threshold** (root of `x³ − 2x² + 9x − 4`) for the
  adiabatic model at β = 1, and the `b = β + 1/β < b_max(x)` boundary at general β. This
  is the closed form behind the `C_adia = 0` dead zone documented under the (r, ξ) sweep —
  the empirical scan there just rediscovers it. **`file_concurrence/concurrence_adia.md`
  §3.3 derives the same cubic independently**, without the `√x` substitution, and its §4
  gives `b_max(x)` in closed form; the two projects agree.
- **The full model's concurrence maximum**: `C = 0.29675` at `x = 0.291`, `β = r = 1`,
  `𝒜 = 2.75` (i.e. **this repo's ξ = 2.75**), 25% above the adiabatic maximum of 0.2376.
  Sweep 2's `ξ_range = (0.1, 1.0)` **excludes it**; sweep 3's `(0.2, 5.0)` brackets it.
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

## `file_concurrence/` — closed-form analysis of the adiabatic concurrence

⚠ **This section describes rev. 4 and the file has moved on.** A session note dated later
cites **rev. 5**, carrying a *universal* `∂C/∂u < 0` proof via `Q > 0 on ℝ × (0,1)` — which
would supersede the rev. 4 story below of a degree-174 discriminant needed on the
`x < x_dagger` sliver. Nothing here has been rewritten from that note alone: **read the .md
before trusting the revision-specific claims in this section.** `adia_derivatives.jl`'s
check (a2) adds a third, independent large-u route (`A₀, A₁ < 0`) that is revision-agnostic.

One file, `concurrence_adia.md` (~1220 lines, 56 KB, rev. 4), a handoff document that
solves the adiabatic model's concurrence in closed form. It is **not** part of
`files_online/` and has no `claude.md` of its own — **this** file governs it. Tracked by
git as of the `.gitignore` exception described under Conventions; watched by the session
hook.

**Conventions match this repo exactly**, which is the opposite of the `files_online/`
situation and worth stating plainly: it defines `x = 4η²/(κ_1κ_2)` and
`β = g_1²κ_2/(g_2²κ_1)` identically (its line 125), and ξ, r, μ, α never appear as
physical parameters at all — because `ρ_adia` does not depend on them (proved under the
(r, ξ) sweep above). Nothing needs reconciling before it becomes a figure. Two local
traps only:

- **⚠ `x*` means the argmax there, not the threshold.** That file uses `x_c` for the
  0.48389 entanglement threshold and `x*` for the 0.14855 argmax. **This file has been
  brought into line with it** — `x_c` is the threshold everywhere below, `x*` the argmax.
  Older prose and any figure captions predating this may still say `x*` for 0.48389.
- `alpha`/`beta` in its §10.3 are the coefficients of `L = alpha(x)·u + beta(x)`, not the
  physical β. Confined to that subsection.

**Everything is done in `u = β + 1/β ≥ 2`**, and that substitution is the whole reason the
file works. It makes the `β ↔ 1/β` symmetry automatic (`u = 2` is exactly β = 1), and
turns the concurrence into

```
C(u,x) = 2[ v·g − sqrt(S) ] / D      v = sqrt(u+2),  D linear in u,  S quadratic in u
```

### Result 1 — existence, PROVEN (§3)

`D > 0`, so `sign C = sign Φ` with `Φ = (u+2)g² − S`, and **`Φ` is a downward parabola in
u**. That single fact carries the argument:

- `Φ(2,x) > 0 ⟺ x³ − 2x² + 9x − 4 < 0`, giving **`x_c = 0.4838882550`** — the same
  threshold `files_online/06` reaches, but derived without the `√x` substitution, which
  makes the x-cubic the cleaner object.
- **No-revival lemma:** the parabola's vertex crosses `u = 2` at `x_v ≈ 0.20813`, and
  **`x_v < x_c` is the load-bearing inequality.** Had it gone the other way there would be
  a band `(x_c, x_v)` where C vanishes at β = 1 but revives at larger asymmetry. That band
  is empty.
- **Theorem.** `x ≥ x_c` ⟹ `C_adia = 0` for every β. `x < x_c` ⟹ `C_adia > 0` exactly on
  `u ∈ [2, u_c(x))`.

### Result 2 — `u_c(x)` in closed form (§4)

**This is the object this repo has nowhere.** The sweeps only ever see the dead zone as
`C_adia = 0` cells on a grid; §4 gives its boundary *in β* analytically. Two equivalent
forms — §4.1 compact in `q1, p1, p2`, §4.2 with a degree-10 radicand (better for root
isolation). Cross-checked against previously stored root-finding values: 577.76 / 18.4304
/ 2.16302 at x = 0.02 / 0.1 / 0.4. Asymptotically `u_c ≈ 1/(4x²)` at small x — the window
widens quadratically as the drive is turned down — and `u_c → 2` as `x → x_c⁻`, i.e. the
window pinches shut at β = 1 exactly at threshold. In β the window is `(1/β₊, β₊)` with
`β₊ = [u_c + sqrt(u_c²−4)]/2`, symmetric on a log axis as the symmetry demands.

### Result 3 — the maximum, PROVEN (§10, new in rev. 4)

> **`dC/du < 0` on the entire physical window, for every `x ∈ (0, x_c)`.**
> Hence **β = 1 is the unique global maximiser**, and the maximum is `C(2,x)` below.

This is the upgrade of what this file previously recorded as sampled. Two elementary
structural facts do ~90% of the work, and they are the transferable part:

- **`F_u < 0` identically** (§10.2). Since `dC/du = F_u/√F − G_u/√G`, this means
  **`G_u ≥ 0` is by itself sufficient** — negative minus non-negative. No comparison of
  magnitudes is ever needed.
- **`L := S'D − 2SD'` is linear in u with strictly positive slope** (§10.3). So it has at
  most one crossing and its sign on the window is decided by its value at `u = 2` alone;
  the right endpoint never needs separate examination.

Only the strip `x < x_dagger ≈ 0.0518` needs the heavy machinery, and even there only the
degree-174 discriminant's *root-count constancy*, settled by one exact Sturm count at
`x = 1/50`. **The cheap structural fact beat the expensive computation** — the planned CAD
existential was never run.

### Result 4 — `C(2,x)` in closed form (§10.7), and a three-way cross-check

> **`C_adia(β=1, x) = [ 2√x (1−x) − x(1+x) ] / [ 2(x²+1) ]`**

Its numerator vanishing recovers the same Cardano cubic as §3.3, and its argmax is
`x* = (t*)²` where `t*` is the root of `R(t) = t⁶ + t⁵ − 3t⁴ − 2t³ − 3t² − t + 1`. **So
`adia_concurrence_max.jl`'s numerical `x* = 0.148550` now has an exact algebraic source.**

Checked against this repo's own stored numbers — it reproduces them to the last digit:

| from `C(2,x)` | value | this repo's number | source |
|---|---|---|---|
| `C(2, 0.12875)` | 0.23604367495075856 | 0.2360436749507588 | run A's stored `C_adia` max |
| `C(2, x*)` | 0.23761539689930783 | 0.2376153969 | `adia_concurrence_max.jl` |
| `C(2, 0.3663)` | 0.11752557550041336 | 0.118 / 0.1175 | the β-table under the (r, ξ) sweep |
| `x*` | 0.1485501998979953 | 0.148550 | `adia_concurrence_max.jl` |
| `x_c` | 0.48388825504305466 | 0.48388826 | this file, and `files_online/06` |

That is **three independent routes agreeing**: hand-plus-sympy algebra, this repo's
numerical `run_adia` + Wootters, and `files_online`'s separate analytics. Reproduce it in
one line if you doubt it:

```bash
python3 -c "import math
C=lambda x:(2*math.sqrt(x)*(1-x)-x*(1+x))/(2*(x*x+1))
print(C(0.12875), C(0.1485501998979953), C(0.3663))"
```

### ⚠ What is NOT verified (its own §10.8 — read before quoting)

- **Sections 2–4 are hand-derived and CAS-unconfirmed.** Its §7 Mathematica block has
  **never been run.** §10 is built on top of §2's element set and inherits any error in it.
- **§10 is verified in one CAS only** (sympy). A second engine is on its own to-do list.
- Its §1.2 domain claim (`x ∈ (0,1)` from the below-threshold condition) is flagged by the
  author as not checked against a specific equation in a specific reference, and its §1.3
  Yu–Eberly X-state citation as lacking a verified equation number. Do not put either in
  front of a supervisor as sourced.

Against that: the closed forms reproduce this repo's independent numerics to 14 digits
(Result 4), which is a stronger check on the element set than anything internal to either
project — a wrong `ρ_adia` would not land on `run_adia`'s answer by accident. Two real
errors were caught during its development by exactly that kind of spot-evaluation (a
`sqrt(S)(2)` off by a factor of x; a wrong coefficient in `W_old`).

### The methodological warning worth importing (§5.2)

**Sign-preserving transformations are not extremum-preserving.** `Φ`, `N` and `C` differ by
strictly positive factors, so they share every *sign* and none of their *critical points*.
Explicit counterexample at x = 0.11: `N` is larger off-centre (0.412056 vs 0.411708) while
`C` is smaller (0.218355 vs 0.231331) — `argmax N ≠ argmax C`, demonstrably. Its §3
existence proof is safe precisely because existence only needs a sign. This repo computes
concurrence as a ratio too, so the same trap is reachable here.

## The entanglement boundary, measured — `adia_boundary.jl`

The numerical check on everything above. It maps `C_adia` over the **(x, u) plane** with
`run_adia` and puts the closed-form `u_c(x)` of §4.1 on top of it. Output is
`Adia_boundary/adia_boundary.pdf` plus six checks on stdout. ~10 s, no threads, no `.jld2`
read or written.

```
VERDICT: closed form and run_adia agree on the boundary to 1.2e-12 (worst, relative)
```

**The load-bearing check is (b), and it is a comparison of two things that share no code.**
For each x column the script brackets the `C > 0 → C = 0` transition by *doubling* u and
then bisects it, never consulting `uc_closed`. Over 181 live columns the bisected boundary
matches §4.1 to **1.2e-12 worst, 2.2e-14 median**. `C_adia` vanishes *linearly* in
`(u_c − u)` — `N = vg − √S` is linear near its root — and `calc_concurrence` clamps at
exactly `0.0`, so the crossing is a clean step and bisection resolves it to nearly the last
bit. That precision is the reason this is worth running rather than eyeballing.

The other five:

- **(a) §4.4's three stored values.** Two match to ~5e-6 relative. **The `x = 0.02` row of
  §4.4 is a last-digit slip in the `.md`**: it prints `577.76`, both routes here give
  `577.751476`, and §4.4's *own* second column (the older root-finding result) says
  `~577.75`. Display value only, nothing depends on it — but do not "correct" this repo to
  match it.
- **(c) the curve meets `u = 2` exactly at `x_c`.** `u_c(x_c) = 2.000000000000`, and
  root-finding `u_c(x) − 2` returns `x_c` to `4e-16` against the cubic. Two independent
  routes to the threshold, and the visual centrepiece of the figure: the coloured wedge
  pinches shut precisely on the `x_c` line.
- **(d) the small-x asymptote** `u_c ≈ 1/(4x²) − 1/x + 4` (§4.5), ratio 0.9999 at x = 0.005
  drifting to 0.989 at x = 0.05 — the expected `O(x)` error, and a straight line of slope
  −2 at the left edge of the log-log figure.
- **(e) no revival.** All 39 grid columns above `x_c` × 180 u-values each come out
  **exactly `0.000e+00`**. This is the half of §3.5 that needed `x_v < x_c`, tested
  directly rather than argued.
- **(f) anchored outside itself.** `C_of(x*, β=1)` matches §10.7's `C(2,x*)` to `1.4e-16`,
  and `C_of(0.12875, β=1)` reproduces run A's stored `0.2360436749507588` to `0.0e+00`.
- **(g) the proof's own geometry**, added when the derivative maps were built. `u*(x)` from
  §10.3 — where `G_u` changes sign, the region-2/region-3 divide, a *rational* function
  with no square root — must cross `u = 2` exactly at `x_dagger` from §10.4. It does:
  `u*(x_dagger) = 2.000000000000`, and root-finding `u*(x) − 2` returns `x_dagger` to
  `3.9e-16`. **Nothing in the `.md` ties `V2`, `V3` and `V` together**, so this checks all
  three transcriptions at once. The figure draws `u*(x)` and the `x_dagger` vertical, which
  shows **region 3 to scale** — the sliver `x < 0.0518, u ∈ [2, u*(x))` that needed the
  degree-174 discriminant is a visibly small corner of the domain, and everywhere else
  `dC/du < 0` falls out of "negative minus non-negative".

**Axes are forced, not chosen.** `u_c ~ 1/(4x²)` as `x → 0`, so `u_c` spans 2 → ~10⁴ and
**the y axis must be log**; `x` cannot start at 0 (every denominator carries `8x²q1`), so
`x_lo` defaults to 0.02 where `u_c = 577.75`. `x_hi` defaults to `x_c + 0.1` on purpose —
the columns past the threshold are what check (e) tests.

Two implementation notes worth keeping:

- **`C_of` is duplicated from `adia_concurrence_max.jl`, not included.** That file has no
  `main()` guard, so including it would run its whole scan — the same reason `observables.jl`
  exists. If the two copies ever disagree, `adia_concurrence_max.jl` is the older and
  better-checked one.
- **This is the one figure whose log axis is vertical.** `plot_map` hardcodes the row
  coordinate onto the *horizontal* log axis; here the log coordinate is `u` and belongs on
  the vertical. So `flatten_grid` is reused for what it is good at — pairing cells with
  coordinates and dropping non-finite ones — and its two outputs are handed to `scatter`
  **swapped**. Dead cells are drawn in grey rather than dropped, deliberately: an absence
  of markers would not distinguish "measured, came out zero" from "never sampled", and the
  measured dead zone is the whole point.

## Inside the region — `adia_derivatives.jl`

Where `adia_boundary.jl` settles the region's *edge*, this maps its *interior*: `∂C/∂u` and
`∂C/∂x` over the same (x, u) plane, from **central finite differences of `run_adia`**,
checked against §10.1. Two figures into `Adia_boundary/`: **`dC_du.pdf`, now a single map
panel** (its band panel was removed — see the panel history below), and **`dC_dx.pdf`, two
panels** — a map, and under it a panel carrying the same claim without colour. ~40 s, no
threads.

```
VERDICT: 0 of 6163 sampled ∂C/∂u are >= 0 (theorem holds numerically);  ridge at u=2 hits x* to 3.6e-11
```

**The two fields are not equally interesting, and the figures differ because of it.**

**`∂C/∂u` is single-signed** — that *is* §10's theorem — so check (a) tests it at every one
of 6163 live points and finds a maximum of `−9.91e-06`. As a picture that is one boolean:
on a diverging colormap half the range goes unused, and on a linear one the panel is a
single flat colour (median `−0.149` against a max of `−9.7e-06`, **5 decades**).

**Two earlier encodings both failed the same way, and the reason is worth keeping.**
`log10|∂C/∂u|` on `:magma` recovered the structure but *deleted the sign* — the panel would
have looked identical had a stray positive value existed, so "it is negative" was true only
because the title said so. `sign(v)·log10(1+|v|/1e-6)` on `:balance` then put every cell on
one half of a symmetric bar, which is better but still asks the reader to notice that the
other half is empty **and** still prints a `|·|` on the colorbar. Neither shows the sign of
a cell; both report a magnitude and assert a sign beside it.

The current default is **`DU_SCALE = :neg_log`: colour = `log10(−∂C/∂u)`, no absolute value
anywhere.** The point is not the label — it is that the transform **does not exist** off the
theorem. A cell with `∂C/∂u ≥ 0` has no `log10(−v)`, so it is split out before the colour
map (`pos`/`neg`), counted on stdout, and over-plotted as red markers whose legend entry is
drawn *only if the set is non-empty* — an always-present "0 cells" entry would be an
assertion, a series that appears from nowhere is evidence. `:signed_log`, `:log` and
`:linear` are kept as comparison branches.

**The colormap is `cgrad(:magma, rev = true)`, and the reversal is not cosmetic.** Plain
`:magma` puts bright at high `log10(−∂C/∂u)` — i.e. **bright = most negative**, which inverts
the ordering of the underlying quantity and makes the panel's most conspicuous feature its
*least* marginal region. A reader scanning for "is this negative?" was pulled to the `u = 2`,
`x → x_c` corner where `∂C/∂u = −0.94` and the answer is least in doubt, while the cells that
actually come closest to zero (`−9.9e-06`, at `u ≈ 565`) sat in near-black obscurity. This
was misread in exactly that way once, as "the maximum of the whole plot is at `u ≈ 2`" — it
is the *minimum*. Reversed, bright = closest to zero, so the eye goes where the claim is
tightest. **The direction is stated on the figure, in the title's second line** —
`bright = closest to zero (weakest) · dark = most negative`. It has moved twice and the
history is the warning: it started in the title, was lost when the title was cut to the
single line `"∂C/∂u < 0 EVERYWHERE"`, survived only in panel 2's legend string "the dark
band of the map", and came back to the title when panel 2 was deleted. **A reversed colormap
whose direction is stated nowhere is worse than an unreversed one**, so if that title line is
ever trimmed, the direction has to land somewhere else in the same edit.

**Figure 1 is ONE PANEL — the map. The `−∂C/∂u` band panel was removed on request.**
What it drew and where that content went:

- its **lower** edge (`min_x`, at `x → x_max(u)`, the death boundary) was the one intrinsic
  curve on it, and is now unillustrated. If any of it is ever wanted back, that is the part
  worth redrawing — **not** the band
- its **upper** edge was `DU[:, 1]`, a slice along `X_LO`, i.e. the window artifact below.
  Deleting it deletes a curve that was being read as a margin
- it never carried the sign claim and **could not**: on a log axis a value with `∂C/∂u ≥ 0`
  is *unplottable* — dropped, not drawn on the wrong side — so the panel could not falsify
  its own headline. Check (a)'s violation **count** carries the sign, and it can fail

**The row reductions are kept even though nothing draws them.** They cross-check check (a)'s
extrema by a different reduction of the same grid (`Δ = 0.0e+00`), and the 175-of-175
assertion is what licenses check (a)'s "the weakest sits at column 1" and all of (a2)'s
framing. They were checks that happened to have a picture, not a picture that printed
numbers — which is why deleting the picture cost nothing. Their stdout labels no longer name
a panel.

### ⚠ `max ∂C/∂u` is a window corner, not a margin — the mistake this section is about

This is the most important thing in this section, and it was wrong here for a long time
while looking like a result.

**`∂C/∂u` is strictly monotone decreasing in x at fixed u.** Asserted on every run now: the
weakest value in a row sits at the row's first live column in **175 of 175 live rows**. So
`max_x ∂C/∂u` is *always* attained at the smallest live x, which is `X_LO` — the edge of the
window. The panel's upper edge is literally `DU[:, 1]`.

**The supremum over the physical region is 0.** `A₀ ~ −√x` (below), so `∂C/∂u → 0` as
`x → 0`, and `x → 0` is *inside* the region at every u. It also → 0 as `u → ∞` via `u^{−3/2}`.
**No grid can ever exhibit a strictly negative bound**, so any number quoted as one is a
coordinate of the window. `X_LO` is a **CLI argument**, and the "max" tracks it:

| `X_LO` | reported value | vs FD floor `1.1e-9` |
|---|---|---|
| 0.02 | −9.58e-06 | 8713× |
| 0.005 | −7.18e-08 | 65× |
| 0.002 | −2.88e-09 | 2.6× |
| 0.001 | −2.54e-10 | **0.2× — below the floor** |

`adia_derivatives.jl 0.001 …` used to print `max ∂C/∂u = −2.54e-10 → uniformly negative`: a
pass declared on a number whose sign the file's own finite differences cannot resolve.
**There is now a guard that fires on exactly that invocation** (`WINDOW TOO NARROW`), and it
has been verified to fire — a guard never shown to fail is not a guard.

**The word "maximum" was itself the bug.** `−9.91e-06` is the *least negative* value, while
the largest magnitude in the same grid is `0.937`, five decades away. Calling the small one
"the max" is the exact wording that got the figure read backwards, twice in one session.
**Check (a) is now a boolean and a count** — `points with ∂C/∂u >= 0: 0` — which is what a
sign claim actually is and which can fail; the extremal values are reported after it, named
"strongest sampled" and "weakest sampled", never "max".

Corollary, and the reason the band panel was ultimately deletable: **its two edges were not
the same kind of object.** The lower one (`min_x`) sits at `x → x_max(u)`, the death
boundary, inside the window since `X_HI > x_c` — **intrinsic**. The upper one (`max_x`) was
the window edge — **an artifact**. A "range" between a real curve and an arbitrary one is
not a range of anything.

### Figure 1's panel history — three panels, then two, then one

Worth keeping in full, because every step was driven by a real defect and two of the fixes
were wrong:

1. **Three panels.** Map, plus a linear untransformed `max_x ∂C/∂u` against a bold zero
   line, plus the full-range band. The linear panel existed so the sign claim had somewhere
   an actual zero could be drawn — a log axis never can.
2. **The middle panel went, on request.** Note what left with it: **there is no panel with a
   zero line drawn on it anywhere in this file.** Before that, an even earlier version had
   plotted only each row's *maximum*, which silently dropped the map's most conspicuous
   feature (the strong band along `u = 2`, `∂C/∂u → −0.94`) and read as the two panels
   disagreeing about the data. **A summary panel that drops the extremum reads as missing
   data** — that is the transferable lesson, and it is why the range panel replaced it.
3. **The y-axis flip: tried, shipped, reverted, all in one session.** The complaint was
   real — at `u ≈ 565` the band pinches to `9.91e-06` at the bottom of the frame, which
   reads as the most extreme place on the plot while being the value closest to zero. But
   `yflip` fixes the *position* reading by breaking the *slope* reading: both edges fall
   3–5 decades in u and are then drawn **rising**, so anyone not registering the flip reads
   growth where there is decay. It was also answering the wrong question, since the quantity
   being repositioned was half artifact. **Content before encoding.**
4. **The band panel went too, on request. Figure 1 is now the map alone.** This is the
   cheapest resolution of everything above: the artifact edge is gone, the un-falsifiable
   sign display is gone, and nothing was lost that the checks did not already carry.

The through-line: **a log axis cannot demonstrate a sign at all.** A value with
`∂C/∂u ≥ 0` is *unplottable* — dropped, not drawn on the wrong side. Every version of this
panel from step 2 onward was displaying a claim it structurally could not contradict. Check
(a)'s violation count is what carries the sign, and unlike any of the panels, it can fail.

**Do not re-add a panel here without saying which of these four it is.**

**Sizing.** Figure 1 is now a single panel at `size = (900, 640)`. **Figure 2 is still two**,
at `size = (900, 980)` with `heights = [0.62, 0.38]`, and it needs the explicit
`left_margin` — the default clips the lower panel's rotated `ylabel` silently, with no
warning and no error. Both keep `left_margin = 7Plots.mm`; on figure 1 it is now
belt-and-braces rather than load-bearing.

**`∂C/∂x` changes sign** (§10.7), range `[−1.46, +2.80]`, and after four attempts at encoding
it the panel now has **no colour scale at all**. `DX_SCALE = :sign_map` is the default:

- **Colour carries the sign, and only the sign** — two flat pale colours, warm for
  `∂C/∂x > 0`, cool for `∂C/∂x < 0`. That is the one thing a colour can say without a
  transform, so it is the only thing this one says.
- **Magnitude comes back as labelled contours**, at `c ∈ {±2, ±1, ±0.5, ±0.1}` (those that
  exist), each drawn in its own colour and annotated **with its own value**. Nothing on the
  figure is a transformed number any more: every number printed on it is a value of `∂C/∂x`.
- **`c = 0` is the ridge**, drawn black as the heaviest member of the same family — so the
  contour set and the payload are one object rather than a curve laid over an unrelated map.

The three colour-mapped scales are kept as branches (`:split_log`, `:signed_sqrt`,
`:linear`), and the sequence is worth knowing because each fixed the previous one's defect
and introduced its own:

| branch | fixed | broke |
|---|---|---|
| `:linear` | — | painted a nonzero strip the colour of zero → a **fake second zero contour** |
| `:signed_sqrt` | the fake contour | did it with an `abs()`: the sign became a factor on a magnitude, so the panel would render identically if a sign were wrong |
| `:split_log` | the `abs()`; `|z| < 1` empty by construction, so nothing *can* render as zero | a colorbar reading `±[1 + log₁₀(±∂C/∂x / 10⁻⁵)]` — three ideas on one axis |
| **`:sign_map`** | the decoding, by having no scale to decode | magnitude is now sampled at 7 levels rather than continuous |

**`:linear` is retained on purpose**: it is the known-bad case the colour-honesty audit must
fire on, and it does — see the audit numbers below. **Do not delete it to tidy the file.**

Two things that survive the removal of the colour scale, because neither is a plotting
parameter:

- **`DX_FLOOR = 1e-5` is measured, not chosen.** Check (b) now records its worst *absolute*
  FD discrepancy (`8.5e-08` for `∂C/∂x` on the default grid); `1e-5` sits `118×` above it,
  and the run says so and complains if the margin drops below 10×. A cell below it has a
  **sign the numerics do not resolve**, so it is drawn grey rather than as either sign — 1 of
  6163 on the default grid, legend entry only if non-empty. It no longer appears anywhere on
  the figure, which was the complaint that ended `:split_log`.
- **The contours cost nothing.** They are interpolated from `DX`, the grid check (b) has
  already scored against §10.1 — not re-scanned. A per-row `dC_dx_fd` scan would be more
  accurate and cost ~17 s (180 rows × 400 points × 2 `run_adia` calls), roughly doubling the
  script; the `c = 0` cross-check below measures what the interpolation costs instead, and it
  is `4.1e-05` worst against a grid spacing of `2.6e-03`.

**Its zero-locus is the payload**: the ridge `x_opt(u)`, the best pump at each arm
asymmetry, bisected per u row. Two things make it evidence rather than decoration:

- **Its bottom endpoint is known in advance.** §10.7 proves `∂C/∂x|_{u=2}` changes sign
  exactly at `x*`, and the ridge lands there to **3.6e-11** — not a fit.
- **It reproduces three stored numbers and extends them to a curve.** `files_online/06` §6
  gives `x_opt` at three asymmetries only; the ridge passes through all of them:

  | β | 1 | 2 | 10 |
  |---|---|---|---|
  | ridge here | 0.148550 | 0.106837 | 0.042379 |
  | stored | 0.1486 | 0.107 | 0.042 |

**How far it extends — check (h), and the figure hides the answer.** The ridge has a
**genuine lower terminus at `(x*, u = 2)`** — a real endpoint, because `u = β + 1/β ≥ 2` is
a hard domain boundary (β = 1), not because the curve stops. It has **no upper terminus**:
it runs to `u → ∞` with `x_opt → 0`, staying strictly inside the entangled region the whole
way. The asymptotics are exact and clean:

```
x_opt · √u  → 1/6        x_max · √u  → 1/2        x_opt / x_max  → 1/3
   (0.16666646 at u = 1e12)   (0.49999950)            (0.33333326)
```

So `x_opt → 1/(6√u)` while the death boundary goes as `1/(2√u)`: **at extreme arm
asymmetry the optimal pump sits at exactly one third of the way to the boundary.**

**In the figure the ridge exits the LEFT edge, not the top** — it crosses `x = x_lo = 0.02`
at `u = 52.6` (β = 52.6), well below the frame top of `u = 664`, so the upper ~44% of the
plotted u range shows no ridge at all. That is framing, not absence. Lowering `x_lo` shows
more of it, at the cost of a much taller u axis (`u_c ~ 1/(4x²)`).

**Figure 2 carries a second panel too, and it needs no colorbar to be believed**: `∂C/∂x`
against `x`, **untransformed, linear, in its own units**, for `u ∈ {2, 2.5, 4, 10.1}`, each
curve crossing a drawn zero line exactly once. That is check (g) as a picture rather than as
a count, and it is the half of the figure that survives any doubt about what a colour means.
It shares the map's `xlims`, so a vertical dropped from a ridge marker above lands on the
crossing below. `u = 2.5` and `10.1` are β = 2 and β = 10 exactly — check (e)'s stored-`x_opt`
rows — and every crossing is printed against `ridge_x(u)`:

```
u = 2    (β = 1)    crossing 0.148550   ridge_x 0.148550   Δ = 2.8e-07
u = 2.5  (β = 2)    crossing 0.106839   ridge_x 0.106837   Δ = 2.2e-06
u = 4              crossing 0.071443   ridge_x 0.071439   Δ = 3.4e-06
u = 10.1 (β = 10)   crossing 0.042383   ridge_x 0.042379   Δ = 3.5e-06
```

Δ is the linear-interpolation error on a 400-point x line, not a disagreement. All four
curves terminate by leaving the entangled region rather than at the frame edge, which the
run states per curve — frame containment, audit 2. **The contour levels are drawn into this
panel as dotted horizontal lines in the same colours**, which is what ties the two panels
together: a level is a curve upstairs and a line downstairs carrying the same number, so
either can be read against the other and neither is a colour scale.

**Three checks specific to the contours**, all printed, all able to fail:

- **Extraction.** The `c = 0` contour, pulled from the grid by the *same* interpolation every
  drawn level goes through, against `ridge_x`'s independent 120-step bisection: worst
  `4.1e-05`, median `1.3e-05`, over 101 rows, against a grid spacing of `2.6e-03`. This is
  the one level whose answer is known to `3.6e-11`, so it validates the extraction — and on a
  square grid nothing else would catch a transposed row index.
- **Crossings per row.** A level crossed twice in one row would be joined into a curve that
  does not exist, so the count is printed and every contour is drawn as **markers, not a
  joined line**. Check (g) only settles the `c = 0` level; monotonicity of `∂C/∂x` in x at
  the other levels is *not* established, so it is measured rather than assumed. Currently 0
  multi-crossing rows at every level.
- **Frame containment.** Six of the seven levels **exit the left edge** — they run to
  `x → 0` where `∂C/∂x → ±∞` — and only `c = −1.0` ends inside the frame. That is stated per
  level on stdout rather than left to be inferred from curves that stop.

### Check (a2) — the large-u law, and the trap in importing it

`∂C/∂u = A₀(x) u^{−3/2} + A₁(x) u^{−2} + O(u^{−5/2})`, with `A₀ = −√x(1−x)³/p₁` and
`A₁ = −2x^{3/2}(1−x)³V₂ / ((1+x)√q₁ p₁²)`, `V₂ = x⁷−7x⁶+19x⁵−21x⁴+7x³+23x²+5x+5`. Source is
a session note comparing `∂C/∂u` at `u = 2` and `u = 500`.

Two things it buys. **`V₂` is §10.3's polynomial, already Sturm-proven positive on (0,1)**,
and `p₁, q₁ > 0` there, so `A₀ < 0` and `A₁ < 0` — **an independent, cheap proof of
`∂C/∂u < 0` at large u**, by a different route from §10 (it says nothing near `u = 2`, which
is what (c) and §10.5 cover). And it **predicts the window corner**: the weakest sampled
value comes out of `A₀(X_LO)·u_c(X_LO)^{−3/2}` with no grid at all — `−9.8733e-06` against a
measured `−9.9117e-06`, **0.39%** — which is what demotes that number from a result to a
coordinate. Verified at a second, unrelated corner too (`X_LO = 0.001`, `u = 122638`: 0.90%).

**⚠ The law is about the UNCLAMPED bracket `C̃`, not about `C`.** The physical concurrence is
`C = 2·max{0, ·}`, so `∂C/∂u ≡ 0` beyond the death boundary `x_d(u)`, while `C̃`'s derivative
there is analytic continuation — a number, but not a slope of anything. `x_d(500) = 0.02144`,
so the note's own `x = 0.05, 0.1486, 0.4838` rows are **all outside**, and its §1 exists to
separate the two objects. **Those three x values were lifted straight into this check on the
first attempt**, into a script whose grid is the clamped field — the one confusion the note
leads with. Every x is now asserted `x < x_d(u)`. (`x_max_of_u(500)` reproduces the note's
60-digit `x_d` to `1.9e-14`, a free cross-check.)

**⚠ The `u^{−1/2}` correction is not observable inside the region.** A *fixed* x leaves the
region as u grows (`x_d ~ 1/(2√u)`), so the scan must hold `σ = x√u` fixed instead — a
constant fraction of the way to the boundary (`σ → 1/2` is the boundary, `σ = 1/6` the
ridge). But then `x ~ σ/√u`, so `A₁/A₀ ~ 10x ~ 10σ/√u` and the `A₁` term is `O(u^{−1})` —
the **same order as the term after it**. Measured at `σ = 0.35`, both 1-term and 2-term
residual ratios converge to 4 over `u × 4`, and the 2-term is no better than the 1-term at
physical x. So the note's clean `u^{−1/2}` scaling is a fixed-x statement, and fixed x at
large u is outside the region — which is why its demonstration used those three x values.
A first attempt at a fixed `x = 0.005` gave a 1-term ratio of **10.86**: there the `A₁` term
(`2.24e-3`) and the next order (`2.0e-3`) are the same size at `u = 500`, so the residual was
a cancellation between comparable terms and its ratio meant nothing.

Other checks: **(b)** FD-vs-§10.1 worst `1.2e-06` / median `7.3e-09` for `∂C/∂u`,
`2.9e-06` / `8.3e-10` for `∂C/∂x` — two routes, no shared code; **(c)** §10.5's closed form
for the `u = 2` row to `~1e-10`; **(f)** `∂C/∂x · sqrt(x(u+2)) → 1` as `x → 0`, the
divergence that forbids a global sign; **(g)** `∂C/∂x` has **exactly one sign change per u
row** — 1 in all 180 — so `C_adia` is unimodal in x at fixed u and the bisected ridge is
the *whole* zero-locus, not one root of several.

**Physical reading of figure 1**, worth having in words: `∂C/∂u|_{u=2}` runs monotonically
**−0.017 at x = 0.02 to −0.94 at x_c**. Since `C''(β=1) = 2·∂C/∂u|_{u=2}`, **the β = 1 peak
sharpens as the entanglement dies** — the state is most fragile to arm asymmetry exactly
where it is weakest.

### The `∂C/∂x` panel used to show a zero contour that does not exist

Worth recording in full, because the figure was wrong in a way that survived being looked
at, and because it is issue 2's family rather than a one-off.

On a **linear** symmetric scale the panel appears to have **two** zero contours: the ridge,
and a second edge running right along the `u = 2` axis. There is only one. `∂C/∂x` has
exactly one sign change per u row — check (g), 180 of 180.

The cause is that **white in a diverging colormap means "small", not "zero" — and how small
is set by an extreme value somewhere else in the panel entirely.** `∂C/∂x → +∞` as `x → 0`
(§10.7), so `|clims|` came out 2.80 from that corner, and along the bottom edge:

| x | 0.1486 | 0.160 | 0.180 | 0.200 | 0.250 |
|---|---|---|---|---|---|
| `∂C/∂x` | **0.0000** | −0.0837 | −0.2152 | −0.3310 | −0.5675 |
| % of range | 0.0% | 3.0% | 7.7% | 11.8% | 20.3% |

Only the first is a zero; the next two render the same white as it, and the strip's far
edge reads as a contour. **`sign(v)·sqrt(|v|)` fixed it at the source** — compressing the
divergent tail and expanding the near-zero range moves `x = 0.18` from 7.7% to 28% of the
bar, and the white collapses onto the actual zero. A 99th-percentile clip was tried first
and reverted: it helps the bulk but saturates a wide band along the death boundary into
flat dark blue, trading one artifact for another.

**That fix was right about the failure and wrong about the method, and it has been replaced
twice since.** `sign(v)·sqrt(|v|)` cures the fake contour by picking an exponent that happens
to work — nothing forces `1/2`, and the `|v|` in it is exactly what figure 1 refuses to
write, so the panel would have rendered identically had a sign been wrong. `:split_log` fixed
that structurally (empty neutral band) but at the cost of a colorbar nobody can read.
`:sign_map`, the shipped default, ends the sequence by **removing the colour scale**: sign in
the colour, magnitude on labelled contours in real units.

The audit numbers, all still re-producible by one edit to `DX_SCALE`:

| `DX_SCALE` | misleading cells | worst value painted as zero |
|---|---|---|
| `:linear` | **fires** (3 of 61 even on a coarse 20×15 grid) | 0.138 |
| `:signed_sqrt` | 0 | 0.000 |
| `:split_log` | 0, and `|z| < 1` provably empty | 0.000 |
| `:sign_map` | N/A — no colour scale to mislead | — |

The thresholds are independent as the convention demands — perceptual is "within the central
5% of the drawn `clims`", physical is "`|∂C/∂x|` above 1% of the panel maximum", the latter
read off the untransformed data the colour map never touches. **The `:linear` branch is kept
in the file for this reason**: a check never shown to fail on a known-bad panel is not a
check, and this one is re-runnable by one edit.

**The general lesson, which outlived three of the four encodings:** a colour scale can report
a magnitude or it can report a sign, and making one axis do both costs either an `abs()` or a
label nobody reads. Where the sign is the payload — as it is here, since the zero-locus *is*
`x_opt(u)` — giving the colour to the sign and the magnitude to labelled contours is cheaper
than any transform, and leaves no number on the figure that has to be decoded.

An earlier draft of this section claimed the full linear range "costs some mid-tone
contrast and hides nothing". **That was wrong** — it hid a nonzero region by painting it
the colour of zero, which is precisely the failure issue 2 describes for
`concurrence_diff`. The general rule: on a diverging map whose data has a divergence or a
long tail, a linear scale will manufacture a fake zero set, and the fix is a transform, not
a clip.

### Two traps, both of which return plausible wrong numbers rather than errors

- **`calc_concurrence` clamps at exactly `0.0`.** A centred difference straddling `u_c`
  differentiates *the clamp*, not the model, and yields a slope that looks fine. Guard: a
  point is evaluated only when the centre **and both legs** are strictly inside the region;
  otherwise `NaN`, which `flatten_grid` drops. Only 6163 of 39600 grid points survive — the
  rest are the dead zone, consistent with the boundary figure's 84%.
  **`u = 2` is a hard edge** (no real β below it), so that row gets a second-order *forward*
  difference instead. It must not be dropped: it is the row §10.5 is about.
- **On a panel with an active `marker_z` colorbar, GR does not give you the colour you ask
  for.** Measured on `dC_dx.pdf`, all through the same PDF pipeline: `lc = :limegreen` and
  `lc = RGB(0,0.7,0)` both rendered **pale blue-grey**; `lc = :magenta` rendered **dark
  red**; a `:lime` marker rendered **cyan**. `:yellow` markers and `lc = :black` survived.
  The first four are all colours present in that panel's own `:balance` colormap, which is
  suggestive, but **the rule is not established** — the surviving yellow argues against any
  simple "everything is quantized" story. What is established: the same keywords render
  correctly on a panel *without* `marker_z`, and Plots reports the attribute as exactly
  what was requested (`linecolor=RGBA(0,0.7,0,1)`), so **nothing errors and nothing warns**.
  Practical rule: on a `marker_z` panel, pick overlay colours **by looking at the output**.
  This is the colour version of the standing advice to eyeball GR rather than assume.

  **`:sign_map` retired this trap on `dC_dx.pdf` by removing `marker_z` from that panel**,
  which also means **the overlay colours there are no longer the ones that were approved.**
  They were chosen by eye *against the mangling* — `:lime` was picked because it rendered
  cyan — so a keyword that now renders truthfully renders *differently*. The ridge was moved
  to `:black` accordingly, and the contour colours are explicit `RGB` triples, but none of it
  has been checked against the output. The trap still stands for `dC_du.pdf` and for the
  three colour-mapped `DX_SCALE` branches, all of which still carry a `marker_z` colorbar.

### `run_adia` has a domain limit, and it is not where you would guess

Mapped by check (h). `steadystate.eigenvector` throws *"Eigenvalue with smallest absolute
value is not zero"* when **large β is combined with large x** — the deep dead zone, far
from any entangled state:

| u ≈ β \ x | 1e-6 | 1e-4 | 1e-2 | 0.1 | 0.3 | ridge sits at |
|---|---|---|---|---|---|---|
| 1e6 | ok | ok | ok | ok | ok | 1.7e-04 |
| 1e8 | ok | ok | **FAIL** | **FAIL** | **FAIL** | 1.7e-05 |
| 1e10 | ok | **FAIL** | **FAIL** | **FAIL** | **FAIL** | 1.7e-06 |

**The ridge is never in the failing region** — `x_opt ~ 1/(6√u)` tracks the survivable band
at every u tested, which is why check (h)'s asymptotics could have used `run_adia` and use
the closed form only because it is exact and free out there. It fails *loudly*, which is
the good kind; the same cannot be said of `bisect`, which returns a plausible number when
handed a bracket with no sign change. That combination cost a debugging cycle: `x_max_of_u`
had its left bracket at `1e-6`, too high to bracket `x_max ~ 1/(2√u)` past `u ~ 1e12`, so it
silently returned `X_C` — six orders too large — and drove the ridge search into `x ≈ 0.48
at β = 1e8`, killing the script. **The bug was the bracket, not the solver.** `x_max_of_u`
now checks for the sign change explicitly and returns `NaN` rather than guessing.

`C_of` is a third copy, after `adia_concurrence_max.jl` and `adia_boundary.jl`, for the
same no-`main()`-guard reason. `sym_clims` is inlined rather than lifted — see issue 2,
which this makes the fourth site for.

## Conventions worth not breaking

- **Grids are `[row, column]` with the row being the horizontal/log axis** — `[N_β, N_x]`
  for sweep 1, `[N_r, N_ξ]` for sweep 2, `[N_ξ, N_x]` for sweep 3. Every report,
  observable matrix and plot assumes it. Sweeps 1 and 3 both carry an x axis and both put
  it in the *column*, but their rows are different coordinates (β vs ξ) — so comparing one
  sweep's grid against another's means naming both indices explicitly, never reusing an
  index. That is exactly the bug `check_against_run_A` shipped with on its first run.
  Getting it backwards cannot be caught by an assertion: `len = 7` makes every grid square.
- **NaN means "no data"**, everywhere. `nanmat()` prefills the state grids; observable
  grids prefill with `NaN`. A failed point becomes a gap in the heatmap instead of
  killing the run. Don't replace with zeros.
- **A figure gets audited like a number.** Colour honesty, frame containment,
  label–encoding agreement — the three checks under "Figure audits" in Plotting, each with
  the value it returned on a version that was actually wrong. Two rules carry the weight:
  **never state a negative about a figure** ("hides nothing", "shows all the structure")
  without a test that could fail, and **confirm a new check fires on the known-bad case
  before trusting it** — a check whose thresholds are not independent passes everything.
  This convention exists because a panel here showed a zero contour that does not exist,
  through several sessions of being looked at.
- **`res.params` is the truth for plot axes**, not `config.jl` — the config may have been
  edited after the run that produced a given `.jld2`.
- **The report block in `run_sweep` is wrapped in `try/catch` on purpose.** It sits before
  the `return`, and a bug in a summary line once discarded a completed sweep. Keep any new
  reporting inside that block.
- `.gitignore` is deny-by-default (`*` then `!*.jl`). Results, figures and `Manifest.toml`
  are not tracked. A new non-`.jl` file needs an explicit `!` line. **Git does not descend
  into an ignored directory**, so the `!*.jl` exception never reaches subdirectories —
  `files_online/` is entirely untracked and is *not* on GitHub, notes and `.py` scripts
  included. **`file_concurrence/` is the one exception that does reach into a
  subdirectory**, and it takes *two* lines to do it:

  ```gitignore
  !file_concurrence/
  !file_concurrence/*.md
  ```

  The directory must be re-included first. A bare `!file_concurrence/concurrence_adia.md`
  is **inert** — git never descends far enough to consider it — and it fails silently, so
  test with `git add --dry-run`, not with `git check-ignore`. (`check-ignore -v` exits 0
  here and prints the *negation* pattern; that is a match report, not an ignore verdict,
  and it reads exactly like a failure.) Everything tracked is the top-level `.jl`,
  `CLAUDE.md`, `Project.toml`, `.gitignore` and `file_concurrence/concurrence_adia.md` —
  20 files, 22 once `adia_boundary.jl` and `adia_derivatives.jl` are added
  (`git ls-files | wc -l` is the answer, not this sentence; it has drifted before).
- **`origin` is `https://github.com/Yash-27/JC_TMS_ADIA.git`** — this working directory
  *is* that repo, so there is nothing to clone into a subfolder. GitHub's tree is
  `.gitignore`, `CLAUDE.md`, `Project.toml`, the top-level `.jl`, and
  `file_concurrence/concurrence_adia.md`. No `.py`, no `files_online/`, no `.jld2` — and
  `origin/main` is an ancestor of local HEAD, so everything there is already here. **The
  `.py` files and the second `claude.md` are in `files_online/`, which is local-only**; if
  you expect to find them on GitHub, they are not there yet.

## Session start — `.claude/session-status.sh`

A `SessionStart` hook in `.claude/settings.json` runs it before anything else, and its
stdout becomes Claude's opening context. Three sections:

1. **GIT** — `git fetch`, then ahead/behind vs upstream (falling back to `origin/main`),
   a specific flag when the working copy's `CLAUDE.md` has drifted from the pushed one,
   and the uncommitted-file list capped at 20.
2. **CONTENT** — `shasum` inventory of the prose directories named in `NOTES_DIRS`
   (currently `files_online file_concurrence`), diffed against
   `.claude/.content-manifest` from the previous session, reported as ADDED / REMOVED /
   MODIFIED. **Adding a directory means editing that one variable and nothing else**: the
   manifest stores full relative paths, so a new root's files simply appear as ADDED on
   the next run. Missing directories are filtered out before the `find` (which errors on
   a missing root) and reported on a `not present:` line; all four combinations of the two
   directories are verified to exit 0.
3. **READ LIST** — names the changed files and instructs Claude to open them, plus the
   reminder that `files_online/claude.md` governs that project and inverts `ξ`, while
   `file_concurrence/` has no `claude.md`, shares this repo's `x` and `β`, and uses `x*`
   for the argmax (the sense this file now uses throughout).

Two reasons both directories are hashed even though only one is git-invisible: this
section reports **content drift since the last session**, which is a different question
from git's "is it committed" — a file edited *and committed* between sessions is clean to
git and still unread by Claude.

**It reports changes rather than pasting content, deliberately.** `files_online/` is
~138 KB (~35k tokens), `file_concurrence/` another ~56 KB (~14k), and `CLAUDE.md` a
further 57 KB the harness already loads every session; dumping all of it unconditionally
would spend most of a context window re-reading unchanged text. Steady-state cost is ~230
tokens. After an edit, Claude is pointed straight at what moved and reads it in full.

**`fetch` only moves remote-tracking refs** — never merges, rebases, or touches the
working tree, so the hook cannot lose an edit; pulling is left to you. It never exits
non-zero, so no network / no remote / neither prose directory still starts a session. The
only file it writes is its own manifest. Verified against all four content paths (first
run, unchanged, modified, added+removed) and against all four directory-presence
combinations. `.claude/` is itself ignored by the deny-by-default rule, so none of this is
tracked or pushed.

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
`cmap = :balance`). That's the pattern to lift into `plotting_functions.jl` when this gets
fixed properly.

**That pattern is right about the colormap and UNAUDITED about the scale.** An earlier
version of this paragraph said the two were "verified to render correctly through GR",
which was an eyeball recorded as a verification — the thing "Figure audits" now forbids.
What is actually known: `sym_clims` puts the neutral midpoint on zero, which is the part
that matters for reading a sign, and that much is structural. What is **not** known is
whether those panels have the fake-zero-contour problem, because `sym_clims` on a *linear*
scale is exactly the configuration that produced one on `∂C/∂x` — a long tail anywhere in
the grid inflates the near-white band until materially nonzero cells are painted as zero.
`ΔC_full` reaches ~0.06 against `D` grids of ~0.16, so a tail is plausible and unmeasured.
**Run audit 1 from "Figure audits" on the five `compare_runs.jl` panels before quoting them
in anything.** No claim either way until then.

**There are now four sites wanting `sym_clims`** — `compare_runs.jl` (which defines it),
`run_plot_r_xi.jl` and `run_plot_x_xi.jl` (via their duplicated `overrides`), and
`adia_derivatives.jl`, which inlines the three lines for its `∂C/∂x` panel. Each new signed
panel raises the cost of not lifting it. The reason it keeps not happening is the same one
in issue 8: it is a refactor of files that currently work, and the repo's practice is to
duplicate rather than risk them.

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

**Fixed in `run_x_xi.jl` only** (sweep 3), exactly that way: a `Base.Semaphore(MAX_BIG)`
acquired around the `run_sim` call for `dim ≥ BIG_DIM` points, in a `try/finally` so a
thrown solve cannot leak a permit and deadlock the tail. The queue preference is kept.
That driver had to fix it — its big points measured **5.04 GiB peak RSS each** at dim 4096,
so on a 16 GiB machine the unenforced version is an OOM rather than a slowdown. **Sweeps 1
and 2 still have the bug**; the fix is 6 lines and can be lifted from `run_x_xi.jl`.

### 8. The scheduler now exists in THREE near-identical copies

`run_sweep_rxi` is a near-copy of `run_sweep`, and `run_sweep_xxi` is a near-copy of
both: same two-pass structure, same atomic big/small work queues, same `take_job`, same
ETA fit, same `try/catch` report block — roughly 200 lines each. Only a few things
genuinely differ per copy: the coordinate inversion, the cost model's gap factor, the
printf labels, and (in sweep 3 only) the semaphore, the `flush(stdout)` calls and the
checkpoint write.

They were written that way deliberately, to keep zero risk to a working sweep, and the
observables were shared instead (`observables.jl`) because those *could* be lifted
without touching the scheduler. But the consequence is real and now threefold: **a fix to
the scheduler in one file reaches neither of the others.** The current concrete evidence:
issue 6 is fixed in `run_x_xi.jl` and still present in the other two, and `flush(stdout)`
plus incremental checkpointing — both responses to recorded data loss — exist only in the
newest copy.

The clean resolution is a `sweep_core.jl` taking `point_params(i,j) -> (κ_1,κ_2,g_1,g_2,η)`
and `cost_factor(i,j)` as closures, with all three drivers calling it. That is a real
refactor of files that currently work, so it should be done deliberately rather than
opportunistically — but note the verification is cheap and already exists: each sweep has a
stored `.jld2` plus its own adiabatic self-test, and sweep 3's `ξ = 1` row must reproduce
run A, so a migrated scheduler that is wrong will say so.

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
