# CLAUDE.md

Julia code for a cavity-QED reservoir-engineering study: two qubits, each dispersively
coupled to its own cavity, with the two cavities linked by a two-mode-squeezing (TMS)
drive. The question the code answers is how well the **adiabatically-eliminated 2-qubit
model** reproduces the **full 4-body model** across a 2D parameter sweep.

## Running

```bash
julia --project=. -t auto run_jc_pump_disp_asy.jl
```

```bash
julia --project=. run_plot.jl
```

```bash
julia --project=. compare_runs.jl
```

`run_jc_pump_disp_asy.jl` has a top-level `begin` block at the bottom that executes the
whole sweep on include. There is no `main()` guard — `include`ing it *runs* it.

Threads matter a lot: the sweep is thread-parallel and `-t auto` (or `-t N`) is the
difference between minutes and hours. `BLAS.set_num_threads(1)` is set deliberately at
include time to stop BLAS from oversubscribing inside the `@threads` loop — don't remove it.

## Pipeline

```
config_jc_pump_disp_asy.jl   parameters, observable selection, output filename
        │
        ├─ truncation.jl     estimate_truncation(κ_1,κ_2,η,g_1,g_2) -> (N_1,N_2)
        ├─ jc_pump_disp_asy.jl  run_sim(...) -> (ρ_full_reduced, ρ_adia)
        │
run_jc_pump_disp_asy.jl      run_sweep + observables + @save to .jld2
        │
        ├─ run_plot.jl        load ONE .jld2, draw heatmaps
        └─ compare_runs.jl    load TWO .jld2, difference them
                 │
        plotting_functions.jl   plot_map / save_fig, included by both
```

Include order is fixed by `run_jc_pump_disp_asy.jl` lines 8–10. `config` has no
dependencies by design and can sit anywhere in that order.

## Dimensionless coordinates

The sweep is over two dimensionless numbers, not raw parameters:

| symbol | definition | meaning | sweep |
|---|---|---|---|
| `x` | `4η² / (κ_1 κ_2)` | drive strength relative to the linear OPO threshold | linear |
| `β` | `(g_1² κ_2) / (g_2² κ_1)` | coupling asymmetry between the two arms | log10 |

`run_sweep` converts back via `η_of(x) = sqrt(x κ_1 κ_2 / 4)` and
`g1_of(β) = g_2 sqrt(β κ_1 / κ_2)`. `κ_1`, `κ_2`, `g_2` stay fixed for a run.

**`x < 1` is a hard constraint.** The adiabatic model has `4η² − κ_1κ_2` in every
denominator, so it diverges at `x = 1`. The current `x_range` tops out at `0.7225 = 0.85²`.
Anything approaching 1 is also where the Liouvillian gap closes and the iterative solve
gets slow (see the cost model).

## Adding an observable

Three edits, nothing else:

1. Write `calc_foo(ρ)` (or `calc_foo(ρ_full, ρ_adia)`) in Section 2 of
   `run_jc_pump_disp_asy.jl`.
2. Add one `Obs(:foo, :single|:compare, calc_foo, "Label")` row to `OBSERVABLE_REGISTRY`.
3. Name `:foo` in `ACTIVE_OBSERVABLES` in the config.

`:single` auto-expands into `foo_full`, `foo_adia`, `foo_diff` (= full − adia).
`:compare` produces one grid under `:foo`. `select_observables` errors loudly on an
unknown name — it is called *before* the sweep so a typo fails in seconds, not hours.

Observables are pure functions of the stored 4×4 matrices. Section 2 never reaches back
into Section 1 except through the returned `res`, so it can be lifted into its own file
and re-run against a saved `.jld2` without repeating the sweep.

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

Two results files are currently on disk, both `len = 7`, same `x_range` and `β_range`:

| file | κ_1 | κ_2 | g_2 |
|---|---|---|---|
| `results_k1_2.0_k2_2.0_g2_1.0.jld2` | 2.0 | 2.0 | 1.0 |
| `results_k1_2.5_k2_1.5_g2_1.0.jld2` | 2.5 | 1.5 | 1.0 |

The config sits at `κ_1 = 2.5, κ_2 = 1.5`, so `run_plot.jl` currently resolves to the
second. `compare_runs.jl` uses both.

Hardcoded on purpose in `plot_map`, because there is only one right answer:

- β horizontal on log10, x vertical linear
- the flattening convention — `vec` is column-major, so the **row** index (β) varies
  fastest; `outer` on β and `inner` on x reproduces that. Swapping them silently
  transposes the figure with no error, and the length assertions won't catch it on a
  square grid.
- NaN points dropped rather than plotted
- colorbar forced on (GR does not infer it from `marker_z`)

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

Figures go to `Full_vs_Adia/` (what `run_plot.jl` passes) and `Full_vs_Full/` (what
`compare_runs.jl` passes), PDF only. `save_fig`'s own default is `figures3/`, and
`figures/`, `figures2/` also exist from earlier runs — five destinations, none tracked by
git.

## Comparing two runs — `compare_runs.jl`

The second entry point. It loads **two** `.jld2` files and compares each model against
*itself* across parameter sets, which the single-run figures structurally cannot do.

Set `RUN_A` and `RUN_B` at the top to the `(κ_1, κ_2, g_2)` triples that name the files.
Four figures, all on the same (β, x) mesh:

| grid key | quantity | signed? |
|---|---|---|
| `:full` | `½‖ρ_full^A − ρ_full^B‖₁` | no — a norm |
| `:adia` | `½‖ρ_adia^A − ρ_adia^B‖₁` | no — a norm |
| `:concurrence_full` | `C(ρ_full^A) − C(ρ_full^B)` | **yes** |
| `:concurrence_adia` | `C(ρ_adia^A) − C(ρ_adia^B)` | **yes** |

This is **not** the in-file `:tracedist` / `:concurrence_diff`, which are full-vs-adia
*within* one run. These never cross the two models; they cross the two *parameter sets*.

`MODELS` declares `:full` and `:adia` once and drives both figure families. The
trace-distance loop uses the whole row; the concurrence loop takes only `key` and `word`
(`conc_key(m) = Symbol(:concurrence_, m.key)` — which works because the sweep's `:single`
expansion builds those names the same way). Drop a row to skip both its panels.

Things that are load-bearing:

- **The two axes hold (β, x) fixed, not the raw parameters.** At matched (β, x) the two
  runs sit at different η and different `g_1`. The figures answer "how much does the
  state depend on the κ asymmetry at fixed dimensionless coordinates", not "at fixed
  drive".
- **Signed differences get `sym_clims` + `:balance`**, not viridis on a data range. This
  is issue 2's prescription, applied here rather than inherited — on an asymmetric range
  a diverging colormap's neutral midpoint lands on some arbitrary nonzero value and the
  sign becomes unreadable.
- **`clims` is computed from the `si_scale`d grid, not the raw one.** `plot_map` applies
  `clims` to what it is handed; limits taken before scaling are off by that power of ten
  and clamp everything to one end colour.
- **Every panel takes its own colour range.** Deliberate — see the finding below for why
  a shared range would hide the main result. Compare magnitudes via the printed summary,
  never by eye across panels.
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
  `D_adia` the correction would be the same fixed local σ_z on qubit 1 in both runs, and
  trace distance is invariant under a unitary applied to both arguments. `ΔC_adia` is
  safer still — concurrence is a local-unitary invariant, so each operand is individually
  unchanged and the argument doesn't need the two runs to share the rotation.

### Finding: the adiabatic steady state is a function of (β, x) alone

For the pair currently configured, `ΔC_adia` is **identically zero to machine precision**
— max |ΔC_adia| ≈ 1.8e-14 against concurrence values of order 0.24, structureless, signs
scattered. The stored maxima agree to 14 digits:

```
concurrence_adia   run A 0.2360436749507588   run B 0.23604367495075904
concurrence_full   run A 0.2866               run B 0.2598
```

So changing κ_1, κ_2 at fixed (β, x) does not move the adiabatic state at all: the
elimination has absorbed the raw rates completely into the two dimensionless coordinates.
The full model has **not** — it shifts by up to 0.060 in concurrence, ~20% of peak.

That gap is the residual cavity dependence the elimination discards, and it is a
statement about the approximation itself rather than about either run. `D_adia` (2e-14)
was already saying this, but it reads as "barely responds"; the concurrence panel makes
it exact.

This is also why the panels are not on a shared colour scale — on a range set by
`ΔC_full`, the entire adiabatic panel would be uniform white, and *identically zero*
would be indistinguishable from *merely small*.

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
  are not tracked. A new non-`.jl` file needs an explicit `!` line.

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
κ_1, κ_2 separately, there is no need to scan for a parameter regime where the
discrepancy shows up; any single (β, x) is representative.

Note also that none of the four `compare_runs.jl` panels can detect this, by construction
— that is the point of the immunity argument in that section. It has to be tested against
the in-file `:tracedist`.

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
- `print("\r\x1B[K")` in the solve-failure catch is a leftover from a removed progress bar.
- `load_results` pulls `status_grid` out of the `.jld2` and **nothing on the single-run
  path ever reads it**. It's the one thing that would let a figure distinguish a failed
  solve from a failed observable — right now both render as an absent marker.
  `compare_runs.jl` does read it (via its own `load_states`, not `load_results`), so the
  pattern to copy exists.
- Two headers name the wrong file: `plotting_functions.jl` calls itself `plotting_lib.jl`
  and says it's included by `plot_try.jl`; `run_plot.jl` calls itself `plot_try.jl`.
  (`plot_try.jl` was a byte-identical copy of `plotting_functions.jl` plus a stale driver,
  and has been deleted.)
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
