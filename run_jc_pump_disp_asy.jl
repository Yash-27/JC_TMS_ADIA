using Base.Threads
using Printf
using LinearAlgebra

include("truncation.jl")
include("jc_pump_disp_asy.jl")

# ==============================================================================
# sweep_driver.jl
#
# Generates a dimensionless (β, x) parameter grid, converts each point to
# physical parameters, calls estimate_truncation to get (N_1, N_2), then
# calls run_sim. Collects the raw (rho_full, rho_adia) results into a grid.
#
# Multithreaded over the grid via @threads, with a progress bar.
#
# No diagnostics (concurrence, trace distance) and no saving yet -- this
# is just "generate parameters -> truncate -> solve", nothing more.
#
# Dimensionless parameters:
#   x = 4η² / (κ_1 κ_2)              -- drive strength relative to threshold
#   β = (g_1² κ_2) / (g_2² κ_1)      -- coupling asymmetry between the arms
# ==============================================================================

# Disable BLAS multithreading to prevent CPU oversubscription during our
# own @threads loop.
BLAS.set_num_threads(1)

# --- Fixed physical parameters (symmetric decay rates, reference coupling) ---
const κ_1 = 2.0
const κ_2 = 2.0
const g_2 = 1.0
const γ_val     = 0.0
const γ_phi_val = 0.0

# --- Dimensionless grid ---
const len = 7
x_vals = collect(range(0.01, 0.7225, length = len))                       # x = 4η²/(κ1κ2)
β_vals = collect(exp10.(range(log10(0.01), log10(100), length = len)))    # β = g1²κ2/(g2²κ1)

N_β = length(β_vals)
N_x = length(x_vals)

# --- Result storage ---
rho_full_grid = Array{Any}(undef, N_β, N_x)
rho_adia_grid = Array{Any}(undef, N_β, N_x)

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

    # Truncation for this point
    N_1, N_2 = estimate_truncation(κ_1, κ_2, η, g_1, g_2)

    # Solve
    rho_full, rho_adia = run_sim(g_1, g_2, κ_1, κ_2, η, γ_val, γ_phi_val, N_1, N_2)

    rho_full_grid[i, j] = rho_full
    rho_adia_grid[i, j] = rho_adia

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

println("\n\nSweep complete. Results stored in rho_full_grid, rho_adia_grid.")