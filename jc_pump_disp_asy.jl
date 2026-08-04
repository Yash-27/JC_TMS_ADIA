using QuantumOptics

# ==============================================================================
# sim_core.jl
#
# Core physics engine only. Two qubits, each coupled to its own cavity, with
# the cavities linked by a two-mode-squeezing drive (η). This is the
# "reservoir engineering" setup where adiabatically eliminating the cavities
# gives an effective dissipative interaction between the two qubits alone.
#
# This file defines the physical system and solves it -- nothing else.
# Truncation choices, diagnostics (concurrence, trace distance), sweeps, and
# saving all belong in separate files that build on top of this.
# ==============================================================================

# ------------------------------------------------------------------------
# run_sim: given physical parameters and a chosen Fock-space truncation
# (N_1, N_2), solves both the full 4-body model and the adiabatically-
# eliminated 2-qubit effective model, and returns their steady-state
# density matrices.
#
# Arguments:
#   g_1, g_2   : qubit-cavity coupling strengths
#   κ_1, κ_2   : cavity decay rates
#   η          : two-mode-squeezing drive strength linking the cavities
#   γ          : intrinsic qubit decay rate
#   γ_phi      : intrinsic qubit dephasing rate
#   N_1, N_2   : Fock-space truncation for cavity 1 and cavity 2
#
# Returns:
#   (rho_ss_full, rho_ss_adia)
#     rho_ss_full : full model's steady state, untraced (cavity_1 ⊗
#                   cavity_2 ⊗ qubit ⊗ qubit) -- kept untraced so
#                   truncation-convergence checks can see artifacts in
#                   the cavity part of the state, not just the reduced
#                   qubit state. Its dimension depends on N_1, N_2, so
#                   comparing two such states at different truncations
#                   requires padding to a common Fock cutoff first.
#     rho_ss_adia : adiabatic (cavity-eliminated) model's 2-qubit
#                   steady state
# ------------------------------------------------------------------------
function run_sim(g_1::Float64, g_2::Float64, κ_1::Float64, κ_2::Float64,
                  η::Float64, γ::Float64, γ_phi::Float64,
                  N_1::Int, N_2::Int)

    # Bases
    cavity_1 = FockBasis(N_1)
    cavity_2 = FockBasis(N_2)
    qubit = SpinBasis(1 // 2)

    basis_full = tensor(cavity_1, cavity_2, qubit, qubit)
    basis_adia = tensor(qubit, qubit)

    # --- Operators (Full Space) ---
    σ_plus  = sigmap(qubit)
    σ_minus = sigmam(qubit)
    σ_z = sigmaz(qubit)

    a_1 = embed(basis_full, 1, destroy(cavity_1))
    a_2 = embed(basis_full, 2, destroy(cavity_2))
    a_dag_1 = embed(basis_full, 1, create(cavity_1))
    a_dag_2 = embed(basis_full, 2, create(cavity_2))

    σ_p_1 = embed(basis_full, 3, σ_plus)
    σ_p_2 = embed(basis_full, 4, σ_plus)
    σ_m_1 = embed(basis_full, 3, σ_minus)
    σ_m_2 = embed(basis_full, 4, σ_minus)
    σ_z_1 = embed(basis_full, 3, σ_z)
    σ_z_2 = embed(basis_full, 4, σ_z)

    h_tms = 1im * (conj(η) * a_1 * a_2 - η * a_dag_1 * a_dag_2)
    h_int_1 = a_1 * σ_p_1 + a_dag_1 * σ_m_1
    h_int_2 = a_2 * σ_p_2 + a_dag_2 * σ_m_2

    H_full = g_1 * h_int_1 + g_2 * h_int_2 + h_tms

    J_full = Operator[
        sqrt(κ_1) * a_1,
        sqrt(κ_2) * a_2,
        sqrt(γ) * σ_m_1,
        sqrt(γ) * σ_m_2,
        sqrt(γ_phi/2) * σ_z_1,
        sqrt(γ_phi/2) * σ_z_2
    ]

    # --- Operators (Adiabatic Space) ---
    σ_p_1_eff = embed(basis_adia, 1, σ_plus)
    σ_m_1_eff = embed(basis_adia, 1, σ_minus)
    σ_z_1_eff = embed(basis_adia, 1, σ_z)

    σ_p_2_eff = embed(basis_adia, 2, σ_plus)
    σ_m_2_eff = embed(basis_adia, 2, σ_minus)
    σ_z_2_eff = embed(basis_adia, 2, σ_z)

    H_adia = 1im * ( ( η * g_1 * g_2 ) / ( η^2 - ( ( κ_1 * κ_2 ) / 4 ) ) ) *
             ( σ_p_1_eff * σ_p_2_eff - σ_m_1_eff * σ_m_2_eff )

    Γ_sqrt_1 = 2im * ( (κ_1 * sqrt(κ_2 ) ) / ( 4η^2 - ( κ_1 * κ_2 ) ) )
    Γ_sqrt_2 = 2im * ( (κ_2 * sqrt(κ_1 ) ) / ( 4η^2 - ( κ_1 * κ_2 ) ) )

    ϵ_1 = ( 2η ) / κ_1
    ϵ_2 = ( 2η ) / κ_2

    J_adia = Operator[
        ( Γ_sqrt_1 ) * ( g_1 * ϵ_1 * σ_p_1_eff - g_2 * σ_m_2_eff ),
        ( Γ_sqrt_2 ) * ( g_2 * ϵ_2 * σ_p_2_eff - g_1 * σ_m_1_eff ),
        sqrt(γ) * σ_m_1_eff,
        sqrt(γ) * σ_m_2_eff,
        sqrt(γ_phi/2) * σ_z_1_eff,
        sqrt(γ_phi/2) * σ_z_2_eff
    ]

    # --- Solve for steady states ---
    # Full space: large/sparse -> iterative solver.
    rho_ss_full = steadystate.iterative(H_full, J_full)

    # Adiabatic space: tiny (4x4) -> exact dense eigenvector solver.
    H_adia_dense = dense(H_adia)
    J_adia_dense = dense.(J_adia)
    rho_ss_adia = steadystate.eigenvector(H_adia_dense, J_adia_dense)

    rho_ss_full_atoms = ptrace(rho_ss_full, (1, 2))

    return rho_ss_full_atoms, rho_ss_adia
end