using NLsolve

# ==============================================================================
# truncation.jl
#
# Estimates mean photon number and variance for the two cavities via a
# self-consistent second-order cumulant expansion (mean-field + second
# moments), including qubit backaction on the cavities. This replaces the
# earlier closed-form (qubit-free, linear-OPO) truncation estimate with a
# solve of the actual coupled steady-state equations.
#
# Returns n_i and var_i for each cavity; deciding how many σ of margin to
# add on top (and how to combine them into an integer Fock-space cutoff
# N_i) is left to the calling script.
# ==============================================================================

# ---- pack / unpack helpers ---------------------------------------------
const NVARS = 32

function unpack(x::AbstractVector)
    z1, z2, n1, n2 = x[1], x[2], x[3], x[4]
    c = Vector{ComplexF64}(undef, 14)
    for k in 1:14
        c[k] = x[5 + 2(k-1)] + im*x[6 + 2(k-1)]
    end
    m,s11,s22,s12,X1m1,X1p1,X2m2,X2p2,X1m2,X1p2,X2m1,X2p1,Qpp,Qpm = c
    return (z1=z1, z2=z2, n1=n1, n2=n2, m=m, s11=s11, s22=s22, s12=s12,
            X1m1=X1m1, X1p1=X1p1, X2m2=X2m2, X2p2=X2p2,
            X1m2=X1m2, X1p2=X1p2, X2m1=X2m1, X2p1=X2p1,
            Qpp=Qpp, Qpm=Qpm)
end

# ---- residual function --------------------------------------------------
function make_F!(η, g1, g2, κ1, κ2)
    function F!(res, x)
        v = unpack(x)
        z1,z2,n1,n2 = v.z1, v.z2, v.n1, v.n2
        m,s11,s22,s12 = v.m, v.s11, v.s22, v.s12
        X1m1,X1p1,X2m2,X2p2 = v.X1m1, v.X1p1, v.X2m2, v.X2p2
        X1m2,X1p2,X2m1,X2p1 = v.X1m2, v.X1p2, v.X2m1, v.X2p1
        Qpp,Qpm = v.Qpp, v.Qpm

        # --- real equations ---
        rn1 = -κ1*n1 + 2η*real(s12) - 2g1*imag(X1p1)
        rn2 = -κ2*n2 + 2η*real(s12) - 2g2*imag(X2p2)
        rz1 = imag(X1p1)
        rz2 = imag(X2p2)

        # --- boson-boson (complex) ---
        rm   = -(κ1+κ2)/2*m + η*(s22+conj(s11)) + im*g1*X2p1 - im*g2*conj(X1p2)
        rs11 = -κ1*s11 + 2η*conj(m) - 2im*g1*X1m1
        rs22 = -κ2*s22 + 2η*m       - 2im*g2*X2m2
        rs12 = -(κ1+κ2)/2*s12 + η*(n1+n2+1) - im*g1*X2m1 - im*g2*X1m2

        # --- boson-qubit, same arm ---
        rX1m1 = -(κ1/2)*X1m1 + η*conj(X2p1) + im*g1*z1*s11
        rX1p1 = -(κ1/2)*X1p1 + η*conj(X2m1) - im*g1*(1-z1)/2 - im*g1*z1*(n1+1)
        rX2m2 = -(κ2/2)*X2m2 + η*conj(X1p2) + im*g2*z2*s22
        rX2p2 = -(κ2/2)*X2p2 + η*conj(X1m2) - im*g2*(1-z2)/2 - im*g2*z2*(n2+1)

        # --- boson-qubit, cross arm ---
        rX1m2 = -(κ1/2)*X1m2 + η*conj(X2p2) - im*g1*conj(Qpp) + im*g2*z2*s12
        rX1p2 = -(κ1/2)*X1p2 + η*conj(X2m2) - im*g1*conj(Qpm) - im*g2*z2*conj(m)
        rX2m1 = -(κ2/2)*X2m1 + η*conj(X1p1) - im*g2*conj(Qpp) + im*g1*z1*s12
        rX2p1 = -(κ2/2)*X2p1 + η*conj(X1m1) - im*g2*Qpm       - im*g1*z1*m

        # --- qubit-qubit ---
        rQpp = -im*g1*z1*conj(X1m2) - im*g2*z2*conj(X2m1)
        rQpm = -im*g1*z1*conj(X1p2) + im*g2*z2*X2p1

        res[1] = rn1; res[2] = rn2; res[3] = rz1; res[4] = rz2
        for (k,c) in enumerate((rm,rs11,rs22,rs12,rX1m1,rX1p1,rX2m2,rX2p2,
                                 rX1m2,rX1p2,rX2m1,rX2p1,rQpp,rQpm))
            res[5+2(k-1)] = real(c)
            res[6+2(k-1)] = imag(c)
        end
        return res
    end
    return F!
end

# ------------------------------------------------------------------------
# solve_steady_state: homotopy-continuation solve, ramping η from 0 up to
# η_target. At η=0 the exact solution is the trivial vacuum/ground state,
# which is used as the starting point; each subsequent step warm-starts
# from the previous converged solution.
# ------------------------------------------------------------------------
function solve_steady_state(η_target, g1, g2, κ1, κ2; nsteps=60, ftol=1e-12)
    x = zeros(NVARS)                       # trivial solution at η = 0
    ηs = range(0.0, η_target, length=nsteps+1)[2:end]

    local sol
    for η in ηs
        F! = make_F!(η, g1, g2, κ1, κ2)
        sol = nlsolve(F!, x; ftol=ftol, iterations=2000, method=:trust_region)
        if !converged(sol)
            @warn "Did not converge at η=$η; continuing anyway with best estimate"
        end
        x = sol.zero                        # warm-start next step
    end
    return x, sol
end

# ------------------------------------------------------------------------
# photon_number_moments: <n_i^2> and Var(n_i) via Gaussian (Wick)
# factorization of the 4th-order moment <(a_i†)^2 a_i^2>, consistent with
# the 2nd-order cumulant closure used to obtain n_i, s_ii above:
#
#   n_i^2  = (a_i†)^2 a_i^2 + n_i                       (operator identity)
#   <(a_i†)^2 a_i^2> ≈ |s_ii|^2 + 2 n_i^2               (Wick/Gaussian, <a_i>=0)
#
#   => <n_i^2>   ≈ 2*n_i^2 + n_i + |s_ii|^2      (Wick gives TWO n_i^2
#                                                  contractions, not one)
#   => Var(n_i)  = <n_i^2> - n_i^2  ≈ n_i^2 + n_i + |s_ii|^2
#                                    = n_i*(n_i+1) + |s_ii|^2
#
# NOTE: this is an approximation (Gaussian closure at 4th order),
# not exact -- consistent with, but distinct from, the 2nd-order
# cumulant truncation used for the moment equations themselves.
# ------------------------------------------------------------------------
function photon_number_moments(v)
    n1sq_mean = 2*v.n1^2 + v.n1 + abs2(v.s11)
    n2sq_mean = 2*v.n2^2 + v.n2 + abs2(v.s22)
    var1      = v.n1*(v.n1+1) + abs2(v.s11)
    var2      = v.n2*(v.n2+1) + abs2(v.s22)
    return (n1sq=n1sq_mean, n2sq=n2sq_mean, var1=var1, var2=var2)
end

# ------------------------------------------------------------------------
# estimate_truncation: convenience wrapper -- solves the steady state and
# returns just (n1, var1, n2, var2). Deciding how much margin (σ) to add
# on top, and how to round into an integer Fock-space cutoff, is left to
# the caller.
# ------------------------------------------------------------------------
function estimate_truncation(κ_1::Float64, κ_2::Float64, η::Float64,
                              g_1::Float64, g_2::Float64)
    x, sol = solve_steady_state(η, g_1, g_2, κ_1, κ_2)
    v = unpack(x)
    pm = photon_number_moments(v)
    return (n1=v.n1, var1=pm.var1, n2=v.n2, var2=pm.var2)
end

########################################################################
# Example usage
########################################################################
if abspath(PROGRAM_FILE) == @__FILE__
    η   = 0.4
    g1  = 1.0
    g2  = 0.6
    κ1  = 2.0
    κ2  = 3.0

    # stability check (linear/no-qubit threshold, for reference only)
    println("Linear OPO threshold η_th = ", sqrt(κ1*κ2)/2,
            "   (your η = $η)")

    x, sol = solve_steady_state(η, g1, g2, κ1, κ2)
    v = unpack(x)
    pm = photon_number_moments(v)

    println("\nConverged: ", converged(sol))
    println("n1 = ", v.n1)
    println("n2 = ", v.n2)
    println("z1 = ", v.z1, "   z2 = ", v.z2)
    println()
    println("<n1^2>   = ", pm.n1sq)
    println("<n2^2>   = ", pm.n2sq)
    println("Var(n1)  = ", pm.var1)
    println("Var(n2)  = ", pm.var2)
    println()
    println("residual norm = ", norm(sol.residual_norm isa Function ? 0.0 : sol.residual_norm))
end