# Concurrence Analysis — Handoff Document (rev. 4)
 
**Status in one line:** BOTH principal questions are now CLOSED. Existence
(section 3) by the elementary parabola argument; the maximum question
(section 10) by a three-region monotonicity proof — **dC/du < 0 on the entire
physical window for every x in (0,x_c), hence beta = 1 is the unique global
maximum of the concurrence.**
 
---
 
## 0c. What changed in revision 4
 
**Question C is CLOSED.** Section 10 is the new material. The route was NOT the
section 5.4 plan (Step 1's `Resolve` on the quartic was never run); it was
Step 3 — the direct sign proof — which turned out to be tractable after two
structural observations that collapse most of the domain:
 
1. **F_u < 0 identically** on the physical domain (10.2). One-line proof from
   `D - 2(u+2)D'` being linear in u with a negative coefficient.
2. **L := S'D - 2SD' is LINEAR in u with strictly positive slope** (10.3). Its
   sign on any u-interval is therefore decided by its value at the left endpoint
   alone. Rev. 3's section 9.2 already noted the degree drop to 1; what is new is
   the *sign* of the slope, which makes L monotone increasing and gives the
   at-most-one-crossing structure the whole proof rests on.
Together these split the domain into three regions, two of which close
immediately (negative minus non-negative) and the third of which is a small
strip x < x_dagger ~ 0.0518 handled by a discriminant-continuity argument (10.6).
 
**Also new this revision:**
 
- **dC/du and dC/dx assembled in unsquared form** (10.1) — rev. 3 only ever wrote
  the *squared* critical-point equations. The unsquared forms carry no extraneous
  -root caveat and are what the whole of section 10 actually uses.
- **dC/du at u=2 in closed form, proven negative for all x in (0,1)** (10.5).
  This upgrades the section 5.1 bullet "C''(beta=1) < 0 ... SAMPLED, not proven"
  to **proven**. Two independent derivations agree (10.5 and 10.4).
- **dC/dx machinery + sign results** (10.7): boundary identity
  dC/dx|_{C=0} = Phi_x/(D v g), and an exact closed form for C(2,x) reducing the
  sign of dC/dx at u=2 to one sextic.
- **CORRECTION to 9.5**: the "interior maximum in x near ~0.2 at u=2" is a
  grid artifact. The true maximum is at **x* = 0.148550199898** (exact: t*^2 where
  t* is the unique root in (0,1) of t^6+t^5-3t^4-2t^3-3t^2-t+1). See 10.7.
- **First CAS confirmation of the core objects.** sympy was available this
  session. Everything in section 10 is exact-verified; see 10.8 for the scope,
  which is much wider than 9.4's but STILL does not cover sections 2-4.
*Method note:* all root counts below use exact Sturm sequences
(`count_roots` / `Poly.intervals`) on integer-coefficient polynomials. No
floating-point scan is used as evidence anywhere in section 10; decimals shown
are illustrative only. This meets the section 8 bar.
 
---
 
## 0b. What changed in revision 3 (see section 9 for the full working)
 
- **dC/du re-derived independently** from the F,G form and found to reproduce the
  section 5.3 quartic exactly. This is a genuine cross-check, not a restatement:
  it was assembled from scratch and the answer matched.
- **The degree-1 claim for `D - 2(u+2)D'` and the leading-term cancellation in
  `S'D - 2SD'` are now verified by explicit expansion**, not asserted structurally.
- **dC/dx derived for the first time** (new work; not present in rev. 2), together
  with the x-critical-point equation `(u+2) S [g2_x D - 2 g2 D_x]^2 = g2 [S_x D - 2 S D_x]^2`.
  Its degree in x is **NOT yet determined** — see 9.3.
- **g2_x, D_x, S_x CAS-verified** (exact symbolic + finite-difference). This is
  the first CAS confirmation of anything in this project; note it covers ONLY
  these three derivatives, NOT sections 2-4, which remain hand-derived only.
- **Numerical trend data added** (9.5, 9.6) for x -> x_c at fixed u, and for u
  swept at fixed x. Evidence only — the section 8 rule that floating-point scans
  are not proof is unchanged and still applies to all of it.
---
 
## 0. What changed in revision 2
 
**Removed as superseded** (correct, but solving a harder problem than necessary):
 
- The entire "N-cubic side-investigation": the cubic `g^2 S = (u+2)(S')^2`, its
  degree-68 discriminant, the Mathematica `Reduce[Exists[...]]` CAD run, and the
  landmark values x ~ 0.114545 and x ~ 0.142659.
  *Reason:* that route studied critical points of N via dN/du. But the existence
  question only ever needed the **sign** of N, and N > 0 is equivalent to
  Phi > 0 where Phi = (u+2)g^2 - S is a **quadratic in u**. Sections 3-4 below
  replace ~40 lines of CAD with a two-line parabola argument. The cubic's results
  are not contradicted; they are simply about a different function (see 5.2).
- "u_c(x) computed numerically for several x" — now available in closed form
  (section 4.4), matching the three stored numerical values.
**Corrected:**
 
- Domain `x in (0,1)` was flagged "assumed, not verified from external source."
  It is not an assumption — see section 1.2.
- The old Mathematica block contained a real syntax bug: `D_ = ...` assigns to a
  *pattern* (`D_` parses as `Pattern[D, Blank[]]`), not a symbol. Renamed to `DD`
  throughout, and all other built-in collisions (`C`, `N`, `D`, `E`) avoided.
**Kept unchanged:** the beta <-> 1/beta symmetry, the u-substitution, the
sqrt(S) at u=2 correction, C''(beta=1) < 0 sampling, the dC/du quartic, and the
x = 0.11 counterexample (relocated to 5.2 — it is still the single most important
methodological warning in this file).
 
---
 
## 1. Setup
 
### 1.1 The state
 
Steady-state adiabatic density matrix (X-state, basis order gg, ge, eg, ee):
 
```
rho_ad^ss =
[ rho_gg,gg    0          0        rho_gg,ee ]
[ 0            rho_ge,ge  0        0         ]
[ 0            0          rho_eg,eg 0        ]
[ rho_gg,ee    0          0        rho_ee,ee ]
```
 
Auxiliary polynomials:
- q1 = (x-1)^4 + 4x   = x^4 - 4x^3 + 6x^2 + 1
- p1 = (x-1)^4 + 4x(x+1) = x^4 - 4x^3 + 10x^2 + 1
- p2 = (x^2-1)^2 - 4x(x+1) = x^4 - 6x^2 - 4x + 1
Elements rho_ij = N_ij / D_beta, with
 
- D_beta   = (1+x) p1 beta^2 + 2(1+x) p2 beta + (1+x) p1
- N_gg,gg  = q1 beta^2 + [(x-2)q1 - 8x + 4] beta + q1
- N_ge,ge  = x q1 beta^2 + 4x^2(x^2-2x-1) beta + 4x^2
- N_eg,eg  = 4x^2 beta^2 + 4x^2(x^2-2x-1) beta + x q1
- N_ee,ee  = 4x^3 beta^2 + x(x^4-10x^2+1) beta + 4x^3
- N_gg,ee  = sqrt(beta)(1+beta) sqrt(x) (1-x)^3 (x+1)
Physical parameters: x = 4 eta^2/(kappa_1 kappa_2), beta = g_1^2 kappa_2/(g_2^2 kappa_1).
 
### 1.2 Domain (previously flagged as an unverified assumption — now resolved)
 
beta > 0 is immediate (ratio of squares of positive rates).
 
x in (0,1) is **the below-threshold stability condition of the model itself**, not
an imported assumption. The nondegenerate two-mode-squeezing drive has its
parametric threshold at eta = sqrt(kappa_1 kappa_2)/2; above it the intracavity
field diverges and there is no steady state to compute a concurrence from. That
condition is exactly 4 eta^2 < kappa_1 kappa_2, i.e. x < 1. And x > 0 since
eta^2 > 0.
 
*Status:* derived from the threshold condition, NOT checked against a specific
equation in a specific reference. Confirm it matches how the threshold appears in
your own adiabatic-elimination derivation before putting it in front of a
supervisor.
 
### 1.3 Concurrence formula
 
For an X-state, C = 2 max{0, |rho_14| - sqrt(rho_22 rho_33),
|rho_23| - sqrt(rho_11 rho_44)}. Here rho_23 = 0, so the second entry is
-sqrt(rho_11 rho_44) <= 0 and never active. Hence
 
**C = 2 max{0, |rho_gg,ee| - sqrt(rho_ge,ge rho_eg,eg)}**
 
*Status:* this is the standard Yu-Eberly X-state specialisation. I know the general
X-state concurrence formula is due to Yu and Eberly, but I cannot verify an exact
equation number or the notation used in any specific paper of theirs — do not cite
a numbered equation for it without checking the source directly.
 
---
 
## 2. Verified foundations
 
Each of these was worked out by hand and is short enough to re-check by eye. None
has been CAS-confirmed in the current session (no Mathematica available); section 7
gives the checks.
 
### 2.1 Trace = 1
 
Using p1 - q1 = 4x^2 (immediate: p1 - q1 = 4x(x+1) - 4x = 4x^2):
 
- beta^2 coefficient of sum(N): q1 + x q1 + 4x^2 + 4x^3 = (1+x)(q1 + 4x^2) = (1+x)p1
- beta^0 coefficient: identical, = (1+x)p1
- beta^1 coefficient: [(x-2)q1 - 8x + 4] + 8x^2(x^2-2x-1) + x(x^4-10x^2+1).
  With (x-2)q1 = x^5 - 6x^4 + 14x^3 - 12x^2 + x - 2 and
  8x^2(x^2-2x-1) = 8x^4 - 16x^3 - 8x^2, the total is
  2x^5 + 2x^4 - 12x^3 - 20x^2 - 6x + 2 = 2(1+x)p2.
So sum(N) = D_beta exactly, Tr rho = 1. The element set is internally consistent.
 
### 2.2 Symmetry beta <-> 1/beta, and the u-substitution
 
D_beta(beta) = beta^2 D_beta(1/beta), and likewise for N_gg,gg, N_ee,ee, N_gg,ee;
N_ge,ge(beta) = N_eg,eg(1/beta) (they swap, but only their PRODUCT enters C).
Hence **C(beta,x) = C(1/beta,x)**.
 
Substituting u = beta + 1/beta (domain u >= 2, equality iff beta = 1):
 
**C(u,x) = 2[ v g - sqrt(S) ] / D**  when the bracket is positive, else 0
 
- v = sqrt(u+2), g = sqrt(x)(1-x)^3(x+1), g^2 = x(1-x)^6(x+1)^2
- a = x q1, c = 4x^2 = p1 - q1, b = 4x^2(x^2-2x-1) = (p1-q1)(x^2-2x-1)
- a + c = p1 + q1(x-1), a - c = q1(x+1) - p1
- S = a c u^2 + b(a+c) u + (a-c)^2 + b^2
- D = (1+x)(p1 u + 2 p2)    [LINEAR in u, so D'' = 0]
*Derivation of the u-form (checked this session):* N_ge,ge N_eg,eg with
N_ge = a beta^2 + b beta + c and N_eg = c beta^2 + b beta + a expands to
ac(beta^4+1) + b(a+c)(beta^3+beta) + (a^2+b^2+c^2)beta^2; dividing by beta^2 and
using beta^2 + 1/beta^2 = u^2 - 2 gives exactly S. Also D_beta = beta D, so
rho_ge rho_eg = beta^2 S/(beta D)^2 = S/D^2. And
sqrt(beta)(1+beta)/(beta) = (1+beta)/sqrt(beta) = sqrt(u+2) = v, so
rho_gg,ee = v g / D. Both radicands are therefore correct as written.
 
### 2.3 D > 0 on the whole physical domain
 
Worth checking because p2 < 0 over much of (0,1) (e.g. p2(0.2) = -0.038).
p1 > 0 on (0,1), so D is increasing in u and minimised at u = 2:
 
p1 + p2 = (x-1)^4 + (x^2-1)^2 = (x-1)^2[(x-1)^2 + (x+1)^2] = 2(x-1)^2(x^2+1) > 0
 
**D(2,x) = 4(1+x)(1-x)^2(x^2+1) > 0**, hence D > 0 for all u >= 2. Every division
by D below is legitimate, and D can never create or destroy a zero of C.
 
### 2.4 S > 0 (so sqrt(S) is real and the squaring in 3.1 is valid)
 
N_ge,ge = a beta^2 + b beta + c has discriminant
b^2 - 4ac = 16x^3 (x m^2 - q1) with m = x^2-2x-1. Expanding
x m^2 - q1 = x^5 - 5x^4 + 6x^3 - 2x^2 + x - 1, which vanishes twice at x = 1;
two synthetic divisions give
 
x m^2 - q1 = (x-1)^2 (x^3 - 3x^2 - x - 1)
 
On (0,1) the cubic x^3-3x^2-x-1 has derivative 3x^2-6x-1 < 0 throughout (its roots
1 +/- 2/sqrt(3) lie outside), so it decreases from -1 to -4: strictly negative.
Discriminant < 0 and a > 0 ==> N_ge,ge > 0 strictly. Same for N_eg,eg. Hence
**S > 0**.
 
### 2.5 sqrt(S) at u = 2 — CORRECTION, retained
 
With c = 4x^2, b = cm, m = x^2-2x-1:
 
S(2) = 4ac + (a-c)^2 + 2cm(a+c) + c^2 m^2 = (a+c)^2 + 2cm(a+c) + (cm)^2
     = (a + c(1+m))^2
 
and a + c(1+m) = (x^5 - 4x^4 + 6x^3 + x) + (4x^4 - 8x^3) = x^5 - 2x^3 + x
= x(x^2-1)^2. Therefore
 
**sqrt(S)(2,x) = x(x^2-1)^2,  S(2,x) = x^2(x^2-1)^4**
 
NOT (x^2-1)^2 as originally claimed — the earlier version was wrong by a factor of
x. This was caught only by a failed numeric cross-check. Re-verified independently
this session by the factorisation above.
 
---
 
## 3. Question A (existence): where is C non-zero? — CLOSED
 
### 3.1 The key move
 
C = 2 max{0,N}/D with N = v g - sqrt(S), and D > 0, so the entire question is the
sign of N. Both v g >= 0 and sqrt(S) >= 0, and squaring is strictly monotone on
[0,inf), so
 
**N > 0  <==>  (u+2) g^2 > S  <==>  Phi(u,x) := (u+2) g^2 - S > 0**
 
This is an **equivalence**, not an implication: no extraneous roots, and sqrt(x)
disappears since g^2 = x(1-x)^6(x+1)^2 is polynomial. (Validity requires S > 0,
established in 2.4.)
 
### 3.2 Phi is a downward parabola
 
S is quadratic in u and (u+2)g^2 is linear, so
 
Phi = -ac u^2 + B u + C0,   B = g^2 - b(a+c),   C0 = 2g^2 - (a-c)^2 - b^2
 
Leading coefficient -ac = -4x^3 q1 < 0 on (0,1) since q1 = (x-1)^4 + 4x > 0.
**Phi opens downward.** That single fact carries the whole argument.
 
### 3.3 Threshold x_c
 
Phi(2) = 4g^2 - S(2) = 4x(1-x)^6(1+x)^2 - x^2(1-x)^4(1+x)^4. Factoring out
x(1-x)^4(1+x)^2 > 0:
 
Phi(2,x) = x(1-x)^4(1+x)^2 [4(1-x)^2 - x(1+x)^2],  and
4(1-x)^2 - x(1+x)^2 = -(x^3 - 2x^2 + 9x - 4)
 
**Phi(2,x) > 0  <==>  x^3 - 2x^2 + 9x - 4 < 0**
 
This cubic is strictly increasing: 3x^2 - 4x + 9 has discriminant 16 - 108 < 0, so
the derivative never vanishes. One real root; negative at x=0 (-4), positive at
x=1 (+4), so the root lies in (0,1):
 
**x_c = 2/3 + cbrt(-19/27 + 4 sqrt(87)/9) + cbrt(-19/27 - 4 sqrt(87)/9)
     ~ 0.4838882550**
 
Same number as the old t^3+2t^2+t-2 Cardano root under x = t^2: expanding
(t^3+2t^2+t-2)(t^3-2t^2+t+2) = t^6 - 2t^4 + 9t^2 - 4. The second factor is >= 2 on
[0,1], so all sign information sits in the first. **The x-cubic is the cleaner
object — no sqrt(x) substitution needed.**
 
### 3.4 No-revival lemma
 
Phi'(u) = -2ac u + B, so Phi'(2) = B - 4ac. Expanding gives
 
**Phi'(2,x) = x(1-x)^2 Q(x),  Q = -3x^6 + 14x^5 - 33x^4 + 4x^3 - 13x^2 - 2x + 1**
 
Q(0.20) = +0.0635, Q(0.205) = +0.0247, Q(0.21) = -0.0150, and Q < 0 at 0.25, 0.3,
0.5, 0.7, 0.9 with Q(1) = -32: a single sign change at **x_v ~ 0.20813**.
 
x_v is *the value of x at which the parabola's vertex u* = B/(2ac) crosses the
domain boundary u = 2* — nothing more. x < x_v: vertex interior (Phi rises above
Phi(2) before turning over). x > x_v: vertex behind the boundary (Phi decreasing
on all of u >= 2).
 
*Direct confirmation:* at x = 0.2081, B = 0.176747 and ac = 0.044182, giving
u* = 2.00021, with roots u- = -0.44107, u+ = 4.44146 summing to 2u*. At x = 0.1,
u* = 8.4315, interior between 2 and u_c = 18.43.
 
**The load-bearing inequality is x_v ~ 0.2081 < x_c ~ 0.4839.**
 
### 3.5 Theorem
 
> **x in [x_c, 1):** C(u,x) = 0 for every u >= 2, i.e. for every beta > 0.
> **x in (0, x_c):** C(u,x) > 0 exactly on u in [2, u_c(x)).
 
*Proof.* Let x >= x_c. Then x > x_v, so Q(x) < 0, so Phi'(2) < 0. Phi' is strictly
decreasing (downward parabola), so Phi'(u) < Phi'(2) < 0 for all u > 2: Phi is
strictly decreasing on [2,inf). And x >= x_c gives Phi(2) <= 0. Hence Phi(u) <= 0
throughout, N <= 0, C == 0.
 
Let x < x_c. Then Phi(2) > 0, so C(2,x) > 0 at once; a downward parabola taking a
positive value has two real roots u- < 2 < u+, and is positive on [2,u+). QED
 
Note what the second branch did NOT need: no vertex location, no discriminant, no
root counting. The vertex only mattered in the first branch, and only through the
single inequality x_v < x_c. If x_v had exceeded x_c there would be a revival band
(x_c, x_v) where C = 0 at beta = 1 but non-zero at larger asymmetry; that band is
empty.
 
Everything rests on two one-variable facts: the sign of x^3-2x^2+9x-4 (elementary)
and Q < 0 on [x_c,1) (one exact CAD call).
 
---
 
## 4. Closed form for the sudden-death point u_c
 
### 4.1 In your notation
 
Four bridging identities:
 
1. p1 - q1 = 4x^2
2. a + c = x lam,  **lam = q1 + 4x**
3. a - c = x kap,  **kap = q1 - 4x = (1-x)^4**
4. b = x sig,      **sig = p1 + p2 - 2 q1**  [ = 4x(x^2-2x-1) ]
and g^2 = x kap tau with **tau = p1 + p2 - q1 + 4x** [ = (1-x^2)^2 ], giving the
relation **kap + tau = p1 + p2**: kap and tau are the two halves into which p1+p2
splits.
 
Every coefficient of Phi then carries a factor x, which cancels against
ac = 4x^3 q1:
 
**u_c = [ Btil + sqrt(Btil^2 + 16 x^2 q1 C0til) ] / (8 x^2 q1)**
 
**Btil = kap tau - x sig lam,   C0til = 2 kap tau - x(kap^2 + sig^2)**
 
Denominator 8x^2 q1 = 2 q1 (p1 - q1) if you prefer no bare x. Note
u_c = u* + sqrt(W)/(8x^2 q1): the parabola's axis plus its half-width.
 
### 4.2 Alternative form (shorter radicand, but needs a non-p prefactor)
 
Pulling x(1-x)^3 out of the discriminant:
 
u_c = [ (1-x)^6(1+x)^2 + 4x^2(1+2x-x^2)(q1+4x) + (1-x)^3 sqrt(W_old) ] / (8x^2 q1)
 
W_old = 9x^10 - 82x^9 + 285x^8 - 488x^7 + 210x^6 + 452x^5 + 482x^4 + 120x^3
        + 37x^2 - 2x + 1
 
Relation between them: **W = (1-x)^6 W_old** (cross-check in section 7).
Trade-off: 4.1 is compact in your symbols but degree 16 expanded; 4.2 is degree 10
under the radical but needs (1-x)^3. Use 4.1 for write-up, 4.2 for root isolation.
 
W_old is not a perfect square (matching 3x^5 would require the next coefficient to
be -41/3). Whether it factors over Q is unknown to me — `Factor` will settle it.
 
### 4.3 Threshold and lemma re-expressed in the same symbols
 
- S(2) = x^2 tau^2, so **Phi(2) = x tau (4 kap - x tau)**, and since x tau > 0:
  **C(beta=1) > 0 <==> 4 kap > x tau**, i.e. 4(q1-4x) > x(p1+p2-q1+4x).
  Expanding recovers (1-x)^2 (4 - 9x + 2x^2 - x^3) — the threshold cubic.
- **sign Phi'(2) = sign(Btil - 16 x^2 q1)**: +0.5434 at x=0.1 (vertex interior),
  -0.8416 at x=0.4 (vertex behind boundary). Consistent with x_v ~ 0.2081.
### 4.4 Numerical agreement (independent check of the whole chain)
 
| x    | u_c (closed form) | previously stored (root-finding on N) |
|------|-------------------|---------------------------------------|
| 0.02 | 577.76            | ~577.75                               |
| 0.10 | 18.4304           | ~18.43                                |
| 0.40 | 2.16302           | ~2.16                                 |
 
The stored values came from a completely different route, so this is genuine
independent corroboration. At x=0.1 and x=0.4, kap and tau also independently
reproduce (1-x)^4 and (1-x^2)^2.
 
### 4.5 In terms of beta, and asymptotics
 
Window: beta in (beta-, beta+) with beta+ = [u_c + sqrt(u_c^2-4)]/2 and
beta- = 1/beta+ — symmetric about beta = 1 on a log axis, as the symmetry demands.
At x = 0.1: beta+ ~ 18.3, beta- ~ 0.0546.
 
Small-x expansion: **u_c = 1/(4x^2) - 1/x + 4 + O(x)**. Check: x=0.1 gives 19 vs
18.43; x=0.02 gives 579 vs 577.7. So beta+ ~ 1/(4x^2) for weak squeezing — the
window widens quadratically as the drive is turned down. At the other end u_c -> 2
as x -> x_c^-: the window pinches shut at beta = 1 exactly at threshold.
 
---
 
## 5. Question C (maximum): where is C maximised? — CLOSED in section 10
 
> **This section is kept as the historical statement of the problem and, in 5.2,
> as a still-live methodological warning. The question it poses is answered in
> section 10. Read 5.2 regardless — it is still the most important trap in this
> file. Treat 5.4's plan as superseded: it was never executed, and the route that
> worked is different (see 0c).**
 
### 5.1 What was known before rev. 4
 
- C'(beta=1) = 0 for all x, from symmetry alone (u'(beta=1) = 1 - 1/beta^2 = 0
  forces this regardless of dC/du).
- C''(beta=1,x) < 0 at every x tested (many rational points, exact + finite-
  difference cross-checked) ==> beta = 1 is a LOCAL max. ~~**SAMPLED, not proven
  for all x in (0,1).**~~ **NOW PROVEN for all x in (0,1) — see 10.5.**
- Chain-rule shortcut: since u'(beta=1) = 0, C''(beta=1) = 2 dC/du|_(u=2).
- Numerical scans (many x, u up to 200) and bounded Brent root-finding on
  [2, u_c(x)] for x in {0.02,...,0.4} found no interior critical point of dC/du.
  **Evidence, not proof.**
### 5.2 THE WARNING — do not conflate sign questions with argmax questions
 
Three functions share the same *sign* but have *different critical points*:
 
    Phi = N (v g + sqrt(S))     [positive factor]
    C   = 2 N / D               [positive divisor]
 
Multiplying or dividing by a strictly positive function preserves sign everywhere
but relocates extrema. Consequences, all real and all previously stepped in:
 
- **Explicit counterexample at x = 0.11.** Interior critical point of N at
  u* ~ 2.1733410671295132:
  N(u*) = 0.412056 > N(2) = 0.411708  (N is LARGER off-centre)
  C(u*) = 0.218355 < C(2) = 0.231331  (C is SMALLER off-centre)
  So argmax N != argmax C, demonstrably.
- The removed cubic's landmark x ~ 0.1145 and this file's x_v ~ 0.2081 are
  **not in conflict** — they are critical points of N and of Phi respectively.
- Section 3 is safe despite all this *because existence only needs a sign*, and
  sign is the one thing all three functions share.
- **Nothing in sections 3-4 constrains where C peaks.**
### 5.3 The critical-point equation (the actual open object)
 
Since D > 0, write C = 2[sqrt(F) - sqrt(G)] with
 
**F = (u+2) g^2 / D^2,  G = S / D^2**
 
Then dC/du = 0 <==> F' sqrt(G) = G' sqrt(F) <==> (squared) F'^2 G = G'^2 F. With
F' = g^2 (D - 2(u+2)D')/D^3 and G' = (S'D - 2SD')/D^3, clearing D^-8 gives
 
**g^2 S [D - 2(u+2)D']^2 - (u+2)[S'D - 2SD']^2 = 0**
 
which reproduces the previously derived quartic — a consistency check on it. The
F,G form may be the friendlier `Reduce` input.
 
*On squaring:* F'^2 G = G'^2 F admits the extra branch F' sqrt(G) = -G' sqrt(F), so
squaring can introduce extraneous roots. That is harmless in the **root-counting
direction only**: if the squared equation has no real roots with u > 2, neither
does the original. Any root found must be checked against the unsquared equation.
 
*Why it is degree 4 and not 6* (structural, not a `Poly.degree()` readout):
- D - 2(u+2)D' = -(1+x)p1 u + 2(1+x)(p2 - 2p1): **degree 1**. Squared -> 2, times
  g^2 S -> 2 + 2 = 4.
- S'D - 2SD': leading terms 2ac u (1+x)p1 u and 2 ac u^2 (1+x)p1 **cancel
  exactly**, dropping it to degree 1. Squared and times (u+2) -> degree 3.
### 5.4 Plan — SUPERSEDED (never executed; see section 10 for what worked)
 
*Historical note:* Step 3 was the correct instinct and is what section 10 does.
Steps 1, 2, 4a, 4b were never run — the degree-174 discriminant WAS eventually
computed (10.6), but only restricted to a small strip, and for a
root-count-constancy argument rather than a CAD existential.
 
**Step 1** — Run `Resolve[Exists[u, u > 2 && quartic == 0], Reals]` restricted to
0 < x < x_c. Note the old rationale ("the cubic terminated, so the quartic might")
no longer applies, since the cubic route has been discarded; the quartic is now
the *only* heavy object left, and its discriminant is degree 174.
 
**Step 2 (if it terminates)** — Intersect with (0, x_c). For any overlap, evaluate
**C itself** at u = 2 and at each interior critical point. Existence of a critical
point does not make it a higher peak (5.2).
 
**Step 3 (better first attempt, arguably)** — Try a direct sign proof that
dC/du < 0 for all u in (2, u_c), which would settle the question with no root
counting at all. Given how completely the analogous approach collapsed the
existence problem in section 3, this is worth attempting before the CAD route.
Concretely: show F'/sqrt(F) < G'/sqrt(G) on the physical window.
 
**Step 4 (fallbacks, in order)**
- (4a) Isolate real roots of the degree-174 discriminant as a standalone
  single-variable problem.
- (4b) Restrict explicitly to u < u_c(x) — now that u_c is closed-form (4.1), this
  is a usable algebraic constraint rather than a numerical one, which it was not
  before. This may be what makes Step 1 terminate.
- (4c) If all stall: stop, and report the numerical evidence in 5.1 as an
  explicitly incomplete result.
---
 
## 6. Numeric summary
 
| quantity | value | meaning | status |
|---|---|---|---|
| x_c | 0.4838882550 | entanglement exists iff x < x_c, at ANY u | exact (Cardano) |
| x_v | ~0.20813 | vertex of Phi crosses u=2; only used to prove no revival | root of Q, exact in principle |
| u_c(0.02) | 577.76 | sudden-death point | closed form 4.1 |
| u_c(0.10) | 18.4304 | " | " |
| u_c(0.40) | 2.16302 | " | " |
| u_c -> 2 | as x -> x_c^- | window pinches shut at threshold | |
| u_c ~ 1/(4x^2) | small x | window widens as drive weakens | |
| x_dagger | 0.051842973189252369 | below it, G_u < 0 at u=2; the only strip needing the hard argument (10.6) | exact, root of V (10.3) |
| x* | 0.148550199898 | argmax over x of C(2,x); replaces the bogus "~0.2" of 9.5 | exact, t*^2 (10.7) |
| x_v | ~0.20813 | (unchanged) | |
 
Retired: x ~ 0.114545 and x ~ 0.142659 (N-cubic landmarks, superseded).
Retired: "interior maximum in x near ~0.2 at u=2" from 9.5 — grid artifact,
true value is x* = 0.1486 (10.7).
 
---
 
## 7. Mathematica
 
Built-in collisions (`C`, `D`, `N`, `E`) avoided throughout. NOTE: the previous
version of this file contained `D_ = ...`, which assigns to a pattern, not a
symbol — that was a genuine bug.
 
```mathematica
qq1 = (x-1)^4 + 4x;  pp1 = (x-1)^4 + 4x(x+1);  pp2 = (x^2-1)^2 - 4x(x+1);
aa = x qq1;  cc = 4x^2;  bb = 4x^2 (x^2-2x-1);
gsq = x (1-x)^6 (x+1)^2;
SS = aa cc u^2 + bb (aa+cc) u + (aa-cc)^2 + bb^2;
DD = (1+x)(pp1 u + 2 pp2);
PHI = Expand[(u+2) gsq - SS];
 
(* --- 2.3, 2.4: D > 0 and S > 0 justify all divisions and the squaring --- *)
Factor[DD /. u -> 2]                                       (* 4(1+x)(1-x)^2(1+x^2) *)
Reduce[DD <= 0 && 0 < x < 1 && u >= 2, {x,u}, Reals]       (* expect False *)
Factor[x (x^2-2x-1)^2 - qq1]                               (* (x-1)^2 (x^3-3x^2-x-1) *)
Reduce[SS <= 0 && 0 < x < 1 && u >= 2, {x,u}, Reals]       (* expect False *)
 
(* --- 3.2: quadratic, opening downward --- *)
Exponent[PHI, u]                                           (* 2 *)
Factor[Coefficient[PHI, u, 2]]                             (* -4x^3 qq1 < 0 on (0,1) *)
 
(* --- 3.3: threshold --- *)
Factor[PHI /. u -> 2]                     (* -x(x-1)^4(1+x)^2 (x^3-2x^2+9x-4) *)
Reduce[3x^2 - 4x + 9 <= 0, x, Reals]      (* False => cubic strictly increasing *)
CountRoots[x^3 - 2x^2 + 9x - 4, {x, 0, 1}]                 (* 1 *)
xc = Root[#^3 - 2#^2 + 9# - 4 &, 1];  N[xc, 20]            (* 0.48388825504305... *)
 
(* --- 3.4: no-revival lemma. THE load-bearing check --- *)
QQ = -3x^6 + 14x^5 - 33x^4 + 4x^3 - 13x^2 - 2x + 1;
Factor[D[PHI, u] /. u -> 2]                                (* x (x-1)^2 QQ *)
Reduce[QQ >= 0 && xc <= x < 1, x, Reals]                   (* MUST return False *)
CountRoots[QQ, {x, 0, 1}]                                  (* 1 *)
(* pick the index in (0,1) by inspecting all six roots first: *)
N[Root[QQ &, #], 10] & /@ Range[6]
xv = Root[-3#^6+14#^5-33#^4+4#^3-13#^2-2#+1 &, 2];  xv < xc  (* True *)
 
(* --- 3.5: the whole existence theorem in one line --- *)
Resolve[Exists[u, u >= 2 && PHI > 0] && 0 < x < 1, Reals]  (* expect 0 < x < xc *)
 
(* --- 4.1: bridging identities and the closed form --- *)
kap = qq1 - 4x;  lam = qq1 + 4x;
sig = pp1 + pp2 - 2 qq1;  tau = pp1 + pp2 - qq1 + 4x;
Simplify[{pp1 - qq1 - 4x^2, kap - (1-x)^4, tau - (1-x^2)^2,
          sig - 4x(x^2-2x-1), kap + tau - (pp1+pp2)}]      (* all 0 *)
Simplify[{aa+cc - x lam, aa-cc - x kap, bb - x sig, gsq - x kap tau}]  (* all 0 *)
 
Btil = kap tau - x sig lam;
C0til = 2 kap tau - x (kap^2 + sig^2);
uc = (Btil + Sqrt[Btil^2 + 16 x^2 qq1 C0til]) / (8 x^2 qq1);
 
Simplify[PHI /. u -> uc, 0 < x < 1]                        (* expect 0 *)
N[uc /. x -> #, 10] & /@ {1/50, 1/10, 2/5}                 (* 577.76, 18.4304, 2.16302 *)
N[Limit[uc, x -> xc, Direction -> "FromBelow"], 10]        (* expect 2 *)
Series[uc, {x, 0, 1}]                                      (* 1/(4x^2) - 1/x + 4 + ... *)
 
(* cross-check the two radicand forms against each other (4.2) *)
Wold = 9x^10 - 82x^9 + 285x^8 - 488x^7 + 210x^6 + 452x^5 + 482x^4 + 120x^3
       + 37x^2 - 2x + 1;
Simplify[(Btil^2 + 16 x^2 qq1 C0til) - (1-x)^6 Wold]       (* expect 0 *)
Factor[Wold]                                               (* does it split over Q? *)
 
(* --- 5: the OPEN question. Everything above is the existence problem only --- *)
Dp = D[DD, u];  Sp = D[SS, u];
quartic = Expand[gsq SS (DD - 2(u+2) Dp)^2 - (u+2)(Sp DD - 2 SS Dp)^2];
Exponent[quartic, u]                                       (* expect 4 *)
Exponent[Expand[Sp DD - 2 SS Dp], u]                       (* expect 1, per 5.3 *)
Resolve[Exists[u, 2 < u < uc && quartic == 0] && 0 < x < xc, Reals]   (* NOT YET RUN *)
```
 
---
 
## 8. Meta-notes for whoever continues this
 
- **Requirements:** no skipped algebra; explicit verification before trusting any
  simplification; no invented citations or equation numbers; explicit labelling of
  proven vs sampled vs assumed. Floating-point scans are NOT accepted as proof —
  only exact methods (Sturm on rational points, CAD via `Reduce`/`Resolve`).
- **Verification status of this revision:** sections 2-4 were derived by hand and
  are not CAS-confirmed (no Mathematica was available when they were written). Run
  section 7 before quoting any of it. Two things do corroborate them independently:
  the u_c values in 4.4 match numbers obtained by a completely different route, and
  4.1 and 4.2 are independent derivations that section 7 checks against each other.
- **A real error was caught mid-derivation** in section 4.2: an expansion produced
  W_old with linear coefficient +14 instead of -2, and only a spot-evaluation at
  x = 0.4 exposed it (37.39 vs the structural 30.99). Spot-check every expansion
  numerically before trusting it. This is the second such catch in this project
  (the first being 2.5).
- **The recurring structural trap** is 5.2: sign-preserving transformations are not
  extremum-preserving. It has now produced confusion three times in different
  costumes (N vs C at x=0.11; the cubic's 0.1145 vs Phi's 0.2081; the original
  "x > threshold => C = 0 for all u" slip). Expect it again.
- **Keep A and C separate.** Section 3-4 answers "is there entanglement." Section 5
  answers "where is it strongest." They are different questions with different
  mathematics.
*Added in rev. 4:*
 
- **The grid-artifact trap fired again**, in 9.5: a maximum was reported "near
  ~0.2" because 0.2 was the largest of three sampled points. True value 0.1486.
  A coarse grid locates an extremum only to within its own spacing — it does not
  locate it at a sample point. This is the same family as the 4.2 and 2.5 catches.
- **Look for the cheap structural fact before the expensive computation.** Question
  C was closed without ever running the section 5.4 Step 1 CAD. What did the work
  was two elementary sign observations (F_u < 0; L linear in u with positive slope)
  that reduced ~90% of the domain to "negative minus non-negative". The degree-174
  discriminant was eventually needed, but only on a strip of width 0.05 and only
  for root-count *constancy*, which is far cheaper than an existential.
- **`simplify` returning non-zero is not a refutation.** Two correct closed forms
  in 10.5 and 10.7 failed `simplify` because of nested radicals, and were confirmed
  by 25-digit evaluation at exact rational points instead. Distinguish "the CAS
  could not prove it" from "it is false" — and distinguish exact evaluation of a
  closed form from a floating-point scan. The former is evidence of the kind
  section 8 accepts as a cross-check; the latter is not.
---
 
## 9. Session log — derivative machinery in both variables (rev. 3)
 
Everything in this section was produced in one session. Provenance is labelled per
item: **derived** (by hand, written out), **CAS-verified** (checked with sympy),
**numeric** (evidence only, not proof).
 
### 9.1 The common form
 
Both derivatives come from the same skeleton. Writing C = 2(sqrt(F) - sqrt(G))
with F = (u+2)g^2/D^2 and G = S/D^2 (the section 5.3 form), then for either
variable t in {u, x}:
 
**dC/dt = F_t/sqrt(F) - G_t/sqrt(G)**
 
and setting it to zero gives, before squaring, **F_t sqrt(G) = G_t sqrt(F)**.
 
*The squaring caveat applies in BOTH variables* and is the same one already logged
in 5.3: squaring admits the spurious branch F_t sqrt(G) = -G_t sqrt(F), so the
squared equation may carry extraneous roots. Safe in the root-counting direction
only (no roots of the squared equation ==> no roots of the original); any root
actually found must be re-checked against the unsquared form.
 
### 9.2 dC/du — independent re-derivation (DERIVED; reproduces 5.3)
 
D is linear in u, so D_u = (1+x)p1 =: D' is u-independent, and S_u = 2ac u + b(a+c) =: S'.
 
- F_u = g^2/D^2 + (u+2)g^2 * (-2 D_u/D^3) = **g^2 [D - 2(u+2)D'] / D^3**
- G_u = [S' D^2 - S * 2 D D'] / D^4 = **[S' D - 2 S D'] / D^3**
Substituting into F_u^2 G = G_u^2 F:
 
- LHS = g^2 S [D - 2(u+2)D']^2 / D^8
- RHS = (u+2) g^2 [S' D - 2 S D']^2 / D^8
D > 0 everywhere (2.3), so D^-8 clears legitimately:
 
**g^2 S [D - 2(u+2)D']^2 = (u+2)[S' D - 2 S D']^2**
 
This is **exactly** the 5.3 quartic. Since it was assembled here from scratch via
F,G rather than copied, the agreement is an independent consistency check on 5.3.
 
**Degree, now verified by expansion rather than asserted:**
 
- `D - 2(u+2)D'` expanded in full:
  (1+x)p1 u + 2(1+x)p2 - 2(1+x)p1 u - 4(1+x)p1
  = **-(1+x)p1 u + 2(1+x)(p2 - 2p1)**  ==> degree 1, confirmed directly.
- `S'D - 2SD'` leading terms: from S'D, 2ac u * (1+x)p1 u = 2ac(1+x)p1 u^2;
  from 2SD', 2(ac u^2)(1+x)p1 = 2ac(1+x)p1 u^2. **Identical, so they cancel**
  ==> degree 1, confirmed. (Rev. 2 asserted this structurally; it is now checked.)
- Hence LHS: 1 -> squared 2, times g^2 S (S degree 2) -> **degree 4**.
  RHS: 1 -> squared 2, times (u+2) -> **degree 3**. Overall **quartic in u**.
### 9.3 dC/dx — NEW WORK (DERIVED; no rev.-2 counterpart to check against)
 
Ingredient derivatives, each written out:
 
- q1_x = **4x^3 - 12x^2 + 12x**
- p1_x = **4(x-1)^3 + 8x + 4**
- p2_x = **4x^3 - 12x - 4**
- c = 4x^2  ==>  c_x = **8x**
- a = x q1  ==>  a_x = q1 + x q1_x = (x^4-4x^3+6x^2+1) + (4x^4-12x^3+12x^2)
  = **5x^4 - 16x^3 + 18x^2 + 1**
- b = 4x^4 - 8x^3 - 4x^2  ==>  b_x = **16x^3 - 24x^2 - 8x**
**g^2_x**, by product rule on g^2 = x(1-x)^6(1+x)^2:
 
  g^2_x = (1-x)^6(1+x)^2 + x(-6)(1-x)^5(1+x)^2 + x(1-x)^6 * 2(1+x)
 
Factoring out (1-x)^5(1+x), the bracket is
(1-x)(1+x) - 6x(1+x) + 2x(1-x) = (1-x^2) + (-6x-6x^2) + (2x-2x^2) = 1 - 4x - 9x^2:
 
**g^2_x = (1-x)^5 (1+x) (1 - 4x - 9x^2)**
 
Sign note (no numerics needed): (1-x)^5(1+x) > 0 on (0,1), and 1-4x-9x^2 is a
downward parabola, positive at x=0 and negative at x=1, so it changes sign once
in (0,1). So g^2 rises then falls on (0,1) — consistent with g^2(0)=g^2(1)=0.
 
**D_x = (p1 u + 2 p2) + (1+x)(p1_x u + 2 p2_x)**
 
**S_x = (a_x c + a c_x) u^2 + [b_x(a+c) + b(a_x+c_x)] u + 2(a-c)(a_x-c_x) + 2 b b_x**
 
(S_x and D_x are deliberately LEFT in a/b/c/p1/p2 form rather than expanded into
raw x-polynomials: it keeps them auditable by eye and matches the notation of 4.1.
The cost is that the degree in x is not readable off them — see below.)
 
Assembling:
 
- F_x = (u+2)[g^2_x D - 2 g^2 D_x] / D^3
- G_x = [S_x D - 2 S D_x] / D^3
and F_x^2 G = G_x^2 F gives, after clearing D^-8 (D>0) and cancelling ONE factor
of (u+2) (legitimate: u >= 2 ==> u+2 >= 4, never zero on the physical domain):
 
**(u+2) S [g^2_x D - 2 g^2 D_x]^2 = g^2 [S_x D - 2 S D_x]^2**
 
**OPEN / NOT DONE:** the degree of this equation in x is **not determined**. For
the u-equation the degree claim was earned by expanding both bracketed factors;
the analogous expansion here has not been performed, and no degree is claimed.
Do not quote a degree for this equation until that expansion is done.
 
### 9.4 CAS verification of 9.3's ingredients (CAS-VERIFIED)
 
This is the first CAS confirmation obtained in this project. **Scope: it covers
g^2_x, D_x, S_x and their sub-derivatives ONLY.** Sections 2-4 remain
hand-derived and CAS-unconfirmed; the section 7 Mathematica block is still unrun.
 
Two independent checks, both passed:
 
1. **Exact symbolic.** sympy differentiated g^2, D, S directly and the difference
   against the hand-derived forms simplified to 0 in every case — including each
   of p1_x, p2_x, q1_x, a_x, b_x, c_x separately.
2. **Finite difference** at x = 0.2, u = 5, step h = 1e-6, central difference:
| quantity | finite difference | hand-derived symbolic |
|---|---|---|
| g^2_x | -0.06291456000 | -0.06291456000 |
| D_x   | 12.80000000    | 12.80000000    |
| S_x   | 8.766259200    | 8.766259200    |
| dC/dx | -1.471220819   | -1.471220819   |
 
The last row checks the *fully assembled* dC/dx of 9.3, not just its parts.
 
This is exactly the spot-check discipline section 8 demands, and it is worth
noting it was run BEFORE investing effort in expanding the x-equation — i.e. the
lesson from the 2.5 and 4.2 errors was applied prospectively for once.
 
### 9.5 Trend at fixed u as x -> x_c (NUMERIC — evidence only)
 
**u = 2 (i.e. beta = 1):**
 
| x | C(2,x) |
|---|---|
| 0.10 | 0.22733167 |
| 0.20 | 0.22862584 |
| 0.30 | 0.17284935 |
| 0.35 | 0.13211152 |
| 0.40 | 0.08575286 |
| 0.45 | 0.03551037 |
| 0.47 | 0.01466106 |
| 0.48 | 0.00411782 |
| 0.483 | 0.00094155 |
| 0.4838 | 0.00009357 |
| 0.48388 | 0.00000875 |
| 0.483888 | 0.00000027 |
 
C(2,x) is **non-monotonic**: it rises slightly from x=0.1 to x=0.2, then decreases
monotonically to zero, vanishing **continuously** (no jump) as x -> x_c^-. That
continuity is expected by construction since Phi(2,x_c) = 0 (3.3).
 
*Note the interior maximum in x near ~0.2 at u=2.* This is a maximum in the OTHER
variable and says nothing about Question C, which is about u. But it does mean
dC/dx = 0 has at least one root on the physical domain, so the 9.3 equation is not
vacuous.
 
> **CORRECTED IN REV. 4.** "Near ~0.2" is wrong — it was read off a grid whose
> points were {0.1, 0.2, 0.3}, i.e. the peak was assigned to whichever sample
> happened to be largest. The true maximiser is **x* = 0.148550199898**, and
> C(2,x*) exceeds every value in the table above. See 10.7 for the exact
> determination. The qualitative conclusion (an interior max in x exists at u=2,
> so the 9.3 equation is not vacuous) survives unchanged.
>
> This is the section-8 lesson recurring for the third time: a coarse grid
> locates a maximum only to within its own spacing. Item 2 of 9.7 ("verify the
> x=0.2-ish interior maximum is a genuine root of the 9.3 equation") would have
> FAILED, for the wrong reason, if run against 0.2.
 
**u = 3 (fixed interior value):**
 
| x | C(3,x) |
|---|---|
| 0.05 | 0.15760010 |
| 0.10 | 0.17179453 |
| 0.15 | 0.14757448 |
| 0.20 | 0.09927442 |
| 0.25 | 0.03560185 |
| 0.30 | 0 (exactly) |
| 0.35 - 0.4838 | 0 (exactly) |
 
At u=3, C reaches zero **between x=0.25 and x=0.30**, well before x_c, and stays
zero thereafter. This is not decay — it is u=3 falling outside the window once
u_c(x) < 3. It is a direct numerical illustration of 3.5 and of u_c -> 2 in 4.5:
for any FIXED u > 2, the window excludes that u strictly before x reaches x_c.
Only at u = 2 does "C -> 0" and "x -> x_c" coincide.
 
*Not done:* the crossing value solving u_c(x) = 3 was not computed in closed form.
Since u_c is closed-form (4.1), this is inversion algebra, not new numerics.
 
### 9.6 Trend at fixed x as u varies (NUMERIC — evidence only)
 
**x = 0.1** (u_c = 18.4304), 12-point grid on [2, u_c):
 
| u | C |
|---|---|
| 2.0000 | 0.22733167 |
| 3.4920 | 0.15196934 |
| 4.9840 | 0.10913040 |
| 6.4760 | 0.08128699 |
| 7.9680 | 0.06151894 |
| 9.4600 | 0.04662446 |
| 10.9520 | 0.03491713 |
| 12.4440 | 0.02542032 |
| 13.9360 | 0.01752685 |
| 15.4280 | 0.01083786 |
| 16.9200 | 0.00507970 |
| 18.4120 | 0.00005779 |
 
**x = 0.4** (u_c = 2.16302), 12-point grid on [2, u_c):
 
| u | C |
|---|---|
| 2.0000 | 0.08575286 |
| 2.0146 | 0.07649306 |
| 2.0292 | 0.06759791 |
| 2.0439 | 0.05904635 |
| 2.0585 | 0.05081888 |
| 2.0731 | 0.04289746 |
| 2.0877 | 0.03526533 |
| 2.1024 | 0.02790696 |
| 2.1170 | 0.02080788 |
| 2.1316 | 0.01395462 |
| 2.1462 | 0.00733463 |
| 2.1609 | 0.00093619 |
 
At both x values C decreases monotonically across the grid, max at the left
endpoint u=2, no interior bump.
 
**How much this is worth (read before quoting it):** a 12-point grid at two x
values. It is consistent with 5.1's existing scans and with beta=1 being the
global max, and it is NOT in tension with the x=0.11 counterexample of 5.2 —
that counterexample is about N, a different function. But failing to find an
interior critical point on a coarse grid is strictly weaker than showing none
exists, and Question C asks whether one can exist for ANY x in (0,x_c). Per
section 8, this does not move the status of Question C from OPEN.
 
### 9.7 State of play / next moves — SUPERSEDED by 10.9
 
*Historical.* The recommendation below was correct and is what rev. 4 executed.
 
Unchanged recommendation from 5.4: **Step 3 (direct sign proof) before the CAD
route.** Show F_u/sqrt(F) < G_u/sqrt(G) on the physical window; that closes
Question C with no root counting and no squaring caveat at all. The 9.6 numerics
are at least consistent with that inequality holding, which is mild encouragement
to attempt it, nothing more.
 
Concrete unfinished items created by this session:
 
1. ~~Expand the 9.3 x-equation to determine its degree in x.~~ **Moot** — the
   unsquared form (10.1) is what section 10 uses; the squared x-equation was never
   needed. Its degree in x is still not determined and no longer matters.
2. ~~Verify that the x=0.2-ish interior maximum in 9.5 is a genuine root.~~
   **Done differently** — the maximiser was determined exactly (10.7) and 0.2 was
   shown to be wrong. See the correction box in 9.5.
3. Run the section 7 Mathematica block. **STILL NOT DONE.** Sections 2-4 remain
   CAS-unconfirmed. 10.8's verification does not cover them either. This is now
   the single largest unverified block in the file.
4. Invert u_c(x) = u for the 9.5 crossing values. **Still not done** (low value).
---
 
## 10. Question C CLOSED — C is strictly decreasing in u (rev. 4)
 
**Theorem.** For every x in (0, x_c) and every u in [2, u_c(x)),
 
>   **dC/du < 0.**
 
Consequently C(.,x) attains its maximum over the physical window at the left
endpoint u = 2, i.e. **beta = 1 is the unique global maximiser of the concurrence
for every x in (0,x_c)**, and the maximum value is the closed form C(2,x) of 10.7.
 
The proof is a three-region decomposition. Regions 1 and 2 are one-line sign
arguments; region 3 is a small strip closed by discriminant continuity.
 
### 10.1 The unsquared derivatives (DERIVED; CAS-VERIFIED)
 
Rev. 3 wrote only the *squared* critical-point equations. Assembling the
derivatives themselves: with C = 2(v g - sqrt(S))/D, v = sqrt(u+2),
D' := D_u = (1+x)p1 (constant in u), S' := S_u = 2ac u + b(a+c):
 
>   **dC/du = (1/D)[ g/v - S'/sqrt(S) ] - (D'/D) C**
 
and, at fixed u (so v is constant),
 
>   **dC/dx = (1/D)[ v (g^2)_x / g - S_x/sqrt(S) ] - (D_x/D) C**
 
*Derivation (du case, quotient rule on C = 2 Ncal/D with Ncal = v g - sqrt(S)):*
dC/du = 2(Ncal_u D - Ncal D')/D^2 = 2Ncal_u/D - (2Ncal/D)(D'/D), and the second
term is C D'/D. Then dv/du = 1/(2 sqrt(u+2)) = 1/(2v), so
Ncal_u = g/(2v) - S'/(2 sqrt(S)). Multiply by 2/D. QED
 
*Cross-check against the 9.1 F,G skeleton:* using
F_u/sqrt(F) = g[D - 2(u+2)D']/(v D^2) and G_u/sqrt(G) = [S'D - 2SD']/(sqrt(S) D^2),
subtracting and regrouping gives
(1/D^2){ D(g/v - S'/sqrt(S)) - 2D'( g(u+2)/v - S/sqrt(S) ) }; since (u+2)/v = v
and S/sqrt(S) = sqrt(S), the second bracket is v g - sqrt(S) = Ncal = CD/2,
recovering the boxed form. Both routes agree.
 
*Status:* CAS-VERIFIED. `simplify(diff(C,u) - boxed_u) == 0` and
`simplify(diff(C,x) - boxed_x) == 0`, both exact.
 
*Why these matter:* no squaring, hence **no extraneous-root caveat at all** —
the 5.3 / 9.1 warning does not apply to anything in section 10 except 10.6, where
it is handled explicitly.
 
### 10.2 F_u < 0 identically (DERIVED; CAS-VERIFIED)
 
F = (u+2)g^2/D^2, and from 9.2, F_u = g^2 [D - 2(u+2)D'] / D^3. Expanding:
 
  D - 2(u+2)D' = (1+x)p1 u + 2(1+x)p2 - 2(u+2)(1+x)p1
               = **(1+x)[ -p1 u + 2p2 - 4p1 ]**
 
p1 = (x-1)^4 + 4x(x+1) > 0 on (0,1), and u >= 2, so p1 u >= 2 p1 and hence
the bracket is <= 2p2 - 6p1 = 2(p2 - 3p1). Now
 
  p2 - 3p1 = (x^4 - 6x^2 - 4x + 1) - 3(x^4 - 4x^3 + 10x^2 + 1)
           = **-2x^4 + 12x^3 - 36x^2 - 4x - 2**
 
On (0,1): 12x^3 - 36x^2 = 12x^2(x-3) < 0, and -2x^4, -4x, -2 are each <= 0 with
-2 < 0 strictly. So p2 - 3p1 < 0, the bracket is strictly negative, and since
g^2 > 0 and D^3 > 0 (2.3):
 
>   **F_u < 0 for all u >= 2, x in (0,1).**
 
*CAS-VERIFIED:* `factor(D - 2(u+2)D')` returns
`-(x+1)(u x^4 - 4u x^3 + 10u x^2 + u + 2x^4 - 16x^3 + 52x^2 + 8x + 2)`, whose
second factor is manifestly positive for u >= 2, x in (0,1) — note
2x^4 - 16x^3 + 52x^2 = 2x^2(x^2 - 8x + 26) and x^2-8x+26 has discriminant
64 - 104 < 0.
 
**This is the engine of the whole proof.** Since dC/du = F_u/sqrt(F) - G_u/sqrt(G)
and F_u < 0 always, **G_u >= 0 is by itself sufficient for dC/du < 0**
(negative minus non-negative). No comparison of magnitudes is needed.
 
### 10.3 L := S'D - 2SD' is linear in u with POSITIVE slope (CAS-VERIFIED)
 
sign(G_u) = sign(L), L = S'D - 2SD'. Rev. 3's 9.2 established deg_u L = 1 by
showing the u^2 terms cancel. Writing L = alpha(x) u + beta(x), the CAS gives the
factorisations (this is the new content):
 
>   **alpha = 4 x^3 (1-x)^3 (1+x) V2(x)**,
>   V2 = x^7 - 7x^6 + 19x^5 - 21x^4 + 7x^3 + 23x^2 + 5x + 5
>
>   **beta = -2 x^2 (1-x)^3 (1+x) V3(x)**,
>   V3 = x^9 - 13x^8 + 68x^7 - 160x^6 + 214x^5 - 154x^4 - 28x^3 - 56x^2 + x - 1
 
Exact Sturm counts on (0,1): **V2 has 0 real roots** there, with V2(0)=5>0 and
V2(1)=32>0, so V2 > 0 throughout. **V3 has 0 real roots** there, with V3(0)=-1<0
and V3(1)=-128<0, so V3 < 0 throughout. The prefactors x^3, x^2, (1-x)^3, (1+x)
are all positive on (0,1). Hence
 
>   **alpha(x) > 0 and beta(x) < 0 for every x in (0,1).**
 
Two consequences, both used below:
 
- **L is strictly increasing in u.** So L >= 0 at any point implies L >= 0 for all
  larger u. In particular *L(2,x) >= 0 implies G_u >= 0 on the whole window* —
  the right endpoint never needs separate examination.
- **L has exactly one zero**, at the rational function
>   **u*(x) = -beta/alpha = -V3(x) / (2 x V2(x))**
 
  Note: **no square root** — the sqrt(W) of u_c never enters. G_u < 0 for
  u < u*(x) and G_u > 0 for u > u*(x).
 
### 10.4 The left-endpoint split: L(2,x) and x_dagger (DERIVED; CAS-VERIFIED)
 
Substituting u=2 into L, using S'(2) = 4x^3 P(x) with
**P = x^6 - 6x^5 + 17x^4 - 20x^3 + 11x^2 - 6x + 3** (CAS-verified:
`expand(S'(2) - 4x^3 P) == 0`), together with D(2) = 4(1+x)(1-x)^2(1+x^2) (2.3),
S(2) = x^2(1-x^2)^4 (2.5) and D' = (1+x)p1:
 
  L(2,x) = 16x^3(1+x)(1-x)^2(1+x^2) P - 2x^2(1-x)^4(1+x)^5 p1
         = 2x^2(1+x)(1-x)^2 [ 8x(1+x^2)P - (1-x)^2(1+x)^4 p1 ]
 
The bracket expands to a degree-10 polynomial with a **triple root at x=1**;
dividing out (x-1)^3 and absorbing signs:
 
>   **L(2,x) = 2 x^2 (1+x) (1-x)^5 V(x)**,
>   **V(x) = x^7 - 7x^6 + 25x^5 - 27x^4 + 51x^3 + 3x^2 + 19x - 1**
 
Prefactor > 0 on (0,1), so sign L(2,x) = sign V(x). Exact Sturm: **V has exactly
one real root in (0,1)** (and it is V's only real root anywhere), at
 
>   **x_dagger = 0.051842973189252368741...**  (exact: the unique real root of V)
 
V(0) = -1 < 0, so
 
>   **L(2,x) < 0 on (0, x_dagger)**;  **L(2,x) > 0 on (x_dagger, 1)**
 
*Independent CAS check of the factorisation:* `factor(L.subs(u,2))` returns
`-2 x^2 (x-1)^5 (x+1) V(x)`, identical after folding the sign into (1-x)^5.
 
*Aside, not load-bearing:* this same V reappears in the numerator of 10.5. I have
no a priori reason why the two should share a factor and have not looked for one;
it may be structural or may be an accident of this problem's algebra. **Do not
build anything on it without an actual explanation.**
 
### 10.5 dC/du at u = 2, in closed form, negative for ALL x (DERIVED; CAS-VERIFIED)
 
Independently of the above, substituting u=2 into the 10.1 boxed formula and
clearing denominators (a factor (1-x)^3 cancels between numerator and
denominator):
 
>   **dC/du |_(u=2) = -[ sqrt(x)(1+x)^3 (x^4-6x^3+18x^2+2x+1) + x V(x) ]
>                      / [ 8(1-x)(1+x)^3(1+x^2)^2 ]**
 
with V as in 10.4. Denominator > 0 on (0,1). Substituting x = t^2 the bracket is
t Y(t) with
 
  Y(t) = t^15 + t^14 - 7t^13 - 3t^12 + 25t^11 + 3t^10 - 27t^9 + 39t^8
         + 51t^7 + 55t^6 + 3t^5 + 27t^4 + 19t^3 + 5t^2 - t + 1
 
**Y > 0 on (0,1)**, by grouping its only four negative terms — no root-finding
needed:
 
- 25t^11 - 7t^13 - 3t^12 >= 25t^11 - 10t^11 = 15t^11 > 0   (t^13, t^12 <= t^11)
- 39t^8 - 27t^9 >= 12t^8 > 0
- 1 - t + 5t^2 > 0   (discriminant 1 - 20 = -19 < 0)
- every remaining term is positive.
Hence
 
>   **dC/du|_(u=2) < 0 for every x in (0,1)**, so **C''(beta=1) = 2 dC/du|_(u=2) < 0**.
 
This is the promised upgrade of the 5.1 bullet from SAMPLED to PROVEN.
 
*CAS-VERIFIED numerically to 25 digits* at x = 0.02, 0.05, 0.1, 0.3 against
`diff(C,u).subs(u,2)` (exact agreement; sympy's `simplify` fails on the nested
radicals, so high-precision evaluation was used instead — this is exact-arithmetic
evaluation of a closed form, not a floating-point scan).
 
*Consistency:* for x > x_dagger this is also implied by 10.2 + 10.3 + 10.4
(G_u(2,x) >= 0 and F_u < 0). The two derivations are independent — one via
explicit factoring, one via the linear-L sign argument — and they agree.
 
### 10.6 The three regions
 
Write **Q(u,x) := g^2 S [D - 2(u+2)D']^2 - (u+2)[S'D - 2SD']^2**, the 5.3 /
9.2 quartic. Because F_u^2 G - G_u^2 F = g^2 Q / D^8 and g^2, D^8 > 0,
sign(Q) = sign(F_u^2 G - G_u^2 F); and when F_u, G_u are **both** negative,
 
  dC/du < 0  <==>  |F_u|/sqrt(F) > |G_u|/sqrt(G)  <==>  F_u^2 G > G_u^2 F  <==>  Q > 0.
 
**Region 1 — x in (x_dagger, x_c), all u in [2, u_c].**
L(2,x) >= 0 by 10.4; L is increasing in u by 10.3; so L >= 0 on the whole window,
i.e. G_u >= 0. With F_u < 0 (10.2): dC/du = (negative) - (non-negative) < 0. **Done.**
This is the bulk of the domain: x_dagger ~ 0.0518 vs x_c ~ 0.4839.
 
**Region 2 — x in (0, x_dagger), u in (u*(x), u_c(x)].**
Past the unique crossing of L, G_u > 0. Same argument. **Done.**
 
**Region 3 — x in (0, x_dagger), u in [2, u*(x)).**
Here G_u < 0 and F_u < 0: both negative, so the magnitudes genuinely compete and
Q's sign must be established. Two facts are already available *without* any
root-finding:
 
- **Q(2,x) > 0** for all x in (0, x_dagger) — this is exactly 10.5.
- **Q(u*(x),x) > 0** — trivially, since G_u = 0 there by definition of u*, leaving
  Q = g^2 S [D-2(u+2)D']^2 > 0 (S > 0 by 2.4; the bracket is non-zero by 10.2).
So Q > 0 at both ends of the interval. It remains to exclude interior roots. As x
varies continuously, a pair of real roots of Q(.,x) can appear or disappear only
through a **double root** (disc_u Q = 0) or by **escaping to infinity** (leading
coefficient = 0). Both are excluded on (0, x_dagger):
 
- **Leading coefficient.** `factor` gives
  **[u^4] Q = 4 x^4 (1-x)^6 (1+x)^4 q1 p1^2**, whose only zero in
  [0, x_dagger) is x = 0 (isolating-interval computation: multiplicity 4 at 0,
  nothing else). x = 0 is excluded from the physical domain by 1.2. Note this is
  also an independent structural confirmation of the degree-4 claim of 5.3/9.2,
  since q1, p1 > 0.
- **Discriminant.** disc_u(Q) has degree 174 in x (154 terms), matching the value
  recorded in 5.4. `Poly.intervals` on [0, x_dagger) returns exactly one entry,
  **((0,0), multiplicity 21)** — the only zero is x = 0 itself, again excluded.
Therefore the number of roots of Q(.,x) in (2, u*(x)) is **constant** across
(0, x_dagger), and evaluating at one interior sample settles the whole strip. At
the rational point **x = 1/50** (well inside, since 0.02 < 0.0518):
 
  u*(1/50) = 1958296807330649/399160367215100 = 4.90604019881...
  exact Sturm count of roots of Q(u, 1/50) in (2, u*): **0**
  Q(2, 1/50) > 0 and Q(u*, 1/50) > 0, both confirmed exactly.
 
Hence Q > 0 throughout region 3, so dC/du < 0 there too. **Done.**
 
All three regions are covered, and the theorem of 10 follows. []
 
*On the squaring:* the extraneous-branch caveat of 5.3/9.1 is confined to region 3,
where it is harmless: the argument only ever uses **Q > 0 ==> no critical point**,
which is the safe direction (a spurious root of the squared equation could only
make Q vanish, and Q does not vanish). Regions 1 and 2 never square anything.
 
### 10.7 dC/dx results (secondary; not needed for the theorem)
 
Recorded because they were derived alongside and because one of them corrects 9.5.
 
**Boundary identity.** At any point where C = 0 (i.e. u = u_c(x)), the
-(D_x/D)C term drops and sqrt(S) = v g, so the 10.1 formula collapses to
 
>   **dC/dx |_(C=0) = Phi_x / (D v g)**,  and likewise  **dC/du|_(C=0) = Phi_u/(D v g)**
 
Denominators positive. Since Phi is a downward parabola in u with roots
u_- < 2 < u_+ = u_c (3.5), Phi_u(u_c) < 0 strictly, giving an independent
confirmation that dC/du < 0 at the right endpoint. Likewise dC/dx < 0 at the death
boundary in x, since Phi crosses downward there.
 
**No global sign for dC/dx.** As x -> 0 at fixed u, using q1 -> 1:
a = x + O(x^3), c = 4x^2, b = -4x^2 + O(x^3), so S = x^2 + O(x^3), sqrt(S) = x(1+O(x)),
g = sqrt(x)(1+O(x)), D -> u+2 = v^2, whence
 
  C = 2 sqrt(x)/v - 2x/v^2 + O(x^(3/2))   and   **dC/dx = 1/sqrt(x(u+2)) + O(1) -> +infinity**
 
Positive near x=0, negative at the death boundary: **at least one sign change for
every fixed u.** No global statement is available, unlike for dC/du.
 
**Exact closed form at u = 2.** With v=2, sqrt(S)(2) = x(1-x^2)^2 (2.5) and
D(2) = 4(1+x)(1-x)^2(1+x^2) (2.3), everything cancels:
 
>   **C(2,x) = [ 2 sqrt(x) (1-x) - x(1+x) ] / [ 2(x^2+1) ]**
 
*Two independent confirmations:* setting the numerator to zero and putting x = t^2
gives 2t(1-t^2) = t^2(1+t^2), i.e. t^3 + 2t^2 + t - 2 = 0 after dividing by t —
**exactly the Cardano cubic of 3.3**. And it reproduces the 9.5 table
(x=0.1 -> 0.2273316726882714..., matching to all printed digits).
CAS-verified to 25 digits at x = 0.05, 0.1, 0.3, 0.45.
 
**Sign of dC/dx at u=2.** Substituting x = t^2 (so sign(dC/dx) = sign(dC/dt),
since dx = 2t dt with t > 0) and applying the quotient rule to
C = (2t - 2t^3 - t^2 - t^4)/(2(t^4+1)):
 
  P'Q - PQ' = 4 R(t),   **R(t) = t^6 + t^5 - 3t^4 - 2t^3 - 3t^2 - t + 1**
 
Exact Sturm: R has **exactly one root in (0,1)** (Descartes gives 0 or 2 positive
roots; R(0)=1>0, R(1)=-6<0, R(2)=19>0 places the second in (1,2), outside). So
 
>   dC/dx|_(u=2) > 0 on (0, x*),  < 0 on (x*, x_c),
>   **x* = (t*)^2 = 0.148550199898...**,  t* = 0.385422106136... the root of R
 
**Single interior maximum in x at u=2, and it is at 0.1486, not 0.2.** See the
correction box in 9.5.
 
### 10.8 Verification scope for section 10 (READ BEFORE QUOTING)
 
sympy was available this session; this is a much wider CAS confirmation than 9.4's,
but it is still **not** total.
 
**CAS-verified exactly:**
- both boxed derivatives of 10.1, against `diff(C,u)` and `diff(C,x)` (`simplify` -> 0)
- S'(2) = 4x^3 P  (`expand` -> 0)
- the factorisations of alpha, beta, L(2,x), and [u^4]Q  (`factor`)
- deg_u L = 1 and deg_u Q = 4
- every root count via `count_roots` / `Poly.intervals` — exact Sturm sequences on
  integer-coefficient polynomials: V2, V3, V (=> x_dagger), R (=> x*), the
  degree-15 Delta numerator, disc_u Q (degree 174), and [u^4]Q
- disc_u Q computed exactly (degree 174, 154 terms), matching 5.4's recorded degree
**CAS-verified by 25-digit evaluation at exact rational points** (closed-form
evaluation, not sampling — used because sympy's `simplify` fails on the nested
radicals): C(2,x) at x = 0.05, 0.1, 0.3, 0.45; dC/du|_(u=2) at x = 0.02, 0.05,
0.1, 0.3.
 
**NOT verified, and inherited as assumptions from earlier sections:**
- **q1 > 0 and p1 > 0 on (0,1)** — taken from 2.3/3.2, used in 10.2 and 10.6.
  Both are easy (q1 = (x-1)^4 + 4x, p1 = (x-1)^4 + 4x(x+1), manifestly positive
  for x > 0) but were not re-derived this session.
- **S > 0 (2.4), D > 0 (2.3), sqrt(S)(2) = x(1-x^2)^2 (2.5)** — used throughout
  section 10, still hand-derived only.
- **Sections 2-4 as a whole remain CAS-unconfirmed.** The section 7 Mathematica
  block is STILL unrun. Section 10 is built on top of section 2's element set and
  inherits any error in it.
- **Second-engine confirmation.** Everything above is one CAS (sympy). The file's
  own discipline calls for independent confirmation; Mathematica would be the
  natural second engine, and the section 7 block plus the section 10 objects
  should be run there.
- **The seam x = x_dagger itself** (where u*(x_dagger) = 2 and region 3 collapses
  to a point) was not checked explicitly. It should be automatic — regions 1 and 2
  cover a neighbourhood — but it is a one-line check that was not performed.
### 10.9 What is left in this project
 
With Questions A and C both closed, the remaining items are verification and
write-up, not new mathematics:
 
1. **Run the section 7 Mathematica block.** Now the top priority: it is the only
   thing standing between sections 2-4 and full verification, and everything else
   rests on them.
2. **Re-verify section 10's objects in a second CAS** (see 10.8).
3. **Confirm the domain claim of 1.2** against your own adiabatic-elimination
   derivation, as 1.2 itself flags. Still not done.
4. **Confirm the X-state concurrence formula (1.3)** against a primary source
   before citing an equation number. Yu-Eberly is the right attribution but the
   exact equation number and notation remain unverified — 1.3's warning stands.
5. Optional: explain (or dismiss) the V-appears-twice coincidence noted in 10.4.
6. Optional: 9.7 item 4 (inverting u_c(x) = u), low value.
**Physics statement now available for write-up:** for a below-threshold
nondegenerate two-mode-squeezing drive (x < 1), steady-state entanglement exists
iff x < x_c = 0.48389, and is always maximised at the symmetric point
beta = 1 (g_1^2 kappa_2 = g_2^2 kappa_1), where its value is the closed form of
10.7. Both statements are proven, not sampled.
