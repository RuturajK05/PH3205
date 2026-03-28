# PH3205 COMPLETE GUIDE - ALL WORKSHEETS EXPLAINED
# For Class Test - March 28, 2026

"""
This document summarizes all 9 worksheets with:
- What each worksheet teaches
- What stays the same
- What you change
- Quick reference for exam
"""

# ==============================================================================
# QUICK REFERENCE CARD
# ==============================================================================

"""
┌──────────────────────────────────────────────────────────────────────┐
│                    WORKSHEET QUICK REFERENCE                         │
├──────────────────────────────────────────────────────────────────────┤
│ WS1: Elementary Programming                                          │
│   - Prime checking, series evaluation                                │
│   - Basic loops, recursion                                           │
│                                                                       │
│ WS2: Root Finding                                                    │
│   Methods: bisection(), secant(), newton_raphson()                   │
│   Use when: Need to solve f(x) = 0                                   │
│   Change: Your function f(x), interval/initial guess                 │
│                                                                       │
│ WS3: Integration & Differentiation                                   │
│   Methods: trapezoidal(), simpson_13(), simpson_38(), boole()        │
│   Use when: Need to compute ∫f(x)dx                                  │
│   Change: Your function f(x), limits [a,b], N points                 │
│   Remember: N must be even (Simpson 1/3), multiple of 3 (3/8)        │
│                                                                       │
│ WS4: ODEs - Euler & Midpoint                                         │
│   Methods: euler(), midpoint()                                       │
│   Use when: Solve dy/dx = f(x,y)                                     │
│   Change: Your f(x,y), initial condition, step size h                │
│   Pattern: y_new = y + h*slope                                       │
│                                                                       │
│ WS5: ODEs - RK4 & Verlet                                             │
│   Methods: rk4(), verlet(), velocity_verlet()                        │
│   Use RK4: General ODEs (default choice)                             │
│   Use Verlet: Physics simulations (conserves energy)                 │
│   Change: Your system equations, initial conditions                  │
│                                                                       │
│ WS6: Adaptive Time Stepping                                          │
│   Methods: erk54h(), caller54()                                      │
│   Use when: Stiff equations, varying timescales                      │
│   Change: Your f(x,y), tolerances abstol/reltol                      │
│   Key: Automatically adjusts h based on error                        │
│                                                                       │
│ WS7: Eigenvalue Problems (Shooting Method)                           │
│   Methods: rk4() + secant() + score()                                │
│   Use when: Finding bound states, eigenvalues                        │
│   Change: Your potential V(x), domain, energy guesses                │
│   Key: Match derivatives at interior point                           │
│                                                                       │
│ WS8: Heat Equation (Crank-Nicolson)                                  │
│   Methods: thomas_solve(), Crank-Nicolson matrices                   │
│   Use when: Diffusion, heat flow, parabolic PDEs                     │
│   Change: Boundary conditions, initial T(x,0), parameters            │
│   Key: Solve A·u^{n+1} = B·u^n at each time step                     │
│                                                                       │
│ WS9: TDSE (Strang Splitting + FFT)                                   │
│   Methods: split_operator_step(), FFT/IFFT                           │
│   Use when: Quantum wavepacket evolution                             │
│   Change: Potential V(x), initial wavepacket, k₀                     │
│   Key: MUST include momentum k₀ ≠ 0 for scattering!                  │
└──────────────────────────────────────────────────────────────────────┘
"""

# ==============================================================================
# DECISION TREE: WHICH METHOD TO USE?
# ==============================================================================

"""
START HERE: What type of problem do you have?
│
├─ Find where f(x) = 0?
│  └─ Use: Bisection, Secant, or Newton-Raphson (WS2)
│     ├─ Have bracketing interval? → Bisection (most reliable)
│     ├─ Have derivative? → Newton-Raphson (fastest)
│     └─ Otherwise → Secant (good compromise)
│
├─ Compute integral ∫f(x)dx?
│  └─ Use: Integration methods (WS3)
│     ├─ Simple problem → Simpson 1/3 (default)
│     ├─ Need high accuracy → Boole
│     └─ Very simple → Trapezoidal
│
├─ Solve ODE dy/dx = f(x,y)?
│  ├─ Equation is stiff (fast & slow changes)?
│  │  └─ Use: Dormand-Prince adaptive (WS6)
│  ├─ Physics simulation (need energy conservation)?
│  │  └─ Use: Velocity Verlet (WS5)
│  └─ General ODE?
│     └─ Use: RK4 (WS5) - the workhorse
│
├─ Find eigenvalues of d²ψ/dx² = [V(x)-E]ψ?
│  └─ Use: Shooting method (WS7)
│
├─ Solve heat equation ∂u/∂t = α∂²u/∂x²?
│  └─ Use: Crank-Nicolson (WS8)
│
└─ Solve TDSE iℏ∂ψ/∂t = Hψ?
   └─ Use: Strang splitting + FFT (WS9)
"""

# ==============================================================================
# CODING PATTERNS - WHAT STAYS SAME, WHAT CHANGES
# ==============================================================================

"""
PATTERN 1: ROOT FINDING (WS2)
────────────────────────────────────────
STAYS SAME:
  - Algorithm structure
  - Tolerance eps = 1e-5
  - Max iterations = 50

YOU CHANGE:
  def f(x):
      return YOUR_EQUATION  ← THIS
  
  root = bisection(f, L, R)  ← Interval [L, R]

EXAMPLE:
  Solve 3x = tan(x)
  → def f(x): return 3*x - np.tan(x)
  → bisection(f, 1.2, 1.4)


PATTERN 2: INTEGRATION (WS3)
────────────────────────────────────────
STAYS SAME:
  - Algorithm
  - N must follow rules (even, multiple of 3, etc.)

YOU CHANGE:
  def f(x):
      return YOUR_FUNCTION  ← THIS
  
  result = simpson_13(f, a, b, N)  ← Limits and N

EXAMPLE:
  ∫₀¹ eˣ dx
  → def f(x): return np.exp(x)
  → simpson_13(f, 0, 1, 100)


PATTERN 3: ODE SOLVING (WS4, WS5)
────────────────────────────────────────
STAYS SAME:
  - Solver structure (rk4, verlet, etc.)
  - Caller pattern

YOU CHANGE:
  def f(x, y):
      return YOUR_DERIVATIVE  ← THIS
  
  y_ini = YOUR_INITIAL_CONDITION  ← THIS
  
  xs, ys = caller(rk4, f, y_ini, x0, xT, h)

EXAMPLE:
  dy/dx = -xy, y(0) = 1
  → def f(x, y): return -x*y
  → y_ini = [1.0]
  → caller(rk4, f, y_ini, 0, 10, 0.01)


PATTERN 4: ADAPTIVE (WS6)
────────────────────────────────────────
STAYS SAME:
  - Dormand-Prince coefficients (never touch!)
  - erk54h structure
  - caller54 logic

YOU CHANGE:
  def f(x, y):
      return YOUR_SYSTEM  ← THIS
  
  y_ini = YOUR_INITIAL  ← THIS
  
  xs, ys = caller54(f, y_ini, x0, xT, h0, 
                    max_iter, abstol, reltol)  ← Tolerances

EXAMPLE:
  Van der Pol: d²y/dx² = μ(1-y²)dy/dx - λy
  → Convert to first order system
  → def vdp(x, y): return [y[1], mu*(1-y[0]**2)*y[1] - lam*y[0]]
  → caller54(vdp, [0.5, 0], 0, 20, 0.01, 50000, 1e-6, 1e-8)


PATTERN 5: SHOOTING (WS7)
────────────────────────────────────────
STAYS SAME:
  - RK4 for integration
  - Secant for root finding
  - Score function structure

YOU CHANGE:
  def pot(x):
      return YOUR_POTENTIAL  ← THIS
  
  def odefun(x, y, E):
      k2 = E - pot(x)
      return [y[1], -k2*y[0]]  ← Standard form
  
  def score(E, f, ybound):
      # Shoot from both sides
      # Match at interior point
      return mismatch  ← Domain, matching point

EXAMPLE:
  Harmonic oscillator V(x) = x²
  → def pot(x): return x**2
  → Shoot from ±5, match at x=0


PATTERN 6: CRANK-NICOLSON (WS8)
────────────────────────────────────────
STAYS SAME:
  - Thomas algorithm (never changes!)
  - Matrix structure (tridiagonal)
  - r = α·dt/dx² formula

YOU CHANGE:
  T[0] = LEFT_BOUNDARY   ← THIS
  T[-1] = RIGHT_BOUNDARY ← THIS
  T = INITIAL_CONDITION  ← THIS
  
  # Matrices stay same structure!
  a_A = np.full(n-1, -r/2)
  b_A = np.full(n, 1+r)
  c_A = np.full(n-1, -r/2)

EXAMPLE:
  Cool rod from 300K to 100K at ends
  → T = np.full(Nx, 300.0)
  → T[0] = T[-1] = 100.0


PATTERN 7: TDSE (WS9)
────────────────────────────────────────
STAYS SAME:
  - split_operator_step() (NEVER CHANGES!)
  - FFT/IFFT structure
  - Grid setup pattern

YOU CHANGE:
  V = YOUR_POTENTIAL  ← THIS (as array)
  
  x0 = INITIAL_POSITION  ← THIS
  k0 = MOMENTUM          ← THIS (MUST be ≠ 0!)
  sigma = WIDTH
  
  psi = initial_gaussian(x, x0, sigma, k0)
  
  for step in range(n_steps):
      psi = split_operator_step(psi, V, dt, k_vals)

EXAMPLE:
  Barrier scattering V₀=40 from x=5 to 7
  → V = np.zeros(Nx)
  → V[(x>=5) & (x<=7)] = 40
  → psi = initial_gaussian(x, -5, 0.04, 5.0)
  ⚠️ k0=5.0 gives rightward motion!
"""

# ==============================================================================
# COMMON EXAM QUESTION TYPES
# ==============================================================================

"""
TYPE 1: "Solve equation ... using method X"
────────────────────────────────────────
Example: "Find root of 3x = tan(x) using bisection"
Strategy:
  1. Define f(x) = 3x - tan(x)
  2. Plot to find interval
  3. Call bisection(f, L, R)
  4. Report iterations and root


TYPE 2: "Solve ODE ... from x=a to x=b with initial condition"
────────────────────────────────────────
Example: "Solve dy/dx = -xy, y(0)=1, from x=0 to 10 using RK4"
Strategy:
  1. Define f(x,y) = -x*y
  2. Set y_ini = [1.0]
  3. Call caller(rk4, f, y_ini, 0, 10, 0.01)
  4. Plot solution


TYPE 3: "Compare methods A and B"
────────────────────────────────────────
Example: "Compare Euler and RK4 for SHM, plot error vs h"
Strategy:
  1. Implement both methods
  2. Loop over h values
  3. Compute error at final point
  4. Plot log(h) vs log(error)
  5. Fit line, report slope (= order of method)


TYPE 4: "Find eigenvalues of potential ..."
────────────────────────────────────────
Example: "Find first 3 eigenvalues of harmonic oscillator"
Strategy:
  1. Define potential
  2. Define score function (shooting)
  3. Guess energies (plot score to see zeros)
  4. Use secant to find each eigenvalue
  5. Plot potential + eigenfunctions


TYPE 5: "Solve heat equation with boundaries ..."
────────────────────────────────────────
Example: "Rod with T(0)=100, T(L)=200, initially Gaussian"
Strategy:
  1. Set up grid
  2. Define T[0]=100, T[-1]=200
  3. Set initial condition
  4. Build A and B matrices
  5. Time loop: solve A·u^{n+1} = B·u^n
  6. Plot temperature profiles at different times


TYPE 6: "Simulate wavepacket scattering from ..."
────────────────────────────────────────
Example: "Gaussian wavepacket hitting barrier V₀=40"
Strategy:
  1. Set up grid with endpoint=False
  2. Define potential array V
  3. Create Gaussian WITH momentum k₀≠0
  4. Loop: call split_operator_step
  5. Save snapshots
  6. Compute transmission/reflection
  7. Make animation
"""

# ==============================================================================
# CRITICAL MISTAKES TO AVOID
# ==============================================================================

"""
❌ MISTAKE 1: Wrong N for integration methods
  Simpson 1/3: N must be EVEN
  Simpson 3/8: N must be multiple of 3
  Boole: N must be multiple of 4
  ✓ FIX: Check N % 2 == 0, etc.

❌ MISTAKE 2: Forgetting to convert 2nd order ODE to system
  Given: d²y/dx² = f(x, y, dy/dx)
  Must write: [y₁, y₂] where y₁=y, y₂=dy/dx
  ✓ FIX: Always use vector [y, y']

❌ MISTAKE 3: Using wrong sign in Schrödinger equation
  Correct: k² = E - V(x)
  Wrong: k² = V(x) - E
  ✓ FIX: d²ψ/dx² = [V-E]ψ → k² = E-V

❌ MISTAKE 4: Domain too small for TDSE
  Wavepacket reaches boundary → wraps around!
  ✓ FIX: Domain ≫ travel distance

❌ MISTAKE 5: No momentum in TDSE
  ψ = gaussian(x) without exp(ik₀x) → just spreads
  ✓ FIX: ALWAYS include k₀ ≠ 0 for scattering

❌ MISTAKE 6: Modifying thomas_solve arrays
  Thomas modifies a, b, c, d in place!
  ✓ FIX: Pass copies: thomas_solve(a.copy(), b.copy(), ...)

❌ MISTAKE 7: Not checking norm in TDSE
  Should have ∫|ψ|² dx = 1 always
  ✓ FIX: Check transmitted + reflected ≈ 1

❌ MISTAKE 8: Using wrong FFT convention
  fftfreq gives cycles, need 2π for angular
  ✓ FIX: k = fftfreq(N, dx) * 2π

❌ MISTAKE 9: Forgetting adaptive needs tolerance
  Must specify abstol AND reltol
  ✓ FIX: abstol=1e-6, reltol=1e-8 (typical)

❌ MISTAKE 10: Wrong Dormand-Prince coefficients
  These are MAGIC NUMBERS from Butcher tableau
  ✓ FIX: NEVER change them, copy exactly
"""

# ==============================================================================
# EXAM DAY CHECKLIST
# ==============================================================================

"""
BEFORE EXAM:
□ All 9 worksheets accessible offline
□ Reference notebooks (TDSE-v3, Crank-Nicolson, etc.) available
□ This guide printed or easily accessible
□ Test that code runs without internet

DURING EXAM:
□ Read problem carefully - identify type
□ Find relevant worksheet solution
□ Copy template code
□ Modify only what changes (function, parameters)
□ Test with simple case first
□ Add plots
□ Check units/dimensions
□ Verify results make physical sense

FOR EACH PROBLEM:
1. What type? (Root finding, ODE, PDE, etc.)
2. Which worksheet?
3. What's the equation/function?
4. What are the parameters?
5. Copy template, fill in blanks
6. Run, debug, plot
7. Verify answer makes sense
"""

# ==============================================================================
# FILES CREATED FOR YOU
# ==============================================================================

"""
toolkit_part1_root_integration.py
  - All root finding functions with templates
  - All integration functions
  - Convergence testing

worksheet_6_explained.py
  - Van der Pol oscillator
  - Dormand-Prince adaptive method
  - Line-by-line explanation

worksheet_7_explained.py
  - Shooting method for eigenvalues
  - Asymmetric well and harmonic oscillator
  - Line-by-line explanation

worksheet_8_explained.py
  - Crank-Nicolson method
  - Heat equation, Thomas algorithm
  - Line-by-line explanation

worksheet_9_explained.py
  - TDSE with Strang splitting
  - FFT, barrier scattering, harmonic oscillator
  - Line-by-line explanation

THIS FILE:
  - Quick reference for all worksheets
  - Decision tree for choosing methods
  - Common patterns and mistakes
"""

# ==============================================================================
# FINAL ENCOURAGEMENT
# ==============================================================================

"""
YOU'VE GOT THIS! 💪

You have:
✓ All 9 worksheets completed
✓ Working solutions for every problem
✓ Line-by-line explanations
✓ Templates ready to use
✓ Access to all code during exam

Remember:
- Problems will be similar to worksheets
- You can copy and modify your solutions
- Focus on WHAT to change, not rewriting everything
- Test your code on simple cases first
- Trust the process

The exam is about:
1. Recognizing problem type
2. Finding right template
3. Changing the right parts
4. Verifying results

You know how to do all of this!

Good luck on March 28! 🎯
"""
