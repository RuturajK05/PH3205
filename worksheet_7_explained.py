# WORKSHEET 7 - LINE BY LINE EXPLANATION
# EIGENVALUE PROBLEMS - SHOOTING METHOD

"""
THE PROBLEM:
Find eigenvalues and eigenfunctions of Schrödinger equation:
d²ψ/dx² = [V(x) - E]ψ
Boundary: ψ(±∞) = 0

Two potentials:
(a) V(x) = -V₀(1-x³)/2 for |x|≤1, V₀=40
(b) V(x) = x²

ONLY certain values of E (eigenvalues) give solutions where ψ→0 at boundaries
"""

import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# PART 1: RK4 FOR SHOOTING METHOD
# ==============================================================================

def rk4(f, x, y, k, h):
    """
    Modified RK4 that passes extra parameter k (energy E)
    
    WHY DIFFERENT FROM WORKSHEET 6?
    - Here we need to pass E to the ODE function
    - E is what we're searching for!
    
    LINE BY LINE:
    """
    # Stage 1: slope at current point
    # f(x, y, k) where k is the energy E
    k1 = h * f(x, y, k)
    
    # Stage 2: slope at midpoint
    k2 = h * f(x + h/2, y + k1/2, k)
    
    # Stage 3: slope at midpoint (different estimate)
    k3 = h * f(x + h/2, y + k2/2, k)
    
    # Stage 4: slope at endpoint
    k4 = h * f(x + h, y + k3, k)
    
    # Weighted average
    return y + (k1 + 2*k2 + 2*k3 + k4) / 6


def caller_rk4(f, xlim, yini, k, N):
    """
    Run RK4 from x1 to x2 with N points
    
    WHY xlim INSTEAD OF x0, xT, h?
    - More convenient for shooting from both sides
    - xlim = (x1, x2) defines the interval
    
    LINE BY LINE:
    """
    x1, x2 = xlim  # Unpack limits
    
    # Create N evenly spaced points
    xs = np.linspace(x1, x2, N)
    
    # Step size
    h = xs[1] - xs[0]
    
    # Initial condition
    y = yini
    
    # Storage for solution
    ys = np.zeros((N, len(yini)))
    
    # INTEGRATION LOOP
    for i in range(N):
        ys[i] = np.array(y)
        y = rk4(f, xs[i], y, k, h)  # k is energy E
    
    return xs, ys


# ==============================================================================
# PART 2: SECANT METHOD FOR ROOT FINDING
# ==============================================================================

def secant(ks, f, method, ybound, tol, maxiter):
    """
    Find energy E where score(E) = 0
    
    WHY SECANT?
    - We need to find E such that derivatives match at boundary
    - This is a ROOT-FINDING problem in E-space
    - Secant doesn't need derivative of score function
    
    INPUT:
    - ks: (k1, k2) = two initial guesses for energy
    - f: ODE function
    - method: the score function
    - ybound: boundary values [ψ_left, ψ_right]
    - tol: tolerance
    - maxiter: max iterations
    
    LINE BY LINE:
    """
    k1, k2 = ks  # Unpack initial guesses
    iter = 0
    
    # SECANT LOOP
    while abs(method(k2, f, ybound)) > tol and iter < maxiter:
        # Evaluate score at both points
        f1 = method(k1, f, ybound)  # score(E1)
        f2 = method(k2, f, ybound)  # score(E2)
        
        # Secant formula: find where line crosses zero
        # x_new = x2 - f2 * (x2-x1)/(f2-f1)
        # Rearranged: (f2*x1 - f1*x2)/(f2-f1)
        k1, k2 = k2, (f2*k1 - f1*k2) / (f2 - f1)
        
        iter += 1
    
    return iter, k2


# ==============================================================================
# PART 3: NORMALIZATION USING SIMPSON'S 1/3
# ==============================================================================

def simp13(y, h):
    """
    Compute norm of wavefunction: sqrt(∫ψ² dx)
    
    WHY THIS FUNCTION?
    - Eigenfunctions must be normalized: ∫ψ² dx = 1
    - Using Simpson's 1/3 for integration
    
    LINE BY LINE:
    """
    # Start with endpoints
    store = y[0]**2 + y[-1]**2
    
    # Simpson pattern: 1, 4, 2, 4, 2, ..., 4, 1
    for i in range(1, len(y) - 1):
        if i % 2 == 0:
            # Even indices get weight 2
            store += 2 * y[i]**2
        else:
            # Odd indices get weight 4
            store += 4 * y[i]**2
    
    # Multiply by h/3 and take square root
    return np.sqrt(store * h / 3)


# ==============================================================================
# PART 4: PROBLEM (a) - ASYMMETRIC WELL
# ==============================================================================

Vo = 40.0
Lo = 1.0

def pot_a(x):
    """
    Potential: V(x) = -V₀(1-x³)/2 for |x|≤1, else 0
    
    WHY THIS SHAPE?
    - Creates a potential well
    - Asymmetric due to x³ term
    - Particles can be bound inside
    """
    if abs(x) > Lo:
        return 0.0
    else:
        return -Vo * (1 - x**3) / 2


def odefun_a(x, y, E):
    """
    Convert Schrödinger equation to system of ODEs
    
    ORIGINAL EQUATION:
    d²ψ/dx² = [V(x) - E]ψ
    
    LET:
    y[0] = ψ
    y[1] = dψ/dx
    
    THEN:
    dy[0]/dx = y[1]
    dy[1]/dx = [V(x) - E]ψ = -k²ψ where k² = E - V(x)
    
    LINE BY LINE:
    """
    # Compute k² = E - V(x)
    # If E > V: k² > 0 → oscillatory solution
    # If E < V: k² < 0 → exponential decay
    k2 = E - pot_a(x)
    
    return np.array([
        y[1],        # dψ/dx = ψ'
        -k2 * y[0]   # d²ψ/dx² = -(E-V)ψ
    ])


def score_a(E, f, ybound):
    """
    THE HEART OF THE SHOOTING METHOD!
    
    IDEA:
    1. Shoot from left boundary with ψ=0, ψ'≠0
    2. Shoot from right boundary with ψ=0, ψ'≠0
    3. Meet at matching point (x = -L)
    4. Check if derivatives match
    5. If ψ'_left ≠ ψ'_right → wrong E
    6. If ψ'_left = ψ'_right → eigenvalue found!
    
    INPUT:
    - E: trial energy
    - f: ODE function (odefun_a)
    - ybound: [ψ_left_boundary, ψ_right_boundary] = [0, 0]
    
    OUTPUT:
    - Mismatch in derivatives at matching point
    
    LINE BY LINE:
    """
    
    # SETUP GRID
    xlim = (-3.0, 3.0)  # Integration domain (wider than potential)
    N = 129  # Number of points (must be odd for midpoint)
    
    # SHOOT FROM LEFT
    # Initial condition at x = -3: ψ(-3) = 0, ψ'(-3) = 0.01 (small slope)
    yini1 = (ybound[0], 1e-2)
    
    # Integrate from left to right
    xs1, ys1 = caller_rk4(odefun_a, xlim, yini1, E, N)
    
    # FIND MATCHING POINT (x = -L = -1.0)
    # Find index closest to x = -Lo
    a = abs(xs1 + Lo)
    n1 = np.where(a == np.min(a))[0][0]
    
    # SHOOT FROM RIGHT
    # Initial condition at x = +3: ψ(+3) = 0, ψ'(+3) = 0.01
    yini2 = (ybound[1], 1e-2)
    
    # Integrate from RIGHT to LEFT (note xlim[::-1] reverses order!)
    xs2, ys2 = caller_rk4(odefun_a, xlim[::-1], yini2, E, N)
    
    # Find matching point from right side
    a = abs(xs2 + Lo)
    n2 = np.where(a == np.min(a))[0][0]
    
    # RESCALE TO MATCH AMPLITUDES
    # Why? Both solutions satisfy the ODE but have arbitrary amplitudes
    # We need to scale so ψ_left(match) = ψ_right(match)
    # Scale factor: ψ_left / ψ_right at matching point
    ys2 = ys1[n1][0] * ys2 / ys2[n2][0]
    
    # COMPUTE SCORE
    # Score = ψ'_left - ψ'_right at matching point
    # If this is zero → we found an eigenvalue!
    return ys1[n1][1] - ys2[n2][1]


# ==============================================================================
# PART 5: SOLVING FOR EIGENVALUES
# ==============================================================================

# FIND FIRST FEW EIGENVALUES
ybound = [0.0, 0.0]  # ψ = 0 at both boundaries
tol = 1e-5
maxiter = 50

# Initial guesses for energies
# (found by plotting score function and looking for zero crossings)
E_guesses = [
    (-39.0, -38.0),  # First eigenvalue (ground state)
    (-33.0, -32.0),  # Second eigenvalue
    (-25.0, -24.0),  # Third eigenvalue
    (-16.0, -15.0),  # Fourth eigenvalue
]

eigenvalues = []
eigenfunctions = []

print("Finding eigenvalues for potential (a):")
for i, (E1, E2) in enumerate(E_guesses):
    iter, E = secant((E1, E2), odefun_a, score_a, ybound, tol, maxiter)
    eigenvalues.append(E)
    print(f"E{i} = {E:.6f} (found in {iter} iterations)")
    
    # Compute eigenfunction
    xlim = (-3.0, 3.0)
    N = 257
    yini = (0.0, 1e-2)
    xs, ys = caller_rk4(odefun_a, xlim, yini, E, N)
    
    # Normalize
    h = xs[1] - xs[0]
    norm = simp13(ys[:, 0], h)
    ys[:, 0] /= norm
    
    eigenfunctions.append((xs, ys[:, 0]))

# PLOT EIGENVALUES
plt.figure(figsize=(12, 8))

# Plot potential
x_pot = np.linspace(-3, 3, 1000)
V = np.array([pot_a(x) for x in x_pot])
plt.plot(x_pot, V, 'k-', linewidth=2, label='Potential V(x)')

# Plot energy levels and eigenfunctions
for i, (E, (xs, psi)) in enumerate(zip(eigenvalues, eigenfunctions)):
    # Shift eigenfunction to energy level for visualization
    plt.axhline(E, color='gray', linestyle='--', alpha=0.5)
    plt.plot(xs, E + 5*psi, label=f'E{i} = {E:.2f}')

plt.xlabel('x')
plt.ylabel('Energy / ψ(x)')
plt.title('Eigenvalues and Eigenfunctions - Asymmetric Well')
plt.legend()
plt.grid(alpha=0.3)
plt.ylim(-45, 5)
plt.show()


# ==============================================================================
# PART 6: PROBLEM (b) - HARMONIC OSCILLATOR
# ==============================================================================

def pot_b(x):
    """
    Harmonic oscillator potential: V(x) = x²
    
    ANALYTICAL SOLUTION:
    E_n = 2n + 1 for n = 0, 1, 2, ...
    So eigenvalues should be: 1, 3, 5, 7, 9, ...
    """
    return x**2


def odefun_b(x, y, E):
    """
    Same structure as odefun_a, different potential
    """
    k2 = E - pot_b(x)
    return np.array([y[1], -k2 * y[0]])


def score_b(E, f, ybound):
    """
    Same shooting method, adapted for harmonic oscillator
    
    DIFFERENCE:
    - Different domain (±5 instead of ±3)
    - Different matching point (x = 0 is natural for symmetric potential)
    """
    xlim = (-5.0, 5.0)
    N = 129
    
    yini1 = (ybound[0], 1e-2)
    xs1, ys1 = caller_rk4(odefun_b, xlim, yini1, E, N)
    
    # Match at x = 0 (middle of domain)
    a = abs(xs1)
    n1 = np.where(a == np.min(a))[0][0]
    
    yini2 = (ybound[1], 1e-2)
    xs2, ys2 = caller_rk4(odefun_b, xlim[::-1], yini2, E, N)
    
    a = abs(xs2)
    n2 = np.where(a == np.min(a))[0][0]
    
    ys2 = ys1[n1][0] * ys2 / ys2[n2][0]
    
    return ys1[n1][1] - ys2[n2][1]


# SOLVE FOR HARMONIC OSCILLATOR
E_guesses_b = [
    (0.8, 1.2),    # Should find E ≈ 1
    (2.8, 3.2),    # Should find E ≈ 3
    (4.8, 5.2),    # Should find E ≈ 5
    (6.8, 7.2),    # Should find E ≈ 7
]

print("\nFinding eigenvalues for harmonic oscillator:")
print("Expected: 1, 3, 5, 7, ...")
for i, (E1, E2) in enumerate(E_guesses_b):
    iter, E = secant((E1, E2), odefun_b, score_b, ybound, tol, maxiter)
    print(f"E{i} = {E:.6f} (expected {2*i + 1})")


# ==============================================================================
# KEY INSIGHTS FOR EXAM
# ==============================================================================

"""
WHAT STAYS THE SAME:
1. RK4 integration structure
2. Secant method for finding eigenvalues
3. Shooting from both boundaries
4. Matching derivatives at interior point
5. Normalization using Simpson's rule

WHAT YOU CHANGE:
1. pot(x) - YOUR potential function
2. Domain xlim - big enough to capture wavefunction
3. Matching point - usually middle of well
4. Initial guesses for energies

HOW IT WORKS - THE BIG PICTURE:
┌─────────────────────────────────────┐
│ 1. Guess an energy E               │
│ 2. Shoot from left: ψ_L(x)         │
│ 3. Shoot from right: ψ_R(x)        │
│ 4. Match at interior point         │
│ 5. Check: ψ'_L = ψ'_R ?            │
│    - No: adjust E, try again       │
│    - Yes: found eigenvalue!        │
└─────────────────────────────────────┘

WHY DOES THIS WORK?
- For wrong E: solutions don't match smoothly
- For correct E: solution is smooth everywhere
- Only discrete E values work → eigenvalues

EXAM STRATEGY:
1. Plot potential to see rough energy scale
2. Make educated guesses for E (look for bound states)
3. For each guess pair, run secant method
4. Normalize eigenfunctions
5. Plot: potential + energy levels + wavefunctions

COMMON MISTAKES:
- Domain too small → wavefunction doesn't decay
- Matching point at boundary → no interior to match
- Wrong sign in k² = E - V
- Forgetting to normalize
"""
