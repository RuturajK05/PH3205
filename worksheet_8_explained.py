# WORKSHEET 8 - LINE BY LINE EXPLANATION
# HEAT EQUATION - CRANK-NICOLSON METHOD

"""
THE PROBLEM:
Heat equation (diffusion): ∂u/∂t = α ∂²u/∂x²

Two scenarios:
(a) Rod with both ends at T₀ = 100K, initially at Tᵢₙᵢₜ = 300K
(b) Left end at 200K, right end at 400K, Gaussian initial condition

METHOD: Crank-Nicolson (implicit, 2nd order accurate, unconditionally stable)
"""

import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# PART 1: THOMAS ALGORITHM (TRIDIAGONAL SOLVER)
# ==============================================================================

def thomas_solve(a, b, c, d):
    """
    Solve tridiagonal system: A·x = d
    
    Matrix A looks like:
    [b0  c0   0   0  ...   0 ]   [x0]   [d0]
    [a0  b1  c1   0  ...   0 ]   [x1]   [d1]
    [ 0  a1  b2  c2  ...   0 ] × [x2] = [d2]
    [ .   .   .   .  ...   . ]   [.]    [.]
    [ 0   0   0   0  ... bn-1]   [xn]   [dn]
    
    WHY TRIDIAGONAL?
    - Heat equation with nearest-neighbor interactions
    - Sparse matrix → special efficient algorithm
    - O(n) instead of O(n³) for general matrix!
    
    THOMAS ALGORITHM:
    1. Forward elimination: eliminate lower diagonal
    2. Back substitution: solve for unknowns
    
    INPUT:
    - a: lower diagonal [a0, a1, ..., a_{n-2}] (length n-1)
    - b: main diagonal [b0, b1, ..., b_{n-1}] (length n)
    - c: upper diagonal [c0, c1, ..., c_{n-2}] (length n-1)
    - d: right hand side [d0, d1, ..., d_{n-1}] (length n)
    
    OUTPUT:
    - x: solution vector
    
    LINE BY LINE:
    """
    n = len(b)  # Number of unknowns
    
    # FORWARD ELIMINATION
    # Eliminate lower diagonal by row operations
    for i in range(1, n):
        # Compute multiplier: how much to subtract row i-1 from row i
        w = a[i-1] / b[i-1]
        
        # Update diagonal: b_i = b_i - w*c_{i-1}
        b[i] = b[i] - w * c[i-1]
        
        # Update RHS: d_i = d_i - w*d_{i-1}
        d[i] = d[i] - w * d[i-1]
    
    # Now matrix is upper triangular!
    
    # BACK SUBSTITUTION
    x = np.zeros(n)
    
    # Last equation: b_{n-1}*x_{n-1} = d_{n-1}
    x[-1] = d[-1] / b[-1]
    
    # Work backwards
    for i in range(n-2, -1, -1):
        # b_i*x_i + c_i*x_{i+1} = d_i
        # x_i = (d_i - c_i*x_{i+1}) / b_i
        x[i] = (d[i] - c[i] * x[i+1]) / b[i]
    
    return x


# ==============================================================================
# PART 2: CRANK-NICOLSON DISCRETIZATION
# ==============================================================================

"""
HEAT EQUATION: ∂u/∂t = α ∂²u/∂x²

CRANK-NICOLSON: Average of implicit and explicit
u^{n+1} - u^n     α   [∂²u^{n+1}   ∂²u^n]
─────────────── = ─── [───────── + ─────]
      Δt          2   [   ∂x²        ∂x² ]

DISCRETIZE SPACE: u_i^n ≈ u(x_i, t_n) where x_i = i·Δx

SECOND DERIVATIVE:
∂²u     u_{i+1} - 2u_i + u_{i-1}
──── ≈ ─────────────────────────
∂x²              Δx²

SUBSTITUTE:
u_i^{n+1} - u_i^n     α   [u_{i+1}^{n+1} - 2u_i^{n+1} + u_{i-1}^{n+1}
───────────────────  = ─── [─────────────────────────────────────────
        Δt            2   [              Δx²

                          u_{i+1}^n - 2u_i^n + u_{i-1}^n]
                        + ─────────────────────────────]
                                      Δx²                ]

DEFINE: r = α·Δt/Δx²

REARRANGE:
-r/2·u_{i-1}^{n+1} + (1+r)·u_i^{n+1} - r/2·u_{i+1}^{n+1} 
    = r/2·u_{i-1}^n + (1-r)·u_i^n + r/2·u_{i+1}^n

LEFT SIDE: Future time (unknown) → Matrix A
RIGHT SIDE: Current time (known) → Matrix B

MATRIX FORM: A·u^{n+1} = B·u^n
"""


# ==============================================================================
# PART 3: PROBLEM (a) - COOLING ROD
# ==============================================================================

print("=" * 60)
print("PROBLEM (a): COOLING ROD")
print("=" * 60)

# PARAMETERS
L = 1.0        # Rod length
Nx = 51        # Number of spatial points
dx = L / (Nx - 1)  # Spatial step
alpha = 1e-4   # Thermal diffusivity
dt = 0.1       # Time step
Nt = 5000      # Number of time steps
r = alpha * dt / dx**2  # Crank-Nicolson parameter

print(f"Grid: {Nx} points, dx = {dx:.4f}")
print(f"Time: {Nt} steps, dt = {dt:.4f}")
print(f"r = α·dt/dx² = {r:.6f}")

# SPATIAL GRID
x = np.linspace(0, L, Nx)

# BOUNDARY CONDITIONS
T0 = 100.0  # Temperature at both ends
Tinit = 300.0  # Initial temperature

# INITIAL CONDITION
T = np.full(Nx, Tinit)  # All points at Tinit
T[0] = T0   # Left boundary
T[-1] = T0  # Right boundary

print(f"\nInitial: T(interior) = {Tinit}K, T(boundaries) = {T0}K")

# SETUP MATRICES
# We only solve for INTERIOR points (indices 1 to Nx-2)
# Boundary points are fixed!

n_interior = Nx - 2  # Number of interior points

# MATRIX A (left side): tridiagonal
# Diagonal pattern: [-r/2, 1+r, -r/2]
a_A = np.full(n_interior - 1, -r/2)  # Lower diagonal
b_A = np.full(n_interior, 1 + r)     # Main diagonal
c_A = np.full(n_interior - 1, -r/2)  # Upper diagonal

# MATRIX B (right side): also tridiagonal
# Diagonal pattern: [r/2, 1-r, r/2]
a_B = np.full(n_interior - 1, r/2)
b_B = np.full(n_interior, 1 - r)
c_B = np.full(n_interior - 1, r/2)

"""
MATRIX A STRUCTURE (for n_interior = 4):
[1+r  -r/2    0     0  ]
[-r/2  1+r  -r/2    0  ]
[ 0   -r/2   1+r  -r/2 ]
[ 0     0   -r/2   1+r ]

MATRIX B STRUCTURE:
[1-r   r/2    0     0  ]
[ r/2  1-r   r/2    0  ]
[ 0    r/2   1-r   r/2 ]
[ 0     0    r/2   1-r ]
"""

# STORAGE FOR TIME SNAPSHOTS
snapshots_time = []
snapshots_temp = []

# TIME-STEPPING LOOP
for n in range(Nt):
    # RHS: B·u^n for interior points
    # This is matrix-vector multiplication
    
    d = np.zeros(n_interior)
    
    # COMPUTE RHS = B·u^n
    for i in range(n_interior):
        # Interior point index in full array
        idx = i + 1
        
        # Multiply row i of B with current temperature
        d[i] = b_B[i] * T[idx]  # Diagonal term
        
        if i > 0:  # Not first interior point
            d[i] += a_B[i-1] * T[idx-1]  # Lower diagonal
        
        if i < n_interior - 1:  # Not last interior point
            d[i] += c_B[i] * T[idx+1]  # Upper diagonal
    
    # BOUNDARY CONDITIONS
    # Add contribution from fixed boundaries
    # First interior point (i=1) neighbors T[0] = T0
    d[0] += (r/2) * T[0]
    
    # Last interior point (i=Nx-2) neighbors T[Nx-1] = T0
    d[-1] += (r/2) * T[-1]
    
    # SOLVE: A·u^{n+1} = d
    T_interior = thomas_solve(a_A.copy(), b_A.copy(), c_A.copy(), d.copy())
    
    # UPDATE TEMPERATURE
    T[1:-1] = T_interior
    # Boundaries stay fixed: T[0] = T0, T[-1] = T0
    
    # SAVE SNAPSHOTS
    if n % 500 == 0:
        snapshots_time.append(n * dt)
        snapshots_temp.append(T.copy())

# FINAL STATE
snapshots_time.append(Nt * dt)
snapshots_temp.append(T.copy())

# PLOT EVOLUTION
plt.figure(figsize=(10, 6))
for i, (t, temp) in enumerate(zip(snapshots_time, snapshots_temp)):
    plt.plot(x, temp, label=f't = {t:.1f}', linewidth=2)

plt.xlabel('Position x')
plt.ylabel('Temperature (K)')
plt.title('Cooling Rod: Both Ends at T₀ = 100K')
plt.legend()
plt.grid(alpha=0.3)
plt.show()

print(f"\nInitial center temperature: {Tinit:.2f} K")
print(f"Final center temperature: {T[Nx//2]:.2f} K")
print(f"Boundary temperature: {T0:.2f} K")
print(f"\nAs t → ∞, the rod equilibrates to T₀ = {T0:.2f} K everywhere.")


# ==============================================================================
# PART 4: PROBLEM (b) - DIFFERENT BOUNDARY TEMPERATURES
# ==============================================================================

print("\n" + "=" * 60)
print("PROBLEM (b): DIFFERENT BOUNDARY TEMPERATURES")
print("=" * 60)

# NEW PARAMETERS
Tleft = 200.0   # Left boundary
Tright = 400.0  # Right boundary
sigma = 0.05

# GAUSSIAN INITIAL CONDITION
# T(x,0) = 300 + exp(-(x-L/2)²/(2σ²))
T_new = 300 + np.exp(-(x - L/2)**2 / (2 * sigma**2))

# SET BOUNDARIES
T_new[0] = Tleft
T_new[-1] = Tright

print(f"Left boundary: {Tleft}K")
print(f"Right boundary: {Tright}K")
print(f"Initial condition: Gaussian peak at center")

# TIME-STEPPING (same as before but different boundary handling)
snapshots_time_b = []
snapshots_temp_b = []

for n in range(Nt):
    # COMPUTE RHS
    d = np.zeros(n_interior)
    
    for i in range(n_interior):
        idx = i + 1
        d[i] = b_B[i] * T_new[idx]
        
        if i > 0:
            d[i] += a_B[i-1] * T_new[idx-1]
        
        if i < n_interior - 1:
            d[i] += c_B[i] * T_new[idx+1]
    
    # DIFFERENT BOUNDARIES!
    # Left: T[0] = Tleft
    d[0] += (r/2) * Tleft
    
    # Right: T[-1] = Tright
    d[-1] += (r/2) * Tright
    
    # SOLVE
    T_interior = thomas_solve(a_A.copy(), b_A.copy(), c_A.copy(), d.copy())
    
    # UPDATE
    T_new[1:-1] = T_interior
    T_new[0] = Tleft   # Keep boundaries fixed
    T_new[-1] = Tright
    
    # SNAPSHOTS
    if n % 500 == 0:
        snapshots_time_b.append(n * dt)
        snapshots_temp_b.append(T_new.copy())

snapshots_time_b.append(Nt * dt)
snapshots_temp_b.append(T_new.copy())

# PLOT
plt.figure(figsize=(12, 8))
for i, (t, temp) in enumerate(zip(snapshots_time_b, snapshots_temp_b)):
    plt.plot(x, temp, label=f't = {t:.1f}', linewidth=2)

# STEADY STATE (linear profile)
T_steady = Tleft + (Tright - Tleft) * x / L
plt.plot(x, T_steady, 'k--', linewidth=2, label='Steady state (linear)')

plt.xlabel('Position x')
plt.ylabel('Temperature (K)')
plt.title('Heating: Left=200K, Right=400K')
plt.legend()
plt.grid(alpha=0.3)
plt.show()

print(f"\nSteady state: Linear temperature profile")
print(f"T(x) = {Tleft} + {Tright-Tleft}·x")


# ==============================================================================
# KEY INSIGHTS FOR EXAM
# ==============================================================================

"""
WHAT STAYS THE SAME:
1. Thomas algorithm - NEVER CHANGES
2. Matrix structure (tridiagonal)
3. r = α·dt/dx² definition
4. Crank-Nicolson weights: -r/2, 1±r, ±r/2

WHAT YOU CHANGE:
1. Boundary conditions (Dirichlet, Neumann, etc.)
2. Initial condition T(x,0)
3. Parameters: L, α, dx, dt
4. Domain and number of points

CRANK-NICOLSON ADVANTAGES:
✓ Unconditionally stable (can use large dt)
✓ 2nd order accurate in time and space
✓ Implicit method (need to solve system)

WHY TRIDIAGONAL?
- Heat equation couples each point to nearest neighbors only
- Sparse matrix → efficient solver
- Thomas algorithm: O(n) complexity

BOUNDARY CONDITIONS:
- Dirichlet: u(0,t) = u₀, u(L,t) = u_L (fixed temperature)
- Neumann: ∂u/∂x = 0 at boundaries (insulated)
- Robin: combination of value and derivative

STEADY STATE:
- When ∂u/∂t = 0, we get ∂²u/∂x² = 0
- Solution: u(x) = Ax + B (linear!)
- Determined by boundary conditions

EXAM STRATEGY:
1. Set up grid (x, t)
2. Compute r = α·dt/dx²
3. Build A and B matrices (just diagonals!)
4. Handle boundaries separately
5. Time loop:
   - Compute RHS = B·u^n + boundary terms
   - Solve A·u^{n+1} = RHS using Thomas
6. Plot snapshots

COMMON MISTAKES:
- Forgetting to add boundary contributions to RHS
- Using wrong signs in matrices
- Not copying arrays in thomas_solve (modifies in place!)
- Confusing interior points vs full grid
"""
