# WORKSHEET 6 - LINE BY LINE EXPLANATION
# Van der Pol Oscillator with Adaptive Time Stepping

"""
THE PROBLEM:
d²y/dx² = μ(1-y²)dy/dx - λy
where μ = 5, λ = 40, from x=0 to x=20
Initial conditions: y(0) = 0.5, y'(0) = 0

This is a STIFF equation - has both fast and slow changes.
"""

# ==============================================================================
# PART 1: IMPORTS AND COEFFICIENTS
# ==============================================================================

import numpy as np
import matplotlib.pyplot as plt

# DORMAND-PRINCE COEFFICIENTS
# These are MAGIC NUMBERS from the Butcher tableau - DON'T CHANGE THEM!
# They define the 7-stage Runge-Kutta method

# a_i values: time fractions for each stage
a2 = 1/5    # Stage 2 evaluates at x + h/5
a3 = 3/10   # Stage 3 evaluates at x + 3h/10
a4 = 4/5    # Stage 4 evaluates at x + 4h/5
a5 = 8/9    # Stage 5 evaluates at x + 8h/9
a6 = 1      # Stage 6 evaluates at x + h
a7 = 1      # Stage 7 evaluates at x + h (for next step)

# b_ij values: how to combine previous k's to get input for next stage
b21 = 1/5
b31 = 3/40;       b32 = 9/40
b41 = 44/45;      b42 = -56/15;      b43 = 32/9
b51 = 19372/6561; b52 = -25360/2187; b53 = 64448/6561; b54 = -212/729
b61 = 9017/3168;  b62 = -355/33;     b63 = 46732/5247; b64 = 49/176; b65 = -5103/18656

# c_i values: weights for 5th order solution (more accurate)
c1 = 35/384; c2 = 0; c3 = 500/1113; c4 = 125/192; c5 = -2187/6784; c6 = 11/84; c7 = 0

# c_i* values: weights for 4th order solution (less accurate)
# The DIFFERENCE between y5 and y4 gives error estimate!
c1s = 5179/57600; c2s = 0; c3s = 7571/16695; c4s = 393/640
c5s = -92097/339200; c6s = 187/2100; c7s = 1/40


# ==============================================================================
# PART 2: REGULAR RK4 METHOD (for comparison)
# ==============================================================================

def rk4(fn, x, y, h):
    """
    LINE BY LINE:
    
    INPUT:
    - fn: the function f(x,y) that defines dy/dx = f(x,y)
    - x: current x value
    - y: current y value (can be array for coupled ODEs)
    - h: step size
    
    OUTPUT:
    - y_new: value at next step
    """
    # k1: slope at beginning of interval
    k1 = h * fn(x, y)
    
    # k2: slope at midpoint, using k1 to estimate y at midpoint
    k2 = h * fn(x + h/2, y + k1/2)
    
    # k3: slope at midpoint again, but using k2's estimate
    k3 = h * fn(x + h/2, y + k2/2)
    
    # k4: slope at end of interval, using k3's estimate
    k4 = h * fn(x + h, y + k3)
    
    # Weighted average: 1*k1 + 2*k2 + 2*k3 + 1*k4, divided by 6
    return y + (k1 + 2*k2 + 2*k3 + k4) / 6


# ==============================================================================
# PART 3: CALLER FOR REGULAR RK4
# ==============================================================================

def caller(my_method, fn, y_ini, x0, xT, h):
    """
    This runs ANY ODE method repeatedly from x0 to xT
    
    LINE BY LINE:
    """
    # Create array of x values from x0 to xT with step h
    xs = np.arange(x0, xT, h)
    
    # How many steps?
    N = len(xs)
    
    # Convert initial y to array (works for single or multiple equations)
    y = np.asarray(y_ini)
    
    # Storage for all y values
    # Shape: (N steps, number of variables)
    ys = np.zeros((N, len(y_ini)))
    
    # MAIN LOOP: iterate through all time steps
    for i in range(N):
        # Store current y values
        ys[i, :] = y
        
        # Take one step using the method
        # my_method could be euler, rk4, etc.
        y = my_method(fn, xs[i], y, h)
    
    return xs, ys


# ==============================================================================
# PART 4: DORMAND-PRINCE SINGLE STEP (THE HEART OF ADAPTIVE METHOD)
# ==============================================================================

def erk54h(f, x, y, h, k7bh):
    """
    Embedded Runge-Kutta 5(4) - one adaptive step
    
    WHY 7 STAGES?
    - Computes k1, k2, k3, k4, k5, k6, k7
    - Uses them to get BOTH 4th order AND 5th order estimates
    - Difference = error estimate
    
    FSAL PROPERTY (First Same As Last):
    - k7 from previous step = k1 for next step
    - Saves one function evaluation!
    
    INPUT:
    - f: function defining dy/dx = f(x,y)
    - x, y: current state
    - h: current step size
    - k7bh: k7 from PREVIOUS step (or zeros for first step)
    
    OUTPUT:
    - y5: 5th order solution (USE THIS - more accurate)
    - err: error estimate |y5 - y4|
    - k7bh: k7/h for NEXT step's FSAL
    """
    
    # STAGE 1: Use k7 from previous step if available
    if sum(k7bh) == 0:
        # First step - compute k1 normally
        k1 = h * f(x, y)
    else:
        # FSAL: reuse k7 from last step
        k1 = k7bh * h
    
    # STAGE 2: Evaluate at x + h/5
    # Input: y + (1/5)*k1
    k2 = h * f(x + a2*h, y + b21*k1)
    
    # STAGE 3: Evaluate at x + 3h/10
    # Input: y + (3/40)*k1 + (9/40)*k2
    k3 = h * f(x + a3*h, y + b31*k1 + b32*k2)
    
    # STAGE 4: Evaluate at x + 4h/5
    # Input: weighted combination of k1, k2, k3
    k4 = h * f(x + a4*h, y + b41*k1 + b42*k2 + b43*k3)
    
    # STAGE 5: Evaluate at x + 8h/9
    k5 = h * f(x + a5*h, y + b51*k1 + b52*k2 + b53*k3 + b54*k4)
    
    # STAGE 6: Evaluate at x + h
    k6 = h * f(x + a6*h, y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5)
    
    # COMPUTE 5TH ORDER SOLUTION
    # This uses weights c1, c2, c3, c4, c5, c6
    # Notice c2=0 and c7=0, so we skip k2 and k7
    y5 = y + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5 + c6*k6
    
    # STAGE 7: Evaluate at the NEW point for FSAL
    # Store k7/h (not k7) to avoid multiplying by h later
    k7bh = f(x + a7*h, y5)
    
    # COMPUTE 4TH ORDER SOLUTION
    # Uses different weights c1s, c2s, etc.
    # Now c7s ≠ 0, so we need k7
    y4 = y + c1s*k1 + c2s*k2 + c3s*k3 + c4s*k4 + c5s*k5 + c6s*k6 + c7s*k7bh*h
    
    # ERROR ESTIMATE
    # |y5 - y4| tells us how much 4th and 5th order solutions differ
    # Small difference = accurate step
    # Large difference = need smaller h
    err = abs(y5 - y4)
    
    return y5, err, k7bh


# ==============================================================================
# PART 5: ADAPTIVE CALLER (THE CONTROL LOOP)
# ==============================================================================

def caller54(fn, y_ini, x0, xT, h0, max_iter, abstol, reltol):
    """
    Adaptive time stepping controller
    
    THE ALGORITHM:
    1. Try a step with current h
    2. Check if error < tolerance
    3. If yes: accept, move forward, maybe increase h
    4. If no: reject, retry with smaller h
    
    INPUT:
    - fn: function f(x,y)
    - y_ini: initial condition
    - x0, xT: start and end x
    - h0: initial step size guess
    - max_iter: safety limit on number of steps
    - abstol: absolute tolerance (for when y ≈ 0)
    - reltol: relative tolerance (scales with |y|)
    
    OUTPUT:
    - xs, ys: arrays of accepted steps
    """
    
    # Initialize
    y = np.asarray(y_ini)
    h = h0  # Current step size
    
    # Storage (will trim later)
    ys = np.zeros((max_iter, len(y_ini)))
    xs = np.zeros((max_iter, 1))
    
    # Store initial condition
    xs[0] = x0
    ys[0, :] = y
    
    x = x0
    i = 0  # Step counter (counts ACCEPTED steps only)
    
    # FSAL: k7 from previous step
    k7bh = np.zeros(len(y_ini))
    
    # MAIN ADAPTIVE LOOP
    while x <= xT and i < max_iter:
        
        # COMPUTE TOLERANCE
        # tol = abstol + reltol * |y|
        # abstol handles y≈0, reltol handles large y
        tol = abstol + reltol * np.linalg.norm(ys[i, :])
        
        # TRY ONE STEP
        y_new, err, k7bh_new = erk54h(fn, xs[i], ys[i, :], h, k7bh)
        
        # Maximum error across all components
        merr = np.max(err)
        
        # Safety: avoid division by zero
        if merr == 0.0:
            merr = tol / 100
        
        # DECISION: ACCEPT OR REJECT?
        if merr <= tol:
            # ✅ ACCEPT THE STEP
            i += 1  # Move to next slot
            xs[i] = xs[i-1] + h  # Update x
            ys[i, :] = y_new  # Store new y
            x = xs[i]
            k7bh = k7bh_new  # Save k7 for FSAL
        
        # ADJUST STEP SIZE
        # Formula: h_new = 0.9 * h * (tol/err)^(1/5)
        # The 0.9 is safety factor
        # The 1/5 is for 5th order method
        h = 0.9 * h * (tol / merr)**0.2
        
        # SAFETY LIMITS on step size changes
        # Don't let h change too drastically
        if h > 2 * h0:
            h = 2 * h0  # Don't grow too fast
        if h < h0 / 100:
            h = h0 / 100  # Don't shrink too small
    
    # TRIM arrays to actual size used
    xs = xs[:i+1]
    ys = ys[:i+1, :]
    
    return xs, ys


# ==============================================================================
# PART 6: PROBLEM SETUP - VAN DER POL
# ==============================================================================

# Convert 2nd order ODE to system of 1st order:
# d²y/dx² = μ(1-y²)dy/dx - λy
#
# Let: y1 = y, y2 = dy/dx
# Then:
# dy1/dx = y2
# dy2/dx = μ(1-y1²)y2 - λy1

def vdp(x, y):
    """
    Van der Pol oscillator as system of ODEs
    
    INPUT:
    - x: independent variable (time)
    - y: array [y1, y2] where y1=y, y2=y'
    
    OUTPUT:
    - [dy1/dx, dy2/dx]
    """
    mu = 5
    lam = 40
    
    dy1 = y[1]  # dy/dx = y'
    dy2 = mu * (1 - y[0]**2) * y[1] - lam * y[0]  # d²y/dx²
    
    return np.array([dy1, dy2])


# Initial conditions
y_ini = [0.5, 0.0]  # y(0) = 0.5, y'(0) = 0
x0 = 0.0
xT = 20.0


# ==============================================================================
# PART 7: SOLVE THE PROBLEM
# ==============================================================================

# (a) Regular RK4 with small fixed h
print("Running RK4 with h = 1e-4...")
xs_rk4, ys_rk4 = caller(rk4, vdp, y_ini, x0, xT, 1e-4)
print(f"RK4 used {len(xs_rk4)} steps")

# Plot
plt.figure(figsize=(10, 5))
plt.plot(xs_rk4, ys_rk4[:, 0], 'b-', linewidth=0.8, label='y(x)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Van der Pol with RK4 (h=1e-4)')
plt.grid(alpha=0.3)
plt.legend()
plt.show()


# (b) Dormand-Prince adaptive
print("\nRunning Dormand-Prince adaptive...")
xs_dp, ys_dp = caller54(vdp, y_ini, x0, xT, h0=0.01, 
                        max_iter=50000, abstol=1e-6, reltol=1e-8)
print(f"Dormand-Prince used {len(xs_dp)} steps")
print(f"Speedup: {len(xs_rk4) / len(xs_dp):.1f}x fewer steps!")

# Plot
plt.figure(figsize=(10, 5))
plt.plot(xs_dp, ys_dp[:, 0], 'r-', linewidth=0.8, label='y(x)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Van der Pol with Dormand-Prince Adaptive')
plt.grid(alpha=0.3)
plt.legend()
plt.show()


# (c) Phase space plot
plt.figure(figsize=(7, 7))
plt.plot(ys_rk4[:, 0], ys_rk4[:, 1], 'b-', linewidth=0.5, 
         label='RK4', alpha=0.7)
plt.plot(ys_dp[:, 0], ys_dp[:, 1], 'r--', linewidth=0.8, 
         label='Dormand-Prince', alpha=0.8)
plt.xlabel('y')
plt.ylabel("y'")
plt.title('Phase Space: Van der Pol Oscillator')
plt.legend()
plt.grid(alpha=0.3)
plt.show()


# ==============================================================================
# KEY INSIGHTS FOR EXAM
# ==============================================================================

"""
WHAT STAYS THE SAME FOR ANY PROBLEM:
1. All the Dormand-Prince coefficients (a, b, c values)
2. The erk54h function structure
3. The caller54 adaptive logic
4. The tolerance formula

WHAT YOU CHANGE:
1. The function vdp(x, y) → YOUR differential equation
2. Initial conditions y_ini
3. x0, xT limits
4. Tolerances (abstol, reltol) if specified

HOW TO RECOGNIZE THIS TYPE OF PROBLEM:
- Keywords: "adaptive", "stiff", "varying timescales"
- Equation has both fast and slow dynamics
- Problem asks to compare with fixed-step method

TYPICAL EXAM QUESTION:
"Solve the equation [...] using adaptive RK45 with abstol=1e-6, reltol=1e-8.
Compare with RK4 using h=1e-4. How many steps does each use?"

YOUR ANSWER TEMPLATE:
1. Convert equation to 1st order system
2. Define function f(x, y)
3. Call caller54() with specified tolerances
4. Call caller() with RK4
5. Count steps, compare
6. Plot both solutions
"""
