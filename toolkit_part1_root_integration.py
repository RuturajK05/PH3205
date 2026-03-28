"""
PH3205 COMPUTATIONAL PHYSICS - MASTER TOOLKIT PART 1
Root Finding & Numerical Integration
For Class Test - March 28, 2026

HOW TO USE THIS FILE:
1. Copy the function you need
2. Define YOUR specific f(x) or problem
3. Call the function with appropriate parameters
4. Plot results

TOPICS COVERED:
- Root Finding (Bisection, Secant, Newton-Raphson)
- Numerical Integration (Trapezoidal, Simpson 1/3, Simpson 3/8, Boole)
- Error Analysis & Convergence Testing
"""

import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# SECTION 1: ROOT FINDING METHODS
# ==============================================================================

# ------------------------------------------------------------------------------
# 1.1 BISECTION METHOD
# ------------------------------------------------------------------------------
def bisection(f, L, R, eps=1e-5, maxiter=50):
    """
    WHAT IT DOES: Finds root of f(x) = 0 in interval [L, R]
    
    WHEN TO USE:
    - Most reliable method
    - When you have a bracketing interval where f(L) and f(R) have different signs
    
    WHAT STAYS SAME:
    - The algorithm structure
    - eps and maxiter (unless problem specifies)
    
    WHAT YOU CHANGE:
    - f: YOUR function
    - L, R: YOUR interval endpoints
    
    EXAMPLE USAGE:
        def my_func(x):
            return 3*x - np.tan(x)
        root = bisection(my_func, 1.2, 1.4)
    
    RETURNS: (iterations_used, root_value)
    """
    # Safety check: proper bracketing
    if f(L) * f(R) > 0:
        print("Error: f(L) and f(R) must have different signs")
        return None, None
    
    # Check if endpoints are roots
    if abs(f(L)) < eps:
        return 0, L
    if abs(f(R)) < eps:
        return 0, R
    
    # Main bisection loop
    for iteration in range(maxiter):
        C = (L + R) / 2  # Midpoint
        
        if abs(f(C)) < eps:  # Found root!
            return iteration + 1, C
        
        # Decide which half contains the root
        if f(L) * f(C) < 0:
            R = C  # Root in left half
        else:
            L = C  # Root in right half
    
    # Max iterations reached
    return maxiter, (L + R) / 2


# ------------------------------------------------------------------------------
# 1.2 SECANT METHOD
# ------------------------------------------------------------------------------
def secant(f, x1, x2, eps=1e-5, maxiter=50):
    """
    WHAT IT DOES: Finds root using two points, approximating derivative
    
    WHEN TO USE:
    - When you DON'T have derivative
    - When you want faster convergence than bisection
    
    WHAT STAYS SAME:
    - The algorithm
    - eps, maxiter
    
    WHAT YOU CHANGE:
    - f: YOUR function
    - x1, x2: YOUR two starting guesses (close to root)
    
    EXAMPLE USAGE:
        def my_func(x):
            return x**4 + 2*x**3/3 - x**2 + 0.125
        root = secant(my_func, 1.2, 1.25)
    
    RETURNS: (iterations_used, root_value)
    """
    # Check if starting points are roots
    if abs(f(x1)) < eps:
        return 0, x1
    if abs(f(x2)) < eps:
        return 0, x2
    
    for iteration in range(maxiter):
        # Compute next point using secant formula
        f1 = f(x1)
        f2 = f(x2)
        
        if abs(f2 - f1) < 1e-12:  # Avoid division by zero
            print("Warning: f(x1) ≈ f(x2), method may fail")
            break
        
        x3 = x2 - f2 * (x2 - x1) / (f2 - f1)
        
        if abs(f(x3)) < eps:  # Found root!
            return iteration + 1, x3
        
        # Update points
        x1, x2 = x2, x3
    
    return maxiter, x2


# ------------------------------------------------------------------------------
# 1.3 NEWTON-RAPHSON METHOD
# ------------------------------------------------------------------------------
def newton_raphson(f, df, x0, eps=1e-5, maxiter=50):
    """
    WHAT IT DOES: Finds root using function AND its derivative
    
    WHEN TO USE:
    - When you HAVE the derivative
    - When you want fastest convergence
    
    WHAT STAYS SAME:
    - The algorithm
    - eps, maxiter
    
    WHAT YOU CHANGE:
    - f: YOUR function
    - df: YOUR derivative function
    - x0: YOUR starting guess
    
    EXAMPLE USAGE:
        def my_func(x):
            return 3*x - np.tan(x)
        def my_deriv(x):
            return 3 - 1/np.cos(x)**2
        root = newton_raphson(my_func, my_deriv, 1.2)
    
    RETURNS: (iterations_used, root_value)
    """
    x = x0
    
    if abs(f(x)) < eps:  # Starting point is root
        return 0, x
    
    for iteration in range(maxiter):
        fx = f(x)
        dfx = df(x)
        
        if abs(dfx) < 1e-12:  # Avoid division by zero
            print("Warning: derivative ≈ 0, method may fail")
            break
        
        x_new = x - fx / dfx
        
        if abs(f(x_new)) < eps:  # Found root!
            return iteration + 1, x_new
        
        x = x_new
    
    return maxiter, x


# ==============================================================================
# SECTION 2: NUMERICAL INTEGRATION
# ==============================================================================

# ------------------------------------------------------------------------------
# 2.1 TRAPEZOIDAL RULE
# ------------------------------------------------------------------------------
def trapezoidal(f, a, b, N):
    """
    WHAT IT DOES: Integrates f(x) from a to b using N+1 points
    
    FORMULA: h/2 * [f0 + 2*f1 + 2*f2 + ... + 2*f_{N-1} + fN]
    
    WHEN TO USE:
    - Simple problems
    - When you need okay accuracy quickly
    
    ERROR: O(h²) where h = (b-a)/N
    
    WHAT STAYS SAME:
    - The algorithm
    
    WHAT YOU CHANGE:
    - f: YOUR function to integrate
    - a, b: YOUR limits
    - N: Number of intervals (try N = 2*m where m = 4, 8, 16, 32, 64)
    
    EXAMPLE USAGE:
        def my_func(x):
            return np.exp(x)
        result = trapezoidal(my_func, 0, 1, 100)
    
    RETURNS: integral_value
    """
    h = (b - a) / N
    x = np.linspace(a, b, N + 1)
    y = f(x)
    
    # Trapezoidal: sum = f0 + 2*f1 + 2*f2 + ... + 2*f_{N-1} + fN
    integral = (y[0] + y[-1]) / 2 + np.sum(y[1:-1])
    integral *= h
    
    return integral


# ------------------------------------------------------------------------------
# 2.2 SIMPSON'S 1/3 RULE
# ------------------------------------------------------------------------------
def simpson_13(f, a, b, N):
    """
    WHAT IT DOES: Integrates using parabolic interpolation
    
    FORMULA: h/3 * [f0 + 4*f1 + 2*f2 + 4*f3 + ... + fN]
    PATTERN: 1, 4, 2, 4, 2, 4, ..., 2, 4, 1
    
    WHEN TO USE:
    - DEFAULT CHOICE for integration
    - Better accuracy than trapezoidal
    
    IMPORTANT: N must be EVEN!
    
    ERROR: O(h⁴)
    
    WHAT YOU CHANGE:
    - f: YOUR function
    - a, b: YOUR limits
    - N: EVEN number (N = 2*m where m = 4, 8, 16, ...)
    
    EXAMPLE USAGE:
        result = simpson_13(lambda x: np.exp(x), 0, 1, 100)
    
    RETURNS: integral_value
    """
    if N % 2 != 0:
        print("Warning: N must be even for Simpson 1/3")
        N = N + 1
    
    h = (b - a) / N
    x = np.linspace(a, b, N + 1)
    y = f(x)
    
    # Simpson 1/3: pattern is 1, 4, 2, 4, 2, ..., 4, 1
    integral = y[0] + y[-1]
    integral += 4 * np.sum(y[1:-1:2])  # Odd indices: 1, 3, 5, ...
    integral += 2 * np.sum(y[2:-1:2])  # Even indices: 2, 4, 6, ...
    integral *= h / 3
    
    return integral


# ------------------------------------------------------------------------------
# 2.3 SIMPSON'S 3/8 RULE
# ------------------------------------------------------------------------------
def simpson_38(f, a, b, N):
    """
    WHAT IT DOES: Uses cubic interpolation
    
    FORMULA: 3h/8 * [f0 + 3*f1 + 3*f2 + 2*f3 + 3*f4 + ...]
    PATTERN: 1, 3, 3, 2, 3, 3, 2, ..., 3, 3, 1
    
    IMPORTANT: N must be multiple of 3!
    
    ERROR: O(h⁴) (same as Simpson 1/3, but slightly different)
    
    WHAT YOU CHANGE:
    - N: Must be 3*m where m = 4, 8, 16, ...
    
    RETURNS: integral_value
    """
    if N % 3 != 0:
        print("Warning: N must be multiple of 3")
        N = N + (3 - N % 3)
    
    h = (b - a) / N
    x = np.linspace(a, b, N + 1)
    y = f(x)
    
    # Simpson 3/8: pattern 1, 3, 3, 2, 3, 3, 2, ..., 1
    integral = y[0] + y[-1]
    
    for i in range(1, N):
        if i % 3 == 0:
            integral += 2 * y[i]
        else:
            integral += 3 * y[i]
    
    integral *= 3 * h / 8
    
    return integral


# ------------------------------------------------------------------------------
# 2.4 BOOLE'S RULE
# ------------------------------------------------------------------------------
def boole(f, a, b, N):
    """
    WHAT IT DOES: Uses 4th order polynomial interpolation
    
    IMPORTANT: N must be multiple of 4!
    
    ERROR: O(h⁶) - most accurate!
    
    WHAT YOU CHANGE:
    - N: Must be 4*m where m = 4, 8, 16, ...
    
    RETURNS: integral_value
    """
    if N % 4 != 0:
        print("Warning: N must be multiple of 4")
        N = N + (4 - N % 4)
    
    h = (b - a) / N
    x = np.linspace(a, b, N + 1)
    y = f(x)
    
    # Boole: pattern 7, 32, 12, 32, 7, 32, 12, ...
    integral = 0
    for i in range(0, N, 4):
        integral += (7*y[i] + 32*y[i+1] + 12*y[i+2] + 32*y[i+3] + 7*y[i+4])
    
    integral *= 2 * h / 45
    
    return integral


# ==============================================================================
# SECTION 3: CONVERGENCE ANALYSIS TOOLS
# ==============================================================================

def convergence_test_integration(f, a, b, exact, method_name='simpson_13'):
    """
    WHAT IT DOES: Tests how error decreases with smaller h
    
    WHEN TO USE:
    - Worksheet questions asking to "compare methods"
    - When you need to verify order of accuracy
    
    WHAT YOU CHANGE:
    - f, a, b, exact: YOUR problem
    - method_name: 'trapezoidal', 'simpson_13', 'simpson_38', 'boole'
    
    EXAMPLE USAGE:
        convergence_test_integration(np.exp, 0, 1, np.e - 1, 'simpson_13')
    
    WHAT IT PLOTS:
    - log(h) vs log(error)
    - Slope of line = order of method
    """
    methods = {
        'trapezoidal': trapezoidal,
        'simpson_13': simpson_13,
        'simpson_38': simpson_38,
        'boole': boole
    }
    
    method = methods[method_name]
    
    m_values = [4, 8, 16, 32, 64]
    h_values = []
    errors = []
    
    for m in m_values:
        if method_name == 'trapezoidal' or method_name == 'simpson_13':
            N = 2 * m
        elif method_name == 'simpson_38':
            N = 3 * m
        else:  # boole
            N = 4 * m
        
        h = (b - a) / N
        result = method(f, a, b, N)
        error = abs(result - exact)
        
        h_values.append(h)
        errors.append(error)
    
    # Plot log-log
    plt.figure(figsize=(8, 6))
    plt.loglog(h_values, errors, 'o-', label=method_name)
    
    # Fit line to get slope
    coeffs = np.polyfit(np.log10(h_values), np.log10(errors), 1)
    slope = coeffs[0]
    
    plt.title(f'{method_name}: Order ≈ {slope:.2f}')
    plt.xlabel('log(h)')
    plt.ylabel('log(error)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.show()
    
    print(f"Method: {method_name}")
    print(f"Estimated order: {slope:.2f}")
    
    return h_values, errors


# ==============================================================================
# SECTION 4: QUICK PROBLEM TEMPLATES
# ==============================================================================

def solve_root_finding_problem():
    """
    TEMPLATE: Copy this to solve root-finding questions
    """
    print("=" * 60)
    print("ROOT FINDING TEMPLATE")
    print("=" * 60)
    
    # STEP 1: Define your function
    def f(x):
        return 3*x - np.tan(x)  # CHANGE THIS
    
    # STEP 2: Plot to see where roots are
    x = np.linspace(0, 1.57, 200)
    plt.plot(x, f(x))
    plt.axhline(0, color='k', linestyle='--', alpha=0.3)
    plt.grid(True, alpha=0.3)
    plt.title('Plot f(x) to find root location')
    plt.show()
    
    # STEP 3: Use bisection
    iter_b, root_b = bisection(f, 1.2, 1.4)
    print(f"\nBisection: {iter_b} iterations, root = {root_b:.6f}")
    
    # STEP 4: Use secant
    iter_s, root_s = secant(f, 1.2, 1.25)
    print(f"Secant: {iter_s} iterations, root = {root_s:.6f}")
    
    # STEP 5: Use Newton (if you have derivative)
    def df(x):
        return 3 - 1/np.cos(x)**2  # CHANGE THIS
    
    iter_n, root_n = newton_raphson(f, df, 1.2)
    print(f"Newton-Raphson: {iter_n} iterations, root = {root_n:.6f}")


def solve_integration_problem():
    """
    TEMPLATE: Copy this to solve integration questions
    """
    print("=" * 60)
    print("INTEGRATION TEMPLATE")
    print("=" * 60)
    
    # STEP 1: Define function and exact answer
    def f(x):
        return np.exp(x)  # CHANGE THIS
    
    a, b = 0, 1  # CHANGE LIMITS
    exact = np.e - 1  # CHANGE EXACT ANSWER
    
    # STEP 2: Compute with different methods
    N = 64  # or 2*32, 3*32, 4*32 depending on method
    
    trap = trapezoidal(f, a, b, N)
    simp13 = simpson_13(f, a, b, N)
    
    print(f"Trapezoidal: {trap:.8f}, Error: {abs(trap - exact):.2e}")
    print(f"Simpson 1/3: {simp13:.8f}, Error: {abs(simp13 - exact):.2e}")
    print(f"Exact: {exact:.8f}")
    
    # STEP 3: Convergence test
    convergence_test_integration(f, a, b, exact, 'simpson_13')


# ==============================================================================
# TESTING SECTION - Uncomment to test
# ==============================================================================
if __name__ == "__main__":
    print("PH3205 Toolkit Part 1 - Loaded Successfully!")
    print("\nAvailable Functions:")
    print("  Root Finding: bisection(), secant(), newton_raphson()")
    print("  Integration: trapezoidal(), simpson_13(), simpson_38(), boole()")
    print("  Analysis: convergence_test_integration()")
    print("\nTemplates: solve_root_finding_problem(), solve_integration_problem()")
    
    # Uncomment to run examples:
    # solve_root_finding_problem()
    # solve_integration_problem()
