# WORKSHEET 9 - LINE BY LINE EXPLANATION
# TIME-DEPENDENT SCHRÖDINGER EQUATION (TDSE)
# STRANG SPLITTING + FFT METHOD

"""
THE PROBLEM:
Solve: iℏ ∂ψ/∂t = [-ℏ²/2m ∂²/∂x² + V(x)]ψ

Two parts:
Q1: FFT of (a) square function, (b) double slit
Q2: TDSE for (a) barrier scattering, (b) harmonic oscillator

KEY INSIGHT:
Hamiltonian = Kinetic + Potential
- Kinetic is diagonal in MOMENTUM space
- Potential is diagonal in POSITION space
- FFT switches between them efficiently!
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.display import HTML

# ==============================================================================
# PART 1: QUESTION 1(a) - FFT OF SQUARE FUNCTION
# ==============================================================================

"""
FOURIER TRANSFORM:
F(k) = ∫ f(x) e^(-ikx) dx

For square function (rect function):
- In real space: sharp edges
- In k-space: sinc function (oscillations)

Why sinc? Sharp edges → high frequencies needed
"""

# SETUP GRID
N = 1024  # Number of points (power of 2 for efficient FFT)
L = 20.0  # Domain length
x = np.linspace(-L/2, L/2, N, endpoint=False)
dx = L / N

# DEFINE SQUARE FUNCTION
width = 2.0
square = np.zeros(N)
mask = np.abs(x) <= width/2  # True where |x| ≤ 1
square[mask] = 1.0

"""
LINE BY LINE - Square function:
- square[mask] = 1.0: Set to 1 inside [-1, 1]
- square is zero elsewhere
- This is the "rect" or "top-hat" function
"""

# COMPUTE FFT
fft_square = np.fft.fft(square)

"""
np.fft.fft() computes:
F[k] = Σ_n f[n] e^(-2πikn/N)

This is DISCRETE Fourier transform
- Fast: O(N log N) instead of O(N²)
- Returns complex numbers
- Output is shifted (need fftshift)
"""

# SHIFT FOR PLOTTING
fft_square_shifted = np.fft.fftshift(fft_square)
freqs = np.fft.fftfreq(N, d=dx)
freqs_shifted = np.fft.fftshift(freqs)

"""
WHY FFTSHIFT?
- FFT outputs frequencies in order: [0, 1, 2, ..., N/2, -N/2, ..., -2, -1]
- fftshift rearranges to: [-N/2, ..., -1, 0, 1, ..., N/2]
- This puts zero frequency at center (better for plotting)

np.fft.fftfreq(N, d=dx):
- Returns frequency array
- Frequencies: k_n = 2πn/L for n = 0, 1, ..., N/2, -N/2, ..., -1
"""

# PLOT
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))

ax1.plot(x, square, 'b-', linewidth=2)
ax1.set_xlabel('x')
ax1.set_ylabel('f(x)')
ax1.set_title('Square Function')
ax1.grid(alpha=0.3)

ax2.plot(freqs_shifted, np.abs(fft_square_shifted), 'r-', linewidth=1)
ax2.set_xlabel('Frequency k')
ax2.set_ylabel('|F(k)|')
ax2.set_title('Magnitude of FFT (Sinc Pattern)')
ax2.set_xlim(-5, 5)
ax2.grid(alpha=0.3)

ax3.plot(freqs_shifted, np.real(fft_square_shifted), 'g-', linewidth=1, label='Real')
ax3.plot(freqs_shifted, np.imag(fft_square_shifted), 'm-', linewidth=1, label='Imag')
ax3.set_xlabel('Frequency k')
ax3.set_ylabel('F(k)')
ax3.set_title('Real and Imaginary Parts')
ax3.set_xlim(-5, 5)
ax3.legend()
ax3.grid(alpha=0.3)

plt.tight_layout()
plt.show()

"""
WHAT TO OBSERVE:
- Magnitude shows sinc pattern: sin(x)/x shape
- First zero at k ≈ π/width (wider function → narrower FFT)
- Real part: even function (symmetric)
- Imaginary part: should be ≈0 (real input → conjugate symmetric FFT)
"""


# ==============================================================================
# PART 2: QUESTION 1(b) - FFT OF DOUBLE SLIT
# ==============================================================================

N = 1024
L = 20.0
x = np.linspace(-L/2, L/2, N, endpoint=False)
dx = L / N

# DEFINE DOUBLE SLIT (two separate rectangles)
slit1_center = -2.0
slit1_width = 0.5
slit2_center = 2.0
slit2_width = 1.0

double_slit = np.zeros(N)

# First slit
mask1 = np.abs(x - slit1_center) <= slit1_width/2
double_slit[mask1] = 1.0

# Second slit
mask2 = np.abs(x - slit2_center) <= slit2_width/2
double_slit[mask2] = 1.0

"""
DOUBLE SLIT:
- Two rectangles at different positions
- Unequal widths (0.5 vs 1.0)
- Separation: 4.0 units

FFT will show:
- Interference fringes (from separation)
- Modulated by sinc envelope (from widths)
"""

# COMPUTE FFT
fft_double = np.fft.fft(double_slit)
fft_double_shifted = np.fft.fftshift(fft_double)
freqs = np.fft.fftfreq(N, d=dx)
freqs_shifted = np.fft.fftshift(freqs)

# PLOT (similar to square function)
# ... [plotting code omitted for brevity]

"""
WHAT TO OBSERVE:
- Rapid oscillations (interference from two slits)
- Frequency of oscillations ∝ slit separation
- Envelope shape from individual slit widths
- This is like optical diffraction pattern!
"""


# ==============================================================================
# PART 3: STRANG SPLITTING - THE CORE ALGORITHM
# ==============================================================================

def split_operator_step(psi, V, dt, k_vals, hbar=1.0, m=1.0):
    """
    ONE TIME STEP OF TDSE USING STRANG SPLITTING
    
    THE EQUATION: iℏ ∂ψ/∂t = Hψ where H = T + V
    
    EXACT SOLUTION: ψ(t+Δt) = e^(-iHΔt/ℏ) ψ(t)
    
    PROBLEM: T and V don't commute → can't split exponential
    
    STRANG SPLITTING (2nd order accurate):
    e^(-i(T+V)Δt/ℏ) ≈ e^(-iVΔt/2ℏ) e^(-iTΔt/ℏ) e^(-iVΔt/2ℏ)
    
    ALGORITHM:
    1. Half-step potential (real space)
    2. FFT to momentum space
    3. Full kinetic step (momentum space)
    4. IFFT back to real space
    5. Half-step potential (real space)
    
    INPUT:
    - psi: wavefunction (complex array)
    - V: potential (real array, same length as psi)
    - dt: time step
    - k_vals: momentum grid (from fftfreq)
    - hbar, m: physical constants (set to 1 typically)
    
    OUTPUT:
    - psi: evolved wavefunction
    
    LINE BY LINE:
    """
    
    # KINETIC ENERGY IN MOMENTUM SPACE
    # T = p²/2m = ℏ²k²/2m
    E_kin = (hbar**2 * k_vals**2) / (2.0 * m)
    
    # STEP 1: HALF POTENTIAL (real space)
    # ψ → e^(-iVΔt/2ℏ) ψ
    # This is pointwise multiplication!
    psi_k = np.fft.fft(psi)
    psi_k *= np.exp(-0.5j * E_kin * dt / hbar)
    psi = np.fft.ifft(psi_k)
    
    """
    WHY THIS ORDER?
    - We want: potential, kinetic, potential
    - But we START with potential!
    - This uses FSAL (First Same As Last) trick
    - Actually: half-pot from previous step + this step = full pot
    """
    
    # STEP 2: FULL POTENTIAL (real space)
    # ψ → e^(-iVΔt/ℏ) ψ
    psi *= np.exp(-1j * V * dt / hbar)
    
    # STEP 3: KINETIC (momentum space)
    # FFT to momentum space
    psi_k = np.fft.fft(psi)
    
    # Apply kinetic propagator
    # ψ_k → e^(-iE_kin·Δt/ℏ) ψ_k
    psi_k *= np.exp(-0.5j * E_kin * dt / hbar)
    
    # IFFT back to real space
    psi = np.fft.ifft(psi_k)
    
    return psi


"""
CRITICAL UNDERSTANDING:

STRANG SPLITTING FORMULA:
ψ^{n+1} = e^{-iVΔt/2} · F^{-1}[e^{-ik²Δt/2} · F[e^{-iVΔt/2} · ψ^n]]

Read from right to left:
1. e^{-iVΔt/2} · ψ^n          : Half potential
2. F[...]                      : FFT to k-space
3. e^{-ik²Δt/2} · (...)        : Full kinetic
4. F^{-1}[...]                 : IFFT to x-space
5. e^{-iVΔt/2} · (...)         : Half potential

WHY IT WORKS:
- Potential diagonal in x: just multiply
- Kinetic diagonal in k: just multiply
- FFT connects x and k: O(N log N)
- No matrix inversion needed!
- Unitary: preserves norm |ψ|²
"""


def initial_gaussian(x, x0, sigma, k0=0.0):
    """
    Create Gaussian wavepacket
    
    FORMULA: ψ(x,0) = (1/σ√(2π))^(1/2) exp[-(x-x₀)²/2σ²] exp(ik₀x)
    
    PARTS:
    - Gaussian envelope: centered at x₀, width σ
    - Momentum: exp(ik₀x) gives mean momentum ℏk₀
    - Normalization: (1/σ√(2π))^(1/2)
    
    CRITICAL: k₀ ≠ 0 FOR SCATTERING!
    - k₀ = 0: stationary packet (just spreads)
    - k₀ ≠ 0: moving packet (can scatter)
    
    LINE BY LINE:
    """
    # Normalization constant
    norm = 1.0 / (sigma * np.sqrt(2 * np.pi))
    
    # Gaussian · phase
    return np.sqrt(norm) * np.exp(-(x - x0)**2 / (2 * sigma**2)) * np.exp(1j * k0 * x)
    
    """
    exp(ik₀x) means:
    - Real part: cos(k₀x) - oscillates in space
    - Imaginary part: sin(k₀x) - oscillates in space
    - Together: traveling wave with wavelength λ = 2π/k₀
    """


# ==============================================================================
# PART 4: QUESTION 2(a) - BARRIER SCATTERING
# ==============================================================================

print("=" * 60)
print("QUESTION 2(a): BARRIER SCATTERING")
print("=" * 60)

# PARAMETERS (ℏ = m = 1 units)
hbar = 1.0
m = 1.0

# SPATIAL GRID
Nx = 512  # Number of points
x_min = -20.0
x_max = 20.0
L = x_max - x_min
dx = L / Nx

x = np.linspace(x_min, x_max, Nx, endpoint=False)

"""
endpoint=False:
- Excludes x_max from grid
- FFT assumes periodic domain: x(N) = x(0) + L
- So we want [x_min, x_min+dx, ..., x_max-dx]
"""

# MOMENTUM GRID
k_vals = np.fft.fftfreq(Nx, d=dx) * 2.0 * np.pi

"""
fftfreq returns: k_n = n/L for n = 0, 1, ..., N/2-1, -N/2, ..., -1
Multiply by 2π to get angular frequency
k_vals[0] = 0 (zero frequency)
k_vals[N/2] = π/dx (Nyquist frequency - highest resolvable)
"""

# POTENTIAL BARRIER
V0 = 40.0
barrier_left = 5.0
barrier_right = 7.0

V_barrier = np.zeros(Nx)
mask = (x >= barrier_left) & (x <= barrier_right)
V_barrier[mask] = V0

"""
RECTANGULAR BARRIER:
- V = 0 for x < 5 or x > 7
- V = 40 for 5 ≤ x ≤ 7
- Width: 2 units
- Height: 40 (energy units)

PHYSICS:
- Incoming wavepacket will partially transmit, partially reflect
- Quantum tunneling through classically forbidden region
"""

# INITIAL WAVEPACKET
x0 = -5.0   # Start left of barrier
sigma = 0.04  # Narrow packet
k0 = 5.0    # CRITICAL: gives momentum to the right!

psi_barrier = initial_gaussian(x, x0, sigma, k0)

"""
k₀ = 5.0 means:
- Wavelength λ = 2π/5 ≈ 1.26
- Momentum p = ℏk₀ = 5
- Kinetic energy E = p²/2m = 12.5
- E < V₀ = 40 → classically forbidden!
- But quantum: can tunnel through
"""

print(f"Wavepacket:")
print(f"  Position: x₀ = {x0}")
print(f"  Width: σ = {sigma}")
print(f"  Momentum: k₀ = {k0}")
print(f"  Energy: E ≈ {k0**2 / 2:.2f}")
print(f"Barrier: V₀ = {V0}, from x = {barrier_left} to {barrier_right}")

# TIME EVOLUTION
dt = 0.01  # Time step
n_steps = 2000  # Total steps
plot_every = 10  # Save every 10th step

# STORAGE
snapshots_barrier = [psi_barrier.copy()]
times_barrier = [0]

print(f"\nTime evolution:")
print(f"  dt = {dt}")
print(f"  Total time: {n_steps * dt}")
print(f"  Saving {n_steps // plot_every} snapshots")

# MAIN EVOLUTION LOOP
for step in range(1, n_steps + 1):
    # ONE STRANG STEP
    psi_barrier = split_operator_step(psi_barrier, V_barrier, dt, k_vals, hbar, m)
    
    # SAVE SNAPSHOT
    if step % plot_every == 0:
        snapshots_barrier.append(psi_barrier.copy())
        times_barrier.append(step * dt)

"""
WHAT'S HAPPENING:
t=0: Wavepacket at x=-5, moving right
t→5: Approaching barrier
t≈7: Hitting barrier - splits!
t>10: Transmitted part continues right
      Reflected part moves left
t→20: Separated completely
"""

# ANALYSIS: TRANSMISSION AND REFLECTION
transmitted = np.sum(np.abs(snapshots_barrier[-1][x > barrier_right])**2) * dx
reflected = np.sum(np.abs(snapshots_barrier[-1][x < barrier_left])**2) * dx

print(f"\nResults:")
print(f"  Transmitted: {transmitted:.4f}")
print(f"  Reflected: {reflected:.4f}")
print(f"  Total: {transmitted + reflected:.4f} (should be ≈1)")

"""
PROBABILITY INTERPRETATION:
- |ψ|² = probability density
- ∫|ψ|² dx = total probability = 1
- Transmitted = probability in region x > 7
- Reflected = probability in region x < 5
"""

# CREATE ANIMATION
fig, ax = plt.subplots(figsize=(12, 6))
line_barrier, = ax.plot(x, np.abs(snapshots_barrier[0])**2, 'b-', linewidth=2)

# Show barrier
ax.fill_between(x, 0, V_barrier/V0 * 0.15, alpha=0.3, color='red', label='Barrier')

ax.set_xlim(x_min, x_max)
ax.set_ylim(0, 0.2)
ax.set_xlabel('x')
ax.set_ylabel('|ψ|²')
ax.legend()
ax.grid(alpha=0.3)

def update_barrier(frame):
    """Update function for animation"""
    prob = np.abs(snapshots_barrier[frame])**2
    line_barrier.set_ydata(prob)
    ax.set_title(f'Barrier Scattering - t = {times_barrier[frame]:.2f}')
    return line_barrier,

# Create animation object
ani_barrier = FuncAnimation(fig, update_barrier, frames=len(snapshots_barrier), 
                           blit=True, interval=50)

# Convert to HTML (for Jupyter)
plt.close()
HTML(ani_barrier.to_jshtml())

"""
ANIMATION SHOWS:
- Blue line: |ψ(x,t)|² - probability density
- Red region: potential barrier
- Watch wavepacket split at barrier!
- Transmitted and reflected parts separate
"""


# ==============================================================================
# PART 5: QUESTION 2(b) - HARMONIC OSCILLATOR
# ==============================================================================

# Similar structure, different potential
V_harmonic = 0.1 * x**2

x0 = -5.0
sigma = 0.04
k0 = 3.0  # Less momentum than barrier case

psi_harmonic = initial_gaussian(x, x0, sigma, k0)

# Time evolution (longer to see oscillations)
dt = 0.01
n_steps = 3000
plot_every = 15

snapshots_harmonic = [psi_harmonic.copy()]
times_harmonic = [0]
expectation_x = [np.sum(x * np.abs(psi_harmonic)**2) * dx]

for step in range(1, n_steps + 1):
    psi_harmonic = split_operator_step(psi_harmonic, V_harmonic, dt, k_vals, hbar, m)
    
    if step % plot_every == 0:
        snapshots_harmonic.append(psi_harmonic.copy())
        times_harmonic.append(step * dt)
        expectation_x.append(np.sum(x * np.abs(psi_harmonic)**2) * dx)

"""
HARMONIC OSCILLATOR:
- Potential: V(x) = 0.1x²
- Parabolic well
- Classical: mass oscillates back and forth
- Quantum: wavepacket oscillates AND spreads

EXPECTATION VALUE:
⟨x⟩ = ∫ x|ψ|² dx
Should oscillate like classical particle!
"""

# Plot ⟨x⟩ vs time
plt.figure(figsize=(10, 5))
plt.plot(times_harmonic, expectation_x, 'b-', linewidth=2)
plt.xlabel('Time')
plt.ylabel('⟨x⟩')
plt.title('Expectation Value of Position vs Time')
plt.grid(alpha=0.3)
plt.show()


# ==============================================================================
# KEY INSIGHTS FOR EXAM
# ==============================================================================

"""
WHAT STAYS THE SAME:
1. split_operator_step() function - NEVER CHANGES
2. FFT/IFFT structure
3. Gaussian initialization
4. Grid setup with endpoint=False
5. k_vals from fftfreq

WHAT YOU CHANGE:
1. Potential V(x) - YOUR problem
2. Initial position x₀
3. Momentum k₀ - MUST BE ≠0 for scattering!
4. Domain size (must be big enough)
5. dt, n_steps

CRITICAL POINTS:

1. MOMENTUM IS ESSENTIAL:
   ψ = gaussian * exp(ik₀x)
   Without k₀: packet just spreads
   With k₀: packet moves and scatters

2. GRID MUST BE LARGE:
   - Wavepacket shouldn't reach boundaries
   - Otherwise: wraps around (periodic!)
   - Rule: domain ≫ distance packet travels

3. NORM CONSERVATION:
   ∫|ψ|² dx = 1 always
   Check: transmitted + reflected ≈ 1

4. FFT CONVENTIONS:
   - fftfreq gives k in cycles
   - Multiply by 2π for angular frequency
   - fftshift for plotting

EXAM STRATEGY:
1. Set up grid with endpoint=False
2. Get k_vals from fftfreq
3. Define potential
4. Create Gaussian with k₀ ≠ 0
5. Loop: call split_operator_step
6. Save snapshots
7. Analyze transmission/reflection
8. Make animation

COMMON MISTAKES:
- Forgetting k₀ (no scattering!)
- Domain too small (wraparound)
- dt too large (numerical instability)
- Using endpoint=True (breaks periodicity)
- Not checking norm conservation
"""
