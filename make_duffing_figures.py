"""Generate the Duffing-oscillator figures used in the README.

Uses the functions in Forced_Oscillator.py with:
    m = 1, k = -1 (double-well), b = 0.3 (damping), alpha = 1 (x^3 coefficient),
    drive amplitude F0 = 0.5, drive frequency omega = 1.2, x(0) = 1, v(0) = 0.

The run length is chosen so that one drive period is exactly 500 time steps,
which keeps the once-per-period sampling for the Poincare section exact.

Saves: poincare.png, bifurcation.png, chaos_map.png (takes about a minute).
"""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import Forced_Oscillator as fo

k, m, b, alpha = -1.0, 1.0, 0.3, 1.0
F0, omega = 0.5, 1.2
x0, v0 = 1.0, 0.0

T = 2 * np.pi / omega            # drive period
steps_per_period = 500

# 1) Poincare section: long run (2000 drive periods) so there are plenty of points
N_long = 2000 * steps_per_period + 1
t_max_long = 2000 * T
dt = t_max_long / (N_long - 1)
t, x, v = fo.duffing(x0, v0, k, m, t_max_long, N_long, b, F0, omega, alpha, dt)
fo.plot_duffing(t, x, v, N_long, omega, t_max_long, dt)

# 2) Bifurcation diagram and 3) Lyapunov exponent map: 100 drive periods per run
N = 100 * steps_per_period + 1
t_max = 100 * T
dt = t_max / (N - 1)
fo.bifurcation(x0, v0, k, m, t_max, N, b, F0, omega, alpha, dt)
fo.chaos_map(x0, v0, k, m, b, F0, omega, alpha, dt)

lam = fo.lyapunov_duffing(x0, v0, k, m, b, F0, omega, alpha, dt)
print(f"Lyapunov exponent at omega = {omega}: {lam:.3f}")
