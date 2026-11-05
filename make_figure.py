"""
Use the integrators from Numerical_Oscillator.py with:
m = 1 kg, k = 1 N/m, x0 = 1 m, v0 = 0 m/s, t_max = 100 s, N = 10001 (dt = 0.01 s).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Numerical_Oscillator import euler, verlet, rk4

x0, v0, k, m, t_max, N = 1.0, 0.0, 1.0, 1.0, 100.0, 10001
t = np.linspace(0, t_max, N)

xe, ve, Ee, E0e = euler(x0, v0, t_max, N, k, m)
xv, vv, Ev, E0v = verlet(x0, v0, k, m, t_max, N)
xr, vr, Er, E0r = rk4(x0, v0, k, m, t_max, N)

err_e = (Ee - E0e) / E0e
err_v = np.abs((Ev - E0v) / E0v)
err_r = np.abs((Er - E0r) / E0r)

print(f"Euler   final relative energy error: {err_e[-1]*100:.1f}%  (max {np.max(np.abs(err_e))*100:.1f}%)")
print(f"Verlet  max |relative energy error|: {np.max(err_v):.2e}")
print(f"RK4     max |relative energy error|: {np.max(err_r):.2e}")
print(f"RK4     final |relative energy error|: {err_r[-1]:.2e}")
print(f"Verlet  final |relative energy error|: {err_v[-1]:.2e}")
print("Euler E(t_max)/E(0) =", Ee[-1] / E0e)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
ax1.plot(xe, ve, color="tab:blue", lw=1, label="Euler")
ax1.plot(xv, vv, color="tab:red", lw=1.5, label="Verlet")
ax1.plot(xr, vr, color="tab:green", lw=1, ls="--", label="RK4")
ax1.set_xlabel("Displacement (m)")
ax1.set_ylabel("Velocity (m/s)")
ax1.set_title("Phase space")
ax1.set_aspect("equal")
ax1.legend(loc="upper right")

ax2.semilogy(t, np.maximum(np.abs(err_e), 1e-16), color="tab:blue", label="Euler")
ax2.semilogy(t, np.maximum(err_v, 1e-16), color="tab:red", label="Verlet")
ax2.semilogy(t, np.maximum(err_r, 1e-16), color="tab:green", label="RK4")
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("|Relative energy error|")
ax2.set_title("Energy conservation")
ax2.legend(loc="lower right")

plt.tight_layout()
plt.savefig("energy_error_comparison.png", dpi=200)
