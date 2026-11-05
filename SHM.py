import numpy as np
import matplotlib.pyplot as plt

m = float(input("What is the mass of the system? "))
k = float(input("What is the spring constant of the system? "))
A = float(input("What is the amplitude of the oscillations? "))
phi = np.deg2rad(float(input("What is the phase constant of the system? ")))
t_min = float(input("What is the starting time?"))
t_max = float(input("What is the end time?"))
delta_t = float(input("What is the time step size?"))

omega = np.sqrt(k/m)

t = np.arange(t_min, t_max, delta_t)
x = A*np.cos(omega*t + phi)
v = -A*omega*np.sin(omega*t + phi)
K = 0.5*m*v**2
U = 0.5*k*x**2
E = K + U

plt.figure()
plt.plot(t, x, label="x(t)", color="red")
plt.plot(t, v, label="v(t)", color="blue")
plt.xlabel("Time (s)")
plt.ylabel("Displacement (m) / Velocity (m/s)")
plt.title("Simple Harmonic Motion")
plt.legend()
plt.show()

plt.figure()
plt.plot(t, K, label="K(t)", color="blue")
plt.plot(t, U, label="U(t)", color="red")
plt.plot(t, E, label="E(t)", color="black")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.title("SHM Energies")
plt.legend()
plt.show()

E_theory = 0.5*k*A**2
print("Theoretical energy:", E_theory)
print("Mean numerical energy:", np.mean(E))
print("Max deviation:", np.max(np.abs(E - E_theory)))