import numpy as np
import matplotlib.pyplot as plt

m = float(input("Mass m: "))
k = float(input("Spring constant k: "))
t_min = float(input("Start time: "))
t_max = float(input("End time: "))
N = int(input("Number of samples: "))
gamma = float(input("Damping constant gamma: "))
x0 = float(input("Initial displacement x0: "))
v0 = float(input("Initial velocity v0: "))


omega0 = np.sqrt(k/m)
alpha = gamma/2
t = np.linspace(t_min, t_max, N)
gamma_critical = 2*omega0
A0 = np.sqrt(x0**2 + (v0/omega0)**2)



def shm_underdamped(t, x0, v0, omega0, alpha):
    omega_d = np.sqrt(omega0**2 - alpha**2)

    C = x0
    D = (v0 + alpha*x0)/omega_d

    x = np.exp(-alpha*t)*(C*np.cos(omega_d*t) + D*np.sin(omega_d*t))
    v = np.exp(-alpha*t)*(-alpha*(C*np.cos(omega_d*t) + D*np.sin(omega_d*t))
        + omega_d*(-C*np.sin(omega_d*t) + D*np.cos(omega_d*t)))

    return x, v

def shm_critical(t, x0, v0, alpha):
    C1 = x0
    C2 = v0 + alpha*x0

    x = (C1 + C2*t)*np.exp(-alpha*t)
    v = C2*np.exp(-alpha*t) - alpha*(C1 + C2*t)*np.exp(-alpha*t)

    return x, v

def shm_overdamped(t, x0, v0, omega0, alpha):
    r1 = -alpha + np.sqrt(alpha**2 - omega0**2)
    r2 = -alpha - np.sqrt(alpha**2 - omega0**2)

    C = (v0 - r2*x0)/(r1 - r2)
    D = (r1*x0 - v0)/(r1 - r2)

    x = C*np.exp(r1*t) + D*np.exp(r2*t)
    v = C*r1*np.exp(r1*t) + D*r2*np.exp(r2*t)

    return x, v


if alpha < omega0:
    print("System is UNDERDAMPED")
    x, v = shm_underdamped(t, x0, v0, omega0, alpha)
    label = "Underdamped"

elif np.isclose(alpha, omega0):
    print("System is CRITICALLY DAMPED")
    x, v = shm_critical(t, x0, v0, alpha)
    label = "Critical damping"

else:
    print("System is OVERDAMPED")
    x, v = shm_overdamped(t, x0, v0, omega0, alpha)
    label = "Overdamped"


plt.figure(figsize=(10,5))
plt.plot(t, x, lw=2)
plt.xlabel("Time (s)")
plt.ylabel("Displacement (m)")
plt.title(f"Damped Harmonic Oscillator ({label})")
plt.grid(True)
plt.show()

plt.figure(figsize=(10,5))
plt.plot(x, v, lw=2, color="blue")
plt.xlabel("Displacement (m)")
plt.ylabel("Velocity (m/s)")
plt.title(f"Phase Space ({label})")
plt.grid(True)
plt.show()

X = x / A0
V = v / (omega0 * A0)
plt.figure(figsize=(6,6))
plt.plot(X, V, lw=2, color="red")
plt.xlabel("Displacement")
plt.ylabel("Velocity")
plt.title(f"Normalised Phase Space ({label})")
plt.grid(True)
plt.show()


print(f"Natural frequency = {omega0:.3f}")
print(f"Critical gamma = {gamma_critical:.3f}")


