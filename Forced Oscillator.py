import numpy as np
import matplotlib.pyplot as plt

def rk4(prev_x, prev_v, k, m, t_max, N, b, F0, omega):
    dt = t_max/(N-1)
    t = np.linspace(0, t_max, N)
    x = np.zeros(N)
    v = np.zeros(N)
    E = np.zeros(N)
    
    x[0] = prev_x
    v[0] = prev_v
    E[0] = 0.5*m*v[0]**2 + 0.5*k*x[0]**2

    for i in range(N-1):
        k1x = v[i]
        k1v = (-b/m)*v[i] - (k/m)*x[i] + (F0/m)*np.cos(omega*t[i])
        k2x = v[i] + (dt/2)*k1v
        k2v = (-b/m)*(v[i] + (dt/2)*k1v) - (k/m)*(x[i] + (dt/2)*k1x) + (F0/m)*np.cos(omega*(t[i] + dt/2))
        k3x = v[i] + (dt/2)*k2v
        k3v = (-b/m)*(v[i] + (dt/2)*k2v) - (k/m)*(x[i] + (dt/2)*k2x) + (F0/m)*np.cos(omega*(t[i] + dt/2))
        k4x = v[i] + dt*k3v
        k4v = (-b/m)*(v[i] + dt*k3v) - (k/m)*(x[i] + dt*k3x) + (F0/m)*np.cos(omega*(t[i] + dt))
        
        x[i+1] = x[i] + dt/6*(k1x + 2*k2x + 2*k3x + k4x)
        v[i+1] = v[i] + dt/6*(k1v + 2*k2v + 2*k3v + k4v)
        
        E[i+1] = 0.5*m*v[i+1]**2 + 0.5*k*x[i+1]**2
    return t, x, v, E

def duffing(prev_x, prev_v, k, m, t_max, N, b, F0, omega, alpha, dt):
    t = np.linspace(0, t_max, N)

    x = np.zeros(N)
    v = np.zeros(N)

    x[0] = prev_x
    v[0] = prev_v

    for i in range(N-1):

        k1x = v[i]
        k1v = (-b/m)*v[i] - (k/m)*x[i] - (alpha/m)*x[i]**3 + (F0/m)*np.cos(omega*t[i])

        k2x = v[i] + (dt/2)*k1v
        k2v = (-b/m)*(v[i] + (dt/2)*k1v) - (k/m)*(x[i] + (dt/2)*k1x) \
              - (alpha/m)*(x[i] + (dt/2)*k1x)**3 + (F0/m)*np.cos(omega*(t[i]+dt/2))

        k3x = v[i] + (dt/2)*k2v
        k3v = (-b/m)*(v[i] + (dt/2)*k2v) - (k/m)*(x[i] + (dt/2)*k2x) \
              - (alpha/m)*(x[i] + (dt/2)*k2x)**3 + (F0/m)*np.cos(omega*(t[i]+dt/2))

        k4x = v[i] + dt*k3v
        k4v = (-b/m)*(v[i] + dt*k3v) - (k/m)*(x[i] + dt*k3x) \
              - (alpha/m)*(x[i] + dt*k3x)**3 + (F0/m)*np.cos(omega*(t[i]+dt))

        x[i+1] = x[i] + dt/6*(k1x + 2*k2x + 2*k3x + k4x)
        v[i+1] = v[i] + dt/6*(k1v + 2*k2v + 2*k3v + k4v)

    return t, x, v

def duffing_step(x, v, t, dt, k, m, b, F0, omega, alpha):

    k1x = v
    k1v = (-b/m)*v - (k/m)*x - (alpha/m)*x**3 + (F0/m)*np.cos(omega*t)

    k2x = v + 0.5*dt*k1v
    k2v = (-b/m)*(v + 0.5*dt*k1v) - (k/m)*(x + 0.5*dt*k1x) \
          - (alpha/m)*(x + 0.5*dt*k1x)**3 + (F0/m)*np.cos(omega*(t+0.5*dt))

    k3x = v + 0.5*dt*k2v
    k3v = (-b/m)*(v + 0.5*dt*k2v) - (k/m)*(x + 0.5*dt*k2x) \
          - (alpha/m)*(x + 0.5*dt*k2x)**3 + (F0/m)*np.cos(omega*(t+0.5*dt))

    k4x = v + dt*k3v
    k4v = (-b/m)*(v + dt*k3v) - (k/m)*(x + dt*k3x) \
          - (alpha/m)*(x + dt*k3x)**3 + (F0/m)*np.cos(omega*(t+dt))

    x_new = x + dt/6*(k1x + 2*k2x + 2*k3x + k4x)
    v_new = v + dt/6*(k1v + 2*k2v + 2*k3v + k4v)

    return x_new, v_new

def lyapunov_duffing(prev_x, prev_v, k, m, b, F0, omega, alpha, dt):
    chunk_steps = 20
    renorm_steps = 2000
    delta_t = chunk_steps*dt
    delta_0 = 10**-8 
    lyap_sum = 0
    t = 0
    x1_current = prev_x
    v1_current = prev_v
    x2_current = prev_x + delta_0
    v2_current = prev_v 
    
    for i in range(renorm_steps):
        for j in range(chunk_steps):
            x1_current, v1_current = duffing_step(x1_current, v1_current, t, dt, k, m, b, F0, omega, alpha)
            x2_current, v2_current = duffing_step(x2_current, v2_current, t, dt, k, m, b, F0, omega, alpha)
            
            t += dt
            
        delta = np.sqrt((x1_current - x2_current)**2 + (v1_current - v2_current)**2)
        if delta == 0:
            delta = 1e-16
        
        lyap_sum += np.log(delta/delta_0)
        
        dx = x2_current - x1_current
        dv = v2_current - v1_current
        
        x2_current = x1_current + delta_0*dx/delta        
        v2_current = v1_current + delta_0*dv/delta        
    total_time = renorm_steps * chunk_steps * dt
    lyap_lambda = lyap_sum/total_time
    return lyap_lambda

def plot_linear(t, x, v, E):
    plt.figure()
    plt.title("Displacement vs time")
    plt.xlabel("Time (s)")
    plt.ylabel("Displacement (m)")
    plt.plot(t, x, color = "blue")
    plt.savefig("displacement.png", dpi=300)
    plt.show(block=False)
    plt.pause(0.1)
        
    plt.figure()
    plt.title("Phase Space")
    plt.xlabel("Displacement (m)")
    plt.ylabel("Velocity (m/s)")
    plt.plot(x, v, color = "red")
    plt.axis("equal")
    plt.savefig("phase_space.png", dpi=300)
    plt.show(block=False)
    plt.pause(0.1)

    plt.figure()
    plt.title("Energy Variation")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (J)")
    plt.plot(t, E, color = "green")
    plt.savefig("energy.png", dpi=300)
    plt.show(block=False)
    plt.pause(0.1)

def amplitudes_linear(prev_x, prev_v, k, m, t_max, N, b, F0, omega_0):
    omega_max = 2*omega_0
    omega_N = 500
    omega_values = np.linspace(0.01, omega_max, omega_N)
    amplitude_values = np.zeros(omega_N)
    
    for i in range(omega_N):
        current_omega = omega_values[i]
        t, x, v, E = rk4(prev_x, prev_v, k, m, t_max, N, b, F0, current_omega)
        cut_index=int(0.7*N)
        steady_state = x[cut_index:]
        amplitude = (np.max(steady_state) - np.min(steady_state))/2
        amplitude_values[i] = amplitude
        
    return amplitude_values, omega_values

def amplitude_plot(omega_0, omega_values, amplitude_values):
    plt.figure()
    plt.title("Amplitude as a function of driving frequency")
    plt.xlabel("Driving frequency (rad/s)")
    plt.ylabel("Amplitude")
    plt.plot(omega_values, amplitude_values, color = "black")
    plt.axvline(omega_0, color = "green", linestyle = "--", label = "Natural frequency")
    plt.legend()
    plt.savefig("amplitude.png", dpi=300)
    plt.show(block=False)
    plt.pause(0.1)

def phase_plot(b, m, k, omega_0, omega_values):
    phi = np.arctan2(b*omega_values,(k - m*omega_values**2))
    plt.figure()
    plt.title("Phase as a function of driving frequency")
    plt.xlabel("Driving frequency (rad/s)")
    plt.ylabel("Phase (rad)")
    plt.plot(omega_values, phi, color = "black")
    plt.axvline(omega_0, color = "green", linestyle = "--", label = "Natural frequency")
    plt.legend()
    plt.savefig("phase.png", dpi=300)
    plt.show(block=False)
    plt.pause(0.1)
    
def plot_duffing(t, x, v, N, omega_0, t_max, dt):
    plt.figure()
    plt.title("Phase Space")
    plt.xlabel("Displacement (m)")
    plt.ylabel("Velocity (m/s)")
    plt.plot(x, v, color = "red")
    plt.axis("equal")
    plt.show(block=False)
    plt.pause(0.1)
    
    x_samples = []
    v_samples = []
    T = (2*np.pi)/omega_0
    steps_per_period = int(T/dt)
    start_index = int(0.7*N)
    for i in range(start_index, N, steps_per_period):
        x_samples.append(x[i])
        v_samples.append(v[i])
        
    plt.figure()
    plt.title("Poincare Section")
    plt.xlabel("Displacement (m)")
    plt.ylabel("Velocity (m/s)")
    plt.scatter(x_samples, v_samples, c=range(len(x_samples)), s=5, cmap="plasma")
    plt.colorbar(label="time")
    plt.axis("equal")
    plt.savefig("poincare.png", dpi=300)
    plt.show(block=False)
    plt.pause(0.1)

def bifurcation(prev_x, prev_v, k, m, t_max, N, b, F0, omega, alpha, dt):
    F_values = np.linspace(0, F0, 200)
    plt.figure()
    x_current = prev_x
    v_current = prev_v
           
    T = 2*np.pi/omega
    steps_per_period = int(T/dt)
        
    for F in F_values:
        t, x, v = duffing(x_current, v_current, k, m, t_max, N, b, F, omega, alpha, dt)
        x_current = x[-1]
        v_current = v[-1]

        for i in range(int(0.7*N), N, steps_per_period):
            plt.plot(F, x[i], "k.", markersize=0.5)
            
    plt.xlabel("Driving amplitude F0")
    plt.ylabel("Displacement (Poincare points)")
    plt.title("Duffing bifurcation plot")
    plt.savefig("bifurcation.png", dpi=300)
    plt.show()

def chaos_map(prev_x, prev_v, k, m, b, F0, omega, alpha, dt):
    omega_max = 2*omega
    omega_N = 75
    omega_values = np.linspace(0.2, omega_max, omega_N)
    lambda_values = np.zeros(len(omega_values))
    
    for i in range(len(omega_values)):
        omega = omega_values[i]
        lambda_values[i] = lyapunov_duffing(prev_x, prev_v, k, m, b, F0, omega, alpha, dt)
        
    plt.figure()
    plt.title("Duffing chaos map")
    plt.xlabel("Driving frequency")
    plt.ylabel("Lyapunov exponent")
    plt.plot(omega_values, lambda_values)
    plt.axhline(0, linestyle="--")
    plt.savefig("chaos_map.png", dpi=300)
    plt.show()
    
def main():
    prev_x = float(input("What is the initial displacement? "))
    prev_v = float(input("What is the initial velocity? "))
    t_max = float(input("What is the total time of the measurements? "))
    N = int(input("How many measurements were taken? "))
    k = float(input("What is the spring constant of the system? "))
    m = float(input("What is the mass of the system? "))
    b = float(input("What is the damping constant of the system? "))
    F0 = float(input("What is the initial force applied? "))
    dt = t_max/(N-1)

    choice = str(input("Which oscillator would you like? Linear (L), Duffing (D) ")).upper()
    if choice == "L":
        omega_0 = np.sqrt(k/m)
        t, x, v, E = rk4(prev_x, prev_v, k, m, t_max, N, b, F0, omega_0)
        plot_linear(t, x, v, E)
        amplitude_values, omega_values = amplitudes_linear(prev_x, prev_v, k, m, t_max, N, b, F0, omega_0)
        amplitude_plot(omega_0, omega_values, amplitude_values)
        phase_plot(b, m, k, omega_0, omega_values)
    else:
        alpha = float(input("What is the coefficient of x^3? "))
        omega = float(input("What is the driving frequency? "))
        t, x, v = duffing(prev_x, prev_v, k, m, t_max, N, b, F0, omega, alpha, dt)
        plot_duffing(t, x, v, N, omega, t_max, dt)
        bifurcation(prev_x, prev_v, k, m, t_max, N, b, F0, omega, alpha, dt)
        chaos_map(prev_x, prev_v, k, m, b, F0, omega, alpha, dt)
        
if __name__ == "__main__":
    main()
