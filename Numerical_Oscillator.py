import numpy as np
import matplotlib.pyplot as plt

def euler(prev_x, prev_v, t_max, N, k, m):
    dt = t_max/(N-1)
    x = np.zeros(N)
    v = np.zeros(N)
    E = np.zeros(N)
    
    x[0] = prev_x
    v[0] = prev_v
    E[0] = 0.5*m*v[0]**2 + 0.5*k*x[0]**2

    for i in range(N-1):
        a = -(k/m)*x[i]
        
        x[i+1] = x[i] + v[i]*dt
        v[i+1] = v[i] + a*dt
    
        E[i+1] = 0.5*m*v[i+1]**2 + 0.5*k*x[i+1]**2
    return x, v, E, E[0]

    
def verlet(prev_x, prev_v, k, m, t_max, N):
    dt = t_max/(N-1)
    x = np.zeros(N)
    v = np.zeros(N)
    a = np.zeros(N)
    E = np.zeros(N)
    
    x[0] = prev_x
    v[0] = prev_v
    a[0] = -(k/m)*x[0]
    E[0] = 0.5*m*v[0]**2 + 0.5*k*x[0]**2

    for i in range(N-1):   
        x[i+1] = x[i] + v[i]*dt + 0.5*a[i]*dt**2
        a[i+1] = -(k/m)*x[i+1]
        v[i+1] = v[i] + 0.5*(a[i] + a[i+1])*dt
        
        E[i+1] = 0.5*m*v[i+1]**2 + 0.5*k*x[i+1]**2

    return x, v, E, E[0]

def rk4(prev_x, prev_v, k, m, t_max, N):
    dt = t_max/(N-1)
    x = np.zeros(N)
    v = np.zeros(N)
    E = np.zeros(N)
    
    x[0] = prev_x
    v[0] = prev_v
    E[0] = 0.5*m*v[0]**2 + 0.5*k*x[0]**2
    
    for i in range(N-1):
        k1x = v[i]
        k1v = (-k/m)*x[i]
        k2x = v[i] + (dt/2)*k1v
        k2v = (-k/m)*(x[i] + (dt/2)*k1x)
        k3x = v[i] + (dt/2)*k2v
        k3v = (-k/m)*(x[i] + (dt/2)*k2x)
        k4x = v[i] + dt*k3v
        k4v = (-k/m)*(x[i] + dt*k3x)
        
        x[i+1] = x[i] + dt/6*(k1x + 2*k2x + 2*k3x + k4x)
        v[i+1] = v[i] + dt/6*(k1v + 2*k2v + 2*k3v + k4v)
        
        E[i+1] = 0.5*m*v[i+1]**2 + 0.5*k*x[i+1]**2
        
    return x, v, E, E[0]
        
def plot(t, x, v, E, E_0):
    plt.figure()
    plt.title("Displacement vs time")
    plt.xlabel("Time (s)")
    plt.ylabel("Displacement (m)")
    plt.plot(t, x, color = "blue")
    plt.show()
        
    plt.figure()
    plt.title("Phase Space")
    plt.xlabel("Displacement (m)")
    plt.ylabel("Velocity (m/s)")
    plt.plot(x, v, color = "red")
    plt.axis("equal")
    plt.show()
    
    plt.figure()
    plt.title("Energy Variation")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (J)")
    plt.plot(t, E, color = "green")
    plt.show()
    
    plt.figure()
    plt.title("Relative Energy Error")
    plt.xlabel("Time (s)")
    plt.ylabel("Fractional Energy Difference (J)")
    plt.plot(t, (E-E_0)/E_0, color = "black")
    plt.show()
    
def comparison(t, x1, v1, E1, E_0_1, x2, v2, E2, E_0_2, x3, v3, E3, E_0_3):
    plt.figure()
    plt.title("Euler Displacement ")
    plt.xlabel("Time (s)")
    plt.ylabel("Displacement (m)")
    plt.plot(t, x1, color = "blue", label = "Euler")
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.title("Displacement Comparison")
    plt.xlabel("Time (s)")
    plt.ylabel("Displacement (m)")
    plt.plot(t, x2, color = "red", label = "Verlet")
    plt.plot(t, x3, color = "green", label = "RK4")
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.title("Euler Phase Space ")
    plt.xlabel("Displacement (m)")
    plt.ylabel("Velocity (m/s)")
    plt.plot(x1, v1, color = "blue", label = "Euler")
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.title("Phase Space Comparison")
    plt.xlabel("Displacement (m)")
    plt.ylabel("Velocity (m/s)")
    plt.plot(x2, v2, color = "red", label = "Verlet")
    plt.plot(x3, v3, color = "green", label = "RK4")
    plt.legend()
    plt.axis("equal")
    plt.show()
    
    plt.figure()
    plt.title("Euler Energy Variation")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (J)")
    plt.plot(t, E1, color = "blue", label = "Euler")
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.title("Energy Variation Comparison")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (J)")
    plt.plot(t, E2, color = "red", label = "Verlet")
    plt.plot(t, E3, color = "green", label = "RK4")
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.title("Euler Energy Error")
    plt.xlabel("Time (s)")
    plt.ylabel("Fractional Energy Difference")
    plt.plot(t, (E1-E_0_1)/E_0_1, color = "blue", label = "Euler")
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.title("Relative Energy Error")
    plt.xlabel("Time (s)")
    plt.ylabel("Fractional Energy Difference")
    plt.plot(t, abs((E2-E_0_2)/E_0_2), color = "red", label = "Verlet")
    plt.plot(t, abs((E3-E_0_3)/E_0_3), color = "green", label = "RK4")
    plt.yscale("log")
    plt.legend()
    plt.show()
    
def main():
    
    prev_x = float(input("What is the initial displacement? "))
    prev_v = float(input("What is the initial velocity? "))
    t_max = float(input("What is the total time of the measurements? "))
    N = int(input("How many measurements were taken? "))
    k = float(input("What is the spring constant of the system? "))
    m = float(input("What is the mass of the system? "))

    t = np.linspace(0, t_max, N)
    
    choice = str(input("Which simulation would you like? Euler (E), Verlet (V), Runge-Kutta 4 (R), comparison (C) ")).upper()
    
    if choice == "E":
        x, v, E, E0 = euler(prev_x, prev_v, t_max, N, k, m)
        plot(t, x, v, E, E0)
        
    elif choice == "V":
        x, v, E, E0 = verlet(prev_x, prev_v, k, m, t_max, N)
        plot(t, x, v, E, E0)
        
    elif choice == "R":
        x, v, E, E0 = rk4(prev_x, prev_v, k, m, t_max, N)
        plot(t, x, v, E, E0)
        
    elif choice == "C":
        x_euler, v_euler, E_euler, E_0_euler = euler(prev_x, prev_v, t_max, N, k, m)
        x_verlet, v_verlet, E_verlet, E_0_verlet = verlet(prev_x, prev_v, k, m, t_max, N)
        x_rk4, v_rk4, E_rk4, E_0_rk4 = rk4(prev_x, prev_v, k, m, t_max, N)
        comparison(t, x_euler, v_euler, E_euler, E_0_euler, x_verlet, v_verlet, E_verlet, E_0_verlet, x_rk4, v_rk4, E_rk4, E_0_rk4)

if __name__ == "__main__":
    main()
