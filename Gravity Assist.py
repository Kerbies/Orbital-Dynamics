import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def main():
    # Masses 
    M_sun = 1
    M_Jupiter = 9.55e-4

    # Gravitational constant
    G = 4*np.pi**2

    # Time parameters (years)
    t0 = 0
    tf = 2.5
    dt = 0.001
    t = np.arange(t0, tf, dt)

    # Path of Jupiter w.r.t. the Sun
    Initial_phase = np.pi/2 + 0.164
    x_Jupiter = 5.2 * np.cos(2 * np.pi * t / 11.86 + Initial_phase) # AU
    y_Jupiter = 5.2 * np.sin(2 * np.pi * t / 11.86 + Initial_phase) # AU

    # Distance to Jupiter
    r_jupiter0 = np.sqrt( (1-x_Jupiter[0])**2 + y_Jupiter[0]**2)

    # Initial conditions for Voyager 2
    x0_Voyager2 = 1        # AU
    y0_Voyager2 = 0

    vx0_Voyager2 = 0
    vy0_Voyager2 = 8.5 # AU/year

    ax0_Voyager2 = G * (-M_sun * x0_Voyager2 - M_Jupiter / r_jupiter0**3 * (x0_Voyager2 - x_Jupiter[0]))
    ay0_Voyager2 = G * (-M_sun * y0_Voyager2 - M_Jupiter / r_jupiter0**3 * (y0_Voyager2 - y_Jupiter[0]))

    # Euler method 
    x_eul, y_eul, vx_eul, vy_eul = three_body_euler(t, dt, x0_Voyager2, y0_Voyager2, vx0_Voyager2, vy0_Voyager2, ax0_Voyager2, ay0_Voyager2)

    # Leapfrog method
    x_leap, y_leap, vx_leap, vy_leap = three_body_leapfrog(t, dt, x0_Voyager2, y0_Voyager2, vx0_Voyager2, vy0_Voyager2, ax0_Voyager2, ay0_Voyager2)

    # Plotting
    plt.figure(figsize=(6, 6))
    plt.plot(x_eul, y_eul, label='Euler', color='red')
    plt.plot(x_leap, y_leap, label='Leapfrog', color='blue')
    plt.plot(x_Jupiter, y_Jupiter, color = "green")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Orbit of the planet')
    plt.xlim(-6, 6)
    plt.ylim(-6, 6)
    plt.legend()
    plt.show()

    v_leap = np.sqrt(vx_leap**2 + vy_leap**2)
    v_eul = np.sqrt(vx_eul**2 + vy_eul**2)
    plt.figure(figsize=(6, 6))
    plt.plot(t, v_leap, label='Euler', color='red')
    plt.plot(t, v_eul, label='Leapfrog', color='blue')
    plt.xlabel('Time (years)')
    plt.ylabel('Speed (AU/year)')
    plt.title('Speed of the Voyager 2')
    plt.legend()
    plt.show()

    def animate(i):
        plt.clf()
        plt.plot(x_leap[:i], y_leap[:i], label='Voyager 2', color='blue')
        plt.plot(x_Jupiter[:i], y_Jupiter[:i], label = "Jupiter", color = "green")
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Orbit of the planet')
        plt.legend()
        plt.xlim(-6, 1.2)
        plt.ylim(-0, 6)

        # Add time subtitle

        t_months = t[i] * 12
        speed_km_s = v_leap[i] * 4.743 # Convert AU/year to km/s

        plt.text(
            0.05, 0.95, f't = {t_months:.2f} months \nSpeed = {speed_km_s:.2f} km/s',
            transform=plt.gca().transAxes,
            fontsize=12,
            verticalalignment='top'
        )

    fig = plt.figure(figsize=(6, 6))
    ani = animation.FuncAnimation(fig, animate, frames=len(t), interval=30)

    print("Saving animation...")
    ani.save('orbit_animation.gif', writer='pillow')
    print("Animation saved as 'orbit_animation.gif'")




# Leapfrog integration
def three_body_leapfrog(t, dt, x0, y0, vx0, vy0, ax0, ay0):

    # Gravitational constant
    G = 4*np.pi**2

    # Masses
    M_sun = 1
    M_Jupiter = 9.55e-4

    # Number of steps 
    N = len(t)

    # Position, velocity and acceleration arrays
    x = np.zeros(N)
    y = np.zeros(N)

    vx = np.zeros(N)
    vy = np.zeros(N)

    ax = np.zeros(N)
    ay = np.zeros(N)

    # Set initial conditions
    x[0] = x0
    y[0] = y0

    vx[0] = vx0
    vy[0] = vy0

    ax[0] = ax0
    ay[0] = ay0

    # Path of Jupiter w.r.t. the Sun
    Initial_phase = np.pi/2 + 0.164
    x_Jupiter = 5.2 * np.cos(2 * np.pi * t / 11.86 + Initial_phase) # AU
    y_Jupiter = 5.2 * np.sin(2 * np.pi * t / 11.86 + Initial_phase) # AU

    # leapfrog integration
    for i in range(1, N):

        # Update velocity at half step
        vx_half = vx[i-1] + ax[i-1] * dt / 2
        vy_half = vy[i-1] + ay[i-1] * dt / 2

        # Update position at full step
        x[i] = x[i-1] + vx_half * dt
        y[i] = y[i-1] + vy_half * dt

        # New distance from the sun
        r_sun = np.sqrt(x[i]**2 + y[i]**2)

        # Distance to Jupiter
        r_Jupiter = np.sqrt((x[i] - x_Jupiter[i])**2 + (y[i] - y_Jupiter[i])**2)

        # Update acceleration at full step
        ax_sun = -G * M_sun * x[i] / r_sun**3
        ay_sun = -G * M_sun * y[i] / r_sun**3

        ax_jupiter = -G * M_Jupiter * (x[i] - x_Jupiter[i]) / r_Jupiter**3
        ay_jupiter = -G * M_Jupiter * (y[i] - y_Jupiter[i]) / r_Jupiter**3
        
        ax[i] = ax_sun + ax_jupiter
        ay[i] = ay_sun + ay_jupiter

        # Update velocity at full step
        vx[i] = vx_half + ax[i] * dt / 2
        vy[i] = vy_half + ay[i] * dt / 2

    return x, y, vx, vy

# Euler integration
def three_body_euler(t, dt, x0, y0, vx0, vy0, ax0, ay0):

    # Gravitational constant
    G = 4*np.pi**2

    # Masses
    M_sun = 1
    M_Jupiter = 9.55e-4

    # Number of steps
    N = len(t)

    # Position, velocity and acceleration arrays
    x = np.zeros(N)
    y = np.zeros(N)

    vx = np.zeros(N)
    vy = np.zeros(N)

    ax = np.zeros(N)
    ay = np.zeros(N)

    # Set initial conditions
    x[0] = x0
    y[0] = y0

    vx[0] = vx0
    vy[0] = vy0

    ax[0] = ax0
    ay[0] = ay0

    # Path of Jupiter w.r.t. the Sun
    Initial_phase = np.pi/2 + 0.164
    x_Jupiter = 5.2 * np.cos(2 * np.pi * t / 11.86 + Initial_phase) # AU
    y_Jupiter = 5.2 * np.sin(2 * np.pi * t / 11.86 + Initial_phase) # AU

    # Euler integration
    for i in range(1, N):

        # Update position at full step
        x[i] = x[i-1] + vx[i-1] * dt
        y[i] = y[i-1] + vy[i-1] * dt

        # Update velocity at full step
        vx[i] = vx[i-1] + ax[i-1] * dt
        vy[i] = vy[i-1] + ay[i-1] * dt

        # New distance from the sun
        r_sun = np.sqrt(x[i]**2 + y[i]**2)

        # Distance to Jupiter
        r_Jupiter = np.sqrt((x[i] - x_Jupiter[i])**2 + (y[i] - y_Jupiter[i])**2)

        # Update acceleration at full step
        ax_sun = -G * M_sun * x[i] / r_sun**3
        ay_sun = -G * M_sun * y[i] / r_sun**3

        ax_jupiter = -G * M_Jupiter * (x[i] - x_Jupiter[i]) / r_Jupiter**3
        ay_jupiter = -G * M_Jupiter * (y[i] - y_Jupiter[i]) / r_Jupiter**3
        
        ax[i] = ax_sun + ax_jupiter
        ay[i] = ay_sun + ay_jupiter

    return x, y, vx, vy

if __name__ == "__main__":
    main()
