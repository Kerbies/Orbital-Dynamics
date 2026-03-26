import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

def main():

    # Gravitational constant
    G = 4*np.pi**2

    # Mass of the sun
    M = 1

    # time
    t0 = 0
    tf = 3
    dt = 0.0001
    t = np.arange(t0, tf, dt)

    # inital position 
    x0 = 1
    y0 = 0

    # inital velocity
    vx0 = 0
    vy0 = 2*np.pi

    # inital Acceleration
    ax0 = -G * M / x0**2 
    ay0 = 0

    # Euler method 
    x_eul, y_eul, vx_eul, vy_eul = euler(t, dt, x0, y0, vx0, vy0, ax0, ay0)

    # Leapfrog method
    x_leap, y_leap, vx_leap, vy_leap = leapfrog(t, dt, x0, y0, vx0, vy0, ax0, ay0)

    # Plotting
    plt.figure(figsize=(6, 6))
    plt.plot(x_eul, y_eul, label='Euler', color='red')
    plt.plot(x_leap, y_leap, label='Leapfrog', color='blue')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Orbit of the planet')
    plt.legend()
    plt.show()

    # Plotting speed
    speed_eul = np.sqrt(vx_eul**2 + vy_eul**2)
    speed_leap = np.sqrt(vx_leap**2 + vy_leap**2)

    plt.figure(figsize=(10, 6))
    plt.plot(t, speed_eul, label='Euler', color='red')
    plt.plot(t, speed_leap, label='Leapfrog', color='blue')
    plt.xlabel('Time (years)')
    plt.ylabel('Speed (AU/year)')
    plt.title('Speed of the planet over time')
    plt.legend()
    plt.show()

    # Calculate Specific Energy
    Energy_eul = 0.5 * speed_eul**2 - G * M / np.sqrt(x_eul**2 + y_eul**2)
    Energy_leap = 0.5 * speed_leap**2 - G * M / np.sqrt(x_leap**2 + y_leap**2)

    # Plot Energy 
    plt.figure(figsize=(10, 6))
    plt.plot(t, Energy_eul, label='Euler', color='red')
    plt.plot(t, Energy_leap, label='Leapfrog', color='blue')
    plt.xlabel('Time (years)')
    plt.ylabel('Specific Energy (AU^2/year^2)')
    plt.title('Specific Energy of the planet over time')
    plt.legend()
    plt.show()

    # Repeat Calculations for smaller inital velocity and for five orbits
    t0 = 0
    tf = 5
    dt = 0.0001
    t = np.arange(t0, tf, dt)

    # inital velocity
    vy0 = vy0 * 0.8

    # Euler method 
    x_eul, y_eul, vx_eul, vy_eul = euler(t, dt, x0, y0, vx0, vy0, ax0, ay0)

    # Leapfrog method
    x_leap, y_leap, vx_leap, vy_leap = leapfrog(t, dt, x0, y0, vx0, vy0, ax0, ay0)

    # Plotting new trajectory 
    plt.figure(figsize=(6, 6))
    plt.plot(x_eul, y_eul, label='Euler', color='red')
    plt.plot(x_leap, y_leap, label='Leapfrog', color='blue')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Orbit of the planet')
    plt.legend()
    plt.show()

    # Calculate Specific Energy

    speed_eul = np.sqrt(vx_eul**2 + vy_eul**2)
    speed_leap = np.sqrt(vx_leap**2 + vy_leap**2)

    Energy_eul = 0.5 * speed_eul**2 - G * M / np.sqrt(x_eul**2 + y_eul**2)
    Energy_leap = 0.5 * speed_leap**2 - G * M / np.sqrt(x_leap**2 + y_leap**2)

    # Plot Energy
    plt.figure(figsize=(10, 6))
    plt.plot(t, Energy_eul, label='Euler', color='red')
    plt.plot(t, Energy_leap, label='Leapfrog', color='blue')
    plt.xlabel('Time (years)')
    plt.ylabel('Specific Energy (AU^2/year^2)')
    plt.title('Specific Energy of the planet over time')
    plt.legend()
    plt.show()

    # Animation
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    line_eul, = ax.plot([], [], 'r-', label='Euler', lw = 1)
    line_leap, = ax.plot([], [], 'b-', label='Leapfrog', lw = 1)
    ax.legend()

    def animate(i):

        line_eul.set_data(x_eul[:i], y_eul[:i])
        line_leap.set_data(x_leap[:i], y_leap[:i])

        return line_eul, line_leap,


    
    #ani = animation.FuncAnimation(fig, animate, frames=len(t), interval=1, blit=True)

    #ani.save('Orbit.mp4', writer='pillow', fps=30)
    #print("Animation saved as Orbit.mp4")

    plt.show()

# Leapfrog integration
def leapfrog(t, dt, x0, y0, vx0, vy0, ax0, ay0):

    # Gravitational constant
    G = 4*np.pi**2

    # Mass of the sun
    M = 1

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

    # leapfrog integration
    for i in range(1, N):

        # Update velocity at half step
        vx_half = vx[i-1] + ax[i-1] * dt / 2
        vy_half = vy[i-1] + ay[i-1] * dt / 2

        # Update position at full step
        x[i] = x[i-1] + vx_half * dt
        y[i] = y[i-1] + vy_half * dt

        # Update acceleration at full step
        r = np.sqrt(x[i]**2 + y[i]**2)
        ax[i] = -G * M * x[i] / r**3
        ay[i] = -G * M * y[i] / r**3

        # Update velocity at full step
        vx[i] = vx_half + ax[i] * dt / 2
        vy[i] = vy_half + ay[i] * dt / 2

    return x, y, vx, vy

# Euler integration
def euler(t, dt, x0, y0, vx0, vy0, ax0, ay0):

    # Gravitational constant
    G = 4*np.pi**2

    # Mass of the sun
    M = 1

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

    # Euler integration
    for i in range(1, N):

        # Update position at full step
        x[i] = x[i-1] + vx[i-1] * dt
        y[i] = y[i-1] + vy[i-1] * dt

        # Update velocity at full step
        vx[i] = vx[i-1] + ax[i-1] * dt
        vy[i] = vy[i-1] + ay[i-1] * dt

        # Update acceleration at full step
        r = np.sqrt(x[i]**2 + y[i]**2)
        
        ax[i] = -G * M * x[i] / r**3
        ay[i] = -G * M * y[i] / r**3

    return x, y, vx, vy


    return


    

if __name__ == "__main__":
    main()