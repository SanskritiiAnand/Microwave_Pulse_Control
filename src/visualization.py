"""
Visualization utilities for pulse-level simulations of superconducting qubits.
Provides functions for visualizing state populations, Bloch-vector components, and trajectories on the Bloch sphere
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt

#State Populations
def plot_state_populations(time: np.ndarray, ground_population: np.ndarray, excited_population: np.ndarray):
    """
    plots the ground and excited stqate populations of a qubit.

    Parameters:
     time: np.ndarray (simulation time points)
     ground_population: np.ndarray (population of the |0> state)
     excited_population: np.ndarray (population of the |1> state)

    Returns:
     matplotlib.axes.Axes (axes containing the population plot)
    """
    fig, ax= plt.subplots()

    ax.plot(time, ground_population, label=r"$P_0$")
    ax.plot(time, excited_population, label=r"$P_1$")

    ax.set_xlabel("time")
    ax.set_ylabel("population")
    ax.set_title("Qubit State Populations")

    ax.set_ylim(0.0, 1.05)
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    return ax

#Bloch Vector Components
def bloch_components(state_trajectory: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    calculate Bloch-vector components from a single-qubit state trajectory.
    For a qubit state |psi> = alpha|0> + beta|1>,
    the Bloch-vector components are:
      x= 2 Re(alpha*beta)
      y= 2 Im(alpha*beta)
      z= |alpha|^2 -|beta|^2

    Parameters:
     state_trajectory: np.ndarray (state trajectory withe shape (2, N), where N is the number of simulation time points)

    Returns:
     tuple[np.ndarray, np.ndarray, np.ndarray] (arrays containing the x, y, z bloch-vector components)
    """
    alpha= state_trajectory[0]
    beta= state_trajectory[1]

    x= 2.0 * np.real(np.conjugate(alpha)*beta)
    y= 2.0 * np.imag(np.conjugate(alpha)*beta)
    z= np.abs(alpha) ** 2 - np.abs(beta) ** 2

    return x,y,z

#Bloch Components Plot
def plot_bloch_components(time: np.ndarray, state_trajectory: np.ndarray):
    """
    plots the three components of the qubit Bloch vector versus time

    Parameters:
     time: np.ndarray (simulation tome points)
     state_trajectory: np.ndarray (state trajectory with shape (2,N))

    Returns:
     matplotlib.axes.Axes (axes containing the bloch-component plot)
    """
    x, y, z= bloch_components(state_trajectory)

    fig, ax= plt.subplots()

    ax.plot(time, x, label= r"$\langle X \rangle$")
    ax.plot(time, y, label= r"$\langle Y \rangle$")
    ax.plot(time, z, label= r"$\langle Z \rangle$")

    ax.set_xlabel("time")
    ax.set_ylabel("expectation value")
    ax.set_title("Qubit Bloch-Vector Components")

    ax.set_ylim(-1.05, 1.05)
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    return ax

#Bloch Sphere Trajectory
def plot_bloch_trajectory(state_trajectory: np.ndarray):
    """
    plots the qubit state trajectory on the Bloch sphere

    Parameters:
     state-trajectory: np.ndarray (state trajectory with shape (2,N))

    Returns:
     matplotlib.axes.Axes (3D axes containing the Bloch-sphere trajectory)
    """
    x, y, z= bloch_components(state_trajectory)

    fig= plt.figure()
    ax= fig.add_subplot(111, projection="3d")

    #Bloch Sphere
    u= np.linspace(0.0, 2.0*np.pi, 100)
    v= np.linspace(0.0, np.pi, 50)

    sphere_x= np.outer(np.cos(u), np.sin(v))
    sphere_y= np.outer(np.sin(u), np.sin(v))
    sphere_z= np.outer(np.ones_like(u), np.cos(v))

    ax.plot_wireframe(sphere_x, sphere_y, sphere_z, alpha=0.15, linewidth=0.5)

    #State Trajectory
    ax.plot(x, y, z, linewidth= 2, label="State Trajectory")

    #Initial and final states
    ax.scatter(x[0], y[0], z[0], s=50, label="Initial State")

    ax.scatter(x[-1], y[-1], z[-1], s=50, label="Final State")

    #Coordinate Axes
    ax.plot([-1,1], [0,0], [0,0], alpha=0.4)
    ax.plot([0,0], [-1,1], [0,0], alpha=0.4)
    ax.plot([0,0], [0,0], [-1,1], alpha=0.4)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    ax.set_title("Qubit State Trajectory on the Bloch Sphere")

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)

    ax.legend()

    fig.tight_layout()

    return ax

def plot_pulse_envelope(time: np.ndarray, envelope: np.ndarray, title: str= "Microwave Pulse Envelope") -> None:
    """
    Plots the time-dependent microwave pulse envelope

    Parameters:
     time: np.ndarray (time coordinates)
     envelope: np.ndarray (Pulse amplitude evaluated at each time)
     title: str, optional
    """
    fig, ax= plt.subplots(figsize=(8,4))

    ax.plot(time, envelope)

    ax.set_title(title)
    ax.set_xlabel("Time")
    ax.set_ylabel("Amplitude")

    ax.grid(True, alpha=0.3)

    plt.tight_layout()

def plot_drag_pulse(time,pulse):
    """
    plots the I and Q quadratures and magnitude of a DRAG pulse

    Parameters:
     time: np.ndarray (time coordinates)
     pulse: np.ndarray (complex-valued DRAG pulse envelope: pulse= I(t) + i Q(t))
    """
    I= np.real(pulse)
    Q= np.imag(pulse)
    magnitude= np.abs(pulse)

    fig, ax= plt.subplots(figsize=(10, 6))

    ax.plot(time, I, label="I(t)")
    ax.plot(time, Q, label="Q(t)")
    ax.plot(time, magnitude, label="|I + iq|")

    ax.set_title("DRAG microwave pulse")
    ax.set_xlabel("Time")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.show()

def plot_rabi_oscillations(time: np.ndarray, amplitudes: np.ndarray, rabi_populations: np.ndarray):
    """
    plots excited state population as a function of time for different microwave drive amplitudes

    Parameters:
     time: np.ndarray (simulation time points)
     amplitudes: np.ndarray (drive amplitudes used in the Rabi sweep)
     rabi_populations: np.ndarray (excited state populations; shape: (number of amplitudes, number of time points))
    """
    fig, ax= plt.subplots(figsize=(10,6))
    for index, amplitude in enumerate(amplitudes):
        ax.plot(time, rabi_populations[index], label=f"A={amplitude:.3f}")

    ax.set_title("Rabi Oscillations")
    ax.set_xlabel("time")
    ax.set_ylabel(r"excited satte population $P_1$")
    ax.set_ylim(0.0, 1.05)
    ax.grid(True, alpha=0.3)

    if len(amplitudes) <= 10:
        ax.legend()

    fig.tight_layout()

    return ax

def plot_transmon_populations(time: np.ndarray, ground_population: np.ndarray, first_excited_population: np.ndarray, second_excited_population: np.ndarray):
    """
    Plots the populations of the three transmon levels.

    Parameters:
     time: np.ndarray (simulation time points)
     ground_population: np.ndarray (population of |0>)
     first_excited_population: np.ndarray (population of |1>)
     second_excited_population: np.ndarray (population of |2>)
    
    Returns:
     matplotlib.axes.Axes (axes contaning the population plot)
    """
    fig, ax= plt.subplots(figsize=(8,5))

    ax.plot(time, ground_population, label=r"$P_0$")
    ax.plot(time, first_excited_population, label=r"$P_1$")
    ax.plot(time, second_excited_population, label=r"$P_2$")

    ax.set_xlabel("time")
    ax.set_ylabel("population")
    ax.set_title("Three-level Transmon State Populations")

    ax.set_ylim(0.0, 1.05)
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    return ax

def plot_comparison(time: np.ndarray, gaussian: tuple[np.ndarray, np.ndarray, np.ndarray], drag: tuple[np.ndarray, np.ndarray, np.ndarray]) -> tuple[plt.axes, plt.axes]:
    """
    Plots a side-by-side population comparison for Gaussian vs DRAG pulses on a 3-level transmon system.

    Parameters:
     time: np.ndarray (Simulation time grid)
     gaussian: tuple (P0, P1, P2 population arrays for the gaussian pulse)
     drag: tuple (P0, P1, P2 population arrays for the DRAG pulse)

    Returns:
     tuple: Matplotlib axes objects for the two subplots
    """
    fig, (ax1, ax2)= plt.subplots(1, 2, figsize=(14,5))

    #Gaussian subplot
    l1= ax1.plot(time, gaussian[0], label=r"$P_0$ (ground)")
    l2= ax1.plot(time, gaussian[1], label=r"$P_1$ (excited)")
    ax1.set_title("3-level Transmon: Gaussian Pulse")
    ax1.set_xlabel("time")
    ax1.set_ylabel("population ($P_0, P_1$)")
    ax1.set_ylim(-0.02, 1.05)
    ax1.grid(True, alpha=0.3)

    #Log secondary axis for Gaussian plot
    ax1_leak= ax1.twinx()
    l3= ax1.plot(time, gaussian[2], label=r"$P_2$ (leakage)", color="red", linestyle="--")
    ax1_leak.set_ylabel("leakage ($P_2$-Log scale)", color="red")
    ax1_leak.set_yscale("log")
    ax1_leak.set_ylim(1e-14, 1.0)
    ax1_leak.tick_params(axis="y", labelcolor="red")

    #combined legend for gaussian
    lines_1= l1 + l2 + l3
    labesl_1= [line.get_label() for line in lines_1]
    ax1.legend(lines_1, labesl_1, loc="center right")

    #DRAG subplot
    l4= ax2.plot(time, drag[0], label=r"$P_0$ (ground)")
    l5= ax2.plot(time, drag[1], label=r"$P_1$ (excited)")
    ax2.set_title("3-level Transmon: DRAG Pulse")
    ax2.set_xlabel("time")
    ax2.set_ylabel("population ($P_0, P_1$)")
    ax2.set_ylim(-0.02, 1.05)
    ax2.grid(True, alpha=0.3)

    #Log secondary axis for DRAG leakage
    ax2_leak= ax2.twinx()
    l6= ax2.plot(time, drag[2], label=r"$P_2$ (leakage)", color="red", linestyle="--")
    ax2_leak.set_ylabel("leakage ($P_2$-Log scale)", color="red")
    ax2_leak.set_yscale("log")
    ax2_leak.set_ylim(1e-14, 1.0)
    ax2_leak.tick_params(axis="y", labelcolor="red")

    #combined legend for DRAG
    lines_2= l4 + l5 + l6
    labels_2= [line.get_label() for line in lines_2]
    ax2.legend(lines_2, labels_2, loc="center right")

    fig.tight_layout()
    return ax1, ax2