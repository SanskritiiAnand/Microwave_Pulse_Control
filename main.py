"""
Entry point for the pulse-level control simulation.
Performs a basic end-to-end simulation of a resonantly driven two-qubit and visualizes the resulting state dynamics
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from src.hamiltonians import total_hamiltonian, total_IQ_hamiltonian
from src.pulse_library import gaussian_pulse, drag_pulse
from src.simulations import simulate_dynamics, state_probabilities, calibrate_drive_amplitude, calibrate_drag_beta, simulate_rabi_oscillation, transmon_state_probabilities, simulate_three_level_gaussian, simulate_three_level_drag
from src.utils import check_normalisation, create_time_grid
from src.visualization import plot_bloch_components, plot_bloch_trajectory, plot_state_populations, plot_pulse_envelope, plot_drag_pulse, plot_rabi_oscillations, plot_transmon_populations, plot_comparison

#Simulation Parameters
qubit_freq= 2.0 * np.pi * 5.0
drive_freq= qubit_freq
drive_amp= 0.1

start_time= 0.0
end_time= 40.0
num_points= 5000

beta= 5.0
pulse_centre= 20.0
pulse_sigma= 4.0

#Main Simulation
def main() -> None:
    """
    Runs the driven-qubit simulation
    """
    #create simulation time grid
    time= create_time_grid(start_time, end_time, num_points)

    #initial state |0> (2-level)
    initial_state= np.array([1.0, 0.0], dtype= complex)

    #Rabi oscillation (2-level)
    rabi_amplitudes= np.array([0.02, 0.04, 0.06, 0.08, 0.10])

    rabi_populations= simulate_rabi_oscillation(amplitudes=rabi_amplitudes, qubit_frequency=qubit_freq, drive_frequency=drive_freq, initial_state=initial_state, time_span=(start_time, end_time), evaluation_times=time)

    #Three-level transmon simulation
    drive_amp= 0.1
    anharmonicity= -2.0 * np.pi * 0.3
    three_level_initial_state= np.array([1.0, 0.0, 0.0], dtype=complex)

    #3-level DRAG beta calibration
    beta_values= np.linspace(-5.0, 5.0, 201)
    amplitudes= np.linspace(0.10, 0.90, 161)

    best_beta, best_amplitude, best_population, beta_results= calibrate_drag_beta(beta_values=beta_values, amplitudes=amplitudes, qubit_frequency=qubit_freq, drive_frequency=drive_freq, center=pulse_centre, sigma=pulse_sigma, initial_state=three_level_initial_state, anharmonicity=anharmonicity, three_level=True, time_span=(start_time, end_time), evaluation_times=time)

    beta= best_beta
    drive_amp= best_amplitude

    #Three-level gaussian simulation
    gaussian_result= simulate_three_level_gaussian(amplitude=drive_amp, qubit_frequency=qubit_freq, drive_frequency=drive_freq, center=pulse_centre, sigma=pulse_sigma, anharmonicity=anharmonicity, initial_state=three_level_initial_state, time_span=(start_time, end_time), evaluation_times=time)
    
    gaussian_state_trajectory= gaussian_result.y
    
    gauss_P0, gauss_P1, gauss_P2= transmon_state_probabilities(gaussian_state_trajectory)

    #Three-level Calibrated DRAG pulse simulation
    drag_result= simulate_three_level_drag(amplitude=drive_amp, qubit_frequency=qubit_freq, drive_frequency=drive_freq, center=pulse_centre, sigma=pulse_sigma, beta=beta, anharmonicity=anharmonicity, initial_state=three_level_initial_state, time_span=(start_time, end_time), evaluation_times=time)

    drag_state_trajectory= drag_result.y

    drag_P0, drag_P1, drag_P2= transmon_state_probabilities(drag_state_trajectory)

    #2-level simulation with calibrated parameters
    def hamiltonian(time_point):
        drag= drag_pulse(time_point, amplitude=drive_amp, center=pulse_centre, sigma=pulse_sigma, beta=beta)
        I= np.real(drag)
        Q= np.imag(drag)

        return total_IQ_hamiltonian(time= time_point, qubit_frequency=qubit_freq, drive_freq=drive_freq, I=I, Q=Q)

    #run simulation
    result= simulate_dynamics(hamiltonian= hamiltonian, initial_state= initial_state, time_span=(start_time, end_time), evaluation_times= time)

    #extract state trajectory
    state_trajectory= result.y

    #calculate state populations
    ground_population, excited_population= state_probabilities(state_trajectory)

    #verify normalisation
    normalized = check_normalisation(state_trajectory, tolerance= 1e-7)

    norms = np.linalg.norm(state_trajectory, axis=0)

    #Results
    print(f"2-level Simulation successful: {result.success}")
    print(f"State remains normalised: {normalized}")
    print(f"Minimum state norm: {np.min(norms):.12f}")
    print(f"Maximum state norm: {np.max(norms):.12f}")
    print(f"Maximum norm error: {np.max(np.abs(norms - 1.0)):.3e}\n")

    print(f"--- 3-Level DRAG Parameter Calibration ---")
    print(f"Best beta: {best_beta:.3f}")
    print(f"Best drive amplitude: {best_amplitude:.6f}")
    print(f"Calibrated 3-level P1 population: {best_population:.6f}")
    print(f"Uncorrected Gaussian P2 Leakage: {gauss_P2[-1]:.6e}")
    print(f"DRAG Corrected P2 Leakage: {drag_P2[-1]:.6e}")


    #visualisations
    pulse_envelope = gaussian_pulse(time, amplitude=drive_amp, center=pulse_centre, sigma=pulse_sigma)
    drag_signals = drag_pulse(time=time, amplitude=drive_amp, center=pulse_centre, sigma=pulse_sigma, beta=beta)

    # Group population arrays into tuples expected by plot_comparison
    gaussian_tuple = (gauss_P0, gauss_P1, gauss_P2)
    drag_tuple = (drag_P0, drag_P1, drag_P2)

    plot_transmon_populations(time, drag_P0, drag_P1, drag_P2)
    plot_state_populations(time, ground_population, excited_population)
    plot_bloch_components(time, state_trajectory)
    plot_bloch_trajectory(state_trajectory)
    plot_rabi_oscillations(time, rabi_amplitudes, rabi_populations)
    plot_pulse_envelope(time, pulse_envelope, title="Gaussian microwave pulse envelope")
    plot_drag_pulse(time, drag_signals)
    plot_comparison(time, gaussian_tuple, drag_tuple)
    

    plt.show()

if __name__ == "__main__":
    main() 