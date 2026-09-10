"""
Numerical simulation utilities for pulse-level quantum control.
Provides functions for solving the time-dependent Schrodinger equation and extracting state populations from the resulting quantum state trajectory
"""
from __future__ import annotations 
from collections.abc import Callable
from scipy.integrate import solve_ivp
import numpy as np
from src.pulse_library import drag_pulse, gaussian_pulse
from src.hamiltonians import total_IQ_hamiltonian, total_transmon_IQ_hamiltonian

HamiltonianFunc = Callable[[float], np.ndarray]

#Schrodinger Equation
def schrodinger_rhs(time: float, state: np.ndarray, hamiltonian: HamiltonianFunc, hbar: float= 1.0) -> np.ndarray:
   """
   Computes the right-hand side of the time dependent Schrodinger equation:
   i hbar d|psi(t)>/dt = H(t)|psi(t)>
   rearranegd as:
    d|psi(t)>/dt = -i/hbar H(t)|psi(t)>.

   Parameters:
    time: float (current simulaion time)
    state: np.ndarray (current quantum state vector)
    hamiltonian: callable (function that takes time as input and returns the Hamiltonian metrix at that time)
    hbar: float, optional (reduced planck's constant defaults to 1.0, corresponding to natural units)

   Return:
    np.ndarray (time derivative of the quantum state)
   """
   H = hamiltonian(time)

   return -1j / hbar * H @ state

#Quantum State Evolution
def simulate_dynamics(hamiltonian: HamiltonianFunc, initial_state: np.ndarray, time_span: tuple[float, float], evaluation_times: np.ndarray, hbar: float= 1.0, method: str = "DOP853", rtol: float= 1e-9, atol: float= 1e-11):
   """
   simulates quantum state evolution under a time-dependent Hamiltonian

   Parameters:
    hamiltonian: callable (function that takes time as input and returns the Hamiltonian metrix at that time)
    initial_state: np.ndarray (initial quantum state vector)
    time_span: tuple[float, float] (start and end times of the simulation)
    evaluation_times: np.ndarray (times at whoch the state should be returned)
    hbar: float, optional (reduced planck's constant, defaults to 1.0)
    method: str, optional (numerical integration method useed by Scipy. defaults to 'DOP853')
    rtol: float, optional (relative tolerance for the numercial solver)
    atol: float, optional (absolute tolerance for the numerical solver)

   Returns:
    scipy.integrate.OdeResult (solution object containing the simulated state trajectory)
   """
   initial_state= np.asarray(initial_state, dtype= complex)
   evaluation_times= np.asarray(evaluation_times, dtype= float)

   result= solve_ivp(fun= lambda time, state: schrodinger_rhs(time, state, hamiltonian, hbar),
                     t_span= time_span,
                     y0= initial_state,
                     t_eval= evaluation_times,
                     method= method,
                     rtol= rtol,
                     atol= atol
                    )
   if not result.success:
      raise RuntimeError(f"Simulation Failed: {result.message}")

   return result

#State Probabilities
def state_probabilities(state_trajectory: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
   """
   Calculates computational-basis state populations.

   For a single-qubit state:
    |psi(t)> = alpha(t)|0> + beta(t)|1>, the populations are:
       P0(t) = |alpha(t)|^2
       P1(t) = |beta(t)|^2     

   Parameters:
    state_trajectory: np.ndarray (array containing the simulated state vectors. expected shape is (2, N), where N: number of time points)

   Returns:
    tuple[np.ndarray, np.ndarray] (ground state and excited state populations)
   """
   gnd_population= np.abs(state_trajectory[0]) ** 2
   exctd_population= np.abs(state_trajectory[1]) ** 2

   return gnd_population, exctd_population

#Three-level Transmon State Probabilities
def transmon_state_probabilities(state_trajectory: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
   """
   Calculates the populations of the three transmon levels.
   For a three-level state:
    |psi(t)> = c0(t)|0> + c1(t)|1> + c2(t)|2>
   the populations are:
    P0(t)= |c0(t|^2
    P1(t)= |c1(t)|^2
    P2(t)= |c2(t)|^2

   Parameters:
    state_trajectory: np.ndarray (simulated state trajectory with shape (3,N))

   Returns:
    tuple[np.ndarray, np.ndarray, np.ndarray] (ground, first excited, and second excited state populations)
   """
   ground_population= np.abs(state_trajectory[0])**2
   first_excited_population= np.abs(state_trajectory[1])**2
   second_excited_population= np.abs(state_trajectory[2])**2

   return(ground_population, first_excited_population, second_excited_population)

#three-level gaussian pulse simulation
def simulate_three_level_gaussian(amplitude: float, qubit_frequency: float, drive_frequency: float, center: float, sigma: float, anharmonicity: float, initial_state: np.ndarray, time_span: tuple, evaluation_times: np.ndarray):
   """
   Simulates a three-level transmon driven by a gaussian microwave pulse 
   The pulse is applied through the in-phase quadrature:
    I(t)= gaussian(t)
    Q(t)= 0
   The three transmon levels |0>, |1>, |2> are included, allowing population leakage into |2>

   Parameters:
    amplitude: float (peak amplitude of the gaussian pulse)
    qubit_frequency: float (|0> -> |1> tranisition frequency)
    drive_frequency: float (microwave drive frequency)
    center: float (center of the gaussian pulse)
    anharmonicity: float (transmon anharmonicity)
    initial_state: np.ndarray (initial three level state vector)
    time_span: tuple (start and end times of the simulation)
    evaluation_times: np.ndarray (times at which the state is evaluated)

   Returns:
    scipy.integrate.OdeResult (simulated three-level state trajectory)
   """
   def hamiltonian(time_point):
      I= gaussian_pulse(time_point, amplitude=amplitude, center=center, sigma=sigma)
      Q= 0.0

      return total_transmon_IQ_hamiltonian(time=time_point, qubit_frequency=qubit_frequency, drive_freq=drive_frequency, I=I, Q=Q, anharmonicity=anharmonicity)

   return simulate_dynamics(hamiltonian=hamiltonian, initial_state=initial_state, time_span=time_span, evaluation_times=evaluation_times)

#Three-level DRAG Pulse Simulation
def simulate_three_level_drag(amplitude: float, qubit_frequency: float, drive_frequency: float, center: float, sigma: float, beta: float, anharmonicity: float, initial_state: np.ndarray, time_span: tuple, evaluation_times: np.ndarray):
   """
   Simulates a three-level transmon driven by a complex DRAG microwave pulse.
   I(t) controls in-phase drive
   Q(t) controls derivative quadrature drive

   Parameters:
    amplitude: float (peak amplitude of the DRAG pulse)
    qubit_frequency: float (|0> -> |1> tranisition frequency)
    drive_frequency: float (microwave drive frequency)
    center: float (center of the DRAG pulse)
    sigma: float (gaussian pulse width)
    beta: float (DRAG correction coefficient)
    anharmonicity: float (transmon anharmonicity)
    initial_state: np.ndarray (initial three level state vector)
    time_span: tuple (start and end times of the simulation)
    evaluation_times: np.ndarray (times at which the state is evaluated)
   
   Returns:
    scipy.integrate.OdeResult (Simulated three level state trajectory)
   """
   def hamiltonian(time_point):
      drag= drag_pulse(time_point, amplitude=amplitude, center=center, sigma=sigma, beta=beta)
      I= np.real(drag)
      Q= np.imag(drag)

      return total_transmon_IQ_hamiltonian(time=time_point, qubit_frequency=qubit_frequency, drive_freq=drive_frequency, I=I, Q=Q, anharmonicity=anharmonicity)

   return simulate_dynamics(hamiltonian=hamiltonian, initial_state=initial_state, time_span=time_span, evaluation_times=evaluation_times)

#pi-Pulse Amplitude Calibration
def calibrate_drive_amplitude(amplitudes: np.ndarray, qubit_frequency: float, drive_frequency: float, center: float, sigma: float, beta: float, initial_state: np.ndarray, time_span: tuple, evaluation_times: np.ndarray, anharmonicity: float=None, three_level: bool=False):
   """
   finds the drive amplitude that maximises the final excited-state population P1
   
   Parameters:
    amplitudes: array-like (drive amplitudes to test)
    qubit_frequency: float (qubit transition frequency)
    drive_frequency: float (microwave drive frequency)
    center: float (pulse center)
    sigma: float (gaussian pulse width)
    beta: float (DRAG correction coefficient)
    initial_state: float (initial qubit state)
    time_span: tuple (simulation start and end times)
    evaluation_times: np.ndarray (times at which the solutions is evaluated)
    anharmonicity: float, optional (transmon anharmonicity requireed if three_level=True)
    three_level: bool, optional (3-level transmon system)

   Returns:
    best_amplitude: float (amplitude giving the highest final P1)
    best_population: float (corresponding final excited-state population)
    populations: np.ndarray (final P1 for every tested amplitude)
   """
   populations= []

   for amplitude in amplitudes:
      if three_level:
         if anharmonicity is None:
            raise ValueError("Anharmonicity parameter must be provided for 3-level simulations")

         result= simulate_three_level_drag(amplitude=amplitude, qubit_frequency=qubit_frequency, drive_frequency=drive_frequency, center=center, sigma=sigma, beta=beta, anharmonicity=anharmonicity, initial_state=initial_state, time_span=time_span, evaluation_times=evaluation_times)

         _, excited_populations, _= transmon_state_probabilities(result.y)

      else:
         def hamiltonian(time_point):
            drag= drag_pulse(time_point, amplitude= amplitude, center= center, sigma= sigma, beta= beta)

            I= np.real(drag)
            Q= np.imag(drag)

            return total_IQ_hamiltonian(time= time_point, qubit_frequency= qubit_frequency, drive_freq= drive_frequency, I= I, Q= Q)

         result= simulate_dynamics(hamiltonian= hamiltonian, initial_state= initial_state, time_span= time_span, evaluation_times= evaluation_times)
         state_trajectory= result.y

         _, excited_populations= state_probabilities(state_trajectory)

      populations.append(excited_populations[-1])

   populations= np.array(populations)

   best_index= np.argmax(populations)
   best_amplitude= amplitudes[best_index]
   best_population= populations[best_index]

   return best_amplitude, best_population, populations

#beta calibration
def calibrate_drag_beta(beta_values: np.ndarray, amplitudes: np.ndarray, qubit_frequency: float, drive_frequency: float, center: float, sigma: float, initial_state: np.ndarray, time_span: tuple, evaluation_times: np.ndarray, anharmonicity: float=None, three_level: bool=False):
   """
   Calibrates the DRAG coefficient beta by first calibrating amplitude at beta=0, and then selecting beta to minimise final P2 leakage.

   Parameters: 
    beta_value: array (DRAG beta values to test)
    amplitudes: array (drive amplitudes used for amplitude calibration at each beta)
    qubit_frequency: float (qubit transition frequency)
    drive_freq: float(microwave drive frequency)
    center: float (pulse center)
    sigma: float (gaussian pulse width)
    intial_state: np.ndarray (initial qubit state)
    time_span: tuple (start and end times of the simulation)
    evaluation_times: np.ndarray (times at which the state is evaluated)
    anharmonicity: float, optional (transmon anharmonicity)
    three_level: bool, optional (3-level transmon system)
   Returns:
    best_beta: float (beta giving the highest final excited-state population)
    best_amplitude: float (best drive amplitude corresponding to best_beta)
    best_population: float (highest final excited-state population)
    bets_results: list (calibration results for every beta value)
   """
   #Evaluate P2 leakage
   beta_results= []

   for beta in beta_values:
      best_amplitude, _, _= calibrate_drive_amplitude(amplitudes=amplitudes, qubit_frequency=qubit_frequency, drive_frequency=drive_frequency, center=center, sigma=sigma, beta=beta, initial_state=initial_state, time_span=time_span, evaluation_times=evaluation_times, anharmonicity=anharmonicity, three_level=three_level)
      
      if three_level:
         result= simulate_three_level_drag(amplitude=best_amplitude, qubit_frequency=qubit_frequency, drive_frequency=drive_frequency, center=center, sigma=sigma, beta=beta, anharmonicity=anharmonicity, initial_state=initial_state, time_span=time_span, evaluation_times=evaluation_times)

         _, P1, P2= transmon_state_probabilities(result.y)

         P1_final= P1[-1]
         P2_final= P2[-1]

      else:
         #fallback for 2-level
         def hamiltonian(time_point):
            drag= drag_pulse(time_point, amplitude=best_amplitude, center=center, sigma=sigma, beta=beta)
            I= np.real(drag)
            Q= np.imag(drag)
            return total_IQ_hamiltonian(time=time_point, qubit_frequency=qubit_frequency, drive_freq=drive_frequency, I=I, Q=Q)

         result = simulate_dynamics(hamiltonian=hamiltonian, initial_state=initial_state, time_span=time_span, evaluation_times=evaluation_times)
         _, P1 = state_probabilities(result.y)
         P1_final = P1[-1]
         P2_final = 0.0

      beta_results.append({"beta": beta, "amplitude": best_amplitude, "population": P1_final, "P2_leakage": P2_final})

   #Select beta that minimizes P2 leakage
   if three_level:
        best_result = min(beta_results, key=lambda result: result["P2_leakage"])
   else:
        best_result = max(beta_results, key=lambda result: result["population"])

   return (best_result["beta"], best_result["amplitude"], best_result["population"], beta_results)

#Rabi Oscillation Simulation
def simulate_rabi_oscillation(amplitudes: np.ndarray, qubit_frequency: float, drive_frequency: float, initial_state: np.ndarray, time_span: tuple, evaluation_times: np.ndarray):
   """
   Simulates Rabi Oscillations for a range of microwave drive amplitudes.
   For each amplitude a constant IQ drive is applied and the resulting excited state population is recorded as a function of time

   Parameters:
    amplitudes: np.ndarray (drive amplitudes to simulate)
    qubit_frequency: float (qubit transition frequency in rad/s)
    drive_frequency: float (microwave drive frequency in rad/s)
    initial_state: np.ndarray ()
    time_span: tuple (start and end times of the simulation)
    evaluation_times: np.ndarray (times at which the state is evaluated)

   Returns:
    np.ndarray (excited-state populations for every amplitude; shape: (number of amplitudes, numbe of time points))
   """
   rabi_populations= []
   for amplitude in amplitudes:
      def hamiltonian(time_point):
         return total_IQ_hamiltonian(time=time_point, qubit_frequency=qubit_frequency, drive_freq=drive_frequency, I=amplitude, Q=0.0)

      result= simulate_dynamics(hamiltonian=hamiltonian, initial_state=initial_state, time_span=time_span, evaluation_times=evaluation_times)
      _, excited_population= state_probabilities(result.y)

      rabi_populations.append(excited_population)

   return np.array(rabi_populations)