"""
Utility functions for constructing Hamiltonians used in pulse-level simulations of superconducting qubits.
"""
from  __future__ import annotations
import numpy as np

#Pauli Operators
sigma_x= np.array([[0,1], [1,0]], dtype=complex)

sigma_y= np.array([[0,-1j], [1j,0]], dtype=complex)

sigma_z= np.array([[1,0], [0,-1]], dtype=complex)

Identity= np.eye(2, dtype=complex)

#Static Hamiltonian
def static_qubit_hamiltonian(qubit_frequency: float, hbar: float= 1.0) -> np.ndarray:
    """
    Constructs the static hamiltonian of a two-level superconducting qubit.
    
    The static Hamiltonian is:
    H0 = (ħ ωq / 2) σz

    Parameters:
    qubit frequency: float (transition frequency of the qubit in rad/s)
    hbar: float, optional (reduced Planck's constant, defaults to 1.0)

    Returns: np.ndarray, 2 x 2 static hamiltonian
    """
    return 0.5 * hbar * qubit_frequency * sigma_z

#Drive Hamiltonian
def drive_hamiltonian(time: float, drive_envelope, drive_freq: float, phase: float= 0.0, hbar: float= 1.0) -> np.ndarray:
    """
    Constructs the time dependent microwave drive hamiltonian
    
    The drive hamiltonian is:
    H_d(t) = ħ A(t) cos(ω_d t + φ) σ_x

    Parameters:
    time: float (simulation time)
    drive_envelope: callable (function returning the pulse envelope A(t))
    drive_freq: float (microwave drive frequency in rad/s)
    phase: float, optional (initial phase of the microwave pulse)
    hbar: float, optional (reeduced Planck's constant)

    Returns:
    np.ndarray ( 2 x 2 drive hamiltonian evaluated at time t)
    """
    amplitude= drive_envelope(time)

    carrier= np.cos(drive_freq * time + phase)

    return hbar * amplitude * carrier * sigma_x

#I/Q representation
def iq_drive_hamiltonian(I: float, Q: float, hbar: float = 1.0) -> np.ndarray:
    """
    Constructs the hamiltonian of an IQ-modulated microwave drive
    The two quadratures independently control rotations about the X and Y axes:
      H_d(t) = ħ [I(t) σx + Q(t) σy]
    
    Parameters:
      I: float (in-phase quadrature amplitude)
      Q: float (quadrature amplitude)
      hbar: float, optional (reduced planck's constant, defaults to 1.0)

    Returns:
      np.ndarray (2 x 2 IQ-modulated drive hamiltonian)
    """
    return hbar * (I*sigma_x + Q*sigma_y)

#Total I/Q hamiltonian
def total_IQ_hamiltonian(time: float, qubit_frequency: float, I: float, Q: float, drive_freq: float, hbar: float = 1.0) -> np.ndarray:
    """
    Constructs the total rotating-frame Hamiltonian for an IQ-controlled qubit
     H(t) = H0 + Hd(t)
    where, H0 = (ħ Δ / 2) σz
           Hd(t) = ħ [I(t) σx + Q(t) σy]
           and, Δ = ωq - ωd

    
    Parameters:
     time: float (simulation time)
     qubit_frequency: float (qubit transition frequnecy in rad/s)
     I: float (in-phase pulse envelope)
     Q: float (quadrature pulse envelope)
     drive_freq: float (microwave carrier frequency in rad/s)
     hbar: float, optional (reduced planck's constant, defaults to 1.0)

    Returns:
     np.ndarray (2 x 2 total hamiltonian)
    """
    detuning= qubit_frequency - drive_freq

    H0= 0.5 * hbar * detuning * sigma_z

    Hd= iq_drive_hamiltonian(I= I, Q= Q, hbar=hbar)

    return H0 + Hd

#Total Hamiltonian
def total_hamiltonian(time: float, qubit_frequency: float, drive_envelope, drive_freq: float, phase: float=0.0, hbar: float=1.0) -> np.ndarray: 
    """
    Constructs the complete driven qubit hamiltonian

    H(t) = H0 + Hd(t)

    Parameters: 
    time: float (simulation time)
    qubit_frequency: float (transition frequency of the qubit)
    drive_envelope: callable (function returning the pulse envelope A(t))
    drive_freq: float (microwave drive frequency)
    phase: float, optional (microwave phase)
    hbar: float, optional (reduced planck's constant)

    Returns:
    np.ndarray (total hamiltonian evaluated at time t)
    """
    return (static_qubit_hamiltonian(qubit_frequency, hbar)
            + drive_hamiltonian(time, drive_envelope, drive_freq, phase, hbar)
           )

#Three level transmon system
transmon_x= np.array([[0.0, 1.0, 0.0],
                     [1.0, 0.0, np.sqrt(2.0)],
                     [0.0, np.sqrt(2.0), 0.0]],
                     dtype=complex
                    )
transmon_y= np.array([[0.0, -1j, 0.0],
                     [1j, 0.0, -1j*np.sqrt(2.0)],
                     [0.0, 1j*np.sqrt(2.0), 0.0]],
                     dtype=complex
                    )

def transmon_static_hamiltonian(qubit_frequency: float, drive_freq: float, anharmonicity: float, hbar: float=1.0) -> np.ndarray:
  """
  constructs the static hamiltonian of a three-level transmon in the rotating frame. (|0>, |1>, |2>)
  In the rotating frame:
   H0= hbar * diag(0, delta, 2*delta + alpha)
   where, delta= w_01 - w_d
          and alpha is the transmon anharmonicity.

  Parameters:
   qubit_frequency: float (|0> -> |1> transition frquency in rad/s)
   drive_freq: float (microwave drive frequency in rad/s)
   anharmonicity: float (transmon anharmonicity in rad/s, typically negative)
   hbar: float, optional (reduced Planck's constant)

  Returns:
   np.ndarray (3x3 static hamiltonian)
  """
  detuning= qubit_frequency - drive_freq

  return hbar * np.diag([0.0, detuning, 2.0*detuning + anharmonicity]).astype(complex)

def transmon_IQ_hamiltonian(I: float, Q: float, hbar: float=1.0) -> np.ndarray:
    """
    Constructs the IQ-modulated drive hamiltonian for a 3-level transmon.
    The drive couples adjacent transmon levels: |0> <-> |1>
                                                |1> <-> |2>
    with the |1> <-> |2> coupling enhanced by sqrt(2).
    The hamiltonian is:
     Hd(t)= hbar [I(t)X + Q(t)Y]
    where X, Y: truncated transmon operators.

    Parameters:
     I: float (in-phase quadrature amplitude)
     Q: float (quadrature amplitude)
     hbar: float, optional (reduced Planck's constant)

    Returns:
     np.ndarray (3x3 IQ drive hamiltonian)
    """
    return hbar * (I * transmon_x + Q * transmon_y)

def total_transmon_IQ_hamiltonian(time: float, qubit_frequency: float, drive_freq: float, I: float, Q: float, anharmonicity: float, hbar: float=1.0) -> np.ndarray:
    """
    constructs the complete rotating-frame hamiltonian for a three-level transmon under IQ control
     H(t)= H0 + Hd(t)
    
     Parameters:
      time: float(current simulation time)
      qubit_frequency: float (|0> -> |1> transition frequency in rad/s)
      drive_freq: float (microwave drive frequency in rad/s)
      I: float (In-phase pulse amplitude)
      Q: float (quadrature pulse amplitude)
      anharmonicity: float (Transmon anharmonicity in rad/s, typically negative)
      hbar: float, optional (reduced Planck's constant)
    
     Returns:
      np.ndarray (3x3 total hamiltonian)
    """
    H0= transmon_static_hamiltonian(qubit_frequency=qubit_frequency, drive_freq=drive_freq, anharmonicity=anharmonicity, hbar=hbar)

    Hd= transmon_IQ_hamiltonian(I=I, Q=Q, hbar=hbar)

    return H0 + Hd