"""
General-purpose utility functions for pulse-level qiskit-control simulations
"""
from __future__ import annotations
import numpy as np

#Time Grid
def create_time_grid(start_time: float, end_time: float, num_points: int) -> np.ndarray:
    """
    creates an evenly spaced simulation time grid

    Parameters:
     start_time: float (starting time of the simulation)
     end_time: float (ending time of the simulation)
     num_points: int (number of time points)

    Returns:
     np.ndarray (evenly spaced time points)
    """
    if end_time <= start_time:
        raise ValueError("\nend time must be grater that start time")

    if num_points < 2:
        raise ValueError("\nnum_points must be at least 2")

    return np.linspace(start_time, end_time, num_points)

#State Normalization
def normalise_state(state: np.ndarray) -> np.ndarray:
    """
    normalise a quantum state vector

    Parameters:
     state: np.ndarray (quantum state vector)

    Returns:
     np.ndarray (normalised quantum state vector)
    """
    state= np.asarray(state, dtype=complex)

    norm= np.linalg.norm(state)

    if np.isclose(norm, 0.0):
        raise ValueError("\n Cannot normalize a zero state vector")

    return state / norm

#State Normalisation Check
def check_normalisation(state_trajectory: np.ndarray, tolerance: float= 1e-8) -> bool:
    """
    Checks whether a quantum state trajectory remains normalised

    Parameters:
     state_trajectory: np.ndarray (state trajectory with shape (dimension, N), where N is the number of time points)
     tolerance: float, optional (maximum allowed deviation from unit norm)

    Returns:
     bool (true if every state in the trajectory is normalised within the specified tolerance)
    """
    norms= np.linalg.norm(state_trajectory, axis=0)

    return bool(np.allclose(norms, 1.0, atol=tolerance, rtol=0.0))