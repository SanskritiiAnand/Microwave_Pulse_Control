"""
utility functions for generating microwave pulse envelopes used in pulse-level simulations of superconducting qubits
"""
from __future__ import annotations
import numpy as np

#Constant(square) Pulse
def constant_pulse(time: float | np.ndarray, amplitude: float) -> float | np.ndarray:
    """
    generates a constant (square) pulse

    Parameters:
    time: float or np.ndarray (time coordinate(s))

    amplitude: float (pulse amplitude)

    Returns:
    float or np.ndarray (constant pulse envelope)
    """
    if np.isscalar(time):
        return amplitude

    return amplitude * np.ones_like(time, dtype=float)

#Gaussian Pulse
def gaussian_pulse(time: float | np.ndarray, amplitude: float, center: float, sigma: float) -> float | np.ndarray:
    """
    generates a gaussian pulse envelope

    Parameters:
    time: float or np.ndarray (time corrdinates)
    amplitude: float (peak amplitude)
    center: float (pulse center)
    sigma: float (standard deviation of the gaussian)

    Returns:
    float or np.ndarray (gaussian pulse envelope)
    """
    return amplitude * np.exp(-((time - center)** 2)/(2 * sigma ** 2))

#DRAG Pulse
def drag_pulse(time: float | np.ndarray, amplitude: float, center: float, sigma: float, beta: float) -> float | np.ndarray:
    """
    Generates a DRAG (Derivative Removal by Adiabatic Gate) pulse

    Returned envelope is complex-valued:
    I(t) + i Q(t)
    where, I(t) = Gaussian pulse
           Q(t) = -β d/dt Gaussian

    Parameters:
    time: float or np.ndarray (time coordinates)
    amplitude: float (peak gaussian amplitude)
    center: float (pulse center)
    sigma: float (gaussian width)
    beta: float (DRAG coefficient)

    Returns:
    float or np.ndarray (complex-valued DRAG pulse envelope)
    """
    gaussian= gaussian_pulse(time, amplitude, center, sigma)

    derivative= -((time - center)/ sigma ** 2) * gaussian

    return gaussian - 1j * beta * derivative
