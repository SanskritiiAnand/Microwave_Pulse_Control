# Pulse-Level Control of a Superconducting Transmon
## Introduction 
Superconducting circuit architectures serve as one of the leading physical platforms for scalable quantum information processing. While higher-level quantum compilers abstract operations into discrete unitaries, such as Rx, Ry, or two-qubit entangling gates, these operations are physically realised via precisely modulated, time-dependent microwave pulses applied to non-linear microwave resonators.

This repository presents a numerical model for pulse-level quantum control, modeling the open and closed-system dynamics of a driven transmon system directly from first principles. Built on NumPy/Scipy backend, this solver models time-dependent Schrödinger dynamics without relying on high-level abstractions like Qiskit Dynamics.

The project investigates two fundamental physical regimes:
* 1. Two-level Qubit Approximation- Establishing baseline Rabi oscillations, unitary state trajectory evolution on the Bloch sphere, and numerical solver convergence under resonant microwave driving.
* 2. Three-level Transmon System- Incorporating weakly anharmonic energy spectra (alpha < 0) to capture non-computational state leakage into the |2> subspace induced by spectral overlap from short control pulses.

To mitigate non-adiabatic transitions, two pulse-shaping paradigms are modeled and benchmarked:
* 1. Standard Gaussian Envelope- Real-valued, amplitude modulated microwave signals.
* 2. Derivative Removal by Adiabatic Gate (DRAG) Schemes- Phase-shifted derivative quadrature $Q(t) \propto -\beta \frac{d}{dt}I(t)$ control signals designed to suppress spectral density at the |1> -> |2> transition frequency.
 
Furthermore, an automated two-dimensional parameter calibration pipeline is implemented to jointly optimise the drive amplitude (A) and DRAG coefficient ($\beta$). The overarching goal of this work is to demonstrate high-fidelity |0> -> |1> state inversion ($\mathcal{F} > 99.95\%$) while suppressing |2> population leakage to near-zero levels (~$\sim 10^{-11}$).
