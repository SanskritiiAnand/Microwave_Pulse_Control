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
 
Furthermore, an automated two-dimensional parameter calibration pipeline is implemented to jointly optimise the drive amplitude (A) and DRAG coefficient ($\beta$). The overarching goal of this work is to demonstrate high-fidelity |0> -> |1> state inversion ($\mathcal{F} > 99.95\%$) while suppressing |2> population leakage to near-zero levels (~10^(-11)).

## Physical Model
### Two-Level System (qubit)
The baseline dynamics consider an ideal two-level quantum system spanned byt he computational basis states {|0>, |1>}. In the laboratory frame, the bare system Hamiltonian H0 is expressed via the Pauli-Z operator $\sigma_z$:
*$`H_0 = \frac{\hbar \omega_q}{2} \sigma_z = \frac{\hbar \omega_q}{2} \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}`$

where $\omega_q$ denotes the funadamental qubit transition frequency. The interaction with a real-valued classical microwave drive field $\Omega(t) = A(t) \cos(\omega_d t + \phi)$ is governed by the time-dependent drive Hamiltonian:
* $$H_{\text{drive}}(t) = \hbar \Omega(t) \sigma_x$$

Transforming to the rotating frame at the drive frequency $\omega_d$ via the unitary operator $U(t) = \exp\left(i \frac{\omega_d t}{2} \sigma_z\right)$ and applying the Rotating Wave Approximation (RWA)- neglecting fast-oscillating terms at $2\omega_d$- yields the time-independent detuned frame HAmiltonian:
* $$H_{\text{rot}} = -\frac{\hbar \Delta}{2} \sigma_z + \frac{\hbar}{2} \left[ I(t)\sigma_x + Q(t)\sigma_y \right]$$
where $\Delta = \omega_q - \omega_d$ represents the drive detuning, while I(t) and Q(t) are the in-phase and quadrature envelope functions respectively.
