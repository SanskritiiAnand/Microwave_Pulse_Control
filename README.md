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
* $`H_0 = \frac{\hbar \omega_q}{2} \sigma_z = \frac{\hbar \omega_q}{2} \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}`$

where $\omega_q$ denotes the funadamental qubit transition frequency. The interaction with a real-valued classical microwave drive field $\Omega(t) = A(t) \cos(\omega_d t + \phi)$ is governed by the time-dependent drive Hamiltonian:
* $$H_{\text{drive}}(t) = \hbar \Omega(t) \sigma_x$$

Transforming to the rotating frame at the drive frequency $\omega_d$ via the unitary operator $U(t) = \exp\left(i \frac{\omega_d t}{2} \sigma_z\right)$ and applying the Rotating Wave Approximation (RWA)- neglecting fast-oscillating terms at $2\omega_d$- yields the time-independent detuned frame HAmiltonian:
* $$H_{\text{rot}} = -\frac{\hbar \Delta}{2} \sigma_z + \frac{\hbar}{2} \left[ I(t)\sigma_x + Q(t)\sigma_y \right]$$
where $\Delta = \omega_q - \omega_d$ represents the drive detuning, while I(t) and Q(t) are the in-phase and quadrature envelope functions respectively.

### Three-Level Transmon System (qutrit)
A physical transmon qubit consists of a Josephson junction shunted by a large capacitor. The weak anharmonicity $\alpha$ results in a non-equidistant energy spectrum, requiring an expansion into the qutrit subspace {|0>, |1>, |2>} to model control-induced leakage.
Using bosonic creation ($a^\dagger$) and annihilation ($a$) operators truncated to d=3 dimensions, the full transmon Hamiltonian in the laboratory frame is given by the Duffing oscillator model:
* $$H_{\text{transmon}} = \hbar \omega_q a^\dagger a + \frac{\hbar \alpha}{2} a^\dagger a^\dagger a a$$

Transforming to the frame rotating with the drive frequency $\omega_d$ via $U(t) = \exp\left(i \omega_d t a^\dagger a\right)$, the static free Hamiltonian becomes:
* $` H_0^{\text{rot}} = \hbar \Delta a^\dagger a + \frac{\hbar \alpha}{2} a^\dagger a^\dagger a (a - 1) = \hbar \begin{pmatrix} 0 & 0 & 0 \\ 0 & \Delta & 0 \\ 0 & 0 & 2\Delta + \alpha \end{pmatrix} `$

Under strict resonance ($\omega_d = \omega_q \implies \Delta = 0$), the static Hamiltonian simplifies to $H_0^{\text{rot}} = \text{diag}(0, 0, \hbar \alpha)$.
#### Drive Coupling Operators
The interaction Hamiltonian in the rotating frame takes the form:
$$H_{\text{drive}}^{\text{rot}}(t) = \frac{\hbar}{2} \left[ I(t) X + Q(t) Y \right]$$

The dimensionless dipole coupling operators X and Y are generalised Pauli matrices over the 3-level Hilbert space, derived from $X = a + a^\dagger$ and $Y = -i(a - a^\dagger)$:

$`X = \begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & \sqrt{2} \\ 0 & \sqrt{2} & 0 \end{pmatrix}, \quad Y = \begin{pmatrix} 0 & -i & 0 \\ i & 0 & -i\sqrt{2} \\ 0 & i\sqrt{2} & 0 \end{pmatrix}`$

The matrix element factor  $\sqrt{2} \approx 1.414$ explicitly scales the dipole transition strength for the |1> <-> |2> channel relative to the fundamental |0> <-> |1> transition, establishing the physical pathway for non-computational state leakage.
