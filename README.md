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

## Pulse Synthesis & Control Protocols
### Constant Envelope (Square Drive)
For baseline characterisation of Rabi dynamics, a step-function drive of constant amplitude A is applied over the interaction window $t \in [t_{\text{start}}, t_{\text{end}}]$:
$$\Omega(t) = A$$

Under resonant driving ($\Delta = 0$), a square envelope generates deterministic Rabi oscillations between states |0> and |1> with an angular frequency $\Omega_{\text{Rabi}} = A$. While analytically simple, square pulses exhibit broad spectral tails in the frequency domain ($\text{sinc}(\omega)$ profiles), inducing severe non-computational excitation in multi-level systems with small anharmonicities.

### Gaussian Pulse Envelope
To suppress high-frequency spectral sidebands, the control amplitude is modulated using a smooth Gaussian envelope G(t):
$$I(t) = G(t) = A \exp\left[ -\frac{(t - t_0)^2}{2\sigma^2} \right]$$
where A represents the peak signal amplitude, $t_0$ is the temporal midpoint of the control frame, and $\sigma$  parameterizes the variance of the pulse.
#### Spectral Overlap & Leakage Mechanism
Because the Fourier transform of a Gaussian pulse in time is itself a Gaussian in frequency:
$$\tilde{G}(\omega) = A \sigma \sqrt{2\pi} \exp\left[ -\frac{\sigma^2 (\omega - \omega_d)^2}{2} \right]$$
short pulse durations (small $\sigma$) lead to significant spectral broadening. If the tail of $\tilde{G}(\omega)$ at the |1> -> |2> transition frequency ($\omega_{12} = \omega_q + \alpha$) contains non-zero spectral energy, off-resonant driving induces state leakage into the |2> manifold.

### Derivative Removal by Adiabatic Gate (DRAG) Framework
To eliminate leakage without increasing total gate duration, the DRAG protocol applies a phase-shifted quadrature correction Q(t) proportional to the time derivative of the primary envelope I(t).
The complex-valued baseband control signal $\Omega(t)$ is defined as:
$$\Omega(t) = I(t) + i Q(t)$$

$`\begin{aligned} I(t) &= G(t) = A \exp\left[ -\frac{(t - t_0)^2}{2\sigma^2} \right] \\ Q(t) &= -\beta \frac{d}{dt} I(t) = \beta \frac{(t - t_0)}{\sigma^2} G(t) \end{aligned}`$
where $\beta$ is a dimensionless scaling factor optimised to cancel non-adiabatic transitions.
#### Physical Mechanism of Suppression
In the frequency domain, the derivative operation translates to multiplication by $i\omega$. The imaginary quadrature component $Q(t)$ creates a destructive interference path precisely at the detuned |1> -> |2> transition frequency $\omega_{12}$:
$$\tilde{\Omega}(\omega) = \tilde{I}(\omega) \left[ 1 + \beta (\omega - \omega_d) \right]$$

By tuning $\beta \approx -\frac{1}{2\alpha}$ (to first order in perturbation theory), $\tilde{\Omega}(\omega_{12})$ drops to zero, effectively spectral-binding the drive signal to the |2> level while maintaining a fast, high-fidelity $\pi$-pulse on the |0> <-> |1> computational transition.

## Total Dynamic Hamiltonian and Operator Representation
### Baseband I/Q Drive Modulation
To model general complex-valued baseband control signals $\Omega(t) = I(t) + i Q(t)$, the total interaction Hamiltonian is decomposed into orthogonal in-phase (I) and quadrature (Q) control channels. In the rotating frame of the drive frequency $\omega_d$, the total time-dependent Hamiltonian H(t) takes the unified form: $$H(t) = H_0^{\text{rot}} + H_{\text{drive}}^{\text{rot}}(t)$$
$$H(t) = H_0^{\text{rot}} + \frac{\hbar}{2} \left[ I(t) \mathcal{O}_x + Q(t) \mathcal{O}_y \right]$$
where, $H_0^{\text{rot}}$ = static drift Hamiltonian
       $\mathcal{O}_x, \mathcal{O}_y$ = dimensionless system coupling operators 

### Explicit Matrix Representations
#### Two-level Subspace (d=2)
For the simplified N=2 qubit system, the coupling operators correspond directly to the Pauli spin matrices ($\mathcal{O}_x = \sigma_x$, $\mathcal{O}_y = \sigma_y$). Under resonant drive conditions ($\Delta = 0$), $H_0^{\text{rot}} = \mathbf{0}_2$, yielding the 2x2 matrix system:
$`H_{d=2}(t) = \frac{\hbar}{2} \begin{pmatrix} 0 & I(t) - i Q(t) \\ I(t) + i Q(t) & 0 \end{pmatrix}`$
#### Three-level Transmon Subspace (d=3)
For the N=3 qutrit system, the static drift matrix under resonance ($\Delta = 0$) encapsulates the transmon anharmonicity $\alpha$:
$`H_0^{\text{rot}} = \hbar \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & \alpha \end{pmatrix}`$
The drive coupling matrices $\mathcal{O}_x = X$ and $\mathcal{O}_y = Y$ are constructed from the truncated harmonic oscillator creation ($a^\dagger$) and annihilation ($a$) operators:
$`X = a + a^\dagger = \begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & \sqrt{2} \\ 0 & \sqrt{2} & 0 \end{pmatrix}, \quad Y = -i(a - a^\dagger) = \begin{pmatrix} 0 & -i & 0 \\ i & 0 & -i\sqrt{2} \\ 0 & i\sqrt{2} & 0 \end{pmatrix}`$
Substituting these operators into the drive equation produces the explicit 3x3 time-dependent matrix evaluated at each integration timestep t:
$`H_{d=3}(t) = \hbar \begin{pmatrix}  0 & \frac{1}{2}\left(I(t) - i Q(t)\right) & 0 \\  \frac{1}{2}\left(I(t) + i Q(t)\right) & 0 & \frac{\sqrt{2}}{2}\left(I(t) - i Q(t)\right) \\  0 & \frac{\sqrt{2}}{2}\left(I(t) + i Q(t)\right) & \alpha  \end{pmatrix}`$

### Unified Numerical Integration Mechanics
This formulation provides a unified computational interface:
* Pulse Modulations: Gaussian controls set $Q(t) = 0$, reducing $H_{d=3}(t)$ to a real-symmetric matrix driving both the |0> <-> |1> and |1> <-> |2> channels simultaneously.
* DRAG Controls: Activating the derivative quadrature $Q(t) = -\beta \frac{d}{dt}I(t)$ introduces imaginary off-diagonal terms. These non-zero imaginary elements generate a phase shift during state evolution, driving destructive interference that cancels population transfer across the upper |1> <-> |2> coupling branch.
