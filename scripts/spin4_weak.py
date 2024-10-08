from utils import define_system_operators, spin_chain_H, stochastic_solve
from qutip import fock, tensor
import numpy as np

# Define parameters for the system
N = 2  # Dimension for spin-1/2
N_spin = 4  # Number of spins
J = np.pi   # Coupling constant
h = 1       # Magnetic field strength
gamma = 0.0001  # Measurement rate
times = np.linspace(0, 50, 1000)  # Time evolution points

# Define system operators
sigmax_list, sigmay_list, sigmaz_list = define_system_operators(N, N_spin)

# Construct the Hamiltonian for the spin chain
H = spin_chain_H(N_spin, J, h, sigmax_list, sigmay_list, sigmaz_list)

# Initial state
rho0 = tensor([fock(2,1)]+[fock(2,0) for _ in range(N_spin-1)])

# Define Lindblad operator and collapse operators
L = np.sqrt(gamma) * sigmaz_list[2]
c_ops = []
sc_ops = [L]
e_ops = [sigmaz_list[0], sigmaz_list[-1]]

# Run the stochastic solver
result = stochastic_solve(H, times, rho0, gamma, L, c_ops, sc_ops, e_ops)

qsave(result, 'spin4_weak')
