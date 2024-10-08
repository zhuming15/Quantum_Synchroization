# utils.py

from qutip import qeye, sigmax, sigmay, sigmaz, tensor, SMESolver
import numpy as np

def define_system_operators(N, N_spin):
    """
    Define the system operators for a spin chain, including identity matrices
    and the Pauli operators for each spin.

    Parameters:
    N (int): Dimension of the spin (for spin-1/2, N = 2).
    N_spin (int): Number of spins in the system.

    Returns:
    tuple: Three lists containing the tensor products of the Pauli-X, Pauli-Y, 
           and Pauli-Z operators for each spin.
    """
    # Create identity matrices and lists for Pauli operators for each spin
    qeye_list = [qeye(N) for _ in range(N_spin)]
    sigmax_list = []
    sigmay_list = []
    sigmaz_list = []

    # Populate the lists with the tensor products of Pauli operators and identity matrices
    for i in range(N_spin):
        temp_x = qeye_list.copy()
        temp_x[i] = sigmax()
        sigmax_list.append(tensor(temp_x))

        temp_y = qeye_list.copy()
        temp_y[i] = sigmay()
        sigmay_list.append(tensor(temp_y))

        temp_z = qeye_list.copy()
        temp_z[i] = sigmaz()
        sigmaz_list.append(tensor(temp_z))

    return sigmax_list, sigmay_list, sigmaz_list


def spin_chain_H(N_spin, J, h, sigmax_list, sigmay_list, sigmaz_list):
    """
    Constructs the Hamiltonian for a spin chain with nearest-neighbor interactions 
    and a magnetic field term.

    Parameters:
    N_spin (int): Number of spins in the chain.
    J (float): Coupling constant for the spin-spin interaction.
    h (float): Magnetic field strength.
    sigmax_list, sigmay_list, sigmaz_list (list): Lists of Pauli operators for each spin.

    Returns:
    Qobj: The Hamiltonian of the spin chain.
    """
    # Interaction Hamiltonian (nearest-neighbor coupling)
    H_int = sum(J / 2 * (sigmax_list[j] * sigmax_list[j + 1] + sigmay_list[j] * sigmay_list[j + 1])
                for j in range(N_spin - 1))

    # Magnetic field Hamiltonian
    H_field = sum(h * sigmaz_list[j] for j in range(N_spin))

    return H_int + H_field


def stochastic_solve(H, times, rho0, gamma, L, c_ops, sc_ops, e_ops, ntraj=500, nsubsteps=50, parallel='parallel'):
    """
    Solve the stochastic master equation for a given Hamiltonian using homodyne detection.

    Parameters:
    H (Qobj): Hamiltonian of the system.
    times (array-like): Array of time points for the evolution.
    rho0 (Qobj): Initial density matrix of the system.
    gamma (float): Measurement rate.
    L (Qobj): Lindblad operator for the system.
    c_ops (list): List of collapse operators.
    sc_ops (list): List of stochastic collapse operators.
    e_ops (list): List of expectation value operators.
    ntraj (int): Number of trajectories for the stochastic simulation (default is 500).
    nsubsteps (int): Number of substeps for each time step in the solver (default is 50).
    parallel (str): Method to compute trajectories, either 'serial' or 'parallel' (default is parallel).

    Returns:
    Result: Result object from the stochastic master equation solver.
    """
    # Solver options configuration
    options = {
        "store_final_state": True,
        "store_states": True,
        "store_measurement": False,
        "normalize_output": True,
        "method": "platen",
        "map": "serial",
        "progress_bar": "text",
        "keep_runs_results": True
    }

    # Create the SMESolver instance
    sol = SMESolver(H, sc_ops=sc_ops, c_ops=c_ops, options=options, heterodyne=False)

    # Configure measurement operators and Wiener process factors
    sol.m_ops = [L + L.dag()]
    sol.dW_factors = [1 / np.sqrt(gamma)]

    # Execute the stochastic master equation solver
    result = sol.run(rho0, tlist=times, ntraj=ntraj, e_ops=e_ops)

    return result
