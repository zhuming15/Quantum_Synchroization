
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from qutip import *
import time

### DEFINE SYSTEM OPERATORS
N = 2 #for spin 1/2
N_spin = 8

qeye_list = [qeye(N) for i in range(N_spin)] #identity matrices
sigmax_list = [] #list of pauli operators for each spin
sigmay_list = []
sigmaz_list = []

for i in range(N_spin): #setup list of annihilation operators
    '''
    E.g. For three oscillators, this builds the list 
    [tensor(a,qeye,qeye), tensor(qeye,a,qeye), tensor(qeye,qeye,a)]
    '''
    
    temp1 = qeye_list.copy()
    temp1[i] = sigmax() 
    sigmax_list.append(tensor(temp1))
    
    temp2 = qeye_list.copy()
    temp2[i] = sigmay() 
    sigmay_list.append(tensor(temp2))
    
    temp3 = qeye_list.copy()
    temp3[i] = sigmaz() 
    sigmaz_list.append(tensor(temp3))
    
def spin_chain_H(N_spin,J,h):
    
    H_int = 0
    for j in range(N_spin-1):
        H_int += J/2 * (sigmax_list[j]*sigmax_list[j+1] + sigmay_list[j]*sigmay_list[j+1])
    
    H_field = 0
    for j in range(N_spin):
        H_field += h * sigmaz_list[j]
        
    return H_int + H_field

def stochastic_solve(H,gamma,times,rho0,L):

    ntraj=200
    #L = np.sqrt(gamma) * sigmaz_list[2] 
    
    #rho0 = tensor([fock(N, 0) for i in range(N_spin)])
    
    c_ops = []
    sc_ops = [L]
    e_ops = [sigmaz_list[0], sigmaz_list[-1]]
    #e_ops = []
    
    '''qutip 4.5
    opt = Options()
    #opt.store_states = False
    result = smesolve(H, rho0, times, c_ops, sc_ops,
                  e_ops, ntraj=ntraj, nsubsteps=50, solver='platen',
                  m_ops=[(L + L.dag())], dW_factors=[1/np.sqrt(gamma)],
                  method='homodyne', store_measurement=True, normalize=True,
                  map_func=parallel_map, options=opt)
    '''
    
    options = {
    "method": "platen",
    "dt": 0.001,
    "store_measurement": True,
    "store_states": True,
    "keep_runs_results": True,
    "map": "parallel",
}

    solver = SMESolver(H, sc_ops=sc_ops, heterodyne=False, options=options)
    solver.m_ops = [(L + L.dag())]
    solver.dW_factors = [np.sqrt(2)]
    result = solver.run(rho0, times, e_ops=e_ops, ntraj=ntraj)
    
    
    
    
    return result 

def DFS(N_spin):
    basis_list = [basis(2,0) for i in range(N_spin)]
    DFS_k = 0
    DFS_l = 0
    
    for i in range(N_spin):
        temp = basis_list.copy()
        temp[i] = np.sqrt(2/(N_spin+1)) * np.sin((i+1)*np.pi/3) * basis(2,1)
        DFS_k += tensor(temp)
        temp[i] = np.sqrt(2/(N_spin+1)) * np.sin(2*(i+1)*np.pi/3) * basis(2,1)
        DFS_l += tensor(temp)
        
    return DFS_k,DFS_l

H = spin_chain_H(N_spin,3,1)
times = np.linspace(0,10,501)
gamma = 0.7
L = np.sqrt(gamma) * sigmaz_list[2]
DFS_k,DFS_l = DFS(N_spin)
rho0 = 1/np.sqrt(2) *(DFS_k + DFS_l)
result = stochastic_solve(H,gamma,times,rho0,L)

plt.figure()
plt.plot(times,result.expect[0][0],label='$<\sigma_1^z>$')
plt.plot(times,result.expect[1][0],label='$<\sigma_5^z>$')
plt.legend()
plt.xlabel('Time')
plt.title('Trajectory level')
plt.savefig('Trajectory level,N={}.png'.format(N_spin), bbox_inches='tight')

q1 = []
q2 = []
plt.figure()
for i in range(len(times)):
    q1.append(DFS_k.dag() * result.states[0][i] * DFS_k)
    q2.append(DFS_l.dag() * result.states[0][i] * DFS_l)
    
    
plt.plot(times,q1,label='$|q1|^2$')
plt.plot(times,q2,label='$|q2|^2$')
plt.legend()
plt.xlabel('Time')
plt.title('Probability')
plt.savefig('Probabillity,N={}.png'.format(N_spin), bbox_inches='tight')

plt.figure()
plt.plot(times,result.expect[0].mean(axis=0),label='$<\sigma_1^z>$')
plt.plot(times,result.expect[1].mean(axis=0),label='$<\sigma_5^z>$')
plt.legend()
plt.xlabel('Time')
plt.title('Ensemble level')
plt.savefig('Ensemble level,N={}.png'.format(N_spin), bbox_inches='tight')