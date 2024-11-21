#!/bin/bash
#PBS -q normal
#PBS -j oe
#PBS -l select=1:ncpus=32:mem=128G
#PBS -l walltime=11:59:00
#PBS -P 11003769
#PBS -N spin_sync

# commands start here

module load miniforge3
source activate qutip_env
cd spin_sync
python spin_chain_N_8.py