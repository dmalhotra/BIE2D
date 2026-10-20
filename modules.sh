# Copyright (C) 2023-2026 The Simons Foundation, Inc. and the BIE2D authors (see NOTICE)
# Authors: Dhairya Malhotra, Dan Fortunato, Mariana Martinez Aguilar
# SPDX-License-Identifier: Apache-2.0

# # GCC
# module --force purge
# module load modules
# module load slurm tmux git
# module load gcc/11.3.0
# module load intel-oneapi-mkl
# #module load lib/openblas
# module load lib/fftw3
# module load matlab
# module load python
# module load cuda/12.0.0
# module load openmpi4 openmpi-opa
# #export OMPI_CC=gcc
# #export OMPI_CXX=g++



# Intel
module --force purge
module load modules
module load intel-oneapi-compilers openmpi4 openmpi-intel intel-oneapi-mkl fftw
