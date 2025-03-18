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
