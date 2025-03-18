#!/bin/bash -l
#SBATCH --partition=ccm
#SBATCH --constraint=icelake
#SBATCH --time=10080

#SBATCH --nodes=16
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16

#SBATCH -J precomp-close-interac
## SBATCH -o /mnt/home/dmalhotra/ceph/precomp-close-interac-%j.out
## SBATCH -e /mnt/home/dmalhotra/ceph/precomp-close-interac-%j.err
## SBATCH -p ccm --time=1440 -C 'skylake&opa' -N14 --ntasks-per-node=8 --exclusive
## SBATCH --mail-type=ALL --mail-user=dhairya.android@gmail.com

module --force purge
module load modules
module load gcc
module load intel-oneapi-mkl
module load openmpi4

source ~/.bashrc.local
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

WORK_DIR=/mnt/home/dmalhotra/sandbox/close-to-touching/BIE2D/
OUT_DIR=${WORK_DIR}

cd ${WORK_DIR}
mkdir -p ${OUT_DIR}
env | grep -i slurm | tee -a ${OUT_DIR}./out.txt

#time make -B DEBUG=0 bin/stokes-mobility
time mpirun --report-bindings -n $((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE})) --map-by slot:pe=$OMP_NUM_THREADS ${WORK_DIR}./bin/stokes-mobility 1e-10 2>&1 | tee -a ${OUT_DIR}./out.txt

