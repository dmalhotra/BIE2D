#!/bin/sh

export OMP_NUM_THREADS=7
export NUM_MPI_PROCS=9

time make DEBUG=0 -B bin/sedimentation-parallel

vis_path="ball_vis_$(date +%Y%m%d_%H%M%S)/"
mkdir -p "$vis_path"
git diff > "$vis_path"/git_diff.txt


# ball of discs (forward-Euler time-stepping)
export OMP_NUM_THREADS=62
export NUM_MPI_PROCS=1
unbuffer time mpirun -np ${NUM_MPI_PROCS} --map-by slot:pe=${OMP_NUM_THREADS} ./bin/sedimentation-parallel \
  --vis-path ./"$vis_path" --ndisc 127 --init ball --eps 1e-1 --gravity -1 \
  --ts-order 1 --dt 0.07 --T 1000 \
  --gmres-tol 1e-11 --tol 1e-12 \
  --icip-type precond \
  --ksprecon \
  | tee -a "$vis_path"/log.txt

# # ball of discs (high order adaptive time-stepping)
# unbuffer time mpirun -np ${NUM_MPI_PROCS} --map-by slot:pe=${OMP_NUM_THREADS} ./bin/sedimentation-parallel \
#   --vis-path ./"$vis_path" --ndisc 127 --init ball --eps 1e-1 --gravity -1 \
#   --ts-order 10 --ts-tol 1e-3 --dt 6 --T 1000 --ts-adap \
#   --gmres-tol 1e-11 --tol 1e-12 \
#   --icip-type precond \
#   --ksprecon \
#   | tee -a "$vis_path"/log.txt


# # ICIP + ksprecond
# unbuffer time mpirun -np ${NUM_MPI_PROCS} --map-by slot:pe=${OMP_NUM_THREADS} ./bin/sedimentation-parallel \
#   --vis-path ./"$vis_path" --ndisc 32 --init chain --eps 1e-4 --gravity 0 \
#   --ts-order 10 --ts-tol 1.5e-8 --dt 6 --T 120 \
#   --gmres-iter 4000 --gmres-tol 1e-10 --tol 1e-11 \
#   --icip-type precond \
#   --ksprecon \
#   | tee -a "$vis_path"/log.txt


# # ICIP, no-ksprecond
# unbuffer time mpirun -np ${NUM_MPI_PROCS} --map-by slot:pe=${OMP_NUM_THREADS} ./bin/sedimentation-parallel \
#   --vis-path ./"$vis_path" --ndisc 32 --init chain --eps 1e-4 --gravity 0 \
#   --ts-order 10 --ts-tol 1.5e-8 --dt 6 --T 120 \
#   --gmres-iter 4000 --gmres-tol 1e-10 --tol 1e-11 \
#   --icip-type precond \
#   | tee -a "$vis_path"/log.txt


# # Adaptive
# unbuffer time mpirun -np ${NUM_MPI_PROCS} --map-by slot:pe=${OMP_NUM_THREADS} ./bin/sedimentation-parallel \
#   --vis-path ./"$vis_path" --ndisc 32 --init chain --eps 1e-4 --gravity 0 \
#   --ts-order 10 --ts-tol 1.5e-8 --dt 6 --T 120 \
#   --gmres-iter 4000 --gmres-tol 1e-10 --tol 1e-11 \
#   --icip-type adaptive \
#   --verbose \
#   | tee -a "$vis_path"/log.txt

