#!/bin/bash

# time make -B ./bin/sedimentation1 && \
# time unbuffer ./bin/sedimentation1 --verbose \
# --ndisc 44 --eps 1e-1 --init test \
# --icip-type precond --ts-order 10 --dt 2 \
# --vis-path ./vis2/ --geom-fname vis2/X_23.geom --start-idx 23 \
# | tee -a ./vis2/out.txt

export OMP_NUM_THREADS=62
export OMP_PROC_BIND=spread

N=127
tol=1e-3
T=1000
dt=3.375080
TRG_DIR=./vis/ball-pinv-precond-ksprecon-N${N}-tol3
gmres_tol=$(awk "BEGIN { print ${tol}/${T}*1e-2 }")

mkdir ${TRG_DIR}
cat $0 | tee -a ${TRG_DIR}/out1.txt

make -B ./bin/sedimentation1 && \
time unbuffer ./bin/sedimentation1 --verbose \
--ndisc ${N} --ts-order 10 --ts-adap --ts-tol ${tol} --dt ${dt} --T ${T} --eps 1e-1 \
--icip-type precond --ksprecon --init ball --vis-path ${TRG_DIR}/ \
--geom-fname vis/ball-pinv-precond-ksprecon-tol3/X_47.geom --start-idx 47 \
| tee -a ${TRG_DIR}/out1.txt

#TRG_DIR="vis/ball-pinv-precond-ksprecon-N127-tol3/"
#for ((i = 0 ; i <= 98 ; i++ )); do
#  time unbuffer ./bin/sedimentation1 --T 0 --vis-path ${TRG_DIR}/tmp/ \
#  --geom-fname ${TRG_DIR}/X_${i}.geom --start-idx ${i}
#done

