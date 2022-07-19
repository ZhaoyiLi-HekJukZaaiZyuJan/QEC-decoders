#!/bin/sh
#SBATCH -o results/%A_%a_terminal.out #STDOUT                 
#SBATCH --ntasks=1  
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2          
#SBATCH --mail-type=END,FAIL    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ladmon@stanford.edu
#SBATCH --array=1-100         # Array range

ml py-tensorflow/2.6.2_py36
module load gcc/10.1.0
# srun ./simulate -s PLANE --plmin 0 --plmax 0.25 --pmin 0 --pmax 0.033  --Np 20 --Npl 20 -n 100 --Lmin 3 -v 1 -N INDEP --Nl INDEP_loss --fname "results/${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID.out" --thread 1
./simulate -s TORUS --pmin 0 --pmax 0.15 --Np 30 -n 10000 --Lmin 3 --Lmax 21 -v 1 --fname 'results/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out'

echo "----------------------------------"
echo "id:" "$SLURM_ARRAY_JOB_ID" "$SLURM_ARRAY_TASK_ID"
echo "cpu per task" "$SLURM_CPUS_PER_TASK"
echo "nodelist" "$SLURM_JOB_NODELIST"
echo "cluster name:" "$SLURM_CLUSTER_NAME"