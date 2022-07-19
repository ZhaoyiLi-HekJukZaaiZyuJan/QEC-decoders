#!/bin/sh
#SBATCH -o /scratch/users/ladmon/3D/results/%A_%a_terminal.out #STDOUT                 
#SBATCH --ntasks=1  
#SBATCH --time=48:00:00  
#SBATCH --cpus-per-task=2           
#SBATCH --mail-type=END,FAIL    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ladmon@stanford.edu
#SBATCH --array=1-100         # Array range

module load gcc #for slac cluster
# module load gcc/10.1.0 #for sherlock cluster
# ./simulate -s PLANE --qmin 0 --qmax 0.25 --pmin 0 --pmax 0.033  --Np 20 --Nq 20 -n 100 --lmin 3 -v 1 -N INDEP --Nl INDEP_loss --fname "results/${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID.out"

### GATE/EM2/EM2_full run
./simulate -s PLANE --pmin 0.01 --pmax 0.01 --qmin 0 --qmax 0 --Np 1 --Nq 1 -n 1000000 --lmin 3 --lmax 19 -v 1 -N INDEP --fname "/scratch/users/ladmon/3D/results/${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID.out"
# ./simulate -s PLANE --qmin 0.05 --qmax 0.05 --pmin 0.004 --pmax 0.008  --Np 10 --Nq 1 -n 500 --lmin 3 --lmax 13 -v 1 --thread 0 --use_env 1 --fname "results/${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID.out"
# ./simulate -s PLANE --qmin 1000 --pmin 0.005 --pmax 0.01  --Np 30 --Nq 1 -n 1000 --lmin 3 --lmax 21 -v 1 -N GATE_biased --fname "/scratch/users/ladmon/3D/results/${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID.out"

echo "----------------------------------"
echo "id:" "$SLURM_ARRAY_JOB_ID"
echo "cpu per task" "$SLURM_CPUS_PER_TASK"
echo "nodelist" "$SLURM_JOB_NODELIST"
echo "cluster name:" "$SLURM_CLUSTER_NAME"
echo "node list:" "$SLURM_JOB_NODELIST"
echo "node list:" "$SLURM_JOB_NODELIST"