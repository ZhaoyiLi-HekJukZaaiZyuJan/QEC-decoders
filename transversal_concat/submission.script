#!/bin/sh
#SBATCH -o /scratch/users/ladmon/422/results/%A_%a_terminal.out #STDOUT  
#SBATCH --time=48:00:00            
#SBATCH --ntasks=1    
#SBATCH --cpus-per-task=2           
#SBATCH --mail-type=END,FAIL    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ladmon@stanford.edu
#SBATCH --array=1-100           # Array range


module load gcc/10.1.0
./simulate -N INDEP -n 10000 --pmin 0.035 --pmax 0.055 --lmin 3 --Np 30 -v 1 --fname "/scratch/users/ladmon/422/results/${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID.out"

echo "----------------------------------"
echo "id:" "$SLURM_ARRAY_JOB_ID"
echo "cpu per task" "$SLURM_CPUS_PER_TASK"
echo "nodelist" "$SLURM_JOB_NODELIST"
echo "cluster name:" "$SLURM_CLUSTER_NAME"
echo "node list:" "$SLURM_JOB_NODELIST"
echo "node list:" "$SLURM_JOB_NODELIST"