#!/bin/sh
#SBATCH -o /scratch/users/ladmon/ML/results/%A_%a_terminal.out #STDOUT                 
#SBATCH --ntasks=1  
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2          
#SBATCH --mail-type=END,FAIL    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ladmon@stanford.edu
#SBATCH --array=1-100         # Array range


export LIBRARY_PATH=${LIBRARY_PATH}:~/libtensorflow2/lib:/moismon:/usr/lib64;\
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/libtensorflow2/lib:/moismon:/usr/lib64
# srun ./simulate -s PLANE --plmin 0 --plmax 0.25 --pmin 0 --pmax 0.033  --Np 20 --Npl 20 -n 100 --Lmin 3 -v 1 -N INDEP --Nl INDEP_loss --fname "results/${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID.out" --thread 1
# srun singularity exec ~/lolcow.sif bash -c "export LIBRARY_PATH=${LIBRARY_PATH}:~/libtensorflow2/lib:/moismon:/usr/lib64;\
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/libtensorflow2/lib:/moismon:/usr/lib64;\
srun ./simulate -s TORUS --pmin 0.135 --pmax 0.18 --Np 25 -n 10000 --Lmin 3 --Lmax 17 -v 1 --fname "/scratch/users/ladmon/ML/results/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" -d /scratch/users/ladmon/ML/ --dir 0  -m 'model_h,L=5(7),layer=3x128,epochs=100000,p=' --decode_with_NN --binary --cutoff 0.500005 
# srun ./simulate -s PLANE --plmin 0.05 --plmax 0.05 --pmin 0.004 --pmax 0.008  --Np 10 --Npl 1 -n 500 --Lmin 3 --Lmax 13 -v 1 --thread 0 --use_env 1 --fname "results/${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID.out"

echo "----------------------------------"
echo "id:" "$SLURM_ARRAY_JOB_ID" "$SLURM_ARRAY_TASK_ID"
echo "cpu per task" "$SLURM_CPUS_PER_TASK"
echo "nodelist" "$SLURM_JOB_NODELIST"
echo "cluster name:" "$SLURM_CLUSTER_NAME"