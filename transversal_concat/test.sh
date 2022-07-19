SLURM_CLUSTER_NAME=1
SLURM_JOB_NODELIST=2
SLURM_JOB_ID=4
SLURM_ARRAY_TASK_ID=45

echo "results/${SLURM_JOB_ID}_$SLURM_ARRAY_TASK_ID.out"

echo "----------------------------------"
echo "id:" "$SLURM_JOB_ID"
echo "command:" "srun ./simulate -s PLANE --pmin 0.01 --pmax 0.05  --Np 20 --Npl 1 -n 100 --Lmin 3 -v 1"
echo "cluster name:" "$SLURM_CLUSTER_NAME"
echo "node list:" "$SLURM_JOB_NODELIST"
