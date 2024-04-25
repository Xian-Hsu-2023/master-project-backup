#!/bin/bash
#SBATCH --job-name=plot
#SBATCH --partition=mem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=26
#SBATCH --mem=64GB
#SBATCH --time=12:00:00

export CORENUM=26
export SCR=/lfs/data/xianhsu/gamer/bin/weinberger/analysis
export EXEC="gasmass_coolingtime_mpi.py"

module purge
module load python
source activate astro

cd $SCR
echo Job begins at `date`
echo Start to run $EXEC
mpiexec -n $CORENUM python $EXEC
# srun -n $CORENUM python $EXEC
echo Run complete at `date`