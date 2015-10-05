#!/bin/bash
#SBATCH -J kfp_emu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=26000
#SBATCH --exclusive
#SBATCH --cpus-per-task=8
#SBATCH --partition=onegbnet
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anatoliy.lane@rutgers.edu

module purge
module load matlab
#module load cplex    (this line is a comment with #)

matlab -nodesktop -nodisplay < script_runner.m > logfile_$SLURM_JOBID
