#!/bin/bash
#SBATCH --job-name photo_sim      # Set a name for your job. This is especially useful if you
#SBATCH --partition quick     # Slurm partition to use: quick, norm, 
#SBATCH -c 1        # Number of tasks to run 3
#SBATCH --ntasks-per-core=1	# Do not use hyperthreading (this flag typically used for parallel jobs)
#SBATCH --time 0-3:59       # Wall time limit in D-HH:MM
#SBATCH --mem 8000     # Memory limit for each tasks (in MB) # 1500
#SBATCH -o /gpfs/gsfs8/users/loewingergc/error/photo_sim.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e /gpfs/gsfs8/users/loewingergc/error/photo_sim.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --array 1-1000
#module load gcc/7.1.0-fasrc01 R/3.5.1-fasrc02
#export PATH=/n/home12/gloewinger/apps/gurobi811/linux64/bin:$PATH
#export LD_LIBRARY_PATH=/n/home12/gloewinger/apps/glpk-4.65/lib:/n/home12/gloewinger/apps/gurobi811/#linux64/lib:$LD_LIBRARY_PATH
#export LIBRARY_PATH=/n/home12/gloewinger/apps/glpk-4.65/lib:$LIBRARY_PATH
#export CPATH=/n/home12/gloewinger/apps/glpk-4.65/include:$PATH
#export PATH=~/mosek/9.1/tools/platform/linux64x86/bin:$PATH
#export R_LIBS_USER=$HOME/apps/R_3.5.1:$R_LIBS_USER
#you can find bash in: /gpfs/gsfs8/users/loewingergc/bash 
module load R/4.2.2

Rscript '/gpfs/gsfs8/users/loewingergc/photometry_fglmm/code/photometry_sim_fLME-Science-lengthenRewardPeriod.R' $1