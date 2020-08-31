#! /bin/bash -l

#SBATCH -J BB7DGPVS
#SBATCH -n 16 # Number of cores
#SBATCH -p MCMC # Partition Used.
#SBATCH -t 10-00:00 # Runtime in D-HH:MM
#SBATCH -o JOB%j.out # File to which STDOUT will be written
#SBATCH -e JOB%j.err # File to which STDERR will be written
#SBATCH --mail-type=FAIL # Valid type values are BEGIN, END, FAIL, REQUEUE, and ALL
#SBATCH --mail-user=feng.li@cufe.edu.cn


# CONFIG_FILE=config.BB7.LU.SPLITTPOISSON.BABA-TEXTS.R

# CONFIG_FILE=config.BB7.LU.SPLITT.SP100-SP600.R
# CONFIG_FILE=config.BB7.LU.GARCH.SP100-SP600.R
# CONFIG_FILE=config.BB7.LU.STOCHVOL.SP100-SP600.R

# CONFIG_FILE=config.BB7.SPLITT.SP100-SP600.R
# CONFIG_FILE=config.BB7.GARCH.SP100-SP600.R
# CONFIG_FILE=config.BB7.STOCHVOL.SP100-SP600.R

# CONFIG_FILE=config.CLAYTON.SPLITT.SP100-SP600.R
# CONFIG_FILE=config.CLAYTON.GARCH.SP100-SP600.R
# CONFIG_FILE=config.CLAYTON.STOCHVOL.SP100-SP600.R

# CONFIG_FILE=config.GUMBEL.SPLITT.SP100-SP600.R
# CONFIG_FILE=config.GUMBEL.GARCH.SP100-SP600.R
# CONFIG_FILE=config.GUMBEL.STOCHVOL.SP100-SP600.R

# CONFIG_FILE=config.MVT2.SPLITT.SP100-SP600.R
# CONFIG_FILE=config.MVT2.GARCH.SP100-SP600.R
# CONFIG_FILE=config.MVT2.STOCHVOL.SP100-SP600.R

# CONFIG_FILE=config.MVT3.DF.SPLITT.SP100-SP400-SP600.R
# CONFIG_FILE=config.MVT3.SPLITT.SP100-SP400-SP600.R

## CONFIG_FILE=config.BB7.LU.SPLITT.DGPDATA.R

# CONFIG_FILE=config.MVGARCH.R

## MPIRUN

## mpirun -np 1  ~/.bin/CplRun 4 ~/code/cdcopula/inst/config/${CONFIG_FILE}
## export OMP_NUM_THREADS=16
## mpirun -n 1 /data1/fli/code/cdcopula/inst/scripts/Test.MPI.R
## srun ~/.bin/CplRun  ~/code/cdcopula/inst/config/${CONFIG_FILE}

mpirun -np 1 hostname

## SEND RESULTS
# cat JOB${SLURM_JOB_ID}.out JOB${SLURM_JOB_ID}.err |\
# /usr/bin/mail -s "JOB${SLURM_JOB_ID}.${SLURM_JOB_NAME}(${CONFIG_FILE})" feng.li@cufe.edu.cn
