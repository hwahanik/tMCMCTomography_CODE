#!/bin/sh

# Grid Engine Options
#$ -N mksamples
#$ -cwd
#$ -l h_rt=48:00:00
#$ -pe openmpi_fillup_mark2 16
#$ -R y

# Initialise the environment
. /etc/profile
. /etc/profile.d/modules.sh

# Use the Intel v12.0 compiler with IntelMPI
module load intel
module load openmpi-intel

# Launch your MPI program using the Ethernet network
mpiexec -np $NSLOTS ./mksamples
