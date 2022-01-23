#!/bin/sh
#SBATCH -A def-vidotto
#SBATCH -n 100
#SBATCH --cpus-per-task=10
#SBATCH --time=10-0:00:00
#SBATCH --job-name=star_vertex_ampls_2
#SBATCH --output=star_vertex_ampls_2.log
#SBATCH --error=star_vertex_ampls_2.err
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=pfrisoni@uwo.ca

echo "Running on: $SLURM_NODELIST"
echo

# start commands

BASEDIR=/home/frisus95/projects/def-vidotto/frisus95/sl2cfoam_next_aggiornata
DATADIR=/home/frisus95/scratch/data_sl2cfoam_next

export LD_LIBRARY_PATH="${BASEDIR}/lib":$LD_LIBRARY_PATH

# number of OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

IMMIRZI=1.2

SHELLS=20

TspinMin=13
TspinMax=16
TCURRENTSPIN=$TspinMin

while [ $TCURRENTSPIN -le $TspinMax ]
do

TJS=$(( TCURRENTSPIN ))

now=$(date)
echo
echo "Starting Lorentzian fulltensor [ TJS = ${TCURRENTSPIN}, shells = ${SHELLS} ]... (now: $now)"

$BASEDIR/bin/vertex-fulltensor -V -h -m 2000 $DATADIR $IMMIRZI $TJS,$TJS,$TJS,$TJS,$TJS,$TJS,$TJS,$TJS,$TJS,$TJS $SHELLS

now=$(date)
echo "... done (now: $now)"
echo

let TCURRENTSPIN=TCURRENTSPIN+1

done

echo
echo "All completed."
