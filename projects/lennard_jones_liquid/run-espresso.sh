#!/bin/bash
# ====================================================
# for running jobs using Espresso
# ====================================================


#PBS -l select=1:ncpus=1:mem=500mb:scratch_local=2gb:cluster=tarkil
#PBS -l walltime=2:00:00
#PBS -N lj-liquid-tomanoem
trap 'clean_scratch' TERM EXIT
trap 'cp -r $SCRATCHDIR/temporary $DATADIR && clean_scratch ' TERM


# DATA
# ====================================================
job_name="lj-liquid"
date=$(date '%Y-%m-%d--%H:%M:%S')

DATADIR="/storage/vestec1-elixir/home/tomanoem"
cp $DATADIR/$job_name.py $SCRATCHDIR || exit 1
cd $SCRATCHDIR || exit 2
WORKDIR="${job_name}_${date}"
mkdir $WORKDIR

wd="$SCRATCHDIR/$WORKDIR"
cp $job_name.py $wd || exit 2
cd $wd


# COMPUTATIONAL PART
# ======================================================
module add espresso_md-4.1.4
#
for i in {1, 0.5, 0.25}
do
    for j in {0, -3, -6}
    do
        mpirun -n 1 pypresso $job_name.py 1 ${i} ${j}
        j=$((j+1))
    done
    i=$((i+1))
    j=1
done


#mpirun -n 1 pypresso $job_name.py

# COPY DATA AND LEAVE
# =======================================================
cd ..
cp -r $WORKDIR $DATADIR || export CLEAN_SCRATCH=false

exit 0
