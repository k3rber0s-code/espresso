#!/bin/bash
# ====================================================
# for running jobs using Espresso
# ====================================================


#PBS -l select=1:ncpus=1:mem=500mb:scratch_local=2gb:cluster=luna
#PBS -l walltime=72:00:00
#PBS -N minimal-polymer-tomanoem
trap 'clean_scratch' TERM EXIT
trap 'cp -r $SCRATCHDIR/temporary $DATADIR && clean_scratch ' TERM


# DATA
# ====================================================
job_name="minimal-polymer"
date="2022-10-21_11-32"

DATADIR="/storage/vestec1-elixir/home/tomanoem/experiments/minimal_polymer"
cp $DATADIR/$job_name.py $SCRATCHDIR || exit 1
cd $SCRATCHDIR || exit 2
WORKDIR="${job_name}_${date}"
mkdir $WORKDIR

wd="$SCRATCHDIR/$WORKDIR"
cp $job_name.py $wd || exit 2
cd $wd


# COMPUTATIONAL PART
# ======================================================
module add debian10-compat
module add espresso_md-4.1.4

for i in 2 5 10 20 50; do
  mpirun -n 1 pypresso $job_name.py --n_beads_per_chain ${i}
done



# COPY DATA AND LEAVE
# =======================================================
cd ..
cp -r $WORKDIR $DATADIR || export CLEAN_SCRATCH=false

exit 0
