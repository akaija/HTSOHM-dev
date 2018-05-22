#!/bin/bash
#
# $Revision: 1.0 $
# $Date:  2017-09-11 $
# $Author: akaija $

#SBATCH -N 1
#SBATCH -t 0-24:00 # Runtime in D-HH:MM
#SBATCH --cpus-per-task=16
#SBATCH --mem=4g

echo JOB_ID: $SLURM_JOB_ID JOB_NAME: $SLURM_JOB_NAME HOSTNAME: $SLURM_SUBMIT_HOST
echo start_time: `date`

# dependencies
module purge
module load python/3.5.1
module load postgresql/9.5.2
source ~/venv/htsohm/bin/activate

cd $SLURM_SUBMIT_DIR
for ((i = 0; i < $SLURM_CPUS_ON_NODE; i++))
do
    echo TASK NO.: $i
    ./hts.py launch_worker ${SLURM_JOB_NAME} >> ${SLURM_JOB_NAME}/output_${SLURM_SUBMIT_HOST}_$$_$i.log 2>&1 &
done

wait
