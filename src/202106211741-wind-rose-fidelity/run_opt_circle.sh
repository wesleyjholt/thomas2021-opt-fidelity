#!/bin/sh

#SBATCH --mem-per-cpu=2056  # memory per CPU core
#SBATCH --time=00:05:00 # walltime
#SBATCH -J 'Wind Rose Study, Horns Rev 1 Wind Farm'
#SBATCH --array=1-2    # job array number corresponds to the layout numbers
#SBATCH --qos=test

# get input arguments
layout_number=$(printf %3s $SLURM_ARRAY_TASK_ID | tr ' ' 0)
ndirs=$1
nspeeds=$2
wake_model=$3
nturbines=$4
diameter_spacing=$5

module load julia
julia run_opt_circle.jl \
    $layout_number \
    $ndirs \
    $nspeeds \
    $wake_model \
    $nturbines \
    $diameter_spacing