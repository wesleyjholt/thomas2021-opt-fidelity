#!/bin/sh

run_on_HPC=0    # 1 if running on high-performance computer, 0 if doing a test on a local machine
test=1          # 0 for no test, 1 for light test, 2 for heavy test

if [ $test = 1 ]; then
    # light test
    nturbines_vec=(9)           # number of turbines in farm
    diameter_spacing_vec=(7)    # approximate spacing between turbines (in rotor diameters)
    ndirs_vec=(10)              # number of directions in wind rose
    nspeeds_vec=(1)             # number of speeds in wind rose
    wake_models=(Gaussian)      # wake model
elif [ $test = 2 ]; then
    # heavy test
    nturbines_vec=(9 16)
    diameter_spacing_vec=(7)
    ndirs_vec=(10 15)
    nspeeds_vec=(1)
    wake_models=(JensenCosine Gaussian)
else
    # no test, all levels included
    nturbines_vec=(9 16 25 38)
    diameter_spacing_vec=(5 7 10)
    ndirs_vec=(8 10 11 12 13 15 16 18 20 24 30 36 45 60 90 180 360)
    nspeeds_vec=(1 20)
    wake_models=(JensenCosine Gaussian)
fi

# maximum number of cores
ntasksmax=250

# run optimizations for each combination of variables
for nturbines in ${nturbines_vec}
do
    for diameter_spacing in ${diameter_spacing_vec}
    do
        for ndirs in ${ndirs_vec[@]}
        do
            for nspeeds in ${nspeeds_vec[@]}
            do
                for wake_model in ${wake_models[@]}
                do
                    if [ $run_on_HPC = 1 ]; then
                        # run on high-performance computer using slurm

                        # decide how many cores to use
                        ntasksdefault=$(expr $ndirs \* $nspeeds / 500)
                        _ntasks=$(($ntasksdefault<$ntasksmax ? $ntasksdefault : $ntasksmax))
                        _ntasks=$(($_ntasks>1 ? $_ntasks : 1))

                        # submit job    
                        echo "Now submitting jobs for $ndirs directions and $nspeeds speeds, using $_ntasks CPUs."
                        sbatch --ntasks=$_ntasks run_opt_circle.sh \
                        $ndirs \
                        $nspeeds \
                        $wake_model \
                        $nturbines \
                        $diameter_spacing

                        sleep 1

                    else
                        # run locally
                        layout_number=001

                        julia run_opt_circle.jl \
                        $layout_number \
                        $ndirs \
                        $nspeeds \
                        $wake_model \
                        $nturbines \
                        $diameter_spacing

                        sleep 1
                    fi
                    
                done
            done
        done
    done
done