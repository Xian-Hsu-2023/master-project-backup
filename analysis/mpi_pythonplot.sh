# warning: please never run this on the login node
# the numbers of the simulation files you want to plot (for this case, it plots from Data_000000 to Data_000066)
export NUM='67'
# the file prefix for your simulation files (for FLASH users, this may be crbub_hdf5_plt_cnt_)
export FNAME='Data_'
# the width for the plot (units: kpc)
export REGION='100 200 500 1000 8000'
# the dirextion (axis) for the plots
export DIR='y'
# the variable(s) for the plots, they will be used respectively
# export VARIABLE='density'
# export VARIABLE='four_velocity_magnitude'
# export VARIABLE='cooling_time'
# export VARIABLE='magnetic_field_magnitude'
# export VARIABLE='temperature'
export VARIABLE='kinetic_energy_density density four_velocity_magnitude cooling_time magnetic_field_magnitude faraday_rotation_y temperature'
# export VARIABLE='lobe_density shock_density'

for var in $VARIABLE;
    do
        for region in $REGION;
            do
                for dir in $DIR;
                    do
                        # run each program with 10 cores. For this case, the actual resource usage = 10 * size of $REGION = 50 cores
                        mpiexec -n 10 python ani_mpi_arg_all.py -n $NUM -f $FNAME -r $region -d $dir -v $var &
                    done
            done
        wait
    done
wait
echo everything completed