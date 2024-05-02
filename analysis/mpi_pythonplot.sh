export NUM='67'
export FNAME='Data_'
export REGION='100 200 500 1000 8000'
export DIR='y'
# export VARIABLE='density'
# export VARIABLE='kinetic_energy_density'
# export VARIABLE='magnetic_energy_density'
# export VARIABLE='four_velocity_magnitude'
# export VARIABLE='cooling_time'
# export VARIABLE='four_velocity_magnitude'
# export VARIABLE='magnetic_field_magnitude'
# export VARIABLE='faraday_rotation_y'
# export VARIABLE='particle_mass ParDens'
# export VARIABLE='temperature'
# export VARIABLE='magnetic_field_magnitude faraday_rotation_y'
export VARIABLE='kinetic_energy_density density four_velocity_magnitude cooling_time magnetic_field_magnitude faraday_rotation_y temperature'
# export VARIABLE='lobe_density shock_density'

for var in $VARIABLE;
    do
        for region in $REGION;
            do
                for dir in $DIR;
                    do
                        mpiexec -n 10 python ani_mpi_arg_all.py -n $NUM -f $FNAME -r $region -d $dir -v $var &
                    done
            done
        wait
    done
wait
echo everything completed