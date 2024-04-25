export NUM='0'
export FNAME='Data_00'
export DIRECTORY='s008_ref4_0.001'
# export LEN=$(wc -w <<< "$NUM")


# python radial_profile_dens_arg.py -n $NUM -f $FNAME -d $DIRECTORY
for dir in $DIRECTORY;
    do
        python radial_profile_dens_arg.py -n $NUM -f $FNAME -d $DIRECTORY &
    done
    wait