export NUM='101'
export FNAME='Data_'
export DIRECTORY='s005_CompareToFLASH_ref4_0.001 s006_CompareToFLASH_ref4_0.1'

for dir in $DIRECTORY;
    do
        # python write_arg_mpi.py -n $NUM -f $FNAME -d $dir &
        mpiexec -n 40 python write_arg_mpi.py -n $NUM -f $FNAME -d $dir 
        # mpiexec -n 20 python write_arg_mpi.py -n $NUM -f $FNAME -d $dir &
    done

wait
echo everything completed