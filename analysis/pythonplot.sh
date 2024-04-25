export NUM='101'
export FNAME='Data_'
export REGION='100 200 500 1000'
export DIR='y'

for region in $REGION;
    do
        for dir in $DIR;
            do
                python ani_arg_dens.py -n $NUM -f $FNAME -r $region -d $dir &
            done
    done
wait
echo everything completed