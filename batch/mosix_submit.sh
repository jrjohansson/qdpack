#!/bin/sh

# setup the execution environment
#TARGETBIN=$2-`date +%s`
#TARGETBIN=$2-`date +%F-%T`
TARGETBIN=$2-`date +%F-%H-%M-%S`
echo "Starting batch jobs for the binary " $TARGETBIN
cp $2 $TARGETBIN
mkdir /home/rob/qdpack-data/$TARGETBIN


#!/bin/sh

N=$3


for i in `seq 1 $N`;
do
	run python $1 $N $i 1 $TARGETBIN &
	#python $1 25 $SGE_TASK_ID $$ $2
	
	sleep 5
done    

wait
