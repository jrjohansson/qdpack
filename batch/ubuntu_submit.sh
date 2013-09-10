#!/bin/sh

#
# time sh batch/ubuntu_submit.sh batch/batchspec.py run_qubits 1
#

# setup the execution environment
TARGETBIN=$2-`date +%F-%H-%M-%S`
echo "Starting batch jobs for the binary " $TARGETBIN
cp $2 $TARGETBIN
mkdir /home/rob/qdpack-data/$TARGETBIN


#!/bin/sh

N=$3


for i in `seq 1 $N`;
do
	python $1 $N $i 1 $TARGETBIN &
	
	sleep 1
done    

wait
