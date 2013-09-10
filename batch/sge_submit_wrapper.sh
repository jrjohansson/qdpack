#!/bin/sh

# setup the execution environment
#TARGETBIN=$2-`date +%s`
#TARGETBIN=$2-`date +%F-%T`
TARGETBIN=$2-`date +%F-%H-%M-%S`
echo "Starting batch jobs for the binary " $TARGETBIN
cp $2 $TARGETBIN
mkdir /home/rob/qdpack-data/$TARGETBIN

qsub batch/qdp_sge.sh $1 $TARGETBIN
