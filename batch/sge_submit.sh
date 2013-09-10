#!/bin/sh

# SGE configuration:
#$ -cwd
#$ -o output/stdout-1.log
#$ -e output/stderr-1.log
#$ -t 1-100

# Setup scientific python environment
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/beowulf/lib/atlas:/beowulf/lib/
export PYTHONPATH=/beowulf/python-modules
export MATPLOTLIBDATA=/beowulf/python-modules/matplotlib/

# setup the execution environment

# Execute the wrapper script
python $1 100 $SGE_TASK_ID $$ $2
