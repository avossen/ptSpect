#!/bin/bash                                                                                                                                    
#use as argument e.g. "sub14Apr/data"

for d in `ls $1*.sh`
do

#while [ `squeue | grep voss | wc -l` -ge 200 ]
while [ `squeue -p shared| grep voss | wc -l` -ge 5000 ]
do 
echo " waiting 20 secs, we have " `squeue -p shared | grep voss | wc -l`  " in queue"
sleep 20;
echo " check again:, we have " `squeue -p shared| grep voss | wc -l`  " in queue"
done
sbatch -p shared $d
#sbatch -n 1  --share $d
done
