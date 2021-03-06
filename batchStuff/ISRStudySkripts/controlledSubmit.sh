#!/bin/bash                

#typical size is 17MB, checking for 10 should be ok
minimumsize=10000000                                                                           
for tune in 13
#for tune in 00 01 10 11 12 13 14
do
for yn in yes no
#for yn in no
do
for spec in  uds charm 
#for spec in charm
do
myDir=genMC_red_Tune$tune\_ISR$yn\_$spec
for d in `ls $myDir/*.sh`;
do
#this has to be updated based on the numbers of '_'
count1="$(cut -d'_' -f9 <<< $d)"
counter="$(cut -d'.' -f1 <<< $count1)"

#echo count1: $count1 counter: $counter

FILE="/group/belle/users/vossen/ptSpect/ISRStudies/$myDir/job_$counter.root"
#echo checking for $FILE
if [ -f $FILE ]; then
#echo exists
actualsize=$(wc -c <"$FILE")
if [ $minimumsize -ge $actualsize ]; then
echo $FILE size: $actualsize
echo file too small
bsub $PWD/$d; 
fi

else
echo $FILE does not exist yet
echo submitting $PWD/$d
bsub $PWD/$d; 

fi

done
done
done
done
