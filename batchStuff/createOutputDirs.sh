#!/bin/bash                                                                                                                                    

#echo $d " counter is " $counter;                                                                                                              
counter=0;
subCounter=0;
dateString=`date +%d%b%Y`


mkdir /group/belle/users/vossen/ptSpect/
mkdir /group/belle/users/vossen/ptSpectOut/

for ex in 07 09 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73
do
#for res in  continuum
for res in  on_resonance
do
for spec in uds charm mixed charged
do

myDir=subMC_ex$ex\_$res\_$spec
echo " dir: $myDir " ;
mkdir /group/belle/users/vossen/ptSpect/$myDir
mkdir /group/belle/users/vossen/ptSpectOut/$myDir


done
done 
done








