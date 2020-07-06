#!/bin/bash                                                                                                                                    

#echo $d " counter is " $counter;                                                                                                              
counter=0;
subCounter=0;
dateString=`date +%d%b%Y`


mkdir /group/belle/users/vossen/ptSpect/ISRStudies/
mkdir /group/belle/users/vossen/ptSpectOut/ISRStudies
for tune in 00 01 10 11 12 13 14
do
for yn in yes no
do
for spec in uds charm 
do

myDir=genMC_Tune$tune\_ISR$yn\_$spec
echo " dir: $myDir " ;
mkdir /group/belle/users/vossen/ptSpect/ISRStudies/$myDir
mkdir /group/belle/users/vossen/ptSpectOut/ISRStudies/$myDir


done
done 
done








