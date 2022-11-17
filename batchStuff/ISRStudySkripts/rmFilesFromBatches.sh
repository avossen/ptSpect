#!/bin/bash                                                                                                                                    
#need to split up in batches 

#for tune in 00 01 10 11 12 13 14
for tune in 13 
do
#for yn in yes 
for yn in yes no
do
for spec in uds charm
do
for batch in  1 2 3 4 5 
do
myDirBatch=genMC_red_Tune$tune\_ISR$yn\_$spec\_$batch
echo " dir: $myDirBatch " ;
rm /group/belle/users/vossen/ptSpect/ISRStudies/$myDirBatch/*.root
done
done
done 
done








