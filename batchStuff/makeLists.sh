#!/bin/bash


for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
for res in continuum on_resonance
do
find  /group/belle/bdata_b/dstprod/dat/e0000$ex/HadronBJ/ -wholename "*$res*.mdst" > lists/data$ex\_$res.list
#find  /group/belle/bdata_b/mcprod/dat/e0000$ex/evtgen/$spec/ -wholename "*$res*.mdst" > lists/mc$ex\_$res\_$spec.list
#find /pic/projects/Belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP$ex/$res/ > data$ex\_$res.list
done
done