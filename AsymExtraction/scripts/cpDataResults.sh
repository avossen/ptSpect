for ex in 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
for res in continuum on_resonance
#for res in on_resonance
do
for pt in kT qT
do
myDir=subData_ex$ex\_$res\_$pt
cp mDstOut/$myDir/Normal_*.root mDstOut/all_$res\_$pt
done
done
done
