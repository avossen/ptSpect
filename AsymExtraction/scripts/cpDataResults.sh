
for ex in 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
for res in continuum 
do
myDir=subData_ex$ex\_$res
cp mDstOut/$myDir/Normal_*.root mDstOut/all
done
done

