
import re
import os
conf=['eeuu','eess','eecc','tautau']

fileDict={'eeuu':{},'eess':{},'eecc':{},'tautau':{}}


for exp in range(31,73,2):
    expStr='e'+str(exp)
    print("looking at exp"+expStr)
    for c in conf:
        filename='nagoya_'+c+'_continuum.txt'
 #       print("opening " + filename)
        f=open(filename)
        fileDict[c][expStr]=[]
        for line in f:
            match=re.search(expStr,line)
            if match:
#                print("found exp" + expStr+ ": " + line)
                fileDict[c][expStr].append(line)


#now everything is loaded, construct scripts

#10 files seem to take about 6-7 cpu minuts
numLines=50
os.system("mkdir /group/belle/users/vossen/ptSpect/nonQQStudies/")
os.system("mkdir /group/belle/users/vossen/ptSpectOut/nonQQStudies/")
#for exp in range(31,73,2):
for exp in range(55,57,2):
    expStr='e'+str(exp)
    for c in conf:
#counter for output job files
        counter =0
#counter for input files. Since there are so many, process n (e.g. n=10) input files per job
        fileCounter=0
#        for f in fileDict[c][expStr]:
        numFiles=len(fileDict[c][expStr])
        strDir="nonQQ_"+c+"_"+expStr
        os.system("mkdir "+strDir)
        os.system("mkdir /group/belle/users/vossen/ptSpect/nonQQStudies/"+strDir)
        os.system("mkdir /group/belle/users/vossen/ptSpectOut/nonQQStudies/"+strDir)

        while fileCounter<numFiles:
            targetShFile=strDir+"/job_"+c+"_"+str(counter)+".sh"
            print("target file: "+targetShFile)
            counter=counter+1
            os.system("cp batchHead1.sh "+targetShFile)
#            echo "#BSUB -o  /group/belle/users/vossen/ptSpectOut/ISRStudies/$myDir/jobId_$subCounter.out" >> $targetShFile
#echo "#BSUB -e  /group/belle/users/vossen/ptSpectOut/ISRStudies/$myDir/jobId_$subCounter.err" >> $targetShFile
#echo "#BSUB -J ISRStudy_$subCounter"  >> $targetShFile 
            os.system("cat batchHead2.sh >> " +targetShFile)
            fo=open(targetShFile,'a')
            fo.write("module put_parameter ptSpect rfname\\/group/belle/users/vossen/ptSpect/nonQQStudies/"+strDir+"/job_"+str(counter)+".root\n")
            fo.close()
            print("put ptspect modlue")
            os.system("cat batchMiddle.sh >> "+targetShFile)
            print("put batchmiddle")
            fo=open(targetShFile,'a')
            for i in range(0,numLines):
                #[:-1] necessary to strip the \n
                fo.write("process_event 172.17.196.20:"+ fileDict[c][expStr][fileCounter][:-1]+" 0\n")
                fileCounter=fileCounter+1
                if(fileCounter>=numFiles):
                    break

            fo.close()
#echo "process_event /btmp/$fWOExt 0" >> $targetShFile
            os.system("cat batchEnd.sh >> "+targetShFile)
            os.system("chmod a+x "+targetShFile)


