#!/bin/bash
rm Local* -rf
timestamp=`date +%Y%m%d%H%M`
echo time stamp: ${timestamp}

star-submit-template -template submitByRun_perFile.xml -entities listOfFiles=/star/u/lazhang/data01/CF/urqmd/K0CF_urqmd/analysis/UrQMD_cascade_3GeV_70M.list,cutFile=/star/u/lazhang/data01/CF/urqmd/analysis_v3.0/cutDefault.txt
jobId2=`condor_q | grep "sched" | head -1 | sed 's/.*sched//g' | sed 's/_.*$//g'`
echo jobId2: ${jobId2}
./monitorAndSubmit.sh result_${timestamp}.root
./rm.sh
