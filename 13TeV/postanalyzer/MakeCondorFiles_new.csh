#!/bin/tcsh
setenv pwd $PWD
#echo ${pwd}
cat>Job_${1}.csh<<EOF
#!/bin/tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh
cd /uscms_data/d3/varun/13TeV/Qstar/CMSSW_7_1_5/src
setenv SCRAM_ARCH slc6_amd64_gcc481
cmsenv
cd ${pwd}
cd \${_CONDOR_SCRATCH_DIR}
./${1}.exe
EOF

chmod 775 ${pwd}/Job_${1}.csh

cat>condor_${1}<<EOF
universe = vanilla
Executable = Job_${1}.csh
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
transfer_input_files = $pwd/${1}.exe,${2}
Output               = ${1}_\$(Cluster)_\$(Process).stdout
Error                = ${1}_\$(Cluster)_\$(Process).stderr
Log                  = ${1}_\$(Cluster)_\$(Process).log
Queue 1
EOF

condor_submit condor_${1}
