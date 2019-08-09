
scram p -n CMSSW_10_2_16_CPH CMSSW CMSSW_10_2_16
cd CMSSW_10_2_16_CPH/src
cmsenv
git cms-init
git cms-merge-topic -u danielwinterbottom:from-CMSSW_10_2_16-mvaDM
git clone https://github.com/mbluj/usercode-WawTools.git WarsawAnalysis -b 10_2_X
