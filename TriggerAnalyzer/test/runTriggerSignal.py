import os

os.system("cd .. && scram b -j 8 && cd test/ && cmsRun triggerSignalAnalysis_cfg.py")
