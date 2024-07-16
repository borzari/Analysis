import os

os.system("cd ../../ && scram b -j 8 && cd - && cmsRun zToLepProbeTrk_cfg.py")
