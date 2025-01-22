import os
from os.path import exists

crabDir = "crab/"
data2223 = ["crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022C","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022D","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022E","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022F","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022G","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C0","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C1","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023D0","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023D1"] # Data: 2022 + 2023
data22 = ["crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022C","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022D","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022E","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022F","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022G"] # Data: 2022
data23 = ["crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C0","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C1","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023D0","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023D1"] # Data: 2023
bgmc2223 = ["crab_allMuonPt_WToLNu_MET105IsoTrk50Sig_2022","crab_allMuonPt_WToLNu_MET105IsoTrk50Sig_2022ext","crab_allMuonPt_WToLNu_MET105IsoTrk50Sig_2022EE","crab_allMuonPt_WToLNu_MET105IsoTrk50Sig_2022EEext","crab_allMuonPt_WToLNu_MET105IsoTrk50Sig_2023","crab_allMuonPt_WToLNu_MET105IsoTrk50Sig_2023BPix"] # MC: 2022 + 2023
bgmc22 = ["crab_allMuonPt_WToLNu_MET105IsoTrk50Sig_2022","crab_allMuonPt_WToLNu_MET105IsoTrk50Sig_2022ext","crab_allMuonPt_WToLNu_MET105IsoTrk50Sig_2022EE","crab_allMuonPt_WToLNu_MET105IsoTrk50Sig_2022EEext"] # MC: 2022
bgmc23 = ["crab_allMuonPt_WToLNu_MET105IsoTrk50Sig_2023","crab_allMuonPt_WToLNu_MET105IsoTrk50Sig_2023BPix"] # MC: 2023
data23vs = ["crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C0v1","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C0v2","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C0v3","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C1v1","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C1v2","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C1v3","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023D0v1","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023D1v1"] # Data: 2023 missing versions
data23vsleft = ["crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C0v1"] # Data: 2023 missing versions

# checkDirs = data2223 + bgmc2223
checkDirs = data23vsleft

if exists("checkingCRAB.txt"): os.system("rm checkingCRAB.txt")

for dir in checkDirs:
    os.system("crab status -d " + crabDir + dir + " >> checkingCRAB.txt")
