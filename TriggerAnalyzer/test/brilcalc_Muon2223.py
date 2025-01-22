import os
from os.path import exists

crabDir = "crab/"
data2223 = ["crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022C","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022D","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022E","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022F","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022G","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C0","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C1","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023D0","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023D1"] # Data: 2022 + 2023
data22 = ["crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022C","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022D","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022E","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022F","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2022G"] # Data: 2022
data23 = ["crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C0","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023C1","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023D0","crab_allMuonPt_Muon_MET105IsoTrk50Sig_2023D1"] # Data: 2023

resubDirs = data2223

jsonFile = "/results/processedLumis.json"

for dir in resubDirs:
    # os.system("brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i " + crabDir + dir + jsonFile)
    print("brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i " + crabDir + dir + jsonFile + " && \\")