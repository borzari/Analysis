# Analysis

This has the BackgroundEstimation, StandardAnalysis and TriggerAnalysis selections from the DisappTrks analysis for Run3. It will have all the needed corrections, histograms and possibly trees to make further needed selections. It also saves the EDM files, given that if new variables from the MiniAOD are needed, they can be easily accessed.

## HelperFunctions

- Implemented a `helperFunctions` class to clean the plugins and to only have to define the distinct functions once. Mostly used for lengthy calculations and procedures
   - It helps to implement other event selections
   - It helps on the overall organization of the code
   - It helps when needing to apply changes to the selection
      - Might find a better way to make explicit changes to the methods if needed
- Implemented a `plotPrintFunctions` class to use for printing cumulative and individual efficiencies and to use to plot histograms when needed
- Implemented a `selectingFunctions` to clean up the selection plugins

## Done for BackgroundEstimation

- Charged leptons
   - Implemented all the leptons tag skims and compared with original
   - Implemented ZtoEleProbeTrk selection and compared with original
   - Implemented Z mass and elec/track charge verification and veto selection in ZtoEleProbetrk
      - Got number OS T&P pairs
      - Got number SS T&P pairs
      - Got number OS T&P pairs passing veto selection
      - Got number SS T&P pairs passing veto selection
   - Implemented ZtoMuProbeTrk
      - Compared working cuts with original analysis and everything is (almost) the same
      - Had some issues with caloGeometry, but solved with using Run3 era in `zToLepProbeTrk_cfg.py`
      - Solved events difference; the METFilters needed to be updated; the implementation of the new filters was propagated to the lepton tag skims
      - Got the correct number of T&P pairs, plotted the number of pairs per event and it seems consistent
   - Implemented tau (e and mu) ZtoLepProbeTrk selections
      - Number of T&P was different because of a mistake in `goodInvMass` method
      - Compared with 100k events and everything worked fine
   - Implemented all fiducial selections using 2018 maps
      - Muon and ECAL will be used for run 3; electron will be changed to the ML approach
   - Implemented leptonTagPt35 leptonTagPt35MetTrig selections
      - Compared with 100k events and everything worked fine
- Implemented crab config file to run BGEst selections

## Ongoing for BackgroundEstimation

- Determining Nctrl, Pveto, Poffline and Ptrigger for all leptons and nLayers
   - Using this repo to calculate, developing a calculation script
   - Estimate BG for all types and compare with data
- Implementing a way to add leptons SFs to have the correct plots; most probably will rely on passing the .root files with the SFs for the distinct lepton cases

## TODO for BackgroundEstimation

- Charged leptons
   - Implement distinct layer signal regions
   - Implement systematic uncertainties
      - Poffline and Ptrigger uncertainties
- Fake tracks
   - Implement ZtoEE, ZtoMuMu, ZtoEEDisTrkNoD0Cut and ZtoMuMuDisTrkNoD0Cut selections and compare with original (Important! The electron pt cuts are higher in ElectronTagSkim than for ZtoEE. You must run the ZtoEE skim over the full nTuples, *not* the ElectronTagSkim, to get the right events)
   - Determine Nctrl and Pfake for all layers and ZtoLL (this point might be done using the original analysis repo)
      - Estimate BG for all types and compare with data
   - Implement systematic uncertainties
      - Control region uncertainty and transfer factor uncertainties

## Done for StandardAnalysis

- Implemented all selections
   - Basic
   - IsolatedTrack
   - CandidateTrack
   - DisappearingTrack

## Ongoing for StandardAnalysis

- Need to use signal to test and compare selections
   - Used DY so far, and, for 100k events, everything is consistent

## TODO for StandardAnalysis

- Implement lifetime reweighting method
- Implement scale factor corrections
   - Data vs MC ISR
   - PU
   - Data vs MC trigger
   - Missing hits
   - MG5
- Implement method to estimate signal significance

## Done for TriggerAnalysis

- Implemented base trigger selection for data, BG MC and signal MC
- Estimated efficiencies for 2022 (incomplete) and 2023 data
- Compared efficiencies with BG MC 2022 (incomplete) and 2023 data
- Estimated efficiencies in signal MC samples
- Removed the muon pt cut to make plots with the full muon pt range
- Finished estimating uncertainties on problematic datasets of 2022
   - The issue was when checking for muons passing the HLT path, where the name of the filter of path `HLT_IsoMu24` changed from 2022C/D (2022 in MC) to 2022E/F/G (2022EE in MC)

## TODO for TriggerAnalysis

- Implement method to extract trigger efficiency scale factors
- Estimate efficiency over all signal samples and select the best for AN
   - This might not be needed, and all data from 2022 and 2023 should be used
