# Analysis

This has the BackgroundEstimation, StandardAnalysis and TriggerAnalysis selections from the DisappTrks analysis for Run3. It will have all the needed corrections, histograms and possibly trees to make further needed selections. It also saves the EDM files, given that if new variables from the MiniAOD are needed, they can be easily accessed.

## Done for BackgroundEstimation

- Implemented all the leptons tag skims and compared with original
- Implemented ZtoEleProbeTrk selection and compared with original

## TODO for BackgroundEstimation

- Charged leptons
   - Implement Z mass and elec/track charge verification and veto selection in ZtoEleProbetrk
      - Get number OS T&P pairs
      - Get number SS T&P pairs
      - Get number OS T&P pairs passing veto selection
      - Get number SS T&P pairs passing veto selection
   - Implement muon and tau (e and mu) ZtoLepProbeTrk selections and compare with original
   - Implement leptonTagPt55 leptonTagPt55MetTrig selections and compare with original
   - Implement distinct layer signal regions
   - Determine Nctrl, Pveto, Poffline and Ptrigger for all leptons and nLayers (this point might be done using the original analysis repo)
      - Estimate BG for all types and compare with data
   - Implement systematic uncertainties
      - Poffline and Ptrigger uncertainties
- Fake tracks
   - Implement ZtoEE, ZtoMuMu, ZtoEEDisTrkNoD0Cut and ZtoMuMuDisTrkNoD0Cut selections and compare with original (Important! The electron pt cuts are higher in ElectronTagSkim than for ZtoEE. You must run the ZtoEE skim over the full nTuples, *not* the ElectronTagSkim, to get the right events)
   - Determine Nctrl and Pfake for all layers and ZtoLL (this point might be done using the original analysis repo)
      - Estimate BG for all types and compare with data
   - Implement systematic uncertainties
      - Control region uncertainty and transfer factor uncertainties


## Done for StandardAnalysis

## TODO for StandardAnalysis

- Implement all selections
   - Basic
   - IsolatedTrack
   - CandidateTrack
   - DisappearingTrack
- Use signal samples to test and compare selections
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

## TODO for TriggerAnalysis

- Finish estimating uncertainties on problematic datasets of 2022
- Implement method to extract trigger efficiency scale factors
- Estimate uncertainty over all signal samples and select the best for AN
