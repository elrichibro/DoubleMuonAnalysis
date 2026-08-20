# DoubleMuonAnalysis
This program uses RDataFrame to analyze a CMS OpenData DoubleMuone dataset and eventually provide the differential cross section of the Z0 boson decaying in mu+ mu- as a function of p_t, a special angular variable phi* and the rapidity.

This is my final project for the Computing Methods for Experimental Physics exam at University of Pisa, started the 15/08/2026.

Usage:
'''bash



'''

This project will be under active development for the August and September month.

Physics logic:  
-Muon selection:  
    -Muon track reconstruction: (flag)  
        -Standalone-muon tracks  
        -Tracker muon tracks (X)  
        -Global muon tracks (X)  
    -Muon identifications: (flag)  
        -Loose muon ID  
        -Medium muon ID  
        -Tight muon ID (X)  
        -Soft muon ID  
        -High momentum muon ID  
    -Muon isolation: (95% efficiency) (cuts)  
        -PF isolation: DeltaR < 0.4 -> R_iso < 0.15  
        -Track based isolation: DeltaR < 0.3 -> R_iso < 0.05  
-Event selection: (at generator level) (cuts)  
    -p_T > 25 GeV  
    -|eta| < 2.4  
    -60 GeV < m_{l_+l_-} < 120 GeV  