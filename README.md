### DoubleMuonAnalysis
This program uses ROOT's RDataFrame to analyze a CMS OpenData DoubleMuone dataset and eventually provide the differential cross section of the Z0 boson decaying in mu+ mu- as a function of p_t, a special angular variable phi* and the rapidity.

This is my final project for the Computing Methods for Experimental Physics exam at University of Pisa, started the 15/08/2026.

---

## General options
The code support two modes of operation and can run in both: data and MonteCarlo independently:
- Selection: Filters the events of interest creating a smaller subdataset containing the information needed for efficiecies calculus.
- Analysis: Calculates the efficiencies from data and simulation subdataset. Gives an estimation of differential cross section. 

---

## Usage

The program is entirely commanded by a JSON configuration file. 
To run the analysis, use the following syntax:

./analyse_mc <path/to/config.json> [options]

# Command options:

- ** -v, --verbose ** : Enable verbose output (prints event loops progress, debug info).
- ** -s, --save ** : Save the output (saves booked objects in a .root file).
- ** -vis, --visualize ** : Enable visualization through TApplication.

---

## JSON file (config.json)

The config.json file controls all the parameters of the analysis, from I/O paths to physics cuts, allowing the modification the analysis without recompiling the project.

# Overview:
'''
{
  "mode": "both" 
  ,
  "io": {
    "tree_data_name":"Events",
    "in_data_file":"../data/dati0.root",
    "tree_mc_name":"Events",
    "in_mc_file":"../data/datiMC0.root",
    "val_file":"../data/validation_muon_run.json", 
    "output_file":"../output/risultati.root"
  },
  "flag": {
    "en_kinematics":true,
    "en_isolation":true,
    "en_mass_window":false,
    "en_tight_muon":true
  },
  "cut": {
    "pt_cut":25.0,
    "eta_cut":2.4,
    "iso_cut":0.15,
    "mass_min":60.0,
    "mass_max":120.0
  }
}
'''
- General:
- Input/Output:
- Flags:
- Cuts:
- Plots:
- Analysis:

---

## Physics logic:

- Muon track reconstruction: (flag)
    - Standalone-muon tracks
    - Tracker muon tracks (X)
    - Global muon tracks (X)
- Muon identification: (flag)
    - Loose muon ID
    - Medium muon ID (?)
    - Tight muon ID (X)
    - Soft muon ID
    - High momentum muon ID
- Muon isolation: (95% efficiency) (cuts)
    - PF isolation: Delta R < 0.4 -> R_iso < 0.15
    - Track based isolation: Delta R < 0.3 -> R_iso < 0.05

---

## Event Selection -> Fiducial Region

- p_T > 25 GeV
- |eta| < 2.4
- 60 GeV < m_{mu+mu-} < 120 GeV


This project will be under active development for the August and September months.