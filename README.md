# DoubleMuonAnalysis
This program uses ROOT's RDataFrame to analyze a CMS OpenData DoubleMuone dataset and eventually provide the differential cross section of the Z0 boson decaying in mu+ mu- as a function of p_t, a special angular variable phi* and the rapidity.

This is my final project for the Computing Methods for Experimental Physics exam at University of Pisa, started the 15/08/2026.

---

# General options
To obtain the differential cross section several steps must be followed to tune the analysis parameters. The program works on two datasets: reconstructed experiment data and simulated MonteCarlo data. It's important to follow the specific combinations of these two datasets acording to the operation modes.
The workflow is:
- Selection of the 'interest' data.
- Computation of the pre-analysis quantities such as Acceptance and Resolution.
- Analysis procedure: unfold, event selection, cross section final calculus. 

---

# Prerequisites and Dependencies
- **C++17 Compiler**
- **CMake**: minimum version required 3.14
- **ROOT**: with the following components:
  - *ROOTDataFrame*
  - *Core*
  - *Hist*
  - *Tree*
  - *RIO*
  - *MathCore*
  - *RooFit*
  - *RooFitCore*
  - *Minuit*
  - *Unfold*
- **nlohmann_json**
- **Doxygen**: only for documentation
- **GoogleTest**: only for testing

---

# Installation and usage
The user must have 'git' packet in his pc. The first step is to copy the repo into the user device.

```bash
$ git clone https://github.com/elrichibro/DoubleMuonAnalysis.git
$ cd DoubleMuonAnalysis
$ mkdir build 
$ cd build
$ cmake ..
$ make -j$(nproc)
```

Then the following command starts the program with the JSON file configuration setup:

```bash
./analyze_mc ../config.json $[options]
```

## Command options:

- **-v, --verbose** : Enable verbose output (prints event loops progress, debug info).
- **-c, --control** : Enables only the verbose of the configuration setup.
- **-vis, --visualize** : Enable visualization through TApplication.

# Data setup:
To run the program, the required datasets must be linked (online data mode) throght a txt file or locally downloaded (local data mode) from the server.

## Online Mode (recommended)
This mode streams data directly from the CERN OpenData servers via network protocols, avoiding large local downloads. Instead of downloading full ROOT files, download the file index (`.txt`):

- **Experiment data**:
  - Go to https://opendata.cern.ch/record/30555, click on **Download index** and save the **.txt** format file.
  - Download the validation JSON file (11.4 KiB) from https://opendata.cern.ch/record/14220.

- **MonteCarlo data**:
  - Go to https://opendata.cern.ch/record/35671, click on **Download index** and save the **.txt** format file. 

Then, on the project directory move the **.txt** files into the data directory and link them throught the JSON configuration file. Example:

```json
  "io": {
    "tree_data_name":"Events",
    "in_data_file":"",
    "in_online_data_file":"../data/CMS_Run2016H_DoubleMuon_NANOAOD_UL2016_MiniAODv2_NanoAODv9-v1_2510000_file_index.txt",
    "tree_mc_name":"Events",
    "in_mc_file":"",
    "in_online_mc_file":"../data/CMS_mc_RunIISummer20UL16NanoAODv9_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_NANOAODSIM_106X_mcRun2_asymptotic_v17-v1_40000_file_index.txt",
    "val_file":"../data/validation_muon_run.json"
  }
```
---
## Local Mode
Use this option if you prefer to store ROOT files directly on your local pc.

- **Experiment data**:
  - Go to https://opendata.cern.ch/record/30555, open the **List of files** and download one **.root** dataset file.
  - Download the luminosity validation JSON file (11.4 KiB) from https://opendata.cern.ch/record/14220.

- **MonteCarlo data**:  
  - Go to https://opendata.cern.ch/record/35671, open the **List of files** and download one **.root** dataset file.

Then, in the project directory move the downloaded **.root** files into the data directory and link them throught the JSON configuration file.Example:

```json
  "io": {
    "tree_data_name":"Events",
    "in_data_file":"../data/dati0.root",
    "in_online_data_file":"",
    "tree_mc_name":"Events",
    "in_mc_file":"../data/dati0MC.root",
    "in_online_mc_file":"",
    "val_file":"../data/validation_muon_run.json"
  }
```
---

# Operation Mode: Selection
The complete dataset (DATA + MC) is at least 150 GB of size. Two input strategies are supported in this code:

1. Download two individual data packages (1 DATA, 1 MC  2GB each) and store them within the 'data' directory.
2. Use the stream of the OpeData server and use the entire dataset directly. Analysis times in this case depends on the enthernet connection. (recommended)

To maximize the computational efficiency of the analysis, the selection mode filters the usless quantities, reducing the dataset size to approx 500 MB. The output path is configurable via the JSON file.

To use Selection operation mode one must fix:

```json
"general": { 
"data_mode": "",
"operation_mode":"Selection",
}
```

- **data_mode**: 'online' or 'local'.
- **operation_mode**: 'Selection' keyword is required.


```json
"selection": {
  "dataset": "",
  "selection_mode": "",
  "save_sel_data":false,
  "save_sel_plots":false,
  "visual_sel":false,
  "o_sel_file_plots":"../output/selection_plots.root",
  "o_sel_file_data":"../output/selection_data.root"
}
```

- **dataset**: 'MC' (RespMatrix/Event selection mode) and 'DATA' (Event selection mode).
- **selection_mode**: 'RespMatrix' or 'Event'
  - **RespMatrix** builds the quantities needed for Response matrix calculus.
  - **Event** selects transverse momentum, rapidity, phi* and invariant mass of the reconstructed Z0.
- **save_sel_data**: saves selected data (Important).
- **save_sel_plots**: saves plots.
- **visual_sel**: books the visualization of saved plots.
- **o_sel_file_plots** and **o_sel_file_data**: are the output file paths.

The following section is used to tune the selected data.

```json
  "flag_RM": {
    "en_kinematics":true,
    "en_isolation":false,
    "en_mass_window":false,
    "en_tight_muon":false
  },
  "cut_RM": {
    "pt_cut":25.0,
    "eta_cut":2.5,
    "mass_min":60.0,
    "mass_max":120.0
  }
```

An example can be:

```json
"general": { 
  "data_mode": "online",
  "operation_mode":"Selection",
},
"selection": {
  "dataset": "MC",
  "selection_mode": "RespMatrix",
  "save_sel_data":true,
  "save_sel_plots":false,
  "visual_sel":false,
  "o_sel_file_plots":"../output/selection_plots.root",
  "o_sel_file_data":"../output/selection_data.root"
}
```

This configuration, as explained, opens a stream to the OpenData server and filters the entire MonteCarlo dataset reducing its size, saving it into '../output/selection_data.root' in a TTree (ttree name hardcoded). 

---

# Operation Mode: Pre-analysis

The following pre-analysis operation modes are available in the framework:

---

## Acceptance operation mode:
Processes the unfiltered Monte Carlo dataset to compute the detector geometrical efficiency. Using generator-level events (Z -> \mu+\mu-), it calculates the fraction of reconstructible events  from all generated ones that pass the matching condition and the fiducial kinematical cuts.

  ```json
  "general": { 
    "data_mode":"",
    "operation_mode":"Acceptance"
  }
  ```

  - **data_mode** : 'online' or 'local' is avaible
  - **operation_mode**: 'Acceptance' is required.

  ```json
  "acceptance":{
    "dataset":""
  }
  ```

  - **dataset**: 'MC' is the only option avaible.

  An example can be:
  
  ```json
  "general": { 
    "data_mode":"online",
    "operation_mode":"Acceptance"
  },
  "acceptance":{
    "dataset":"MC"
  }
  ```

---

## Resolution operation mode: 
Uses the selected Response Matrix dataset (local or online) to determine the experimental resolution. Using the matched flag, that is true only when the two generated muons match the reconstructed muons and satisfy all kinematic selection. Using fine binning at the generator level, it extracts the relative reconstructed distributions for each bin to study how detector resolution varies(p_t, |y| and phi*).
  
  ```json
  "general": { 
    "data_mode":"",
    "operation_mode":"Resolution"
  }
  ```

  - **data_mode** : 'online' or 'local' is avaible
  - **operation_mode**: 'Resolution' is required.

  ```json
  "resolution":{
    "quantity":"phis",
    "gen_bins":10,
    "min":0.1,
    "max":1
  }
  ```

  - **quantity**: 'pt', 'y' and 'phis' quantities are avaible.
  - **gen_bins**: Is the number of generated bins for resolution study. The value must be > 1.
  - **min**: Minimum bin value. Must be less than 'max'.
  - **max**: Maximum generated bin value. Must be grater then 'min'.

  An exaple:
  
  ```json
  "general": { 
    "data_mode":"online",
    "operation_mode":"Resolution"
  },
  "resolution":{
    "quantity":"pt",
    "gen_bins":1000,
    "min":0.1,
    "max":200
  }
  ```
  This setup creates a generated grid of 1000 bins and plots the Standard deviation as a function of the generated central value.

---

## Unfolding -> Control Histograms Mode:
Computes efficiency, purity, and stability for each bins configuration. This mode is usefull to select a propper binning setup for unfolding procedure.

  ```json
    "general": { 
      "data_mode":"",
      "operation_mode":"Analysis"
    }
  ```
  
  - **data_mode** : 'online' or 'local' is avaible
  - **operation_mode**: 'Analysis' is required.

  ```json
  "analysis": {
    "analysis_mode":"Unfold",
    "o_fit_file":"../output/risultati_fit.root"
  }
  ```
  
  - **analysis_mode** : 'Unfold' must be setted.
  - **o_fit_file**: output analysis file.
  
  ```json
  "unfold":{
    "unfold_quantity": "pt",
    "check_plot":false,
    "use_custom_bins":true,
  }
  ```

  - **unfold_quantity** : 'pt', 'y' and 'phis' quantities are avaible.
  - **check_plot**: 'true' enables the control plots option.
  - **use_custom_bins**: 'true' uses the custom bins settable directly by JSON file input. And 'false' sets the bins from CreateBins method.
  

---

## Bins
The bins are settable from:
  
```json
"pt_bins":{
  "reco_bins":40,
  "gen_bins":15,
  "min":0.1,
  "max":100,
  "distribution":"log",
  "reco_vec": [1, 2, 3],
  "gen_vec": [1.5, 2.5]
},
```

Non custom setup:

- **reco_bins**: Sets the reconstructed bins number.
- **gen_bins**: Sets the generated bins number.
- **min**: Sets the minimum bin value.
- **max**: Sets the maximum bin value.
- **distribution**: Distribution of the bins along the axis.

Custom setup:

- **reco_vec**: Bins vector of reconstructed quantity.
- **gen_vec**: Bins vector of generated quantity.

---

# Operation mode: Analysis
The Analysis mode offers two sub-modes: **Unfold** and **Event**.

- **Unfold**: Calculates the response matrix, and unfold the input signal yield to extract the final differential cross section.
- **Event**: Performs signal fit on the data to obtain signal yields before feeding it into the unfolding procedure.

```json
  "analysis": {
    "analysis_mode":"",
    "o_fit_file":"../output/risultati_fit.root"
  }
```
- **analysis_mode**: 'Unfold' or 'Event'
- **o_fit_file**: analysis output file path. (currently not used)

## Event:
This options is used to check the fit procedure, here the background is subtrancted from total yield. 

Then, the signal yield is used to build the input histogram tofor the unfolding density matrix, the initializing the 
unfolding procedure. 

The fit procedure is based on RooFit workframe and this steps can be slow depending on the minimization(hardcoded) options.

```json
  "event": {
    "event_quantity":"",
    "o_event_file":"../output/event_fit.root",
    "save_fit_plots":false
  }
```

- **event_quantity**: Kinematical quantity to fit: 'pt', 'y' or 'phis'.
- **o_event_file**: Output file path.
- **save_fit_plots**: Flag that enables saving plots. Should be set to 'false' if unfold analysis option is used.

## Unfold: 
This step is crucial for the analysis results. A good pre-analysis phase for bining tuning should me done. 
```json
"unfold":{
  "unfold_quantity": "",
  "check_plot":false,
  "use_custom_bins":true,
  "closure_test":false,
  "bkg_subtraction":false,
}
```

- **unfold_quantity**: Unfolded quantity, can be 'pt', 'y' or 'phis'.
- **check_plot**: Only used in pre-analysis step.
- **use_custom_bins**: Enables custom binning (recommended).
- **closure_test**: Enables the closure test that fix, as unfold input, the MonteCarlo reconstructed muon signal (used for Response matrix generation). So the perfect match between MC truth-level and Unfolding output is expected.
- **bkg_subtraction**: Enables the bkg subtraction throught fit procedure.

### Regularization
The LScan method is used withing the unfolding procedure to find the best regularization parameter value(tau). The setup options for the minimization procedure are:

```json
"l_scan":{
  "n_iter": 100,
  "tau_min":0.0,
  "tau_max":0.0
}
```

- **n_iter**: Number of minimization iterations.
- **tau_min**: Minimum tau value. 
- **tau_max**: Maximum tau value.

If tau_min == tau_max then the algorithm sets the optimal search range. 

An example of final unfolding setup is:
```json
"analysis": {
  "analysis_mode":"Unfold",
  "o_fit_file":"../output/risultati_fit.root"
},

"unfold":{
  "unfold_quantity": "y",
  "check_plot":false,
  "use_custom_bins":true,
  "closure_test":false,
  "bkg_subtraction":true,

  "l_scan":{
    "n_iter": 100,
    "tau_min":0.0,
    "tau_max":0.0
  },
  
  "y_bins":{
    "reco_bins":30,
    "gen_bins":15,
    "min":0.001,
    "max":2.4,
    "distribution":"linear",
    "reco_vec":[0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4],
    "gen_vec":[0.001, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4]
  }
}
```
---

This project will be under active development for the August and September months.