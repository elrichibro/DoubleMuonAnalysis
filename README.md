### DoubleMuonAnalysis
This program uses ROOT's RDataFrame to analyze a CMS OpenData DoubleMuone dataset and eventually provide the differential cross section of the Z0 boson decaying in mu+ mu- as a function of p_t, a special angular variable phi* and the rapidity.

This is my final project for the Computing Methods for Experimental Physics exam at University of Pisa, started the 15/08/2026.

---

## General options
The code support three modes of operation and can run in both: data and MonteCarlo independently:
- Selection: Filters the events of interest creating a smaller subdataset containing the information needed for efficiecies calculus.
- Template: Creates an intermediate data state for better analysis optimization ad for binning separation. The template phase can create both types of data: binned(TH3D) or unbinned(RVecs).
- Analysis: Reads the template input choosing the operation type of datasets(histo(binned) or data(unbinned)). Then starts the Fit process (the fit model is hardcoded) which reads the input parameters values and starts the minimization. The task of this phase is to obtain the scale factors for cross-section calculus. 

---

## Usage

The program is entirely commanded by a JSON configuration file. 
To run the analysis, use the following syntax:

./analyse_mc <path/to/config.json> [options]

# Command options:

- **-v, --verbose** : Enable verbose output (prints event loops progress, debug info).
- **-c, --control** : Enables only the verbose of the configuration settup.
- **-vis, --visualize** : Enable visualization through TApplication.

---

# Run Commands

cmake ..
make -j(n proc)

./analyse_mc <path/to/config.json> [options]

---

# Analysis procedure
If the output directory is empty follow this steps:
(All the JSON parameters are Case-sensitive)

- Data Selection:
  - The **Selection** can be done to **MC** or **DATA** datasets.
  - The Selection must be done with a selection option that can be **TagAndProbe** or **RespoMatrix**.
  - For each selection option there are specific physical cuts. The cuts flas(for implementation) and values are specified in the JSON file.
  - The important JSON parameters are:
    ```json
    "general": { 
    "dataset": "MC",
    "operation_mode":"Selection",
    }
    ```
    - "dataset" : MC or DATA. Applies the selection procedure to this sample.
    - "operation_mode" : Selection -> needed
    ```json
    "selection": {
    "selection_mode": "TagAndProbe",
    "save_sel_plots":true,
    "save_sel_data":true,
    "visual_sel":false,
    "o_sel_file_plots":"../output/risultati_plots.root",
    "o_sel_file_data":"../output/risultati_data.root"
    }
    ```
    - "selection_mode" : **TagAndProbe(TP)** or **RespMatrix(RM)**.
    - "save_sel_data" : true -> For Selection step.
    - "save_sel_plots" or "visual_sel" : optional.
    - "o_sel_file_data": file path of Selection data output.
    ```json
    "flag_TP": {
      "en_kinematics":true,
      "en_mass_window":true
    }

    "cut_TP": {
      "pt_cut":25.0,
      "eta_cut":2.4,
      "iso_cut":0.15,
      "mass_min":60.0,
      "mass_max":120.0
    }
    ```
    - "flag_TP": Are the cuts flags that manage the implementation.
    - "cut_TP": Are the cuts values 

- Template creation:
  - The Template step creates a new subsample from Selection sample. This step its necessary if you want to change the binning of the 3D histogram.
  - One can create binned or unbinned dataset.
  - The important JSON parameters are:
    ```json
    "general": {
      "operation_mode": "Template"
    }
    ```
    - operations_mode : **Template**
    ```json
    "template": {
      "bins_settup":"0_Set_Pt4_Eta4_Mll70",
      "template_type":"HISTO_DATA",
      "o_template_file_data":"../output/risultati_template.root",
      "pt_bins": [25.0, 30.0, 35.0, 40.0, 50.0],
      "eta_bins": [-2.4, -0.9, 0.0, 0.9, 2.4],
      "mll_bins": 70.0
    },
    ```
    - "bins_settup": is the name used to save the data.
    - "template_type": is the binned/unbinned mode. Options are **HISTO**, **DATA** or **HISTO_DATA**.
    - "o_template_file_data": output file path for template subsample.
    - "pt_bins", "eta_bins": are the bins of TH3D and of the unbinned subsamples.
    - "mll_bins": invariant mass bins. Only needed for **HISTO** template mode.

- Analysis:
  - The Template and Analysis steps are strongly connected.
  - The fit phase has several fit options depending on convergence of MINUIT minimization. 
  - The fit model of both: Signal and Background is Hardcoded in the FitFunction(change name..)
  - The user can change dynamicaly the data type(binned/unbinned) ad the initial fit parameters(central value and physical intervals)
  - Important JSON variables:
  ```json
  "general" {
    "operation_mode": "Analysis"
  }
  ```
  - "operation_mode": Analysis is fundamental.
  ```json
  "analysis": {
    "pre_fit":true,
    "o_fit_file":"../output/risultati_fit.root",
    "sample_pass_data":"histo",
    "sample_pass_mc":"histo",
    "sample_fail_data":"histo",
    "sample_fail_mc":"histo",
    "params": {
      "efficiency":[0.9, 0.4, 1.0],
      "n_tot":[0,0,0],
      "mu":[-0.08, -3.0, 3.0],
      "sigma":[0.1, 0.001, 2.0],
      "lambda_pass":[-0.04, -3.0, 0.0],
      "lambda_fail":[-0.04, -3.0, 0.0]
    }
  }
  ```
  - "pre_fit": sets the prefit option that fits only the signal model and obtains the seed values of the mean/sigma (see the model) of the gaussin and use them into the simultaneous fit.
  - "o_fit_file": file path for fit results.
  - "sample_pass_data", "sample_pass_mc", "sample_fail_data", "sample_fail_mc": Are the binned/unbinned options for the samples -> Data(sample: PASS), MonteCarlo(sample: PASS, for FFTConv), Data (sample: FAIL), MonteCarlo (sample: FAIL, for FFTConv).
  - "params": Are the initial values for minimization process (central value, min value, max value).

## JSON file (config.json)

The config.json file controls all the parameters of the analysis, from I/O paths to physics cuts, allowing the modification the analysis without recompiling the project.

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