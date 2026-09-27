### DoubleMuonAnalysis
This program uses ROOT's RDataFrame to analyze a CMS OpenData DoubleMuone dataset and eventually provide the differential cross section of the Z0 boson decaying in mu+ mu- as a function of p_t, a special angular variable phi* and the rapidity.

This is my final project for the Computing Methods for Experimental Physics exam at University of Pisa, started the 15/08/2026.

---

## General options
To obtain the differential cross section several steps must be followed to tune the analysis parameters. The program relies on two datasets: reconstructed experiment data and simulated MonteCarlo data. It's important to follow the specific combinations of these two datasets acording to the operation modes.
The workflow is:
- Selection of the 'interest' data.
- Computation of the pre-analysis quantities such as Acceptance and Resolution.
- Analysis procedure: unfold, event selection, cross section final calculus. 

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

# Operation Mode: Selection
The complete dataset (DATA + MC) is at least 150 GB of size. Two input strategies are supported in this code:

1. Download the two individual data packages (1 DATA, 1 MC  2GB each) and store them within the 'data' directory.
2. Use the stream of the OpeData server and use the dataset directly. Analysis times in this case depends on the enthernet connection.

To maximize the computational efficiency of the analysis, the selection mode filters the usless quantities, reducing the dataset size to approx 500 MB. The output path is configurable via the JSON file.

To use Selection operation mode one must fix:

```json
"general": { 
"data_mode": "",
"operation_mode":"Selection",
}
```
- data_mode to 'online' or 'local' string.
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
---
# Operation Mode: Pre-analysis

The following pre-analysis operation modes are available in the framework:

- Acceptance operation mode: Processes the unfiltered Monte Carlo dataset to compute the detector geometrical efficiency. Using generator-level events (Z -> \mu+\mu-), it calculates the fraction of reconstructible events  from all generated ones that pass the matching condition and the fiducial kinematical cuts.

- Resolution Mode: Uses the selected Response Matrix dataset (local or online) to determine the experimental resolution. Using the matched flag, that is true only when the two generated muons match the reconstructed muons and satisfy all kinematic selection. Using fine binning at the generator level, it extracts the relative reconstructed distributions for each bin to study how detector resolution varies(p_t, |y| and phi*).

- Unfolding -> Control Histograms Mode: Computes efficiency, purity, and stability for each bins configuration. This mode is usefull to select a propper binning settup for unfolding procedure.


# Run Commands

cmake ..
make -j(n proc)

./analyse_mc <path/to/config.json> [options]

---

# Analysis procedure
If the output directory is empty follow this steps:
(All the JSON parameters are Case-sensitive)

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