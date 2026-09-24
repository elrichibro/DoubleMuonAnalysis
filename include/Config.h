#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <unordered_map>
#include <vector>
#include <utility>// std::pair

#include <cstdint>// uint.._t

// ------------------------------------------------------------------------------------------------------------------------------------
// General Configuration struct
// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Stores general settup status
struct general_config {
    std::string operation_mode = "";// Operation mode: Selection, Analysis, Acceptance, Resolution
    int verbose = 0;// General verbose flag
    bool visualize = false;// General visualization flag
};

/// @brief Stores input/output file paths
struct io_config {
    std::string tree_data_name = "";// Data tree name
    std::string in_data_file = "";// Input data file path
    std::string tree_mc_name = "";// MonteCarlo tree name
    std::string in_mc_file = "";// Input MonteCarlo file path
    std::string val_file = "";// Validation json file path
};

// ------------------------------------------------------------------------------------------------------------------------------------

// ----------------
// Modality options
// ----------------

/// @brief Acceptance settup struct
struct acceptance_config {
    std::string dataset = "";// Dataset: MC
};

/// @brief Resolution settup struct
struct resolution_config {
    std::string quantity = "";// Resolution quantity
    int gen_bins = 0;// Number of fine bins for reconstruction process
    double min = 0;// Min bin value
    double max = 0;// Max bin value
};

/// @brief Selection mode settup
struct selection_config {
    std::string dataset = "";// Dataset: MC, DATA
    std::string selection_mode = "";// TagAndProbe or RespMatrix
    
    bool save_sel_plots = false;// Flag for saving selection plots
    bool save_sel_data = false;// Flag for saving selection data
    bool visual_sel = false;// Flag for visualization of selection subdataset
    
    std::string o_sel_file_plots = "";// Output file path for plots
    std::string o_sel_file_data = "";// Output file path for data
};

/// @brief Analysis mode settup config
struct analysis_config {
    std::string analysis_mode = "";// Analysis mode: Unfold, Event
    std::string o_fit_file = "";// Analysis output file path
};

/// @brief Bins struct for Unfolded(Analysis) procedure
struct bins {
    int reco_bins = 50;// Bins of Reconstructed MC/DATA quantity
    int gen_bins = 30;// Bins of Generated MC quantity
    double min = 0.0f;// Min bin value
    double max = 200.0f;// Max bin value
    std::string distribution = "";// Distribution option for CreateBins function: linear, log.
    std::vector<double> reco_vec;// Vector of reconstructed bins
    std::vector<double> gen_vec;// Vector of generated bins
};

/// @brief LScan parameters struct
struct l_scan {
    int n_iter = 100;// Iterations number
    float tau_min = 0.0f;// Min tau value
    float tau_max = 100.0f;// Max tau value
};

struct unfold_config {
    std::string unfold_quantity = "";
    bool check_plot = false;
    bool use_custom_bins = false;   
    bool closure_test = false;
    bool bkg_subtraction; 
    
    l_scan scan;
    bins pt_bins;
    bins y_bins;
    bins phis_bins;    
};

struct event_config {
    std::string event_quantity = "";  
    std::string o_event_file = "";
    bool save_fit_plots = false;
};

// ------------
// Filters Info
// ------------

/// @brief Flags for enabling/disablig specifics selections cuts
struct flags_config {
    bool en_kinematics = false;// Enables/Disables the kinematic cut, composed by transverse momentum and pseudorapidity cuts
    bool en_isolation = false;// Enables/Disables isolation selection
    bool en_mass_window = false;// Enables/Disables the fiducial invariant mass region 
    bool en_tight_muon = false;// Enables/Disables the tight muon ID flag
};

/// @brief Stores the values used for the kinematic cuts
struct cuts_config {
    float pt_cut = 0.0f;// Minimum transeverse momentum value
    float eta_cut = 3.0f;// Pseudorapidity acceptance
    float iso_cut = 0.3f;// Isolation cut value
    float mass_min = 0.0f;// Invariant mass Min value
    float mass_max = 200.0f;// Invariant mass Max value
};

// ----------
// Plots Info
// ----------

/// @brief Plot settup struct.
struct plot_config {
    std::string title_axis = "";// Plot title.
    float axis_min = 0.0f;// Plot min axis value.
    float axis_max = 100.f;// Plot max axis value.
    int nbins = 50;// Plot number of bins.
};

// Canvas settup struct.
struct canvas_config {
    int width = 800;
    int height = 600; 
};

/// @brief General config struct.
struct config_struct {
    general_config general;
    io_config io;

    acceptance_config acceptance;
    resolution_config resolution;
    
    selection_config selection;
    analysis_config analysis;

    unfold_config unfold;
    event_config event;
    
    flags_config flag_ES;
    cuts_config cut_ES;

    flags_config flag_RM;
    cuts_config cut_RM;
    
    plot_config pt_plot;
    plot_config eta_plot;
    plot_config mll_plot;
    
    canvas_config canvas;
};

// ------------------------------------------------------------------------------------------------------------------------------------
// Configuration functions
// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Sets the json file information into the system struct.
/// @param value Reference to the configuration struct.
/// @param json_path Json file path.
/// @return Returns 0 on success, -1 on failure.
int Configure(config_struct& value, const std::string& json_path);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Prints on terminal the information contained in the configuration struct.
/// @param value Constant reference to the configuration struct.
void Verbose_config(const config_struct& value);

// ------------------------------------------------------------------------------------------------------------------------------------

/* Alias to define the validation map type.
   Validation architecture maps the run number and the luminosity blocks associated.
   The use of the unordered map ensures an average time of the order of O(1) and a worst case of O(N).
   The luminosity blocks are contained by a dynamic std::vector of (std::pair)s.
*/
using validation_type = std::unordered_map<std::uint32_t, std::vector<std::pair<std::uint16_t, std::uint16_t>>>;

/// @brief Loads the json validation file info into the validation container.
/// @param json_path Validation json file path.
/// @return Returns the validation map.
validation_type Validation_load(const std::string& json_path);

#endif