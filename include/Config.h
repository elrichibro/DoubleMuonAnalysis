#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <unordered_map>
#include <vector>
#include <utility>// std::pair

#include <cstdint>// uint.._t


/// @brief Stores input/output file paths.
struct io_config {
    std::string tree_data_name = "";// Data tree name.
    std::string in_data_file = "";// Input data file path.
    std::string tree_mc_name = "";// MonteCarlo tree name.
    std::string in_mc_file = "";// Input MonteCarlo file path.
    std::string val_file = "";// Validation json file path.
    std::string output_file = "";// Output file path.
};

/// @brief Flags for enabling/disablig specifics selections cuts.
struct flags_config {
    bool en_kinematics = false;// Enables/Disables the kinematic cut, composed by transverse momentum and pseudorapidity cuts.  
    bool en_isolation = false;// Enables/Disables isolation selection.
    bool en_mass_window = false;// Enables/Disables the fiducial invariant mass region. 
    bool en_tight_muon = false;// Enables/Disables the tight muon ID flag
};

/// @brief Stores the values used for the kinematic cuts.
struct cuts_config {
    float pt_cut = 0.0f;// Minimum transeverse momentum value.
    float eta_cut = 3.0f;// Pseudorapidity acceptance.
    float iso_cut = 0.3f;// Isolation cut value.
    float mass_min = 0.0f;// Invariant mass Min value.
    float mass_max = 200.0f;// Invariant mass Max value
};

/// @brief Stores the other data structs.
struct config_struct {
    std::string mode = "data";// Analysis Operation Mode.
    io_config io;
    flags_config flag_TP;
    cuts_config cut_TP;
    flags_config flag_RM;
    cuts_config cut_RM;
};

/**
 * @brief Sets the json file information into the system struct.
 * @param value Reference to the configuration struct. 
 * @param json_path Json file path.
 * @return Returns 0 on success, -1 on failure.
*/
int Configure(config_struct& value, const std::string& json_path);

/**
 * @brief Prints on terminal the information contained in the struct.
 * @param value Constant reference to the configuration struct.
*/
void Verbose_config(const config_struct& value);

/* Alias to define the validation map type.
   Validation architecture should map the run number and the luminosity blocks associated.
   The use of the unordered map ensures an average time in the order of O(1) and a worst case of O(N).
   The luminosity blocks are contained by a dynamic std::vector of std::pair s.
*/
using validation_type = std::unordered_map<std::uint32_t, std::vector<std::pair<std::uint16_t, std::uint16_t>>>;

/**
 * @brief Loads the json validation file info into the validation container.
 * @param json_path Validation json file path.
 * @return Returns the validation map.
*/
validation_type Validation_load(const std::string& json_path);

#endif