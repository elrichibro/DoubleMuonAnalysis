#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <unordered_map>
#include <vector>
#include <utility>// std::pair

#include <cstdint>// uint.._t


struct io_config {
    std::string tree_name = "";
    std::string input_file = "";
    std::string validation_file = "";
};

struct flags_config {
    bool en_kinematics = false;
    bool en_isolation = false;
    bool en_opposite_charge = false;
    bool en_mass_window = false;
    bool en_tight_muon = false;
};

struct cuts_config {
    float pt_cut = 25.0f;
    float eta_cut = 2.4f;
    float iso_cut = 0.15f;
    float mass_min = 60.0f;
    float mass_max = 120.0f;
};

struct config_struct {
    io_config io;
    flags_config flag;
    cuts_config cut;
};

int Configure(config_struct& value, std::string json_path);

void Verbose_config(const config_struct& value);

using validation_type = std::unordered_map<std::uint32_t, std::vector<std::pair<std::uint16_t, std::uint16_t>>>;

validation_type Validation_load(const std::string& json_path);

#endif