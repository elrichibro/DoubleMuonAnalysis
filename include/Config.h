#ifndef CONFIG_H
#define CONFIG_H

#include <string>

struct io_config {
    std::string tree_name = "Events";
    std::string input_file = "../data/dati0.root";
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

#endif