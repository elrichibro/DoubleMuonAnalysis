#include "Config.h"
#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>

int Configure(config_struct& value, std::string json_path) {
    std::ifstream file(json_path);
    
    if (!file.is_open()) {
        std::cout << "ERROR: Cannot open: " << json_path << ", exiting.\n" << std::endl;
        return 1;
    }

    try {
        nlohmann::json j;
        file >> j;

        if (j.contains("io")) {
            if (j["io"].contains("tree_name")) {
                value.io.tree_name = j["io"]["tree_name"];
            }
            if (j["io"].contains("input_file")) {
                value.io.input_file = j["io"]["input_file"];
            }
        }

        if (j.contains("flag")) {
            if (j["flag"].contains("en_kinematics")) {
                value.flag.en_kinematics = j["flag"]["en_kinematics"];
            }
            if (j["flag"].contains("en_isolation")) {
                value.flag.en_isolation = j["flag"]["en_isolation"];
            }
            if (j["flag"].contains("en_opposite_charge")) {
                value.flag.en_opposite_charge = j["flag"]["en_opposite_charge"];
            }
            if (j["flag"].contains("en_mass_window")) {
                value.flag.en_mass_window = j["flag"]["en_mass_window"];
            }
            if (j["flag"].contains("en_tight_muon")) {
                value.flag.en_tight_muon = j["flag"]["en_tight_muon"];
            }
        }

        if (j.contains("cut")) {
            if (j["cut"].contains("pt_cut")) {
                value.cut.pt_cut = j["cut"]["pt_cut"];
            }
            if (j["cut"].contains("eta_cut")) {
                value.cut.eta_cut = j["cut"]["eta_cut"];
            }
            if (j["cut"].contains("iso_cut")) {
                value.cut.iso_cut = j["cut"]["iso_cut"];
            }
            if (j["cut"].contains("mass_min")) {
                value.cut.mass_min = j["cut"]["mass_min"];
            }
            if (j["cut"].contains("mass_max")) {
                value.cut.mass_max = j["cut"]["mass_max"];
            }
        }
    } catch (const std::exception& except) {
        std::cout << "ERROR: " << except.what() << std::endl;
    }

    return 0;
}

void Verbose_config(const config_struct& value) {
    std::cout << "Init Verbose -> Loaded Configuration:" << std::endl;
    std::cout << "Input file: " << value.io.input_file << " | Tree: " << value.io.tree_name << std::endl;
    
    std::cout << "Flags:" << std::endl;
    std::cout << "    Kinematics: " << value.flag.en_kinematics << std::endl;
    std::cout << "    Isolation: " << value.flag.en_isolation << std::endl;
    std::cout << "    Opposite charge: " << value.flag.en_opposite_charge << std::endl;
    std::cout << "    Mass window: " << value.flag.en_mass_window << std::endl;
    std::cout << "    Tight muon: " << value.flag.en_tight_muon << std::endl;

    
    std::cout << "Cuts:" << std::endl;
    std::cout << "    p_T cut: " << value.cut.pt_cut << std::endl;
    std::cout << "    Eta cut: " << value.cut.eta_cut << std::endl;
    std::cout << "    Isolation cut: " << value.cut.iso_cut << std::endl;
    std::cout << "    Max mass: " << value.cut.mass_max << std::endl;
    std::cout << "    Min mass: " << value.cut.mass_min << std::endl;

    std::cout << "End of Verbose" << std::endl;
}