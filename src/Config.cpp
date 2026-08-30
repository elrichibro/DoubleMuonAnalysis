#include "Config.h"
#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>

int Configure(config_struct& value, const std::string& json_path) {
    std::ifstream file(json_path);
    
    if (!file.is_open()) {
        std::cout << "ERROR: Cannot open: " << json_path << ", exiting.\n" << std::endl;
        return 1;
    }

    try {
        nlohmann::json json_obj;
        file >> json_obj;

        value.mode = json_obj.value("mode", value.mode);

        if (json_obj.contains("io")) {
            value.io.tree_data_name = json_obj["io"].value("tree_data_name", value.io.tree_data_name);
            value.io.in_data_file = json_obj["io"].value("in_data_file", value.io.in_data_file);
            value.io.tree_mc_name = json_obj["io"].value("tree_mc_name", value.io.tree_mc_name);
            value.io.in_mc_file = json_obj["io"].value("in_mc_file", value.io.in_mc_file);
            value.io.val_file = json_obj["io"].value("val_file", value.io.val_file);
            value.io.output_file = json_obj["io"].value("output_file", value.io.output_file);
        }

        if (json_obj.contains("flag_TP")) {
            value.flag_TP.en_kinematics = json_obj["flag_TP"].value("en_kinematics", value.flag_TP.en_kinematics);
            value.flag_TP.en_isolation = json_obj["flag_TP"].value("en_isolation", value.flag_TP.en_isolation);
            value.flag_TP.en_mass_window = json_obj["flag_TP"].value("en_mass_window", value.flag_TP.en_mass_window);
            value.flag_TP.en_tight_muon = json_obj["flag_TP"].value("en_tight_muon", value.flag_TP.en_tight_muon);
        }

        if (json_obj.contains("cut_TP")) {
            value.cut_TP.pt_cut = json_obj["cut_TP"].value("pt_cut", value.cut_TP.pt_cut);
            value.cut_TP.eta_cut = json_obj["cut_TP"].value("eta_cut", value.cut_TP.eta_cut);
            value.cut_TP.iso_cut = json_obj["cut_TP"].value("iso_cut", value.cut_TP.iso_cut);
            value.cut_TP.mass_min = json_obj["cut_TP"].value("mass_min", value.cut_TP.mass_min);
            value.cut_TP.mass_max = json_obj["cut_TP"].value("mass_max", value.cut_TP.mass_max);
        }

        if (json_obj.contains("flag_RM")) {
            value.flag_RM.en_kinematics = json_obj["flag_RM"].value("en_kinematics", value.flag_RM.en_kinematics);
            value.flag_RM.en_isolation = json_obj["flag_RM"].value("en_isolation", value.flag_RM.en_isolation);
            value.flag_RM.en_mass_window = json_obj["flag_RM"].value("en_mass_window", value.flag_RM.en_mass_window);
            value.flag_RM.en_tight_muon = json_obj["flag_RM"].value("en_tight_muon", value.flag_RM.en_tight_muon);
        }

        if (json_obj.contains("cut_RM")) {
            value.cut_RM.pt_cut = json_obj["cut_RM"].value("pt_cut", value.cut_RM.pt_cut);
            value.cut_RM.eta_cut = json_obj["cut_RM"].value("eta_cut", value.cut_RM.eta_cut);
            value.cut_RM.iso_cut = json_obj["cut_RM"].value("iso_cut", value.cut_RM.iso_cut);
            value.cut_RM.mass_min = json_obj["cut_RM"].value("mass_min", value.cut_RM.mass_min);
            value.cut_RM.mass_max = json_obj["cut_RM"].value("mass_max", value.cut_RM.mass_max);
        }
    } catch (const std::exception& except) {
        std::cout << "ERROR: " << except.what() << std::endl;
        return -1;
    }

    return 0;
}

void Verbose_config(const config_struct& value) {
    std::cout << "Init Verbose -> Loaded Configuration:" << std::endl;
    
    std::cout << "" << std::endl;
    
    std::cout << "Input data file: " << value.io.in_data_file << std::endl;
    std::cout << "Data Tree: " << value.io.tree_data_name << std::endl;
    std::cout << "Input MC file: " << value.io.in_mc_file << std::endl;
    std::cout << "MC Tree: " << value.io.tree_mc_name << std::endl;
    
    std::cout << "" << std::endl;

    std::cout << "Validation file: " << value.io.val_file << std::endl;
    std::cout << "Output file: " << value.io.output_file << std::endl;

    std::cout << "" << std::endl;

    std::cout << "TagAndProbe flags:" << std::endl;
    std::cout << "    Kinematics: " << value.flag_TP.en_kinematics << std::endl;
    std::cout << "    Isolation: " << value.flag_TP.en_isolation << std::endl;
    std::cout << "    Mass window: " << value.flag_TP.en_mass_window << std::endl;
    std::cout << "    Tight muon: " << value.flag_TP.en_tight_muon << std::endl;

    std::cout << "" << std::endl;

    std::cout << "TagAndProbe cuts:" << std::endl;
    std::cout << "    p_T cut: " << value.cut_TP.pt_cut << std::endl;
    std::cout << "    Eta cut: " << value.cut_TP.eta_cut << std::endl;
    std::cout << "    Isolation cut: " << value.cut_TP.iso_cut << std::endl;
    std::cout << "    Max mass: " << value.cut_TP.mass_max << std::endl;
    std::cout << "    Min mass: " << value.cut_TP.mass_min << std::endl;

    std::cout << "" << std::endl;

    std::cout << "RespMatrix flags:" << std::endl;
    std::cout << "    Kinematics: " << value.flag_RM.en_kinematics << std::endl;
    std::cout << "    Isolation: " << value.flag_RM.en_isolation << std::endl;
    std::cout << "    Mass window: " << value.flag_RM.en_mass_window << std::endl;
    std::cout << "    Tight muon: " << value.flag_RM.en_tight_muon << std::endl;

    std::cout << "" << std::endl;

    std::cout << "RespMatrix cuts:" << std::endl;
    std::cout << "    p_T cut: " << value.cut_RM.pt_cut << std::endl;
    std::cout << "    Eta cut: " << value.cut_RM.eta_cut << std::endl;
    std::cout << "    Isolation cut: " << value.cut_RM.iso_cut << std::endl;
    std::cout << "    Max mass: " << value.cut_RM.mass_max << std::endl;
    std::cout << "    Min mass: " << value.cut_RM.mass_min << std::endl;

    std::cout << "" << std::endl;

    std::cout << "End of Verbose" << std::endl;
}



validation_type Validation_load(const std::string& json_path) {    
    validation_type validation_map;
    
    std::ifstream file(json_path);
    if (!file.is_open()) {
        throw std::runtime_error("Error: cannot open " + json_path);
    }

    nlohmann::json json_data = nlohmann::json::parse(file);// JSON CORE !!!
    
    for (const auto& [run_number, lum_block] : json_data.items()) {
        std::uint32_t run = static_cast<std::uint32_t>(std::stoul(run_number));

        std::vector<std::pair<std::uint16_t, std::uint16_t>> blocks;
        blocks.reserve(lum_block.size());

        for (const auto& iter : lum_block) {
            std::uint16_t start = static_cast<std::uint16_t>(iter[0].get<unsigned int>());
            std::uint16_t end   = static_cast<std::uint16_t>(iter[1].get<unsigned int>());
            blocks.emplace_back(start, end);
        }
        validation_map[run] = std::move(blocks);
    }
    return validation_map;
}