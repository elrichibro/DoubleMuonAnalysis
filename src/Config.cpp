#include "Config.h"
#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>

// ------------------------------------------------------------------------------------------------------------------------------------

int Configure(config_struct& value, const std::string& json_path) {
    std::ifstream file(json_path);
    
    if (!file.is_open()) {
        std::cout << "ERROR: Cannot open: " << json_path << ", exiting.\n" << std::endl;
        return 1;
    }

    try {
        nlohmann::json json_obj;
        file >> json_obj;

        if (json_obj.contains("general")) {
            const auto& j = json_obj["general"];

            value.general.dataset = j.value("dataset", value.general.dataset);
            value.general.operation_mode = j.value("operation_mode", value.general.operation_mode);
            value.general.analysis_mode = j.value("analysis_mode", value.general.analysis_mode);
            value.general.verbose = j.value("verbose", value.general.verbose);
            value.general.visualize = j.value("visualize", value.general.visualize);
            value.general.save_sel_plots = j.value("save_sel_plots", value.general.save_sel_plots);
            value.general.save_sel_data = j.value("save_sel_data", value.general.save_sel_data);
        }

        if (json_obj.contains("io")) {
            const auto& j = json_obj["io"];

            value.io.tree_data_name = j.value("tree_data_name", value.io.tree_data_name);
            value.io.in_data_file = j.value("in_data_file", value.io.in_data_file);
            value.io.tree_mc_name = j.value("tree_mc_name", value.io.tree_mc_name);
            value.io.in_mc_file = j.value("in_mc_file", value.io.in_mc_file);
            value.io.val_file = j.value("val_file", value.io.val_file);
            value.io.o_file_plots = j.value("o_file_plots", value.io.o_file_plots);
            value.io.o_file_data = j.value("o_file_data", value.io.o_file_data);
        }

        if (json_obj.contains("flag_TP")) {
            const auto& j = json_obj["flag_TP"];

            value.flag_TP.en_kinematics = j.value("en_kinematics", value.flag_TP.en_kinematics);
            value.flag_TP.en_isolation = j.value("en_isolation", value.flag_TP.en_isolation);
            value.flag_TP.en_mass_window = j.value("en_mass_window", value.flag_TP.en_mass_window);
            value.flag_TP.en_tight_muon = j.value("en_tight_muon", value.flag_TP.en_tight_muon);
        }

        if (json_obj.contains("cut_TP")) {
            const auto& j = json_obj["cut_TP"];

            value.cut_TP.pt_cut = j.value("pt_cut", value.cut_TP.pt_cut);
            value.cut_TP.eta_cut = j.value("eta_cut", value.cut_TP.eta_cut);
            value.cut_TP.iso_cut = j.value("iso_cut", value.cut_TP.iso_cut);
            value.cut_TP.mass_min = j.value("mass_min", value.cut_TP.mass_min);
            value.cut_TP.mass_max = j.value("mass_max", value.cut_TP.mass_max);
        }

        if (json_obj.contains("flag_RM")) {
            const auto& j = json_obj["flag_RM"];

            value.flag_RM.en_kinematics = j.value("en_kinematics", value.flag_RM.en_kinematics);
            value.flag_RM.en_isolation = j.value("en_isolation", value.flag_RM.en_isolation);
            value.flag_RM.en_mass_window = j.value("en_mass_window", value.flag_RM.en_mass_window);
            value.flag_RM.en_tight_muon = j.value("en_tight_muon", value.flag_RM.en_tight_muon);
        }

        if (json_obj.contains("cut_RM")) {
            const auto& j = json_obj["cut_RM"];

            value.cut_RM.pt_cut = j.value("pt_cut", value.cut_RM.pt_cut);
            value.cut_RM.eta_cut = j.value("eta_cut", value.cut_RM.eta_cut);
            value.cut_RM.iso_cut = j.value("iso_cut", value.cut_RM.iso_cut);
            value.cut_RM.mass_min = j.value("mass_min", value.cut_RM.mass_min);
            value.cut_RM.mass_max = j.value("mass_max", value.cut_RM.mass_max);
        }

        if (json_obj.contains("plot")) {
            const auto& j = json_obj["plot"];

            if (j.contains("pt")) {
                const auto& j_p = json_obj["plot"]["pt"];
                
                value.pt_plot.title_axis = j_p.value("title_axis", value.pt_plot.title_axis);                
                value.pt_plot.axis_max = j_p.value("axis_max", value.pt_plot.axis_max);
                value.pt_plot.axis_min = j_p.value("axis_min", value.pt_plot.axis_min);
                value.pt_plot.nbins = j_p.value("nbins", value.pt_plot.nbins);
            }

            if (j.contains("eta")) {
                const auto& j_e = json_obj["plot"]["eta"];
                
                value.eta_plot.title_axis = j_e.value("title_axis", value.eta_plot.title_axis);                
                value.eta_plot.axis_max = j_e.value("axis_max", value.eta_plot.axis_max);
                value.eta_plot.axis_min = j_e.value("axis_min", value.eta_plot.axis_min);
                value.eta_plot.nbins = j_e.value("nbins", value.eta_plot.nbins);
            }

            if (j.contains("mll")) {
                const auto& j_m = json_obj["plot"]["mll"];
                
                value.mll_plot.title_axis = j_m.value("title_axis", value.mll_plot.title_axis);                
                value.mll_plot.axis_max = j_m.value("axis_max", value.mll_plot.axis_max);
                value.mll_plot.axis_min = j_m.value("axis_min", value.mll_plot.axis_min);
                value.mll_plot.nbins = j_m.value("nbins", value.mll_plot.nbins);
            }

            if (j.contains("canvas")) {
                const auto& j_c = json_obj["plot"]["canvas"];

                value.canvas.width  = j_c.value("width", value.canvas.width);
                value.canvas.height = j_c.value("height", value.canvas.height);
            }
        }
        if (json_obj.contains("analysis")) {
            const auto& j = json_obj["analysis"];

            if (j.contains("pt_bins")) {
                value.analysis.pt_bins = j["pt_bins"].get<std::vector<float>>();
            }
            if (j.contains("pt_bins")) {
                value.analysis.eta_bins = j["eta_bins"].get<std::vector<float>>();
            }
        }
    } catch (const std::exception& except) {
        std::cout << "ERROR: " << except.what() << std::endl;
        return -1;
    }

    return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------

void Verbose_config(const config_struct& value) {
    std::cout << "Init Verbose -> Loaded Configuration:" << std::endl;
    
    std::cout << "" << std::endl;

    std::cout << "General settings:" << std::endl;
    std::cout << "    Dataset used: " << value.general.dataset << std::endl;
    std::cout << "    Operation mode: " << value.general.operation_mode << std::endl;
    std::cout << "    Analysis mode: " << value.general.analysis_mode << std::endl;
    std::cout << "    Verbose mode: " << value.general.verbose << std::endl;
    std::cout << "    Visualize flag: " << value.general.visualize << std::endl;
    std::cout << "    Save selection plots flag: " << value.general.save_sel_plots << std::endl;
    std::cout << "    Save selection data flag: " << value.general.save_sel_data << std::endl;

    std::cout << "" << std::endl;    
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "I/O settings:" << std::endl;
    std::cout << "    Input data file: " << value.io.in_data_file << std::endl;
    std::cout << "    Data Tree: " << value.io.tree_data_name << std::endl;
    std::cout << "    Input MC file: " << value.io.in_mc_file << std::endl;
    std::cout << "    MC Tree: " << value.io.tree_mc_name << std::endl;
    
    std::cout << "" << std::endl;

    std::cout << "    Validation file: " << value.io.val_file << std::endl;
    std::cout << "    Output file selected plots: " << value.io.o_file_plots << std::endl;
    std::cout << "    Output file selected data: " << value.io.o_file_data << std::endl;

    std::cout << "" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
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
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Plot settings:" << std::endl;

    std::cout << "" << std::endl;

    std::cout << "    Pt plot: " << std::endl;
    std::cout << "        Axis title: " << value.pt_plot.title_axis << std::endl;
    std::cout << "        Axis min: " << value.pt_plot.axis_min << std::endl;
    std::cout << "        Axis max: " << value.pt_plot.axis_max << std::endl;
    std::cout << "        Number of bins: " << value.pt_plot.nbins << std::endl;

    std::cout << "" << std::endl;

    std::cout << "    Eta plot: " << std::endl;
    std::cout << "        Axis title: " << value.eta_plot.title_axis << std::endl;
    std::cout << "        Axis min: " << value.eta_plot.axis_min << std::endl;
    std::cout << "        Axis max: " << value.eta_plot.axis_max << std::endl;
    std::cout << "        Number of bins: " << value.eta_plot.nbins << std::endl;

    std::cout << "" << std::endl;

    std::cout << "    Invariant mass plot: " << std::endl;
    std::cout << "        Axis title: " << value.mll_plot.title_axis << std::endl;
    std::cout << "        Axis min: " << value.mll_plot.axis_min << std::endl;
    std::cout << "        Axis max: " << value.mll_plot.axis_max << std::endl;
    std::cout << "        Number of bins: " << value.mll_plot.nbins << std::endl;

    std::cout << "" << std::endl;

    std::cout << "    Canvas settings: " << std::endl;
    std::cout << "        Width: " << value.canvas.width << std::endl;
    std::cout << "        Height: " << value.canvas.height << std::endl;
  
    std::cout << "" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Analysis settup:" << std::endl;
    std::cout << "    Pt bins intervals: ";
    for (auto it : value.analysis.pt_bins) {
        std::cout << it << " ";
    }

    std::cout << std::endl;

    std::cout << "    Eta bins intervals: ";
    for (auto it : value.analysis.eta_bins) {
        std::cout << it << " ";
    }

    std::cout << "" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "" << std::endl;
    
    std::cout << "End of Verbose" << std::endl;
    std::cout << "" << std::endl;
}

// ------------------------------------------------------------------------------------------------------------------------------------

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