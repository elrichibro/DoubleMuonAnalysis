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

            value.general.operation_mode = j.value("operation_mode", value.general.operation_mode);
            value.general.verbose = j.value("verbose", value.general.verbose);
            value.general.visualize = j.value("visualize", value.general.visualize);
        }

        if (json_obj.contains("io")) {
            const auto& j = json_obj["io"];

            value.io.tree_data_name = j.value("tree_data_name", value.io.tree_data_name);
            value.io.in_data_file = j.value("in_data_file", value.io.in_data_file);
            value.io.tree_mc_name = j.value("tree_mc_name", value.io.tree_mc_name);
            value.io.in_mc_file = j.value("in_mc_file", value.io.in_mc_file);
            value.io.val_file = j.value("val_file", value.io.val_file);
        }

        if (json_obj.contains("acceptance")) {
            const auto& j = json_obj["acceptance"];
            value.acceptance.dataset = j.value("dataset", value.acceptance.dataset);
        }

        if (json_obj.contains("selection")) {
            const auto& j = json_obj["selection"];

            value.selection.dataset = j.value("dataset", value.selection.dataset);
            value.selection.selection_mode = j.value("selection_mode", value.selection.selection_mode);
            value.selection.save_sel_plots = j.value("save_sel_plots", value.selection.save_sel_plots);
            value.selection.save_sel_data = j.value("save_sel_data", value.selection.save_sel_data);
            value.selection.visual_sel = j.value("visual_sel", value.selection.visual_sel);
            value.selection.o_sel_file_plots = j.value("o_sel_file_plots", value.selection.o_sel_file_plots);
            value.selection.o_sel_file_data = j.value("o_sel_file_data", value.selection.o_sel_file_data); 
        }

        if (json_obj.contains("template")) {
            const auto& j = json_obj["template"];
            
            value.templ.dataset = j.value("dataset", value.templ.dataset);
            value.templ.template_type = j.value("template_type", value.templ.template_type);
            value.templ.bins_settup = j.value("bins_settup", value.templ.bins_settup);
            value.templ.o_template_file_data = j.value("o_template_file_data", value.templ.o_template_file_data);

            if (j.contains("pt_bins")) {
                value.templ.pt_bins = j["pt_bins"].get<std::vector<float>>();
            }
            if (j.contains("eta_bins")) {
                value.templ.eta_bins = j["eta_bins"].get<std::vector<float>>();
            }
            value.templ.mll_bins = j.value("mll_bins", value.templ.mll_bins);
        }

        if (json_obj.contains("analysis")) {
            const auto& j = json_obj["analysis"];

            value.analysis.analysis_mode = j.value("analysis_mode", value.analysis.analysis_mode);            
            value.analysis.o_fit_file = j.value("o_fit_file", value.analysis.o_fit_file);
            value.analysis.bins_settup = j.value("bins_settup", value.analysis.bins_settup);
            value.analysis.pre_fit = j.value("pre_fit", value.analysis.pre_fit);

            value.analysis.sample_pass_data = j.value("sample_pass_data", value.analysis.sample_pass_data);
            value.analysis.sample_pass_mc = j.value("sample_pass_mc", value.analysis.sample_pass_mc);
            value.analysis.sample_fail_data = j.value("sample_fail_data", value.analysis.sample_fail_data);
            value.analysis.sample_fail_mc = j.value("sample_fail_mc", value.analysis.sample_fail_mc);
            
            if (j.contains("params")) {
                const auto& j_p = j["params"];
                value.analysis.params.efficiency = j_p["efficiency"].get<std::vector<double>>();
                value.analysis.params.n_tot = j_p["n_tot"].get<std::vector<double>>();
                value.analysis.params.mu = j_p["mu"].get<std::vector<double>>();
                value.analysis.params.sigma = j_p["sigma"].get<std::vector<double>>();
                value.analysis.params.lambda_pass = j_p["lambda_pass"].get<std::vector<double>>();
                value.analysis.params.lambda_fail = j_p["lambda_fail"].get<std::vector<double>>();
            }
        }

        if (json_obj.contains("unfold")) {
            const auto& j = json_obj["unfold"];

            value.unfold.unfold_quantity = j.value("unfold_quantity", value.unfold.unfold_quantity);
            value.unfold.check_plot = j.value("check_plot", value.unfold.check_plot);

            if (j.contains("l_scan")) {
                const auto& j_l = j["l_scan"];

                value.unfold.scan.n_iter = j_l.value("n_iter", value.unfold.scan.n_iter);
                value.unfold.scan.tau_min = j_l.value("tau_min", value.unfold.scan.tau_min);
                value.unfold.scan.tau_max = j_l.value("tau_max", value.unfold.scan.tau_max);
            }

            value.unfold.use_bins = j.value("use_bins", value.unfold.use_bins);

            if (j.contains("pt_bins")) {
                const auto& j_p = j["pt_bins"];
                
                value.unfold.pt_bins.reco_bins = j_p.value("reco_bins", value.unfold.pt_bins.reco_bins);
                value.unfold.pt_bins.gen_bins = j_p.value("gen_bins", value.unfold.pt_bins.gen_bins);
                value.unfold.pt_bins.min = j_p.value("min", value.unfold.pt_bins.min);
                value.unfold.pt_bins.max = j_p.value("max", value.unfold.pt_bins.max);
                value.unfold.pt_bins.distribution = j_p.value("distribution", value.unfold.pt_bins.distribution);
                value.unfold.pt_bins.split = j_p.value("split", value.unfold.pt_bins.split);
                value.unfold.pt_bins.reco_vec = j_p["reco_vec"].get<std::vector<double>>();
                value.unfold.pt_bins.gen_vec = j_p["gen_vec"].get<std::vector<double>>();
            }

            if (j.contains("y_bins")) {
                const auto& j_y = j["y_bins"];
                
                value.unfold.y_bins.reco_bins = j_y.value("reco_bins", value.unfold.y_bins.reco_bins);
                value.unfold.y_bins.gen_bins = j_y.value("gen_bins", value.unfold.y_bins.gen_bins);
                value.unfold.y_bins.min = j_y.value("min", value.unfold.y_bins.min);
                value.unfold.y_bins.max = j_y.value("max", value.unfold.y_bins.max);
                value.unfold.y_bins.distribution = j_y.value("distribution", value.unfold.y_bins.distribution);
                value.unfold.y_bins.split = j_y.value("split", value.unfold.y_bins.split);
                value.unfold.y_bins.reco_vec = j_y["reco_vec"].get<std::vector<double>>();
                value.unfold.y_bins.gen_vec = j_y["gen_vec"].get<std::vector<double>>();
            }

            if (j.contains("phis_bins")) {
                const auto& j_phis = j["phis_bins"];
                
                value.unfold.phis_bins.reco_bins = j_phis.value("reco_bins", value.unfold.phis_bins.reco_bins);
                value.unfold.phis_bins.gen_bins = j_phis.value("gen_bins", value.unfold.phis_bins.gen_bins);
                value.unfold.phis_bins.min = j_phis.value("min", value.unfold.phis_bins.min);
                value.unfold.phis_bins.max = j_phis.value("max", value.unfold.phis_bins.max);
                value.unfold.phis_bins.distribution = j_phis.value("distribution", value.unfold.phis_bins.distribution);
                value.unfold.phis_bins.split = j_phis.value("split", value.unfold.phis_bins.split);
                value.unfold.phis_bins.reco_vec = j_phis["reco_vec"].get<std::vector<double>>();
                value.unfold.phis_bins.gen_vec = j_phis["gen_vec"].get<std::vector<double>>();
            }
        }

        if (json_obj.contains("flag_ES")) {
            const auto& j = json_obj["flag_ES"];

            value.flag_ES.en_kinematics = j.value("en_kinematics", value.flag_ES.en_kinematics);
            value.flag_ES.en_isolation = j.value("en_isolation", value.flag_ES.en_isolation);
            value.flag_ES.en_mass_window = j.value("en_mass_window", value.flag_ES.en_mass_window);
            value.flag_ES.en_tight_muon = j.value("en_tight_muon", value.flag_ES.en_tight_muon);
        }

        if (json_obj.contains("cut_ES")) {
            const auto& j = json_obj["cut_ES"];

            value.cut_ES.pt_cut = j.value("pt_cut", value.cut_ES.pt_cut);
            value.cut_ES.eta_cut = j.value("eta_cut", value.cut_ES.eta_cut);
            value.cut_ES.iso_cut = j.value("iso_cut", value.cut_ES.iso_cut);
            value.cut_ES.mass_min = j.value("mass_min", value.cut_ES.mass_min);
            value.cut_ES.mass_max = j.value("mass_max", value.cut_ES.mass_max);
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
    std::cout << "    Operation mode: " << value.general.operation_mode << std::endl;
    std::cout << "    Verbose mode: " << value.general.verbose << std::endl;
    std::cout << "    Visualize flag: " << value.general.visualize << std::endl;

    std::cout << "" << std::endl;    
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "I/O settings:" << std::endl;
    std::cout << "    Input data file: " << value.io.in_data_file << std::endl;
    std::cout << "    Data Tree: " << value.io.tree_data_name << std::endl;
    std::cout << "    Input MC file: " << value.io.in_mc_file << std::endl;
    std::cout << "    MC Tree: " << value.io.tree_mc_name << std::endl;
    std::cout << "    Validation file: " << value.io.val_file << std::endl;
    
    std::cout << "" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Acceptance settup:" << std::endl;
    std::cout << "    Dataset: " << value.acceptance.dataset << std::endl;

    std::cout << "" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Selection settup:" << std::endl;
    std::cout << "    Dataset: " << value.selection.dataset << std::endl;
    std::cout << "    Mode: " << value.selection.selection_mode << std::endl;
    std::cout << "    Save plots: " << value.selection.save_sel_plots << std::endl;
    std::cout << "    Save data: " << value.selection.save_sel_data << std::endl;
    std::cout << "    Visualize: " << value.selection.visual_sel << std::endl;
    std::cout << "    Output plots file path: " << value.selection.o_sel_file_plots << std::endl;
    std::cout << "    Output data file path: " << value.selection.o_sel_file_data << std::endl;

    std::cout << "" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Template options: " << std::endl;
    std::cout << "    Dataset: " << value.templ.dataset << std::endl;
    std::cout << "    Template output type: " << value.templ.template_type << std::endl;    
    std::cout << "    Bins settup: " << value.templ.bins_settup << std::endl;
    std::cout << "    Output template file path: " << value.templ.o_template_file_data << std::endl;

    std::cout << "    Pt bins intervals: ";
    for (auto it : value.templ.pt_bins) {
        std::cout << it << " ";
    }
    std::cout << std::endl;
    std::cout << "    Eta bins intervals: ";
    for (auto it : value.templ.eta_bins) {
        std::cout << it << " ";
    }
    std::cout << std::endl;
    std::cout << "    Invariant mass bins number: " << value.templ.mll_bins << std::endl;
    
    std::cout << "" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "" << std::endl;
    
    std::cout << "Analysis settup:" << std::endl;

    std::cout << "" << std::endl;

    std::cout << "    Analysis Mode: " << value.analysis.analysis_mode << std::endl;
    std::cout << "    Fit results file path: " << value.analysis.o_fit_file << std::endl;
    std::cout << "    Bins settup: " << value.analysis.bins_settup << std::endl;
    std::cout << "    Pre-Fit flag option: " << value.analysis.pre_fit << std::endl;

    std::cout << "" << std::endl;

    std::cout << "    Sample Data Pass: " << value.analysis.sample_pass_data << std::endl;
    std::cout << "    Sample MC Pass: " << value.analysis.sample_pass_mc << std::endl;
    std::cout << "    Sample Data Fail: " << value.analysis.sample_fail_data << std::endl;
    std::cout << "    Sample MC Fail: " << value.analysis.sample_fail_mc << std::endl;

    std::cout << "" << std::endl;

    std::cout << "    Fit parameters: " << std::endl;
    std::cout << "        Efficiency: ";
    for (auto it : value.analysis.params.efficiency) {
        std::cout << " / " << it;
    }
    std::cout << " /" << std::endl;
    
    std::cout << "        Signal Events: ";
    for (auto it : value.analysis.params.n_tot) {
        std::cout << " / " << it;
    }
    std::cout << " /" << std::endl;    
    
    std::cout << "        Gaussian mean: ";
    for (auto it : value.analysis.params.mu) {
        std::cout << " / " << it;
    }
    std::cout << " /" << std::endl;    
    
    std::cout << "        Gaussian sigma: ";
    for (auto it : value.analysis.params.sigma) {
        std::cout << " / " << it;
    }
    std::cout << " /" << std::endl;    
    
    std::cout << "        Lambda pass: ";
    for (auto it : value.analysis.params.lambda_pass) {
        std::cout << " / " << it;
    }
    std::cout << " /" << std::endl;    
    
    std::cout << "        Lambda fail: ";
    for (auto it : value.analysis.params.lambda_fail) {
        std::cout << " / " << it;
    }
    std::cout << " /" << std::endl;

    std::cout << "" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Unfold settup: " << std::endl;

    std::cout << "" << std::endl;

    std::cout << "    Unfold Quantity: " << value.unfold.unfold_quantity << std::endl;
    std::cout << "    Check Plots flag: " << value.unfold.check_plot << std::endl;

    std::cout << "" << std::endl;

    std::cout << "    L-Scan: " << std::endl;
    std::cout << "        Iterations number: " << value.unfold.scan.n_iter << std::endl;
    std::cout << "        Tau min value: " << value.unfold.scan.tau_min << std::endl;
    std::cout << "        Tau max value: " << value.unfold.scan.tau_max << std::endl;

    std::cout << "" << std::endl;
    std::cout << "    Use bins option: " << value.unfold.use_bins << std::endl;
    std::cout << "" << std::endl;

    std::cout << "    Pt bins: " << std::endl;
    std::cout << "        Reconstructed bins: " << value.unfold.pt_bins.reco_bins << std::endl;
    std::cout << "        Generated bins: " << value.unfold.pt_bins.gen_bins << std::endl;
    std::cout << "        Min bin: " << value.unfold.pt_bins.min << std::endl;
    std::cout << "        Max bin: " << value.unfold.pt_bins.max << std::endl;
    std::cout << "        Distribution of bins: " << value.unfold.pt_bins.distribution << std::endl;
    std::cout << "        Split option: " << value.unfold.pt_bins.split << std::endl;
    
    std::cout << "        Custom Reco bins: ";
    for (auto it : value.unfold.pt_bins.reco_vec) {
        std::cout << it << " ";
    }
    std::cout << "" << std::endl;

    std::cout << "        Custom Gen bins: ";
    for (auto it : value.unfold.pt_bins.gen_vec) {
        std::cout << it << " ";
    }
    std::cout << "" << std::endl;

    std::cout << "" << std::endl;

    std::cout << "    Rapidity bins: " << std::endl;
    std::cout << "        Reconstructed bins: " << value.unfold.y_bins.reco_bins << std::endl;
    std::cout << "        Generated bins: " << value.unfold.y_bins.gen_bins << std::endl;
    std::cout << "        Min bin: " << value.unfold.y_bins.min << std::endl;
    std::cout << "        Max bin: " << value.unfold.y_bins.max << std::endl;
    std::cout << "        Distribution of bins: " << value.unfold.y_bins.distribution << std::endl;
    std::cout << "        Split option: " << value.unfold.y_bins.split << std::endl;
    std::cout << "        Custom Reco bins: ";
    for (auto it : value.unfold.y_bins.reco_vec) {
        std::cout << it << " ";
    }
    std::cout << "" << std::endl;

    std::cout << "        Custom Gen bins: ";
    for (auto it : value.unfold.y_bins.gen_vec) {
        std::cout << it << " ";
    }
    std::cout << "" << std::endl;

    std::cout << "" << std::endl;

    std::cout << "    Phi Star bins: " << std::endl;
    std::cout << "        Reconstructed bins: " << value.unfold.phis_bins.reco_bins << std::endl;
    std::cout << "        Generated bins: " << value.unfold.phis_bins.gen_bins << std::endl;
    std::cout << "        Min bin: " << value.unfold.phis_bins.min << std::endl;
    std::cout << "        Max bin: " << value.unfold.phis_bins.max << std::endl;
    std::cout << "        Distribution of bins: " << value.unfold.phis_bins.distribution << std::endl;
    std::cout << "        Split option: " << value.unfold.phis_bins.split << std::endl;
    std::cout << "        Custom Reco bins: ";
    for (auto it : value.unfold.phis_bins.reco_vec) {
        std::cout << it << " ";
    }
    std::cout << "" << std::endl;

    std::cout << "        Custom Gen bins: ";
    for (auto it : value.unfold.phis_bins.gen_vec) {
        std::cout << it << " ";
    }
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Event Selection flags:" << std::endl;
    std::cout << "    Kinematics: " << value.flag_ES.en_kinematics << std::endl;
    std::cout << "    Isolation: " << value.flag_ES.en_isolation << std::endl;
    std::cout << "    Mass window: " << value.flag_ES.en_mass_window << std::endl;
    std::cout << "    Tight muon: " << value.flag_ES.en_tight_muon << std::endl;

    std::cout << "" << std::endl;

    std::cout << "Event Selection cuts:" << std::endl;
    std::cout << "    p_T cut: " << value.cut_ES.pt_cut << std::endl;
    std::cout << "    Eta cut: " << value.cut_ES.eta_cut << std::endl;
    std::cout << "    Isolation cut: " << value.cut_ES.iso_cut << std::endl;
    std::cout << "    Max mass: " << value.cut_ES.mass_max << std::endl;
    std::cout << "    Min mass: " << value.cut_ES.mass_min << std::endl;

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
    
    std::cout << "End of Verbose" << std::endl;
    std::cout << "" << std::endl;
}

// ------------------------------------------------------------------------------------------------------------------------------------

validation_type Validation_load(const std::string& json_path) {    
    validation_type validation_map;
    
    // Opening stream
    std::ifstream file(json_path);
    if (!file.is_open()) {
        throw std::runtime_error("Error: cannot open " + json_path);
    }

    // Using JSON lib
    nlohmann::json json_data = nlohmann::json::parse(file);// JSON CORE !!!
    
    // Loop on JSON object
    for (const auto& [run_number, lum_block] : json_data.items()) {
        std::uint32_t run = static_cast<std::uint32_t>(std::stoul(run_number));

        std::vector<std::pair<std::uint16_t, std::uint16_t>> blocks;
        blocks.reserve(lum_block.size());

        // Loop on Luminosity block vector
        for (const auto& iter : lum_block) {
            std::uint16_t start = static_cast<std::uint16_t>(iter[0].get<unsigned int>());
            std::uint16_t end   = static_cast<std::uint16_t>(iter[1].get<unsigned int>());
            blocks.emplace_back(start, end);
        }
        validation_map[run] = std::move(blocks);
    }
    return validation_map;
}