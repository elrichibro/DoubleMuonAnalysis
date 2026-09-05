#include <ROOT/RDataFrame.hxx>

#include <iostream>
#include <string>

#include <TApplication.h>
#include <TCanvas.h>
#include "TEfficiency.h"

#include "Config.h"
#include "Filters.h"
#include "Utils.h"
#include "Checks.h"
#include "Manager.h"

int main(int argc, char* argv[]) {

    config_struct cfg;
    std::string json_path = "";
    
    int verbose = 0;
    bool visualize = false;
    bool save_plots = false;
    bool save_data = false;
    
    // ---------
    // INTERFACE
    // ---------
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg.rfind(".json") != std::string::npos) {
            json_path = arg;
        } else if (arg == "--verbose" || arg == "-v") {
            verbose = 1;
        } else if ((arg == "--visualize") || (arg == "-vis")) {
            visualize = true;
        } else {
            std::cout << "ERROR: invalid input command, please try again, exinting.\n" << std::endl;
            return 1;
        }
    }

    // ------------------------------------------------------------------------------------------------------------------------------------
    // JSON CONFIGURATION
    // ------------------------------------------------------------------------------------------------------------------------------------

    if (Configure(cfg, json_path) != 0) {
        std::cout << "Configurations fails, check needed, exiting." << std::endl;
        return 1;
    }

    verbose = cfg.general.verbose;
    visualize = cfg.general.visualize;

    std::string dataset_tree = "";
    std::string dataset_file = "";

    if (cfg.general.dataset == "DATA") {
        dataset_tree = cfg.io.tree_data_name;
        dataset_file = cfg.io.in_data_file;
    } else if (cfg.general.dataset == "MC") {
        dataset_tree = cfg.io.tree_mc_name;
        dataset_file = cfg.io.in_mc_file;
    }

    // Verbose JSON configuration
    if (verbose) {
        Verbose_config(cfg);
    }

    const flags_config flags_TP = cfg.flag_TP;
    const cuts_config cuts_TP = cfg.cut_TP;

    const flags_config flags_RM = cfg.flag_TP;
    const cuts_config cuts_RM = cfg.cut_TP;
    
    // ------------------------------------------------------------------------------------------------------------------------------------
    // VISUALIZATION OPTION
    // ------------------------------------------------------------------------------------------------------------------------------------
    
    TApplication* app = nullptr;
    if (visualize) {
        std::cout << "TApplication inializating ..." << std::endl;
        app = new TApplication("app", &argc, argv);
    }

    if (cfg.general.operation_mode.find("Selection") != std::string::npos) {
        try {
            ROOT::EnableImplicitMT();// MultiThread option: ON

            ROOT::RDataFrame data_frame(dataset_tree, dataset_file);
            
            if (verbose){ 
                std::cout << "RDataFrame object created, unpacking " << dataset_tree << " tree from " << 
                dataset_file << " file, starting selection ..." << std::endl;
            }

            /*
            // ------------------------------------------------------------------------------------------------------------------------------------
            // Invariant Mass
            // ------------------------------------------------------------------------------------------------------------------------------------

            auto node_InvMass_bFSR = CalculateInvMass(data_frame, "bFSR", 1);
            auto node_InvMass_aFSR = CalculateInvMass(data_frame, "aFSR", 2);
            
            auto report_bFSR = node_InvMass_bFSR.Report();
            auto report_aFSR = node_InvMass_aFSR.Report();

            auto h_mass_ll_bFSR = node_InvMass_bFSR.Histo1D({"m_ll_bFSR", "Massa invariante dileptoni Before FSR; m_{#mu^{+}#mu^{-}}; Events", 100, 60, 120}, "InvMass_bFSR");        
            auto h_mass_ll_aFSR = node_InvMass_aFSR.Histo1D({"m_ll_aFSR", "Massa invariante dileptoni After FSR; m_{#mu^{+}#mu^{-}}; Events", 100, 60, 120}, "InvMass_aFSR");
            */
            
            // ------------------------------------------------------------------------------------------------------------------------------------
            // RespMatrix
            // ------------------------------------------------------------------------------------------------------------------------------------
            
            ROOT::RDF::RNode node_RM = data_frame;
            
            if (cfg.general.analysis_mode.find("RespMatrix") != std::string::npos) {
                node_RM = node_RM
                    .Define("RespMatrix_mask", [flags_RM, cuts_RM](const ROOT::RVec<float>& pt_gen, const ROOT::RVec<float>& eta_gen, 
                    const ROOT::RVec<float>& pt_rec, const ROOT::RVec<float>& eta_rec, const ROOT::RVec<UChar_t>& rec_flav,
                    const ROOT::RVec<Int_t>& rec_gen_idx, const ROOT::RVec<Int_t>& gen_status, const ROOT::RVec<Int_t>& gen_pdg_id) {
                        
                        MuonKinematics_RM kin{pt_gen, eta_gen, pt_rec, eta_rec};
                        MuonFlags_RM val{rec_flav, rec_gen_idx, gen_status, gen_pdg_id};
                        
                        return CalculateRespMatrix(kin, val, flags_RM, cuts_RM);
                    
                    }, {"GenPart_pt", "GenPart_eta", "Muon_pt", "Muon_eta", "Muon_genPartFlav", "Muon_genPartIdx", "GenPart_status", "GenPart_pdgId"});

                node_RM = node_RM
                    .Define("Gen_Pt", [](const ResultsRespMatrix& res) { return res.pt_gen_RM; }, {"RespMatrix_mask"})
                    .Define("Gen_Eta", [](const ResultsRespMatrix& res) { return res.eta_gen_RM; }, {"RespMatrix_mask"})
                    .Define("Rec_Pt", [](const ResultsRespMatrix& res) { return res.pt_rec_RM; }, {"RespMatrix_mask"})
                    .Define("Rec_Eta", [](const ResultsRespMatrix& res) { return res.eta_rec_RM; }, {"RespMatrix_mask"});
            }
            
            // ------------------------------------------------------------------------------------------------------------------------------------
            // Tag and Probe
            // ------------------------------------------------------------------------------------------------------------------------------------
            
            ROOT::RDF::RNode node_TP = data_frame;
            
            if (cfg.general.analysis_mode.find("TagAndProbe") != std::string::npos) {

                if (cfg.general.dataset == "MC") {    
                    
                    node_TP = node_TP
                        .Define("TP_Result",
                        [flags_TP, cuts_TP](const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi,
                        const ROOT::RVec<float>& mass, const ROOT::RVec<bool>& tag, const ROOT::RVec<bool>& probe, const ROOT::RVec<bool>& glob, 
                        const ROOT::RVec<float>& iso, const ROOT::RVec<UChar_t>& rec_flav, const ROOT::RVec<Int_t>& rec_gen_idx, 
                        const ROOT::RVec<Int_t>& gen_status, const ROOT::RVec<Int_t>& gen_pdg_id, const ROOT::RVec<float> gen_eta, const ROOT::RVec<float> gen_phi) {
                            
                            MuonKinematics_TP kin{pt, eta, phi, mass};
                            MuonFlags_TP flags{tag, probe, glob, iso};
                            MuonFlags_RM val{rec_flav, rec_gen_idx, gen_status, gen_pdg_id};

                            
                            return CalculateTagAndProbe_MC(kin, flags, flags_TP, cuts_TP, val, gen_eta, gen_phi);

                        }, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_tightId", "Muon_isStandalone", "Muon_isGlobal", "Muon_pfRelIso04_all",
                            "Muon_genPartFlav", "Muon_genPartIdx", "GenPart_status", "GenPart_pdgId", "GenPart_eta", "GenPart_phi"});
                
                } else if (cfg.general.dataset == "DATA") {

                    node_TP = node_TP
                        .Define("TP_Result",
                        [flags_TP, cuts_TP](const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi,
                        const ROOT::RVec<float>& mass, const ROOT::RVec<bool>& tag, const ROOT::RVec<bool>& probe, const ROOT::RVec<bool>& glob, const ROOT::RVec<float>& iso) {
                            
                            MuonKinematics_TP kin{pt, eta, phi, mass};
                            MuonFlags_TP flags{tag, probe, glob, iso};
                            
                            return CalculateTagAndProbe_DATA(kin, flags, flags_TP, cuts_TP);

                        }, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_tightId", "Muon_isStandalone", "Muon_isGlobal", "Muon_pfRelIso04_all"});
                }
            }

            node_TP = node_TP
                .Define(cfg.general.dataset + "_Probe_Pt", [](const ResultsTagAndProbe& res) { return res.pt; }, {"TP_Result"})
                .Define(cfg.general.dataset + "_Probe_Eta", [](const ResultsTagAndProbe& res) { return res.eta; }, {"TP_Result"})
                .Define(cfg.general.dataset + "_Mll", [](const ResultsTagAndProbe& res) { return res.mll; }, {"TP_Result"})
                
                //.Define(cfg.general.dataset + "_Tag_Pt", [](const ResultsTagAndProbe& res) { return res.tag_pt_pass; }, {"TP_Result"})
                //.Define(cfg.general.dataset + "_Tag_Eta", [](const ResultsTagAndProbe& res) { return res.tag_eta_pass; }, {"TP_Result"})
                .Define(cfg.general.dataset + "_Mask_Pass", [](const ResultsTagAndProbe& res) { return res.mask_pass; }, {"TP_Result"});
            
            // ------------------------------------------------------------------------------------------------------------------------------------
            // Output Manager
            // ------------------------------------------------------------------------------------------------------------------------------------

            OutputManager manager(cfg);

            if (cfg.general.analysis_mode == "TagAndProbe") {
                manager.BookAnalysis(node_TP, cfg);
            } else if (cfg.general.analysis_mode == "ResponseMatrix") {
                manager.BookAnalysis(node_RM, cfg);
            }

            manager.Run();

            if (visualize && app != nullptr) {
                std::cout << "Initializing visualization ..." << std::endl;
                app->Run();
                
                delete app; 
            } else {
                std::cout << "No visualization booked.\n" << std::endl;
            }

        } catch (const std::exception& except) {
            std::cerr << "Error nature: " << except.what() << std::endl;
            return 1;
        }
    }

    if (cfg.general.operation_mode.find("Analysis") != std::string::npos) {
        try {        
            ROOT::EnableImplicitMT();

            std::string tree = cfg.general.dataset + cfg.general.analysis_mode + "_Tree";
            ROOT::RDataFrame data_frame(tree, cfg.io.o_file_data);
            
            if (verbose){ 
                std::cout << "RDataFrame object created, unpacking " << cfg.general.analysis_mode
                << " tree from " <<  cfg.io.o_file_data << " file, starting analysis ..." << std::endl;
            }
        
            OutputManager Amanager(cfg);

            Amanager.BookAnalysis(data_frame, cfg);

            Amanager.Run();
            
            if (visualize && app != nullptr) {
                std::cout << "Initializing visualization ..." << std::endl;
                app->Run();
                
                delete app; 
            } else {
                std::cout << "No visualization booked.\n" << std::endl;
            }
            
            if (verbose) {
                std::cout << "Division applied" << std::endl;
            }
            
        } catch (const std::exception& except) {
            std::cerr << "Error nature: " << except.what() << std::endl;
            return 1;
        }
    }

    return 0;
}