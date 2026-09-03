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

    //save_plots = cfg.general.save_sel_plots;
    //save_data = cfg.general.save_sel_data;

    std::string dataset_tree = "";
    std::string dataset_file = "";

    if (cfg.general.dataset == "data") {
        std::string dataset_tree = cfg.io.tree_data_name;
        std::string dataset_file = cfg.io.in_data_file;
    } else if (cfg.general.dataset == "mc") {
        std::string dataset_tree = cfg.io.tree_mc_name;
        std::string dataset_file = cfg.io.in_mc_file;
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

            ROOT::RDataFrame data_frame(cfg.io.tree_mc_name, cfg.io.in_mc_file);
            
            if (verbose){ 
                std::cout << "RDataFrame object created, unpacking " << cfg.io.tree_mc_name 
                << " tree from " << cfg.io.in_mc_file << " file, starting selection ..." << std::endl;
            }

            // ------------------------------------------------------------------------------------------------------------------------------------
            // Invariant Mass
            // ------------------------------------------------------------------------------------------------------------------------------------

            auto node_InvMass_bFSR = CalculateInvMass(data_frame, "bFSR", 1);
            auto node_InvMass_aFSR = CalculateInvMass(data_frame, "aFSR", 2);
            
            auto report_bFSR = node_InvMass_bFSR.Report();
            auto report_aFSR = node_InvMass_aFSR.Report();

            auto h_mass_ll_bFSR = node_InvMass_bFSR.Histo1D({"m_ll_bFSR", "Massa invariante dileptoni Before FSR; m_{#mu^{+}#mu^{-}}; Events", 100, 60, 120}, "InvMass_bFSR");        
            auto h_mass_ll_aFSR = node_InvMass_aFSR.Histo1D({"m_ll_aFSR", "Massa invariante dileptoni After FSR; m_{#mu^{+}#mu^{-}}; Events", 100, 60, 120}, "InvMass_aFSR");

            // ------------------------------------------------------------------------------------------------------------------------------------
            // Efficiency
            // ------------------------------------------------------------------------------------------------------------------------------------

            ROOT::RDF::RNode node_RM = data_frame.Define("RespMatrix_mask",
                [flags_RM, cuts_RM](const ROOT::RVec<float>& pt_gen, const ROOT::RVec<float>& eta_gen, const ROOT::RVec<float>& pt_rec, const ROOT::RVec<float>& eta_rec,
                const ROOT::RVec<UChar_t>& rec_flav, const ROOT::RVec<Int_t>& rec_gen_idx, const ROOT::RVec<Int_t>& gen_status, const ROOT::RVec<Int_t>& gen_pdg_id) 
                {
                    MuonKinematics_RM kin{pt_gen, eta_gen, pt_rec, eta_rec};
                    MuonFlags_RM val{rec_flav, rec_gen_idx, gen_status, gen_pdg_id};
                    
                    return CalculateRespMatrix(kin, val, flags_RM, cuts_RM);
                }, {"GenPart_pt", "GenPart_eta", "Muon_pt", "Muon_eta", "Muon_genPartFlav", "Muon_genPartIdx", "GenPart_status", "GenPart_pdgId"});

            node_RM = node_RM
                .Define("Gen_Pt", [](const ResultsRespMatrix& res) { return res.pt_gen_RM; }, {"RespMatrix_mask"})
                .Define("Gen_Eta", [](const ResultsRespMatrix& res) { return res.eta_gen_RM; }, {"RespMatrix_mask"})
                .Define("Rec_Pt", [](const ResultsRespMatrix& res) { return res.pt_rec_RM; }, {"RespMatrix_mask"})
                .Define("Rec_Eta", [](const ResultsRespMatrix& res) { return res.eta_rec_RM; }, {"RespMatrix_mask"});

            /*
            auto node_GEN_event = node_InvMass_aFSR
                .Filter([](const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, float mass) {
                    return (ROOT::VecOps::All(pt > 25.0f) && ROOT::VecOps::All(abs(eta) < 2.4f) && (mass > 60.0f) && (mass < 120.0f));
                }, {"good_Pt_aFSR", "good_Eta_aFSR", "m_ll_aFSR"}, "Generator Event");

            auto h_pt_gen = node_GEN_event.Histo1D({"h_pt_gen", "Generator; p_{T} [GeV]; Entries", 100, 25, 100}, "good_Pt_aFSR");
            auto h_eta_gen = node_GEN_event.Histo1D({"h_eta_gen", "Generator; #eta; Entries", 100, -2.4, 2.4}, "good_Eta_aFSR");

            ROOT::RDF::RNode node_REC_event = node_recMC(node_GEN_event);
            
            node_REC_event = node_REC_event
                .Filter([](const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, float mass) {
                    return (ROOT::VecOps::All(pt > 25.0f) && ROOT::VecOps::All(abs(eta) < 2.4f) && (mass > 60.0f) && (mass < 120.0f));
                }, {"good_Muon_pt", "good_Muon_eta", "Muon_inv_mass"}, "4. Reconstructed Event");
                
            auto h_pt_rec = node_REC_event.Histo1D({"h_pt_rec", "Reconstructed; p_{T} [GeV]; Entries", 100, 25, 100}, "good_Pt_aFSR");
            auto h_eta_rec = node_REC_event.Histo1D({"h_eta_rec", "Reconstructed; #eta; Entries", 100, -2.4, 2.4}, "good_Eta_aFSR");


            auto eff_pt = std::make_unique<TEfficiency>(*h_pt_rec, *h_pt_gen);
            auto eff_eta = std::make_unique<TEfficiency>(*h_eta_rec, *h_eta_gen);

            eff_pt->SetTitle("p_{T} Efficiency; p_{T} [GeV]; #epsilon_{p_{T}}");
            eff_eta->SetTitle("#eta Efficiency; #eta; #epsilon_{#eta}");

            TCanvas c1("c1", "Efficiency Canvas", 800, 600);
            c1.SetGrid();
            eff_pt->Draw("AP");

            TCanvas c2("c2", "Efficiency Canvas", 800, 600);
            c2.SetGrid();
            eff_eta->Draw("AP");
            */
            
            // ------------------------------------------------------------------------------------------------------------------------------------
            // Tag and Probe
            // ------------------------------------------------------------------------------------------------------------------------------------

            
            ROOT::RDF::RNode node_TP = data_frame.Define("TP_Result",
                [flags_TP, cuts_TP](const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi, const ROOT::RVec<float>& mass,
                const ROOT::RVec<bool>& tag, const ROOT::RVec<bool>& probe, const ROOT::RVec<bool>& glob, const ROOT::RVec<float>& iso) 
                {
                    MuonKinematics_TP kin{pt, eta, phi, mass};
                    MuonFlags_TP flags{tag, probe, glob, iso};
                    
                    return CalculateTagAndProbe(kin, flags, flags_TP, cuts_TP);
                }, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_tightId", "Muon_isStandalone", "Muon_isGlobal", "Muon_pfRelIso04_all"});

            node_TP = node_TP
                .Define("Probe_Pt_All", [](const ResultsTagAndProbe& res) { return res.pt_all; }, {"TP_Result"})
                .Define("Probe_Pt_Pass", [](const ResultsTagAndProbe& res) { return res.pt_pass; }, {"TP_Result"})
                .Define("Probe_Eta_All", [](const ResultsTagAndProbe& res) { return res.eta_all; }, {"TP_Result"})
                .Define("Probe_Eta_Pass", [](const ResultsTagAndProbe& res) { return res.eta_pass; }, {"TP_Result"})
                .Define("Probe_Mll_All", [](const ResultsTagAndProbe& res) { return res.mll_all; }, {"TP_Result"})
                .Define("Probe_Mll_Pass", [](const ResultsTagAndProbe& res) { return res.mll_pass; }, {"TP_Result"})
                .Define("Tag_Pt_All", [](const ResultsTagAndProbe& res) { return res.tag_pt_all; }, {"TP_Result"})
                .Define("Tag_Pt_Pass", [](const ResultsTagAndProbe& res) { return res.tag_pt_pass; }, {"TP_Result"})
                .Define("Tag_Eta_All", [](const ResultsTagAndProbe& res) { return res.tag_eta_all; }, {"TP_Result"})
                .Define("Tag_Eta_Pass", [](const ResultsTagAndProbe& res) { return res.tag_eta_pass; }, {"TP_Result"});

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

            ROOT::RDataFrame data_frame(cfg.general.analysis_mode + "_Tree", cfg.io.o_file_data);
            
            if (verbose){ 
                std::cout << "RDataFrame object created, unpacking " << cfg.general.analysis_mode
                << " tree from " <<  cfg.io.o_file_data << " file, starting analysis ..." << std::endl;
            }
        
            // HARDCODED

            std::vector<float> pt_bins = {25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0};
            std::vector<float> eta_bins = {-2.4, -1.8, -1.2, -0.6, 0.0, 0.6, 1.2, 1.8, 2.4};
            
            ROOT::RDF::RNode node_pt_pass = ApplyPt_DataDivision(data_frame, pt_bins, "Probe_Pt_Pass");
            ROOT::RDF::RNode node_pt_all = ApplyPt_DataDivision(data_frame, pt_bins, "Probe_Pt_All");

            ROOT::RDF::RNode node_eta_pass = ApplyEta_DataDivision(data_frame, eta_bins, "Probe_Eta_Pass");
            ROOT::RDF::RNode node_eta_all = ApplyEta_DataDivision(data_frame, eta_bins, "Probe_Eta_All");

        } catch (const std::exception& except) {
            std::cerr << "Error nature: " << except.what() << std::endl;
            return 1;
        }
    }

    return 0;
}