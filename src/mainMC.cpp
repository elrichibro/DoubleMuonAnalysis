#include <ROOT/RDataFrame.hxx>

#include <iostream>
#include <string>

#include <TApplication.h>
#include <TCanvas.h>

#include "Config.h"
#include "Filters.h"
#include "Utils.h"

int main(int argc, char* argv[]) {

    config_struct cfg;
    std::string json_path = "";
    int verbose = 0;
    int visualize = 0;
    
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
            visualize = 1;
        } else {
            std::cout << "ERROR: invalid input command, please try again, exinting.\n" << std::endl;
            return 1;
        }
    }

    // ------------------
    // JSON CONFIGURATION
    // ------------------

    if (Configure(cfg, json_path) != 0) {
        std::cout << "Configurations fails, check needed, exiting." << std::endl;
        return 1;
    }

    if (verbose) {
        Verbose_config(cfg);
    }

    // --------------------
    // VISUALIZATION OPTION
    // --------------------
    
    TApplication* app = nullptr;
    if (visualize) {
        std::cout << "Visualization inializating ..." << std::endl;
        app = new TApplication("app", &argc, argv);
    }

    try {
        ROOT::EnableImplicitMT();// MultiThread option: ON

        ROOT::RDataFrame data_frame(cfg.io.tree_mc_name, cfg.io.in_mc_file);
        if (verbose){ 
            std::cout << "RDataFrame object created, unpaching" << cfg.io.tree_mc_name 
            << " from " << cfg.io.in_mc_file << " file, starting analysis ..." << std::endl;
        }

        // -----
        // NODES
        // -----
        
        auto node_isEvent = data_frame
        .Filter(is_MC_Event, {"GenPart_pdgId", "GenPart_statusFlags"}, "1. True Event");

        auto node_analysis = node_isEvent
        .Define("Mu_mask", is_MC_Muon, {"GenPart_pdgId", "GenPart_statusFlags"})
        .Define("AntiMu_mask", is_MC_AntiMuon, {"GenPart_pdgId", "GenPart_statusFlags"})

        .Define("event_Mu_pt",  "GenPart_pt[Mu_mask]")
        .Define("event_Mu_eta", "GenPart_eta[Mu_mask]")
        .Define("event_AntiMu_pt",   "GenPart_pt[AntiMu_mask]")
        .Define("event_AntiMu_eta",  "GenPart_eta[AntiMu_mask]")

        .Define("good_Pt", "GenPart_pt[(Mu_mask || AntiMu_mask)]")
        .Define("good_Eta", "GenPart_eta[(Mu_mask || AntiMu_mask)]")
        .Define("good_Phi", "GenPart_phi[(Mu_mask || AntiMu_mask)]")
        .Define("good_Mass", "GenPart_mass[(Mu_mask || AntiMu_mask)]");

        ROOT::RDF::RNode node_inv_mass = node_analysis    
            .Define("m_ll", mass_leptons<float>, {"good_Pt", "good_Eta", "good_Phi", "good_Mass"});
        
        if (cfg.flag.en_mass_window) {
            node_inv_mass = node_inv_mass.Filter([&cfg] (const float mass) { return (mass > cfg.cut.mass_min) && (mass < cfg.cut.mass_max); },
            {"m_ll"}, "2. Z0 range selection");
        }
        
        // ----------
        // HISTOGRAMS
        // ----------

        auto h_mu_pt  = node_analysis.Histo1D({"h_mu_pt", "p_{T}(#mu^{-}); p_{T} [GeV]; Events", 100, 0, 100}, "event_Mu_pt");
        auto h_mu_eta = node_analysis.Histo1D({"h_mu_eta", "#eta(#mu^{-}); #eta; Events", 100, -2.5, 2.5}, "event_Mu_eta");
        auto h_antimu_pt = node_analysis.Histo1D({"h_antimu_pt", "p_{T}(#mu^{+}); p_{T} [GeV]; Events", 100, 0, 100}, "event_AntiMu_pt");
        auto h_antimu_eta = node_analysis.Histo1D({"h_antimu_eta", "#eta(#mu^{+}); #eta; Events", 100, -2.5, 2.5}, "event_AntiMu_eta");
        
        auto h_mass_ll = node_inv_mass.Histo1D({"h_mass_ll", "Massa invariante dileptoni; m_{(#mu^{+})(#mu^{-})}; Events", 100, 0, 200}, "m_ll");

        if (visualize) {
            auto canvas = new TCanvas("c_all", "Muon Kinematics", 1200, 800);
            canvas->Divide(2, 2);

            canvas->cd(1); h_mu_pt->Draw("E HIST");
            canvas->cd(2); h_mu_eta->Draw("E HIST");
            canvas->cd(3); h_antimu_pt->Draw("E HIST");
            canvas->cd(4); h_antimu_eta->Draw("E HIST");

            auto canvas2 = new TCanvas("c", "Zmass", 800, 600);
            h_mass_ll->Draw("E HIST");

            canvas->Connect("Closed()", "TApplication", app, "Terminate()");
            canvas->Update();

            app->Run();
        } else {
            std::cout << "No visualization booked.\n" << std::endl;
        }

    } catch (const std::exception& except) {
        std::cerr << "Error nature: " << except.what() << std::endl;
        return 1;
    }

    return 0;
}