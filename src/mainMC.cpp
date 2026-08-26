#include <ROOT/RDataFrame.hxx>

#include <iostream>
#include <string>

#include <TApplication.h>
#include <TCanvas.h>

#include "Config.h"
#include "Filters.h"

int main(int argc, char* argv[]) {

    // Variables needed
    config_struct cfg;
    std::string json_path = "";
    int verbose = 0;
    int visualize = 0;
    
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
/*
    if (verbose) {
        Verbose_config(cfg);
    }
*/
    if (Configure(cfg, json_path) != 0) {
        std::cout << "Configurations fails, check needed, exiting." << std::endl;
        return 1;
    }

    if (verbose) {
        Verbose_config(cfg);
    }
    
    TApplication* app = nullptr;
    if (visualize) {
        std::cout << "Visualization inializating ..." << std::endl;
        app = new TApplication("app", &argc, argv);
    }

    try {
        // Trying to Enable Implicit MT
        //ROOT::EnableImplicitMT();

        if (verbose){ 
            std::cout << "RDataFrame object created, starting analysis ..." << std::endl;
        }

        ROOT::RDataFrame data_frame(cfg.io.tree_mc_name, cfg.io.in_mc_file);

        auto node_isEvent = data_frame
        .Filter(is_MC_Event, {"GenPart_pdgId", "GenPart_statusFlags"}, "1. True Event");

        auto node_analysis = node_isEvent
        .Define("Mu_mask", is_MC_Muon, {"GenPart_pdgId", "GenPart_statusFlags"})
        .Define("AntiMu_mask", is_MC_AntiMuon, {"GenPart_pdgId", "GenPart_status"})

        .Define("event_Mu_pt",  "GenPart_pt[Mu_mask]")
        .Define("event_Mu_eta", "GenPart_eta[Mu_mask]")
        .Define("event_AntiMu_pt",   "GenPart_pt[AntiMu_mask]")
        .Define("event_AntiMu_eta",  "GenPart_eta[AntiMu_mask]");

        
        auto h_mu_pt  = node_analysis.Histo1D({"h_mu_pt",  "p_{T}(#mu^{-});p_{T} [GeV];Events", 100, 0, 100}, "event_Mu_pt");
        auto h_mu_eta = node_analysis.Histo1D({"h_mu_eta", "#eta(#mu^{-});#eta;Events",          50, -2.5, 2.5}, "event_Mu_eta");
        auto h_antimu_pt = node_analysis.Histo1D({"h_antimu_pt",   "p_{T}(#mu^{+});p_{T} [GeV];Events", 100, 0, 100}, "event_AntiMu_pt");
        auto h_antimu_eta = node_analysis.Histo1D({"h_antimu_eta",  "#eta(#mu^{+});#eta;Events",          50, -2.5, 2.5}, "event_AntiMu_eta");          
        
        if (visualize) {
            TCanvas canvas1("c1", "Muon pt", 800, 600);
            h_mu_pt->Draw("E_HIST");

            TCanvas canvas2("c2", "Muon eta", 800, 600);
            h_mu_eta->Draw("E_HIST");

            TCanvas canvas3("c3", "AntiMuon pt", 800, 600);
            h_antimu_pt->Draw("E_HIST");

            TCanvas canvas4("c4", "AntiMuon eta", 800, 600);
            h_antimu_eta->Draw("E_HIST");
            
            canvas1.Connect("Closed()", "TApplication", app, "Terminate()");
            canvas1.Update();
            
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