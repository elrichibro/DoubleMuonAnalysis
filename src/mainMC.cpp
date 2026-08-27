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

        /*
        ROOT::RDF::RNode node_InvMass_bFSR = data_frame;
        node_InvMass_bFSR = InvMass_BeforeFSR(node_InvMass_bFSR, cfg);

        ROOT::RDF::RNode node_InvMass_aFSR = data_frame;
        node_InvMass_aFSR = InvMass_AfterFSR(node_InvMass_aFSR, cfg);
        
        auto report_bFSR = node_InvMass_bFSR.Report();
        auto report_aFSR = node_InvMass_aFSR.Report();
        */

        auto node_InvMass_bFSR = InvMass(data_frame, cfg, "bFSR", 1);
        auto node_InvMass_aFSR = InvMass(data_frame, cfg, "aFSR", 2);
        
        auto report_bFSR = node_InvMass_bFSR.Report();
        auto report_aFSR = node_InvMass_aFSR.Report();
        
        // ----------
        // HISTOGRAMS
        // ----------

        auto h_mass_ll_bFSR = node_InvMass_bFSR.Histo1D({"m_ll_bFSR", "Massa invariante dileptoni Before FSR; m_{#mu^{+}#mu^{-}}; Events", 100, 60, 120}, "m_ll_bFSR");        
        auto h_mass_ll_aFSR = node_InvMass_aFSR.Histo1D({"m_ll_aFSR", "Massa invariante dileptoni After FSR; m_{#mu^{+}#mu^{-}}; Events", 100, 60, 120}, "m_ll_aFSR");

        report_bFSR->Print();
        report_aFSR->Print();

        if (visualize) {
            auto canvas = new TCanvas("c_masses", "Invariant Mass Comparison", 1200, 600);
            canvas->Divide(2, 1);

            canvas->cd(1);
            h_mass_ll_bFSR->Draw("E HIST");

            canvas->cd(2);
            h_mass_ll_aFSR->Draw("E HIST");

            canvas->Update();
            canvas->Connect("Closed()", "TApplication", app, "Terminate()");

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