#include <ROOT/RDataFrame.hxx>

#include <iostream>
#include <string>

#include <chrono>

#include <TApplication.h>
#include <TCanvas.h>

#include "Config.h"
#include "Utils.h"
#include "Filters.h"

int main(int argc, char* argv[]) {

    // Variables needed
    config_struct cfg;
    std::string json_path = "";
    int verbose = 0;
    int visualize = 0;
    int time = 0;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg.rfind(".json") != std::string::npos) {
            json_path = arg;
        } else if (arg == "--verbose" || arg == "-v") {
            verbose = 1;
        } else if ((arg == "--visualize") || (arg == "-vis")) {
            visualize = 1;
        } else if ((arg == "--time") || (arg == "-t")) {
            time = 1;
        } else {
            std::cout << "ERROR: invalid input command, please try again, exinting.\n" << std::endl;
            return 1;
        }
    }

    if (verbose) {
        Verbose_config(cfg);
    }

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
        ROOT::EnableImplicitMT();

        if (verbose){ 
            std::cout << "RDataFrame object created, starting analysis ..." << std::endl;
        }

        ROOT::RDataFrame data_frame(cfg.io.tree_data_name, cfg.io.in_data_file);

        auto validation_map = Validation_load(cfg.io.val_file);

        // Validation run filter
        auto node_valid_data = data_frame.Filter(Validation_filter(validation_map), {"run" ,"luminosityBlock"}, "1. JSON Validation");

        auto node1_good_muon = node_valid_data
            // Defining a bool variable -> GoodMuon: Kinematics + Identification + Isolation
            .Define("GoodMuon", GoodMuon_filter(cfg), {"Muon_pt", "Muon_eta", "Muon_tightId", "Muon_pfRelIso04_all"})
            .Filter("Sum(GoodMuon) == 2", "2. Good muon selection");
        
        auto node2_good_muon = node1_good_muon
            .Define("good_Charge", "Muon_charge[GoodMuon]")
            .Filter([&cfg] (const ROOT::RVec<int>& charge) { return charge[0] != charge[1]; }, {"good_Charge"}, "3. Opposite charge muons");

        // Subtracting Multiboson background
        auto node3_good_muon = node2_good_muon
            .Define("LooseMuon", "Muon_pt > 10.0 && abs(Muon_eta) < 2.4 && Muon_looseId")
            .Define("LooseElectron", "Electron_pt > 10.0 && abs(Electron_eta) < 2.5 && Electron_cutBased == 2")
            .Filter("Sum(LooseMuon) + Sum(LooseElectron) == 2", "4. Multiboson background");

        // Creating good variables needed for the invariant mass calculus
        auto node_good_variables = node3_good_muon
            .Define("good_Pt",   "Muon_pt[GoodMuon]")
            .Define("good_Eta",  "Muon_eta[GoodMuon]")
            .Define("good_Phi",  "Muon_phi[GoodMuon]")
            .Define("good_Mass", "Muon_mass[GoodMuon]");
            
        ROOT::RDF::RNode node_inv_mass = node_good_variables    
            .Define("m_ll", mass_leptons<float>, {"good_Pt", "good_Eta", "good_Phi", "good_Mass"});
        
        if (cfg.flag.en_mass_window) {
            node_inv_mass = node_inv_mass.Filter([&cfg] (const float mass) { return (mass > cfg.cut.mass_min) && (mass < cfg.cut.mass_max); },
            {"m_ll"}, "5. Z0 range selection");
        }

        auto cut_report = node_inv_mass.Report();

        // Phi* definition -> NOT in Z range
        auto node_phi_star = node_good_variables.Define("phi_star", phi_star<float>, {"good_Eta", "good_Phi"});
        
        // Histogram booked
        auto histo_m_ll = node_inv_mass.Histo1D({"histo_m_ll", "Invariant mass; m_{#mu+#mu-} [GeV]; Events", 100, 60.0, 120.0}, "m_ll");
        auto histo_phis = node_phi_star.Histo1D({"histo_phis", "Angular variable; #phi* [rad]; Events", 100, 0, 6}, "phi_star");
        
        if (verbose){ 
            std::cout << "Initializing Event Loop" << std::endl;
        }
        
        std::chrono::high_resolution_clock::time_point start;
        if (time) {
            start = std::chrono::high_resolution_clock::now();
        }
        
        cut_report->Print();        
        
        if (time) {
            auto end = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double> elapsed = end - start;
            std::cout << "Tempo di esecuzione: " << elapsed.count() << " s" << std::endl;
        }

        if (visualize) {
            TCanvas canvas("c1", "Massa Z", 800, 600);

            histo_m_ll->Draw("E_HIST");
            canvas.SetLogy();

            TCanvas canvas2("c2", "Phi", 800, 600);
            histo_phis->Draw("E_HIST");
            
            canvas.Connect("Closed()", "TApplication", app, "Terminate()");
            canvas.Update();
            
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