#include <ROOT/RDataFrame.hxx>

#include <iostream>
#include <string>

#include <chrono>

#include <TApplication.h>
#include <TCanvas.h>

#include "Config.h"
#include "Utils.h"

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

    if (Configure(cfg, json_path) != 0) {
        std::cout << "Configurations fails, check needed, exiting." << std::endl;
        return 1;
    }

    if (verbose) {
        Verbose_config(cfg);
    }

    
    TApplication* app = nullptr;
    if (visualize) {
        app = new TApplication("app", &argc, argv);
    }
    
    try {
        // Trying to Enable Implicit MT
        ROOT::EnableImplicitMT();
        ROOT::RDataFrame data_frame(cfg.io.tree_name, cfg.io.input_file);


        auto validation_map = Validation_load(cfg.io.validation_file);

        // Validation filter
        auto validation_selection = [validation_map] (const UInt_t run, const UInt_t lum_block) {
            
            thread_local static UInt_t last_run = 0;
            thread_local static UInt_t last_lumi = 0;
            thread_local static bool last_decision = false;

            if (run == last_run && lum_block == last_lumi) {
                return last_decision;
            }
            
            last_run = run;
            last_lumi = lum_block;

            auto it = validation_map.find(run);
            
            if (it != validation_map.end()) {
                for (const auto& blocks : it->second) {
                    if (lum_block >= blocks.first && lum_block <= blocks.second) {
                        last_decision = true;
                        return true;
                    }
                }
            }
            last_decision = false;
            return false;
        };

        // Good Muon Filter
        auto good_selection = [cfg] (const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, 
        const ROOT::RVec<bool>& tight_id, const ROOT::RVec<float>& iso) {
            ROOT::RVec<bool> mask(pt.size(), true);
            
            if (cfg.flag.en_tight_muon) {
                mask = mask && tight_id;
            }
            if (cfg.flag.en_kinematics) {
                mask = mask && (pt > cfg.cut.pt_cut) && (abs(eta) < cfg.cut.eta_cut);
            }
            if (cfg.flag.en_isolation) {
                mask = mask && (iso < cfg.cut.iso_cut);
            }
            
            return mask;
        };

        ROOT::RDF::RNode mass_data = data_frame
            .Filter(validation_selection, {"run", "luminosityBlock"})
            // Defining a bool variable -> GoodMuon: Kinematics + Identification + Isolation
            .Define("GoodMuon", good_selection, {"Muon_pt", "Muon_eta", "Muon_tightId", "Muon_pfRelIso04_all"})
            
            // Minumum number of Good Muons
            .Filter("Sum(GoodMuon) >= 2")

            // Creating good variables needed for the invariant mass calculus
            .Define("good_Pt",   "Muon_pt[GoodMuon]")
            .Define("good_Eta",  "Muon_eta[GoodMuon]")
            .Define("good_Phi",  "Muon_phi[GoodMuon]")
            .Define("good_Mass", "Muon_mass[GoodMuon]")
            .Define("good_Charge", "Muon_charge[GoodMuon]");

        // Opposite Charge filter
        if (cfg.flag.en_opposite_charge) {
            mass_data = mass_data.Filter([cfg] (const ROOT::RVec<int>& charge) { return charge[0] != charge[1]; }, {"good_Charge"});
        }

        // Invarian Mass
        mass_data = mass_data.Define("m_ll", mass_leptons<float>, {"good_Pt", "good_Eta", "good_Phi", "good_Mass"});
        
        if (cfg.flag.en_mass_window) {
            mass_data = mass_data.Filter([cfg] (const float mass) { return (mass > cfg.cut.mass_min) && (mass < cfg.cut.mass_max); },
            {"m_ll"});
        }

        ROOT::RDF::RNode phis_data = mass_data.Define("phi_star", phi_star<float>, {"good_Eta", "good_Phi"});

        auto histo_m_ll = mass_data.Histo1D({"histo_m_ll", "Invariant mass; m_{#mu+#mu-} [GeV]; Events", 100, 60.0, 120.0}, "m_ll");
        auto histo_phis = phis_data.Histo1D({"histo_phis", "Angular variable; #phi* [rad]; Events", 100, 0, 6}, "phi_star");
        //auto start = std::chrono::high_resolution_clock::now();
        mass_data.Report()->Print();

        if (visualize) {
            TCanvas canvas("c1", "Massa Z", 800, 600);
            histo_m_ll->Draw("E_HIST");

            TCanvas canvas2("c2", "Phi", 800, 600);
            histo_phis->Draw("");
            
            canvas.Connect("Closed()", "TApplication", app, "Terminate()");
            canvas.Update();
            
            app->Run(); 
        } else {
            std::cout << "No visualization booked.\n" << std::endl;
        }

        //auto end = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed = end - start;
        //std::cout << "Tempo di esecuzione: " << elapsed.count() << " s" << std::endl;

    } catch (const std::exception& except) {
        std::cerr << "Error nature: " << except.what() << std::endl;
        return 1;
    }

    return 0;
}