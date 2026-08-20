#include <ROOT/RDataFrame.hxx>

#include <iostream>
#include <string>

#include <chrono>

#include <TCanvas.h>

#include "Config.h"

// Funzione C++ per il calcolo della massa invariante
float mass_leptons(const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi, const ROOT::RVec<float>& mass) 
{
    return ROOT::VecOps::InvariantMass( ROOT::RVec<float>{pt[0], pt[1]}, ROOT::RVec<float>{eta[0], eta[1]}, ROOT::RVec<float>{phi[0], phi[1]},
        ROOT::RVec<float>{mass[0], mass[1]}
    );
}

int main(int argc, char* argv[]) {

    std::string json_path = "";
    int verbose = 0;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg.rfind(".json") != std::string::npos) {
            json_path = arg;
        } else if (arg == "--verbose" || arg == "-v") {
            verbose = 1;
        } else {
            std::cout << "ERROR: invalid input command, please try again, exinting.\n" << std::endl;
            return 1;
        }
    }

    config_struct cfg;
    Configure(cfg, json_path);

    try {
        // Trying to Enable Implicit MT
        ROOT::EnableImplicitMT();
        ROOT::RDataFrame data_frame(cfg.io.tree_name, cfg.io.input_file);

        auto mass_data = data_frame
            // Defining a bool variable -> GoodMuon that is true when the condition succed: Kinematics + Identification + Isolation
            .Define("GoodMuon", "(Muon_pt > 25) && (abs(Muon_eta) < 2.4) && (Muon_tightId != 0) && (Muon_pfRelIso04_all < 0.15)")
            
            // Minumum number of Good Muons
            .Filter("Sum(GoodMuon) >= 2")

            // Creating good variables needed for the invariant mass calculus
            .Define("good_Pt",   "Muon_pt[GoodMuon]")
            .Define("good_Eta",  "Muon_eta[GoodMuon]")
            .Define("good_Phi",  "Muon_phi[GoodMuon]")
            .Define("good_Mass", "Muon_mass[GoodMuon]")
            
            // Opposit charge filter
            .Define("good_Charge", "Muon_charge[GoodMuon]")
            .Filter("good_Charge[0] != good_Charge[1]")
            
            // Using external function -> mass_leptons
            .Define("m_ll", mass_leptons, {"good_Pt", "good_Eta", "good_Phi", "good_Mass"})

            // Invariant Mass filter in the Z region
            .Filter("m_ll > 60 && m_ll < 120", "Z mass window");

        auto histo_m_ll = mass_data.Histo1D({"histo_m_ll", "Invariant mass;m_{#mu+#mu-} [GeV];Events", 100, 60.0, 120.0}, "m_ll");
        
        auto start = std::chrono::high_resolution_clock::now();
        mass_data.Report()->Print();

        // Check out -> loop initialized
        histo_m_ll->Draw("HIST");
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Tempo di esecuzione: " << elapsed.count() << " s" << std::endl;

    } catch (const std::exception& except) {
        std::cerr << "Error nature: " << except.what() << std::endl;
        return 1;
    }

    return 0;
}