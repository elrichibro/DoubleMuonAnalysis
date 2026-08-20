#include <ROOT/RDataFrame.hxx>
#include <iostream>
#include <string>
#include <TCanvas.h>

// Funzione C++ per il calcolo della massa invariante
float mass_leptons(const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi, const ROOT::RVec<float>& mass) 
{
    return ROOT::VecOps::InvariantMass( ROOT::RVec<float>{pt[0], pt[1]}, ROOT::RVec<float>{eta[0], eta[1]}, ROOT::RVec<float>{phi[0], phi[1]},
        ROOT::RVec<float>{mass[0], mass[1]}
    );
}

int main(int argc, char** argv) {
    std::string tree_name = "Events";
    std::string file_name = "../data/dati0.root";

    try {
        // Trying to Enable Implicit MT
        ROOT::EnableImplicitMT();
        ROOT::RDataFrame data_frame(tree_name, file_name);

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
        
        mass_data.Report()->Print();

        // Check out -> loop initialized
        histo_m_ll->Draw("HIST");
        
    } catch (const std::exception& except) {
        std::cerr << "Error nature: " << except.what() << std::endl;
        return 1;
    }

    return 0;
}