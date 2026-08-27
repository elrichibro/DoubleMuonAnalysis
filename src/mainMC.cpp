#include <ROOT/RDataFrame.hxx>

#include <iostream>
#include <string>

#include <TApplication.h>
#include <TCanvas.h>
#include "TEfficiency.h"

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

        // --------------
        // Invariant Mass
        // --------------

        auto node_InvMass_bFSR = InvMass(data_frame, "bFSR", 1);
        auto node_InvMass_aFSR = InvMass(data_frame, "aFSR", 2);
        
        auto report_bFSR = node_InvMass_bFSR.Report();
        auto report_aFSR = node_InvMass_aFSR.Report();

        auto h_mass_ll_bFSR = node_InvMass_bFSR.Histo1D({"m_ll_bFSR", "Massa invariante dileptoni Before FSR; m_{#mu^{+}#mu^{-}}; Events", 100, 60, 120}, "m_ll_bFSR");        
        auto h_mass_ll_aFSR = node_InvMass_aFSR.Histo1D({"m_ll_aFSR", "Massa invariante dileptoni After FSR; m_{#mu^{+}#mu^{-}}; Events", 100, 60, 120}, "m_ll_aFSR");

        // ----------
        // Efficiency
        // ----------

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

        // ----------
        // HISTOGRAMS
        // ----------
        
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