#include "AnalysisTools.h"
#include <string>

// ------------------------------------------------------------------------------------------------------------------------------------
// MC Template Maker
// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RResultPtr<TH3D> TemplateMaker_MC(ROOT::RDF::RNode node, const config_struct& cfg, const bool mask) {
    ROOT::RDF::RNode node_hist = node;

    std::vector<float> pt_bins = cfg.analysis.pt_bins;
    std::vector<float> eta_bins = cfg.analysis.eta_bins;
    std::vector<float> mll_bins(cfg.analysis.mll_bins);

    float step = (120.0 - 60.0) / (cfg.analysis.mll_bins - 1.0);
    for (int i = 0; i < cfg.analysis.mll_bins; i++) {
        mll_bins[i] = 60.0 + (i * step);
    }

    std::string name = (mask) ? "h3_pass" : "h3_fail";
    ROOT::RDF::TH3DModel model(name.c_str(), "3D Histogram MC", eta_bins.size() - 1, eta_bins.data(), pt_bins.size() - 1, pt_bins.data(), 
    mll_bins.size() - 1, mll_bins.data() );
    
    if (mask) {
        node_hist = node_hist
            .Define("MC_Probe_Pt_Pass", "MC_Probe_Pt[MC_Mask_Pass]")
            .Define("MC_Probe_Eta_Pass", "MC_Probe_Eta[MC_Mask_Pass]")
            .Define("MC_Mll_Pass", "MC_Mll[MC_Mask_Pass]");

        auto h3 = node_hist.Histo3D(model, "MC_Probe_Eta_Pass", "MC_Probe_Pt_Pass", "MC_Mll_Pass");
        return h3;
    } else {
        node_hist = node_hist
            .Define("MC_Probe_Pt_Fail", "MC_Probe_Pt[!MC_Mask_Pass]")
            .Define("MC_Probe_Eta_Fail", "MC_Probe_Eta[!MC_Mask_Pass]")
            .Define("MC_Mll_Fail", "MC_Mll[!MC_Mask_Pass]");

        auto h3 = node_hist.Histo3D(model, "MC_Probe_Eta_Fail", "MC_Probe_Pt_Fail", "MC_Mll_Fail");
        return h3;
    }
}


int LoadMCTemplate(const config_struct& cfg, std::vector<MC_Template_RooF>& container) {
    std::unique_ptr<TFile> file(TFile::Open(cfg.io.o_file_template.c_str(), "READ"));
    
    if (!file || file->IsZombie()) {
        std::cout << "ERROR: Cannot open the template output file: " << cfg.io.o_file_template << std::endl;
        return 1;
    }

    auto h3_pass = file->Get<TH3D>("h3_pass");
    auto h3_fail = file->Get<TH3D>("h3_fail");

    if ((!h3_pass) && (!h3_fail)) {
        std::cout << "ERROR: null histogram pointer, exiting" << std::endl;
        return 1;
    }

    std::vector<float> pt_bins = cfg.analysis.pt_bins;
    std::vector<float> eta_bins = cfg.analysis.eta_bins;
    
    const int n_pt_bins = pt_bins.size() - 1;
    const int n_eta_bins = eta_bins.size() - 1;

    const int tot_bins = n_pt_bins * n_eta_bins;
    
    container.clear();
    container.reserve(tot_bins);

    for (int i = 0; i < n_pt_bins; i++) {
        int bin_pt = (i + 1);

        for (int j = 0; j < n_eta_bins; j++) {
            int bin_eta = (j + 1);
    
            std::string name_pass = "h_mll_pass_eta" + std::to_string(bin_eta) + "_pt" + std::to_string(bin_pt);
            std::string name_fail = "h_mll_fail_eta" + std::to_string(bin_eta) + "_pt" + std::to_string(bin_pt);
            
            TH1D* h1_pass = h3_pass->ProjectionZ(name_pass.c_str(), bin_eta, bin_eta, bin_pt, bin_pt);
            TH1D* h1_fail = h3_fail->ProjectionZ(name_fail.c_str(), bin_eta, bin_eta, bin_pt, bin_pt);
        
            if (h1_pass) {
                h1_pass->SetDirectory(nullptr);
            }
            if (h1_fail) {
                h1_fail->SetDirectory(nullptr);
            }
            
            container.emplace_back(MC_Template_RooF{
                bin_pt, 
                bin_eta, 
                pt_bins[i], 
                pt_bins[i + 1], 
                eta_bins[j], 
                eta_bins[j + 1],
        
                h1_pass,
                h1_fail
            });
        }
    }

    file->Close();
    return 0;
}