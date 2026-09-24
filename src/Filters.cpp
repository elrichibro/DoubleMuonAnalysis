#include "Filters.h"
#include "Utils.h"
#include "Checks.h"

#include <Rtypes.h>

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RNode ApplyValidationFilter(ROOT::RDF::RNode node, const validation_type& val_map, const std::string& run_name, const std::string& block_name) {
    ROOT::RDF::RNode node_validation = node;

    node_validation = node_validation
        .Filter([&val_map](const UInt_t run, const UInt_t lum_block) {
            // Auto-update variables -> applied to multithread operation mode.
            thread_local static UInt_t last_run = 0;
            thread_local static UInt_t last_lum_block = 0;
            thread_local static bool last_decision = false;// For fast loop

            // Fast loop
            if ((run == last_run) && (lum_block == last_lum_block)) {
                return last_decision;
            }

            // Updating values
            last_run = run;
            last_lum_block = lum_block;

            // Checks the run input within the validation_map
            auto it = val_map.find(run);

            // if succeds
            if (it != val_map.end()) {
                // Checking blocks
                for (const auto& block : it->second) {
                    // Luminosity block within the range
                    if ((lum_block >= block.first) && (lum_block <= block.second)) {
                        last_decision = true;
                        return true;
                    }
                }
            }
            // else option 
            last_decision = false;
            return false;

        }, {run_name , block_name}, "JSON Validation");
    
    // Filtered node
    return node_validation;
}

// ------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CalculateAcceptance(ROOT::RDF::RNode node, const std::string& tag, int FSR) {
    std::vector<float> results;

    // Selection of true event: Z0 -> mu+ mu-
    ROOT::RDF::RNode node_acc = node;

    node_acc = node_acc
        .Define("Z0_Event", [FSR](const ROOT::RVec<Int_t>& pdg, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother) {
            return is_MC_Event(pdg, flags, mother, FSR);
        }, {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"})
        .Filter("Z0_Event", "Is_Z0_Event");

    // All Z0-> mu+ mu- generated events
    auto tot_generated_events = node_acc.Count();

    // Building muon masks for kinematical selection
    node_acc = node_acc
        .Define("Mu_mask_" + tag, (FSR == 1) ? is_MC_Muon_bFSR : is_MC_Muon_aFSR, {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"})
        .Define("AMu_mask_" + tag, (FSR == 1) ? is_MC_AntiMuon_bFSR : is_MC_AntiMuon_aFSR, {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"})

        .Define("Mu_pt", "GenPart_pt[Mu_mask_" + tag + "]")
        .Define("Mu_eta", "GenPart_eta[Mu_mask_" + tag + "]")
        .Define("Mu_phi", "GenPart_phi[Mu_mask_" + tag + "]")
        .Define("Mu_mass", "GenPart_mass[Mu_mask_" + tag + "]")
        
        .Define("AMu_pt", "GenPart_pt[AMu_mask_" + tag + "]")
        .Define("AMu_eta", "GenPart_eta[AMu_mask_" + tag + "]")
        .Define("AMu_phi", "GenPart_phi[AMu_mask_" + tag + "]")
        .Define("AMu_mass", "GenPart_mass[AMu_mask_" + tag + "]")

        // Physical cuts
        .Filter("(Mu_pt[0] > 25) && (AMu_pt[0] > 25)", "Pt_cut")
        .Filter("(abs(Mu_eta[0]) < 2.4) && (abs(AMu_eta[0]) < 2.4)", "Eta_cut")

        // Invariant mass calculus
        .Define("Z0_InvMass", [](const ROOT::RVec<float>& mu_pt, const ROOT::RVec<float>& anti_mu_pt, const ROOT::RVec<float>& mu_eta, 
        const ROOT::RVec<float>& anti_mu_eta, const ROOT::RVec<float>& mu_phi, const ROOT::RVec<float>& anti_mu_phi, 
        const ROOT::RVec<float>& mu_mass, const ROOT::RVec<float>& anti_mu_mass) -> float {
                
            return CalculateInvariantMass_Pair<float>(mu_pt[0], anti_mu_pt[0], mu_eta[0], anti_mu_eta[0], mu_phi[0], anti_mu_phi[0], 
                mu_mass[0], anti_mu_mass[0]);

        }, {"Mu_pt", "AMu_pt", "Mu_eta", "AMu_eta", "Mu_phi", "AMu_phi", "Mu_mass", "AMu_mass"})

        // Invariant mass filter -> fiducial region
        .Filter("(Z0_InvMass > 60.0) && (Z0_InvMass < 120.0)", "Mass_cut");

        // Accepted events
        auto acc_events = node_acc.Count();
        
        auto report_node = node_acc.Report();
        report_node->Print();

        float num_tot = static_cast<float>(tot_generated_events.GetValue());
        float num_acc = static_cast<float>(acc_events.GetValue());

        float acc = (num_tot > 0) ? (num_acc / num_tot) : 0.0;
        results.push_back(acc);

        float acc_err = sqrt((acc * (1 - acc)) / num_tot);
        results.push_back(acc_err);

    return results;
}

// ------------------------------------------------------------------------------------------------------------------------------------

ResolutionResults CalculateResolution(ROOT::RDF::RNode node, const config_struct& cfg) {
    ROOT::RDF::RNode node_resolution = node;
    
    // Resolution quantity
    std::string tag = cfg.resolution.quantity;
    
    // Columns in RDF node
    std::string gen_col = "";
    std::string reco_col = "";

    if (tag == "pt") {
        reco_col = "Rec_Pt";
        gen_col = "Gen_Pt";
    } else if (tag == "y") {
        reco_col = "Rec_Y";
        gen_col = "Gen_Y";
    } else if (tag == "phis") {
        reco_col = "Rec_Phis";
        gen_col = "Gen_Phis";
    }

    // Vector of bins
    std::vector<double> vector_bins = CreateBins(cfg.resolution.gen_bins, cfg.resolution.min, cfg.resolution.max, "linear");
    int n_bins = vector_bins.size() - 1;

    // Initializing struct
    ResolutionResults result;
    result.mean.resize(n_bins);
    result.sigma.resize(n_bins);
    result.events.resize(n_bins);

    node_resolution = node_resolution.Filter("Matched == true");

    auto gen_ptr = node_resolution.Take<float>(gen_col);
    auto reco_ptr = node_resolution.Take<float>(reco_col);

    const auto& gen_vec = *gen_ptr;
    const auto& reco_vec = *reco_ptr;

    std::vector<std::vector<double>> bin_buffers(n_bins);

    // For each event -> unpacking reconstructed quantity and saving it into bin vector.
    size_t n_entries = gen_vec.size();
    for (size_t i = 0; i < n_entries; i++) {

        auto upper_idx = std::upper_bound(vector_bins.begin(), vector_bins.end(), gen_vec[i]);
        int bin_idx = std::distance(vector_bins.begin(), upper_idx) - 1;

        if ((bin_idx >= 0) && (bin_idx < n_bins)) {
            bin_buffers[bin_idx].push_back(reco_vec[i]);
        }
    }

    // For each bin content
    for (int j = 0; j < n_bins; j++) {
        const auto& bin = bin_buffers[j];
        
        int entries = bin.size();
        result.events[j] = entries;

        if (entries == 0) {
            continue;
        }
        
        // Mean calculus
        double sum = 0.0;
        for (int i = 0; i < entries; i++) {
            sum += bin[i];
        }

        double mean = sum / entries;
        result.mean[j] = mean;

        // Standard deviation calculus
        if (entries > 1) {
            double sq_sum = 0;
            for (int i = 0; i < entries; i++) {
                double diff = bin[i] - mean;   
                sq_sum += diff * diff;
            }
            result.sigma[j] = std::sqrt(sq_sum / (entries - 1));
        }
    }
    return result;
}

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RNode EventSelection(ROOT::RDF::RNode node, const config_struct& cfg) {
    ROOT::RDF::RNode node_event = node;

    // Defining Event cuts from JSON info
    std::string good_muon = "true";
    if (cfg.flag_ES.en_kinematics) {
        good_muon += " && (Muon_pt > " + std::to_string(cfg.cut_ES.pt_cut) + " && abs(Muon_eta) < " + std::to_string(cfg.cut_ES.eta_cut)
         + ")";
    }

    if (cfg.flag_ES.en_tight_muon) {
        good_muon += " && (Muon_tightId == true)";
    }

    std::string event_cut = "GoodEvent";
    if (cfg.flag_ES.en_mass_window) {
        event_cut += " && (InvariantMass > " + std::to_string(cfg.cut_ES.mass_min) + " && InvariantMass < " 
        + std::to_string(cfg.cut_ES.mass_max) + ")";
    }

    node_event = node_event
        .Define("GoodMuon", good_muon)

        .Define("GM_Pt", "Muon_pt[GoodMuon]")
        .Define("GM_Eta", "Muon_eta[GoodMuon]")
        .Define("GM_Phi", "Muon_phi[GoodMuon]")
        .Define("GM_Mass", "Muon_mass[GoodMuon]")
        .Define("GM_Charge", "Muon_charge[GoodMuon]")

        .Define("EventPair", "Sum(GoodMuon) == 2")
        .Define("GoodEvent", "EventPair && (GM_Charge[0] != GM_Charge[1])")

        .Define("InvariantMass", [] (const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi,
         const ROOT::RVec<float>& mass) -> float {
            
            return static_cast<float>(ROOT::VecOps::InvariantMass(pt, eta, phi, mass));
        
        }, {"GM_Pt", "GM_Eta", "GM_Phi", "GM_Mass"})

        .Filter(event_cut, "InvMass selection -> Good Event")

        .Define("Pt_Z", [] (const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi,
         const ROOT::RVec<float>& mass) {
            
            float pt_Z = CalculatePtZ0<float>(pt, eta, phi, mass);
            return pt_Z;
        
        }, {"GM_Pt", "GM_Eta", "GM_Phi", "GM_Mass"})
        
        .Define("Y_Z", [] (const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi,
         const ROOT::RVec<float>& mass) {
            
            float y_Z = CalculateRapidityZ0<float>(pt, eta, phi, mass);
            return y_Z;
        
        }, {"GM_Pt", "GM_Eta", "GM_Phi", "GM_Mass"})
        
        .Define("Phis_Z", [] (const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi) {
            
            float phis_Z = CalculatePhiStar<float>(eta, phi);
            return phis_Z;
        
        }, {"GM_Eta", "GM_Phi"});

    return node_event;
}

// ------------------------------------------------------------------------------------------------------------------------------------

EventHisto BuildEventHisto(ROOT::RDF::RNode node, const config_struct& cfg) {
    ROOT::RDF::RNode node_event = node;
    EventHisto histo;
      
    std::vector<double> bins_pt, bins_y, bins_phis;
    int n_pt, n_y, n_phis;

    // Choosing binning options
    if (cfg.unfold.use_custom_bins == true) {
        bins_pt = cfg.unfold.pt_bins.reco_vec;
        bins_y = cfg.unfold.y_bins.reco_vec;
        bins_phis = cfg.unfold.phis_bins.reco_vec;

    } else {
        const auto& pt = cfg.unfold.pt_bins;
        const auto& y = cfg.unfold.y_bins;
        const auto& phis = cfg.unfold.phis_bins;

        bins_pt = CreateBins(pt.reco_bins, pt.min, pt.max, pt.distribution);
        bins_y = CreateBins(y.reco_bins, y.min, y.max, y.distribution);            
        bins_phis = CreateBins(phis.reco_bins, phis.min, phis.max, phis.distribution);
    }

    // Mll bins - HARDCODED
    int n_mll = 70;
    
    std::vector<double> mll_bins(n_mll + 1);
    double step = (120.0 - 60.0) / n_mll;

    for (int i = 0; i < n_mll + 1; i++) {
        mll_bins[i] = 60.0 + (i * step);
    }

    n_pt = static_cast<int>(bins_pt.size()) - 1;
    n_y = static_cast<int>(bins_y.size()) - 1;
    n_phis = static_cast<int>(bins_phis.size()) - 1;

    ROOT::RDF::TH2DModel model_mll_pt("h2_mll_pt", "Pt vs Mll", n_pt, bins_pt.data(), n_mll, mll_bins.data());
    ROOT::RDF::TH2DModel model_mll_y("h2_mll_y", "Y vs Mll", n_y, bins_y.data(), n_mll, mll_bins.data());
    ROOT::RDF::TH2DModel model_mll_phis("h2_mll_phis", "Phis vs Mll", n_phis, bins_phis.data(), n_mll, mll_bins.data());

    histo.h1_mll = node_event.Histo1D({"hInvMass_fit","", 70, 60.0, 120.0}, "InvariantMass");
        
    histo.h1_pt = node_event.Histo1D({"hPt_event", "", n_pt, bins_pt.data()}, "Pt_Z");
    histo.h1_y = node_event.Histo1D({"hRapidity_event", "", n_y, bins_y.data()}, "Y_Z");
    histo.h1_phis = node_event.Histo1D({"hPhis_event", "", n_phis, bins_phis.data()}, "Phis_Z");

    histo.h2_mll_pt = node_event.Histo2D(model_mll_pt, "Pt_Z", "InvariantMass");
    histo.h2_mll_y = node_event.Histo2D(model_mll_y, "Y_Z", "InvariantMass");
    histo.h2_mll_phis = node_event.Histo2D(model_mll_phis, "Phis_Z", "InvariantMass");

    return histo;
}