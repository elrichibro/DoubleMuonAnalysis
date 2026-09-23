#include "Filters.h"
#include "Utils.h"
#include "Checks.h"

#include <Rtypes.h>

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RNode ApplyValidationFilter(ROOT::RDF::RNode node, const validation_type& val_map, const std::string& run_name, const std::string& block_name) {
    ROOT::RDF::RNode node_validation = node
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
    
    return node_validation;
}

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RNode ApplyKinematicalBinDivision(ROOT::RDF::RNode node, const config_struct& cfg, const std::string& pt_col, const std::string& eta_col, 
    const std::string& mll_col, std::vector<std::string>& columns_name,  const int dataset, const int sample) {

    ROOT::RDF::RNode node_kin_cut = node;

    const auto& eta_bins = cfg.templ.eta_bins;
    const auto& pt_bins  = cfg.templ.pt_bins;

    for (int i = 0; i < eta_bins.size() - 1; i++) {
        float eta_min = eta_bins.at(i);
        float eta_max = eta_bins.at(i + 1);

        for (int j = 0; j < pt_bins.size() - 1; j++) {
            float pt_min = pt_bins.at(j);
            float pt_max = pt_bins.at(j + 1);

            std::string column_name = "Mll_Eta" + std::to_string(i + 1) + "_Pt" + std::to_string(j + 1);

            if (dataset == 1) {
                column_name = "DATA_" + column_name;
            } else if (dataset == 2) {
                column_name = "MC_" + column_name;
            } else {
                std::cout << "ERROR: invalid identifier DATA/MC, exiting..." << std::endl;
            }

            if (sample == 1) {
                column_name = column_name + "_Pass";
            } else if (sample == 2) {
                column_name = column_name + "_Fail";
            } else {
                std::cout << "ERROR: invalid identifier pass/fail, exiting..." << std::endl;
            }

            columns_name.push_back(column_name);

            node_kin_cut = node_kin_cut
                .Define(column_name, [eta_min, eta_max, pt_min, pt_max](const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, 
                    const ROOT::RVec<float>& mll) {
                        auto mask = ((pt > pt_min) && (pt <= pt_max) && (eta > eta_min) && (eta <= eta_max));
                    
                    return mll[mask];
                }, {pt_col, eta_col, mll_col});

            std::cout << "Filter applied: " << column_name << std::endl; 
        }
    }

    return node_kin_cut;
}



ROOT::RDF::RNode CalculateTagAndProbeWrapper(ROOT::RDF::RNode node, const validation_type& validation_map, const selection_config& selection, 
    const flags_config& flags_TP, const cuts_config& cuts_TP) {

    ROOT::RDF::RNode node_TP = node;
    
    if (selection.dataset == "MC") {
        node_TP = node_TP
            .Define("TP_Result",
            [flags_TP, cuts_TP] (const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi,
            const ROOT::RVec<float>& mass, const ROOT::RVec<int>& charge,  const ROOT::RVec<bool>& tag, const ROOT::RVec<bool>& probe, const ROOT::RVec<bool>& glob, 
            const ROOT::RVec<float>& iso, const ROOT::RVec<bool>& hlt, const ROOT::RVec<bool>& reco_tight, const ROOT::RVec<UChar_t>& reco_flav_gen,
            const ROOT::RVec<Int_t>& reco_idx_gen, const ROOT::RVec<Int_t>& gen_status, const ROOT::RVec<Int_t>& gen_pdg_idx, 
            const ROOT::RVec<float> gen_eta, const ROOT::RVec<float> gen_phi) {
                
                MuonKinematics_TP kin{pt, eta, phi, mass, charge};
                MuonFlags_TP flags{tag, probe, glob, iso, hlt};
                MuonFlags val{reco_tight, reco_flav_gen, reco_idx_gen, gen_status, gen_pdg_idx};

                
                return CalculateTagAndProbe_MC(kin, flags, flags_TP, cuts_TP, val, gen_eta, gen_phi);

            }, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge", "Muon_tightId", "Muon_isStandalone", "Muon_isGlobal", "Muon_pfRelIso04_all",
                "HLT_Mu17", "Muon_tightId", "Muon_genPartFlav", "Muon_genPartIdx", "GenPart_status", "GenPart_pdgId", "GenPart_eta", "GenPart_phi"});
    
    } else if (selection.dataset == "DATA") {
        node_TP = ApplyValidationFilter(node_TP, validation_map, "run", "luminosityBlock");
        
        node_TP = node_TP
            .Define("TP_Result",
            [flags_TP, cuts_TP] (const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi,
            const ROOT::RVec<float>& mass, const ROOT::RVec<int>& charge,  const ROOT::RVec<bool>& tag, const ROOT::RVec<bool>& probe, 
            const ROOT::RVec<bool>& glob, const ROOT::RVec<float>& iso, const ROOT::RVec<bool>& hlt) {
                
                MuonKinematics_TP kin{pt, eta, phi, mass, charge};
                MuonFlags_TP flags{tag, probe, glob, iso, hlt};
                
                return CalculateTagAndProbe_DATA(kin, flags, flags_TP, cuts_TP);

            }, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge", "Muon_tightId", "Muon_isStandalone", "Muon_isGlobal", "Muon_pfRelIso04_all", "HLT_Mu17"});
    }


    node_TP = node_TP
        .Define(selection.dataset + "_Probe_Pt", [](const ResultsTagAndProbe& res) { return res.pt; }, {"TP_Result"})
        .Define(selection.dataset + "_Probe_Eta", [](const ResultsTagAndProbe& res) { return res.eta; }, {"TP_Result"})
        .Define(selection.dataset + "_Mll", [](const ResultsTagAndProbe& res) { return res.mll; }, {"TP_Result"})
        
        //.Define(selection.dataset + "_Tag_Pt", [](const ResultsTagAndProbe& res) { return res.tag_pt_pass; }, {"TP_Result"})
        //.Define(selection.dataset + "_Tag_Eta", [](const ResultsTagAndProbe& res) { return res.tag_eta_pass; }, {"TP_Result"})
        .Define(selection.dataset + "_Mask_Pass", [](const ResultsTagAndProbe& res) { return res.mask_pass; }, {"TP_Result"});    

    return node_TP;
}

// ------------------------------------------------------------------------------------------------------------------------------------

ResultsTagAndProbe CalculateTagAndProbe_DATA(const MuonKinematics_TP& kin, const MuonFlags_TP& flags, const flags_config cfg_f, 
const cuts_config cfg_c) {
    // Results container
    ResultsTagAndProbe results;

    // Number of particles in the event
    const int n_muons = kin.pt.size();

    results.mask_pass.reserve(n_muons);
    results.pt.reserve(n_muons);
    results.eta.reserve(n_muons);
    results.mll.reserve(n_muons);
    
    //results.tag_pt.reserve(n_muons);
    //results.tag_eta.reserve(n_muons);

    // Fast Exit -> BKG
    int n_tight = ROOT::VecOps::Sum(flags.tight);

    if (n_tight > 2) {
        return results;
    }

    if (n_tight == 2 && (ROOT::VecOps::Sum(flags.stand) > 2)) {
        return results;
    }

    for(int i = 0; i < n_muons; i++) {
        // Fast exit
        if ((!flags.tight[i]) && (!flags.hlt[i])) {
            continue;
        }

        for (int j = 0; j < n_muons; j++) {
            // Fast exit -> same particle
            if (i == j || !flags.stand[j]) {
                continue;
            }

            if (kin.charge[i] == kin.charge[j]) {
                continue;
            }

            // Loop into good probe muons
            const bool pass = ((kin.pt[j] > cfg_c.pt_cut) && (std::abs(kin.eta[j]) < cfg_c.eta_cut)) || (!(cfg_f.en_kinematics));
            if (!pass) {
                continue;
            }
            
            // Invariant mass (tag + probe)
            float mass = CalculateInvariantMass_Pair<float>(kin.pt[i], kin.pt[j], kin.eta[i], kin.eta[j], kin.phi[i], kin.phi[j],
            kin.mass[i], kin.mass[j]);
            
            // Invariant mass range
            if (((mass > cfg_c.mass_min) && (mass < cfg_c.mass_max)) || (!cfg_f.en_mass_window)) {
                // All probes (passed + failed)
                    
                results.pt.push_back(kin.pt[j]);
                results.eta.push_back(kin.eta[j]);   
                results.mll.push_back(mass);
                
                //results.tag_pt.push_back(kin.pt[i]);
                //results.tag_eta.push_back(kin.eta[i]);         
                
                if ((flags.global[j]) && (flags.iso[j] < cfg_c.iso_cut)) {
                    // Passed probes
                    results.mask_pass.push_back(true);
                } else {
                    // Failed probes
                    results.mask_pass.push_back(false);
                }
            }
        }
    }
    return results;
}

// ------------------------------------------------------------------------------------------------------------------------------------

ResultsTagAndProbe CalculateTagAndProbe_MC(const MuonKinematics_TP& kin, const MuonFlags_TP& flags, const flags_config cfg_f, 
const cuts_config cfg_c, const MuonFlags& DeltaR_flags, const ROOT::RVec<float> gen_eta, const ROOT::RVec<float> gen_phi) {
    // Results container
    ResultsTagAndProbe results;

    // Number of particles in the event
    const int n_muons = kin.pt.size();

    results.mask_pass.reserve(n_muons);
    results.pt.reserve(n_muons);
    results.eta.reserve(n_muons);
    results.mll.reserve(n_muons);
    
    //results.tag_pt.reserve(n_muons);
    //results.tag_eta.reserve(n_muons);
    
    // Fast Exit -> BKG
    int n_tight = ROOT::VecOps::Sum(flags.tight);

    if (n_tight > 2) {
        return results;
    }

    if (n_tight == 2 && (ROOT::VecOps::Sum(flags.stand) > 2)) {
        return results;
    }

    for(int i = 0; i < n_muons; i++) {
        // Fast exit
        if ((!flags.tight[i]) && (!flags.hlt[i])) {
            continue;
        }
        
        if (DeltaR_flags.reco_flav_gen[i] != 1) {
            continue;
        }

        int l = DeltaR_flags.reco_idx_gen[i];
        
        // Not valid index -> fast exit
        if ((l < 0) || (l >= gen_eta.size())) {
            continue;
        }

        const bool pass_gen_tag = ((std::abs(DeltaR_flags.gen_pdg_idx[l]) == 13) && (DeltaR_flags.gen_status[l] == 1));
        if (!pass_gen_tag) {
            continue;
        }

        // DeltaR filter -> Only for MC
        float DeltaR_tag = ROOT::VecOps::DeltaR(kin.eta[i], gen_eta[l], kin.phi[i], gen_phi[l]);
        
        if (DeltaR_tag >= 0.3) {
            continue;
        }
        
        for (int j = 0; j < n_muons; j++) {
            // Fast exit -> same particle
            if (i == j || !flags.stand[j]) {
                continue;
            }
            // Opposite charge
            if (kin.charge[i] == kin.charge[j]) {
                continue;
            }
            
            // Matching with the generated muon
            if (DeltaR_flags.reco_flav_gen[j] != 1) {
                continue;
            }
        
            // Relative GenPart index for this reconstructed muon "j".
            int k = DeltaR_flags.reco_idx_gen[j];
        
            if ((k < 0) || (k >= gen_eta.size())) {
                continue;
            }

            const bool pass_gen_probe = ((std::abs(DeltaR_flags.gen_pdg_idx[k]) == 13) && (DeltaR_flags.gen_status[k] == 1)); 
            if (!pass_gen_probe) {
                continue;
            }

            // DeltaR filter.
            float DeltaR_probe = ROOT::VecOps::DeltaR(kin.eta[j], gen_eta[k], kin.phi[j], gen_phi[k]);
            if (DeltaR_probe >= 0.3) {
                continue;
            }
            
            const bool pass = ((kin.pt[j] > cfg_c.pt_cut) && (std::abs(kin.eta[j]) < cfg_c.eta_cut)) || (!(cfg_f.en_kinematics));
            if (!pass) {
                continue;
            }
            
            // Invariant mass (tag + probe)
            float mass = CalculateInvariantMass_Pair<float>(kin.pt[i], kin.pt[j], kin.eta[i], kin.eta[j], kin.phi[i], kin.phi[j],
            kin.mass[i], kin.mass[j]);

            // Invariant mass range
            if (((mass > cfg_c.mass_min) && (mass < cfg_c.mass_max)) || (!cfg_f.en_mass_window)) {
                // All probes (passed + failed)

                results.pt.push_back(kin.pt[j]);
                results.eta.push_back(kin.eta[j]);
                results.mll.push_back(mass);
                
                //results.tag_pt.push_back(kin.pt[i]);
                //results.tag_eta.push_back(kin.eta[i]); 
                
                if ((flags.global[j]) && (flags.iso[j] < 0.15)) {
                    // Passed probes
                    results.mask_pass.push_back(true);
                } else {
                    // Failed probes
                    results.mask_pass.push_back(false);
                }
            }
        }
    }
    return results;
}


// ------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CalculateAcceptance(ROOT::RDF::RNode node, const std::string& tag, int FSR) {
    std::vector<float> results;

    // Selection of true event: Z0 -> mu+ mu-
    ROOT::RDF::RNode node_acc = node
        .Define("Z0_Event", [FSR](const ROOT::RVec<Int_t>& pdg, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother) {
            return is_MC_Event(pdg, flags, mother, FSR);
        }, {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"})
        .Filter("Z0_Event", "Is_Z0_Event");

    // All Z0-> mu+ mu- events
    auto tot_gen_events = node_acc.Count();

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
        .Define("Z0_InvMass", [](const ROOT::RVec<float>& mu_pt, const ROOT::RVec<float>& amu_pt, const ROOT::RVec<float>& mu_eta, 
        const ROOT::RVec<float>& amu_eta, const ROOT::RVec<float>& mu_phi, const ROOT::RVec<float>& amu_phi, const ROOT::RVec<float>& mu_mass, 
        const ROOT::RVec<float>& amu_mass) {
                
            return CalculateInvariantMass_Pair<float>(mu_pt[0], amu_pt[0], mu_eta[0], amu_eta[0], mu_phi[0], amu_phi[0], 
                mu_mass[0], amu_mass[0]);
        }, {"Mu_pt", "AMu_pt", "Mu_eta", "AMu_eta", "Mu_phi", "AMu_phi", "Mu_mass", "AMu_mass"})

        .Filter("(Z0_InvMass > 60.0) && (Z0_InvMass < 120.0)", "Mass_cut");

        // Accepted events
        auto acc_events = node_acc.Count();
        
        auto report_node = node_acc.Report();
        report_node->Print();

        float num_tot = static_cast<float>(tot_gen_events.GetValue());
        float num_acc = static_cast<float>(acc_events.GetValue());

        float acc = (num_tot > 0) ? (num_acc / num_tot) : 0.0;
        results.push_back(acc);

        float acc_err = sqrt((acc * (1 - acc)) / num_tot);
        results.push_back(acc_err);

    return results;
}

ResolutionResults CalculateResolution(ROOT::RDF::RNode node, const config_struct& cfg) {
    ROOT::RDF::RNode node_resolution = node;
    
    std::string tag = cfg.resolution.quantity;
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

    std::vector<double> vector_bins = CreateBins(cfg.resolution.gen_bins, cfg.resolution.min, cfg.resolution.max, "linear");
    int n_bins = vector_bins.size() - 1;

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

    size_t n_entries = gen_vec.size();
    
    for (size_t i = 0; i < n_entries; i++) {

        auto upper_idx = std::upper_bound(vector_bins.begin(), vector_bins.end(), gen_vec[i]);
        int bin_idx = std::distance(vector_bins.begin(), upper_idx) - 1;

        if ((bin_idx >= 0) && (bin_idx < n_bins)) {
            bin_buffers[bin_idx].push_back(reco_vec[i]);
        }
    }

    for (int j = 0; j < n_bins; j++) {
        const auto& bin = bin_buffers[j];
        
        int entries = bin.size();
        result.events[j] = entries;

        if (entries == 0) {
            continue;
        }
        
        double sum = 0.0;
        
        for (double value : bin) {
            sum += value;
        }

        double mean = sum / entries;
        result.mean[j] = mean;

        if (entries > 1) {
            double sq_sum = 0;
            
            for (double value : bin) {
                double diff = value - mean;
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

    std::string good_muon = "true";

    if (cfg.flag_ES.en_kinematics) {
        good_muon += " && (Muon_pt > " + std::to_string(cfg.cut_ES.pt_cut) + " && abs(Muon_eta) < " + std::to_string(cfg.cut_ES.eta_cut)
         + ")";
    }

    if (cfg.flag_ES.en_tight_muon) {
        good_muon += " && (Muon_tightId == true)";
    }

    node_event = node_event
        .Define("GoodMuon", good_muon)

        .Define("GM_Pt", "Muon_pt[GoodMuon]")
        .Define("GM_Eta", "Muon_eta[GoodMuon]")
        .Define("GM_Phi", "Muon_phi[GoodMuon]")
        .Define("GM_Mass", "Muon_mass[GoodMuon]")
        .Define("GM_Charge", "Muon_charge[GoodMuon]");

    node_event = node_event
        .Define("EventPair", "Sum(GoodMuon) == 2")
        .Define("GoodEvent", "EventPair && (GM_Charge[0] != GM_Charge[1])");

    node_event = node_event
        .Define("InvariantMass", [] (const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi,
         const ROOT::RVec<float>& mass) -> float {
            
            return static_cast<float>(ROOT::VecOps::InvariantMass(pt, eta, phi, mass));
        
        }, {"GM_Pt", "GM_Eta", "GM_Phi", "GM_Mass"});

    std::string event_cut = "GoodEvent";
    if (cfg.flag_ES.en_mass_window) {
        event_cut += " && (InvariantMass > " + std::to_string(cfg.cut_ES.mass_min) + " && InvariantMass < " 
        + std::to_string(cfg.cut_ES.mass_max) + ")";
    }

    node_event = node_event
        .Filter(event_cut, "InvMass selection -> Good Event");

            node_event = node_event
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

EventHisto BuildEventHisto(ROOT::RDF::RNode node, const config_struct& cfg) {
    ROOT::RDF::RNode node_event = node;
    EventHisto histo;
      
    std::vector<double> bins_pt, bins_y, bins_phis;
    int n_pt, n_y, n_phis;

    if (cfg.unfold.use_custom_bins == true) {
        bins_pt = cfg.unfold.pt_bins.reco_vec;
        bins_y = cfg.unfold.y_bins.reco_vec;
        bins_phis = cfg.unfold.phis_bins.reco_vec;

    } else {
        const auto& pt = cfg.unfold.pt_bins;
        const auto& y = cfg.unfold.y_bins;
        const auto& phis = cfg.unfold.phis_bins;

        bins_pt = CreateBins(pt.reco_bins, pt.min, pt.max, pt.distribution, pt.split);
        bins_y = CreateBins(y.reco_bins, y.min, y.max, y.distribution, y.split);            
        bins_phis = CreateBins(phis.reco_bins, phis.min, phis.max, phis.distribution, phis.split);
    }

    // Mll bins
    int n_mll = 70;
    
    std::vector<double> mll_bins(n_mll + 1);
    double step = (120.0 - 60.0) / n_mll;

    for (int i = 0; i <= n_mll; i++) {
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