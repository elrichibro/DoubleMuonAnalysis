#include "Filters.h"

#include "Utils.h"

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

ROOT::RDF::RNode ApplyKinMuonFilter(ROOT::RDF::RNode node, const std::string& mask_name, const std::string& pt_column, const std::string& eta_column,
 float pt_cut, float eta_cut) {
    
    ROOT::RDF::RNode node_kin_cut = node
        .Define(mask_name, [pt_cut, eta_cut](const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta) {
            
            // Physical cut
            return ((pt > pt_cut) && (abs(eta) < eta_cut));

        }, {pt_column, eta_column});

    return node_kin_cut;
}

// ------------------------------------------------------------------------------------------------------------------------------------
// TagAndProbe selection - Data
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
     
    for(int i = 0; i < n_muons; i++) {
        // Fast exit
        if (!flags.tight[i]) {
            continue;
        }

        for (int j = 0; j < n_muons; j++) {
            // Fast exit
            if (i == j || !flags.stand[j]) {
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
                    
                    results.pt.push_back(kin.pt[j]);
                    results.eta.push_back(kin.eta[j]);   
                    
                    results.mll.push_back(mass);
                    
                    //results.tag_pt.push_back(kin.pt[i]);
                    //results.tag_eta.push_back(kin.eta[i]);         
                
                // All probes (passed + failed)
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
// TagAndProbe selection - MonteCarlo
// ------------------------------------------------------------------------------------------------------------------------------------

ResultsTagAndProbe CalculateTagAndProbe_MC(const MuonKinematics_TP& kin, const MuonFlags_TP& flags, const flags_config cfg_f, 
const cuts_config cfg_c, const MuonFlags_RM& DeltaR_flags, const ROOT::RVec<float> gen_eta, const ROOT::RVec<float> gen_phi) {
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
    
    for(int i = 0; i < n_muons; i++) {
        // Fast exit
        if (!flags.tight[i]) {
            continue;
        }
        
        if (DeltaR_flags.gen_flav_rec[i] != 1) {
            continue;
        }

        int l = DeltaR_flags.pair_idx_rec[i];
        
        // Not valid index -> fast exit
        if ((l < 0) || (l >= gen_eta.size())) {
            continue;
        }

        const bool pass_gen_tag = ((std::abs(DeltaR_flags.pdg_id_gen[l]) == 13) && (DeltaR_flags.status_gen[l] == 1));
        if (!pass_gen_tag) {
            continue;
        }

        float DeltaR_tag = ROOT::VecOps::DeltaR(kin.eta[i], gen_eta[l], kin.phi[i], gen_phi[l]);
        
        if (DeltaR_tag >= 0.3) {
            continue;
        }
        

        for (int j = 0; j < n_muons; j++) {
            // Fast exit
            if (i == j || !flags.stand[j]) {
                continue;
            }
            // Loop into good probe muons

            
            if (DeltaR_flags.gen_flav_rec[j] != 1) {
                continue;
            }
        
            // Relative GenPart index for this reconstructed muon "j".
            int k = DeltaR_flags.pair_idx_rec[j];
        
            if ((k < 0) || (k >= gen_eta.size())) {
                continue;
            }

            const bool pass_gen_probe = ((std::abs(DeltaR_flags.pdg_id_gen[k]) == 13) && (DeltaR_flags.status_gen[k] == 1)); 
            if (!pass_gen_probe) {
                continue;
            }

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
// Response Matrix calculus
// ------------------------------------------------------------------------------------------------------------------------------------

ResultsRespMatrix CalculateRespMatrix(const MuonKinematics_RM& kin, const MuonFlags_RM& flags, const flags_config cfg_f, 
const cuts_config cfg_c) {
    
    ResultsRespMatrix results;

    // Number reconstructed muons
    const int n_muons_rec = kin.pt_rec.size();
    
    // Loop on reconstructed muons.
    for (int i = 0; i < n_muons_rec; i++) {
        
        // gen_flav_rec == 1 -> GenPart muon is : prompt muon. (!= -> fast exit)
        if (flags.gen_flav_rec[i] != 1) {
            continue;
        }
        
        // Relative GenPart index for this reconstructed muon "i".
        int j = flags.pair_idx_rec[i];
        
        // Not valid index -> fast exit
        if ((j < 0) || (j >= kin.pt_gen.size())) {
            continue;
        }

        // (Pdg index == 13/-13) + (Stable status of GenPart) 
        if ((std::abs(flags.pdg_id_gen[j]) == 13) && (flags.status_gen[j] == 1)) {
            
            const bool pass = ((kin.pt_gen[j] > cfg_c.pt_cut) && (kin.pt_rec[i] > cfg_c.pt_cut) && 
            (std::abs(kin.eta_gen[j]) < cfg_c.eta_cut) && (std::abs(kin.eta_rec[i]) < cfg_c.eta_cut)) || (!(cfg_f.en_kinematics));
            
            if (pass) {    
                results.pt_gen_RM.push_back(kin.pt_gen[j]);
                results.pt_rec_RM.push_back(kin.pt_rec[i]);

                results.eta_gen_RM.push_back(kin.eta_gen[j]);                
                results.eta_rec_RM.push_back(kin.eta_rec[i]);
            }
        }
    }
    return results;
}

// ------------------------------------------------------------------------------------------------------------------------------------
// DataDivision masks maker
// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RNode ApplyPt_DataDivision(ROOT::RDF::RNode node, const std::vector<float>& pt_bins, const std::string& column_name) {
    ROOT::RDF::RNode div_pt_node = node;
    
    for (int i = 0; i < pt_bins.size() - 1; i++) {
        std::string name_mask = column_name + "_DeltaMask" + std::to_string(i);
        std::string name_column = column_name + "_Delta" + std::to_string(i);

        float min_pt = pt_bins[i];
        float max_pt = pt_bins[i + 1];
        
        div_pt_node = div_pt_node
            .Define(name_mask, [min_pt, max_pt](const ROOT::RVec<float>& pt) {
                return ((pt >= min_pt) && (pt < max_pt)); 
            }, {column_name});
    
        div_pt_node = div_pt_node
            .Define(name_column, [](const ROOT::RVec<float>& pt, const ROOT::RVec<bool>& mask) {
                return pt[mask];
            }, {column_name, name_mask} );
    }

    return div_pt_node;
}

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RNode ApplyEta_DataDivision(ROOT::RDF::RNode node, const std::vector<float>& eta_bins, const std::string& column_name) {
    ROOT::RDF::RNode div_eta_node = node;
    
    for(int i = 0; i < eta_bins.size() - 1; i++) {
        std::string name_mask = column_name + "_DeltaMask" + std::to_string(i);
        std::string name_column = column_name + "_Delta" + std::to_string(i);

        float min_eta = eta_bins[i];
        float max_eta = eta_bins[i + 1];
        
        div_eta_node = div_eta_node
            .Define(name_mask, [min_eta, max_eta](const ROOT::RVec<float>& eta) {
                return ((eta >= min_eta) && (eta < max_eta)); 
            }, {column_name});
    
        div_eta_node = div_eta_node
            .Define(name_column, [](const ROOT::RVec<float>& eta, const ROOT::RVec<bool>& mask) {
                return eta[mask];
            }, {column_name, name_mask} );
    }

    return div_eta_node;
}

// ------------------------------------------------------------------------------------------------------------------------------------
// MC Template Maker
// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RResultPtr<TH3D> TemplateMaker_MC(ROOT::RDF::RNode node, const config_struct& cfg, const bool mask) {
    ROOT::RDF::RNode node_hist = node;

    std::vector<float> pt_bins = cfg.analysis.pt_bins;
    std::vector<float> eta_bins = cfg.analysis.eta_bins;
    std::vector<float> mll_bins(100);

    float step = (120.0 - 60.0) / 99.0;
    for (int i = 0; i < 100; i++) {
        mll_bins[i] = 60.0 + (i * step);
    }

    ROOT::RDF::TH3DModel model("h3", "3D Histogram MC", eta_bins.size() - 1, eta_bins.data(), 
    pt_bins.size() - 1, pt_bins.data(), mll_bins.size() - 1, mll_bins.data() );

    
    if (mask) {
        node_hist = node_hist
            .Define("MC_Probe_Pt_Pass", "MC_Pt[MC_Mask_Pass]")
            .Define("MC_Probe_Eta_Pass", "MC_Eta[MC_Mask_Pass]")
            .Define("MC_Mll_Pass", "MC_Eta[MC_Mask_Pass]");

        auto h3 = node_hist.Histo3D(model, "MC_Probe_Eta_Pass", "MC_Probe_Pt_Pass", "MC_Mll_Pass");
        return h3;
    } else {
        node_hist = node_hist
            .Define("MC_Probe_Pt_Fail", "MC_Pt[!MC_Mask_Pass]")
            .Define("MC_Probe_Eta_Fail", "MC_Eta[!MC_Mask_Pass]")
            .Define("MC_Mll_Fail", "MC_Eta[!MC_Mask_Pass]");

        auto h3 = node_hist.Histo3D(model, "MC_Probe_Eta_Fail", "MC_Probe_Pt_Fail", "MC_Mll_Fail");
        return h3;
    }
}