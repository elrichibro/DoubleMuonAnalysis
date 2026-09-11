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
            // Fast exit -> same particle
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
            
            // Matching with the generated muon
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
            
            // Kinematical cut
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

std::vector<float> CalculateAcceptance(ROOT::RDF::RNode node, const std::string& tag, int FSR) {
    std::vector<float> results;
    
    ROOT::RDF::RNode node_acc = node
        // True Event Filter: Z0 -> mu+ mu-
        .Define("Z0_Event", [FSR](const ROOT::RVec<Int_t>& pdg, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother) {
            return is_MC_Event(pdg, flags, mother, FSR);
        }, {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"})
        .Filter("Z0_Event", "Is_Z0_mumu");

    auto tot_gen_events = node_acc.Count();

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

        .Filter("(Mu_pt[0] > 25) && (AMu_pt[0] > 25)", "Pt_cut")
        .Filter("(abs(Mu_eta[0]) < 2.4) && (abs(AMu_eta[0]) < 2.4)", "Eta_cut")

        .Define("Z0_InvMass", [](const ROOT::RVec<float>& mu_pt, const ROOT::RVec<float>& amu_pt, const ROOT::RVec<float>& mu_eta, 
        const ROOT::RVec<float>& amu_eta, const ROOT::RVec<float>& mu_phi, const ROOT::RVec<float>& amu_phi, const ROOT::RVec<float>& mu_mass, 
        const ROOT::RVec<float>& amu_mass) {
                
            return CalculateInvariantMass_Pair<float>(mu_pt[0], amu_pt[0], mu_eta[0], amu_eta[0], mu_phi[0], amu_phi[0], 
                mu_mass[0], amu_mass[0]);
        }, {"Mu_pt", "AMu_pt", "Mu_eta", "AMu_eta", "Mu_phi", "AMu_phi", "Mu_mass", "AMu_mass"})

        .Filter("Z0_InvMass > 60.0 && Z0_InvMass < 120.0", "Mass_cut");

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