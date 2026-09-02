#include "Filters.h"

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

ResultsTagAndProbe CalculateTagAndProbe(const MuonKinematics_TP& kin, const MuonFlags_TP& flags, const flags_config cfg_f, 
const cuts_config cfg_c) {
    // Results container
    ResultsTagAndProbe results;

    // Number of particles in the event
    const unsigned int n_muons = kin.pt.size();

    results.pt_all.reserve(n_muons);
    results.pt_pass.reserve(n_muons);

    results.eta_all.reserve(n_muons);
    results.eta_pass.reserve(n_muons);

    results.mll_pass.reserve(n_muons);
    results.mll_all.reserve(n_muons);
    
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
            float mass = ROOT::VecOps::InvariantMass(ROOT::RVec<float>{kin.pt[i], kin.pt[j]}, 
            ROOT::RVec<float>{kin.eta[i], kin.eta[j]}, ROOT::RVec<float>{kin.phi[i], kin.phi[j]},
            ROOT::RVec<float>{kin.mass[i], kin.mass[j]});
            
            // Invariant mass range
            if (((mass > cfg_c.mass_min) && (mass < cfg_c.mass_max)) || (!cfg_f.en_mass_window)) {
                
                // All probes (passed + failed)
                results.pt_all.push_back(kin.pt[j]);
                results.eta_all.push_back(kin.eta[j]);
                results.mll_all.push_back(mass);

                if ((flags.global[j]) && (flags.iso[j] < 0.15)) {
                    
                    // Passed probes
                    results.pt_pass.push_back(kin.pt[j]);
                    results.eta_pass.push_back(kin.eta[j]);   
                    results.mll_pass.push_back(mass);
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
    for(int i = 0; i < n_muons_rec; i++) {
        
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
