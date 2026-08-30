#include "Filters.h"
#include "Utils.h"

// #################################################################

ROOT::RDF::RNode ApplyValidationFilter(ROOT::RDF::RNode node, const validation_type& val_map, const std::string run_name, const std::string block_name) {
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

// #################################################################

ROOT::RVec<bool> GoodMuon_filter::operator()(const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<bool>& tight_id) const {
    
    ROOT::RVec<bool> mask(pt.size(), true);
    
    if (cfg_struct.flag.en_tight_muon) {
        mask = mask && tight_id;
    }
    
    if (cfg_struct.flag.en_kinematics) {
        mask = mask && (pt > cfg_struct.cut.pt_cut) && (abs(eta) < cfg_struct.cut.eta_cut);
    }
    
    return mask;
}

// #################################################################

ROOT::RDF::RNode node_Kin_cut(const std::string mask_name, ROOT::RDF::RNode node, const std::string pt, 
    const std::string eta, float pt_cut, float eta_cut) {
    
    ROOT::RDF::RNode node_kin_cut = node
    .Define(mask_name, [pt_cut, eta_cut](const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta) {
        
        // Physical cut
        return ((pt > pt_cut) && (abs(eta) < eta_cut));

    }, {pt, eta});

    return node_kin_cut;
}

// #################################################################

ResultsTagAndProbe TagAndProbe(const MuonKinematics& kin, const MuonFlags& flags) {
    // Results container
    ResultsTagAndProbe results;
    
    // Number of particles in the event
    const unsigned int n_muons = kin.pt.size();
    
    for(int i = 0; i < n_muons; i++) {
        // Fast exit
        if (!flags.tag[i]) {
            continue;
        
        // Good condition
        } else if (flags.tag[i]) {

            for (int j = 0; j < n_muons; j++) {
                // Fast exit
                if (i == j || !flags.probe[j]) {
                    continue;

                // Loop into good probe muons
                } else if ((j != i) && flags.probe[j]) {
                    // Invariant mass (tag + probe)
                    float mass = ROOT::VecOps::InvariantMass(ROOT::RVec<float>{kin.pt[i], kin.pt[j]}, 
                        ROOT::RVec<float>{kin.eta[i], kin.eta[j]}, ROOT::RVec<float>{kin.phi[i], kin.phi[j]},
                        ROOT::RVec<float>{kin.mass[i], kin.mass[j]});
                    
                    // Invariant mass range
                    if((mass > 60) && (mass < 120)) {
                        
                        // All probes (passed + failed)
                        results.pt_all.push_back(kin.pt[j]);
                        results.eta_all.push_back(kin.eta[j]);

                        if ((flags.global[j]) && (flags.iso[j] < 0.15)) {
                            // Passed probes
                            results.pt_pass.push_back(kin.pt[j]);
                            results.eta_pass.push_back(kin.eta[j]);
                        } else {
                            continue;
                        }
                    } else {
                        continue;
                    }
                }
            }
        }
    }
    
    return results;
}