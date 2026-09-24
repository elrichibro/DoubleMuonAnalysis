#include "ControlCheck.h"

#include "Utils.h"


ROOT::RDF::RNode ApplyValidationFilter(ROOT::RDF::RNode node, const validation_type& val_map, const std::string& run_name, 
const std::string& block_name) {
    
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

/* Flags summary in MonteCarlo "Generated particles" collection:
    0 : isPrompt                -> SELECTED
    1 : isDecayedLeptonHadron
    2 : isTauDecayProduct
    3 : isPromptTauDecayProduct
    4 : isDirectTauDecayProduct
    5 : isDirectPromptTauDecayProduct
    6 : isDirectHadronDecayProduct
    7 : isHardProcess           -> SELECTED
    8 : fromHardProcess         -> SELECTED
    9 : isHardProcessTauDecayProduct
    10 : isDirectHardProcessTauDecayProduct
    11 : fromHardProcessBeforeFSR
    12 : isFirstCopy
    13 : isLastCopy             -> SELECTED
    14 : isLastCopyBeforeFSR


    Muon After FSR
    Bitwise mask2: 0, 7, 8 -> 2^8 + 2^7 + 2^0 = 256 + 128 + 1 = 385
    0 x ( 0 0 0 1 )( 1 0 0 0 )( 0 0 0 1 ) = 0x181

    Muon/Z0 Before FSR
    Bitwise mask1: 0, 8, 13 -> 2^13 + 2^8 + 2^0 = 8192 + 256 + 1 = 8449
    0 x ( 0 0 1 0 )( 0 0 0 1 )( 0 0 0 0 )( 0 0 0 1 ) = 0x2101
*/

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RNode CalculateInvMass(ROOT::RDF::RNode node, const std::string& tag, int FSR) {
    // Column names construction
    const std::string mask = "good_Mask_" + tag;
    const std::string pt = "good_Pt_" + tag;
    const std::string eta = "good_Eta_" + tag;
    const std::string phi = "good_Phi_" + tag;
    const std::string mass = "good_Mass_" + tag;
    const std::string m_ll = "InvMass_" + tag;

    // Node definition
    ROOT::RDF::RNode node_FSR = node
        // True Event Filter: Z0 -> mu+ mu-
        .Filter([FSR](const ROOT::RVec<Int_t>& pdg, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother) {
            return is_MC_Event(pdg, flags, mother, FSR);
            },
        {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"}, "1. True Event" + tag)
        
        // Mu / AntiMu maks
        .Define("Mu_mask_" + tag, (FSR == 1) ? is_MC_Muon_bFSR : is_MC_Muon_aFSR, {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"})
        .Define("AMu_mask_" + tag, (FSR == 1) ? is_MC_AntiMuon_bFSR : is_MC_AntiMuon_aFSR, {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"})

        // Merge of both masks
        .Define(mask, "Mu_mask_" + tag + " || AMu_mask_" + tag)

        // Selection of Mu/AntiMu pt, eta, phi, mass
        .Define(pt, "GenPart_pt[" + mask + "]")
        .Define(eta, "GenPart_eta[" + mask + "]")
        .Define(phi, "GenPart_phi[" + mask + "]")
        .Define(mass, "GenPart_mass[" + mask + "]")

        // Invariant Mass of the muon pair
        .Define(m_ll, CalculateInvariantMass<float>, {pt, eta, phi, mass});

    return node_FSR;
}

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RNode node_recMC(ROOT::RDF::RNode node) {
    
    ROOT::RDF::RNode node_rec = node
        .Define("Tight_muon", "Muon_tightId == true")
        .Filter("Sum(Tight_muon) == 2", "1. Tight muon selection")
        
        .Define("good_Muon_pt", "Muon_pt[Tight_muon]")
        .Define("good_Muon_eta", "Muon_eta[Tight_muon]")
        .Define("good_Muon_phi", "Muon_phi[Tight_muon]")
        .Define("good_Muon_mass", "Muon_mass[Tight_muon]")
        .Define("tight_muon_Charge", "Muon_charge[Tight_muon]")

        .Filter([](const ROOT::RVec<int>& charge) {
            return charge[0] != charge[1];
        }, {"tight_muon_Charge"}, "2. Dimuon")
        
        .Define("LooseMuon", "Muon_pt > 10.0 && abs(Muon_eta) < 2.4 && Muon_looseId")
        .Define("LooseElectron", "Electron_pt > 10.0 && abs(Electron_eta) < 2.5 && Electron_cutBased == 2")
        .Filter("Sum(LooseMuon) + Sum(LooseElectron) == 2", "4. Multiboson background")

        .Define("Muon_inv_mass", CalculateInvariantMass<float>, {"good_Muon_pt", "good_Muon_eta", "good_Muon_phi", "good_Muon_mass"});
    return node_rec;
}

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RVec<bool> is_MC_Z0(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags) {
    return ((pdgId == 23) && ((flags & 0x2101) == 0x2101));
}

// ------------------------------------------------------------------------------------------------------------------------------------

int get_MC_Z0_idx(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags) {
    auto z_mask = is_MC_Z0(pdgId, flags);
    
    if (!ROOT::VecOps::Any(z_mask)) {
        return -1;
    }
    
    return static_cast<int>(ROOT::VecOps::ArgMax(z_mask));
}

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RVec<bool> is_MC_Muon_bFSR(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother_id) {
    int z_idx = get_MC_Z0_idx(pdgId, flags);
    
    if (z_idx < 0) {
        return ROOT::RVec<bool>(pdgId.size(), false);
    }
    
    return ((pdgId == 13) && ((flags & 0x181) == 0x181) && (mother_id == z_idx));
}

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RVec<bool> is_MC_AntiMuon_bFSR(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother_id) {
    int z_idx = get_MC_Z0_idx(pdgId, flags);
    
    if (z_idx < 0) {
        return ROOT::RVec<bool>(pdgId.size(), false);
    }
    
    return ((pdgId == -13) && ((flags & 0x181) == 0x181) && (mother_id == z_idx));
}

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RVec<bool> is_MC_Muon_aFSR(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>&) {
    return ((pdgId == 13) && ((flags & 0x2101) == 0x2101));
}

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RVec<bool> is_MC_AntiMuon_aFSR(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>&) {
    return ((pdgId == -13) && ((flags & 0x2101) == 0x2101));
}

// ------------------------------------------------------------------------------------------------------------------------------------

bool is_MC_Event(const ROOT::RVec<int>& pdgId, const ROOT::RVec<int>& flags, const ROOT::RVec<int>& mother_id, const int FSR) {
    
    const int num_parts = pdgId.size();
    int count_Z0 = 0;
    int count_Mu = 0;
    int count_aMu = 0;

    const uint16_t mask_aFSR = 0x2101;
    const uint16_t mask_bFSR = 0x181;

    // Before FSR status
    if (FSR == 1) {
        int z_idx = -1;

        for (int i = 0; i < num_parts; i++) {
            if ((pdgId[i] == 23) && ((flags[i] & mask_aFSR) == mask_aFSR)) {
                count_Z0++;
                z_idx = i;
            }
        }

        if (count_Z0 != 1 || z_idx < 0) {
            return false;
        }

        for (int i = 0; i < num_parts; i++) {
            if (((flags[i] & mask_bFSR) == mask_bFSR) && (mother_id[i] == z_idx)) {
                if (pdgId[i] == 13) {
                    count_Mu++;
                } else if (pdgId[i] == -13) {
                    count_aMu++;
                }
            }
        }

        return ((count_Mu == 1) && (count_aMu == 1));

    // After FSR
    } else if (FSR == 2) {
        for (int i = 0; i < num_parts; i++) {
            if ((flags[i] & mask_aFSR) == mask_aFSR) {
                if (pdgId[i] == 23) {
                    count_Z0++;
                } else if (pdgId[i] == 13) {
                    count_Mu++;
                } else if (pdgId[i] == -13) {
                    count_aMu++;
                }
            }
        }
        return ((count_Z0 == 1) && (count_Mu == 1) && (count_aMu == 1));
    }
    return false;
}