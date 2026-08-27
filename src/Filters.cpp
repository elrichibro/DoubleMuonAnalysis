#include "Filters.h"
#include "Utils.h"

bool Validation_filter::operator()(const UInt_t run, const UInt_t lum_block) const {       
    thread_local static UInt_t last_run = 0;
    thread_local static UInt_t last_lum_block = 0;
    thread_local static bool last_decision = false;

    if ((run == last_run) && (lum_block == last_lum_block)) {
        return last_decision;
    }

    last_run = run;
    last_lum_block = lum_block;

    auto it = val_map.find(run);

    if (it != val_map.end()) {
        for (const auto& block : it->second) {
            if ((lum_block >= block.first) && (lum_block <= block.second)) {
                last_decision = true;
                return true;
            }
        }
    }

    last_decision = false;
    
    return false;
}

ROOT::RVec<bool> GoodMuon_filter::operator()(const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<bool>& tight_id, 
    const ROOT::RVec<float>& iso) const {
    
    ROOT::RVec<bool> mask(pt.size(), true);
    
    if (cfg_struct.flag.en_tight_muon) {
        mask = mask && tight_id;
    }
    
    if (cfg_struct.flag.en_kinematics) {
        mask = mask && (pt > cfg_struct.cut.pt_cut) && (abs(eta) < cfg_struct.cut.eta_cut);
    }
    
    if (cfg_struct.flag.en_isolation) {
        mask = mask && (iso < cfg_struct.cut.iso_cut);
    }
    
    return mask;
}


/* Flags summary:
    0 : isPrompt                -> SELECTED
    1 : isDecayedLeptonHadron
    2 : isTauDecayProduct
    3 : isPromptTauDecayProduct
    4 : isDirectTauDecayProduct
    5 : isDirectPromptTauDecayProduct
    6 : isDirectHadronDecayProduct
    7 : isHardProcess
    8 : fromHardProcess         -> SELECTED
    9 : isHardProcessTauDecayProduct
    10 : isDirectHardProcessTauDecayProduct
    11 : fromHardProcessBeforeFSR
    12 : isFirstCopy
    13 : isLastCopy             -> SELECTED
    14 : isLastCopyBeforeFSR


    Bitwie mask: 0, 8, 13 -> 2^13 + 2^8 + 2^0 = 8192 + 256 + 1 = 8449           
    0 x ( 0 0 1 0 )( 0 0 0 1 )( 0 0 0 0 )( 0 0 0 1 ) = 0x2101
*/

ROOT::RDF::RNode InvMass(ROOT::RDF::RNode node, const config_struct& cfg, const std::string& tag, int FSR) {
    // String construction
    const std::string mask = "good_mask_" + tag;
    const std::string pt = "good_Pt_" + tag;
    const std::string eta = "good_Eta_" + tag;
    const std::string phi = "good_Phi_" + tag;
    const std::string mass = "good_Mass_" + tag;
    const std::string mll = "m_ll_" + tag;

    // Node definition
    ROOT::RDF::RNode node_fsr = node
        .Filter([FSR](const ROOT::RVec<Int_t>& pdg, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother) {
            return is_MC_Event(pdg, flags, mother, FSR);
            },
        {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"}, "1. True Event" + tag)
        
        .Define("Mu_mask_" + tag, (FSR == 1) ? is_MC_Muon_bFSR : is_MC_Muon_aFSR, {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"})
        .Define("AMu_mask_" + tag, (FSR == 1) ? is_MC_AntiMuon_bFSR : is_MC_AntiMuon_aFSR, {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"})

        .Define(mask, "Mu_mask_" + tag + " || AMu_mask_" + tag)

        .Define(pt, "GenPart_pt[" + mask + "]")
        .Define(eta, "GenPart_eta[" + mask + "]")
        .Define(phi, "GenPart_phi[" + mask + "]")
        .Define(mass, "GenPart_mass[" + mask + "]")

        .Define(mll, mass_leptons<float>, {pt, eta, phi, mass});

    // Optional mass cut
    if (cfg.flag.en_mass_window) {
        node_fsr = node_fsr.Filter([&cfg] (const float mass) { return (mass > cfg.cut.mass_min) && (mass < cfg.cut.mass_max); },
        {mll}, "2. Z0 range selection" + tag);
    }

    return node_fsr;
}


ROOT::RDF::RNode InvMass_BeforeFSR(ROOT::RDF::RNode node, const config_struct& cfg) {
    auto node_b = node
        .Filter([](const ROOT::RVec<Int_t>& pdg, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother) {
            return is_MC_Event(pdg, flags, mother, 1);
            },
        {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"}, "1.b True Event")
        
        .Define("Mu_mask_bFSR", is_MC_Muon_bFSR, {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"})
        .Define("AMu_mask_bFSR", is_MC_AntiMuon_bFSR, {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"})

        .Define("good_mask_bFSR", "Mu_mask_bFSR || AMu_mask_bFSR")

        .Define("good_Pt_bFSR", "GenPart_pt[good_mask_bFSR]")
        .Define("good_Eta_bFSR", "GenPart_eta[good_mask_bFSR]")
        .Define("good_Phi_bFSR", "GenPart_phi[good_mask_bFSR]")
        .Define("good_Mass_bFSR", "GenPart_mass[good_mask_bFSR]");

    ROOT::RDF::RNode node_inv_mass_bFSR = node_b
    .Define("m_ll_bFSR", mass_leptons<float>, {"good_Pt_bFSR", "good_Eta_bFSR", "good_Phi_bFSR", "good_Mass_bFSR"});

    if (cfg.flag.en_mass_window) {
        node_inv_mass_bFSR = node_inv_mass_bFSR.Filter([&cfg] (const float mass) { return (mass > cfg.cut.mass_min) && (mass < cfg.cut.mass_max); },
        {"m_ll_bFSR"}, "2.b Z0 range selection");
    }

    return node_inv_mass_bFSR;
}


ROOT::RDF::RNode InvMass_AfterFSR(ROOT::RDF::RNode node, const config_struct& cfg) {

    auto node_a = node
        .Filter([](const ROOT::RVec<Int_t>& pdg, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother) {
            return is_MC_Event(pdg, flags, mother, 2);
        }, {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_genPartIdxMother"}, "1.a True Event")      
        
        .Define("Mu_mask_aFSR", is_MC_Muon_aFSR, {"GenPart_pdgId", "GenPart_statusFlags"})
        .Define("AMu_mask_aFSR", is_MC_AntiMuon_aFSR, {"GenPart_pdgId", "GenPart_statusFlags"})

        .Define("good_mask_aFSR", "Mu_mask_aFSR || AMu_mask_aFSR")

        .Define("good_Pt_aFSR", "GenPart_pt[good_mask_aFSR]")
        .Define("good_Eta_aFSR", "GenPart_eta[good_mask_aFSR]")
        .Define("good_Phi_aFSR", "GenPart_phi[good_mask_aFSR]")
        .Define("good_Mass_aFSR", "GenPart_mass[good_mask_aFSR]");

    ROOT::RDF::RNode node_InvMass_aFSR = node_a
    .Define("m_ll_aFSR", mass_leptons<float>, {"good_Pt_aFSR", "good_Eta_aFSR", "good_Phi_aFSR", "good_Mass_aFSR"});

    if (cfg.flag.en_mass_window) {
        node_InvMass_aFSR = node_InvMass_aFSR.Filter([&cfg] (const float mass) { return (mass > cfg.cut.mass_min) && (mass < cfg.cut.mass_max); },
        {"m_ll_aFSR"}, "2.a Z0 range selection");
    }

    return node_InvMass_aFSR;
}


ROOT::RVec<bool> is_MC_Z0(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags) {
    return ((pdgId == 23) && ((flags & 0x2101) == 0x2101));
}

int get_MC_Z0_idx(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags) {
    auto z_mask = is_MC_Z0(pdgId, flags);
    
    if (!ROOT::VecOps::Any(z_mask)) {
        return -1;
    }
    
    return static_cast<int>(ROOT::VecOps::ArgMax(z_mask));
}

ROOT::RVec<bool> is_MC_Muon_bFSR(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother_id) {
    int z_idx = get_MC_Z0_idx(pdgId, flags);
    
    if (z_idx < 0) {
        return ROOT::RVec<bool>(pdgId.size(), false);
    }
    
    return ((pdgId == 13) && ((flags & 0x181) == 0x181) && (mother_id == z_idx));
}

ROOT::RVec<bool> is_MC_AntiMuon_bFSR(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother_id) {
    int z_idx = get_MC_Z0_idx(pdgId, flags);
    
    if (z_idx < 0) {
        return ROOT::RVec<bool>(pdgId.size(), false);
    }
    
    return ((pdgId == -13) && ((flags & 0x181) == 0x181) && (mother_id == z_idx));
}

ROOT::RVec<bool> is_MC_Muon_aFSR(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother_id) {
    return ((pdgId == 13) && ((flags & 0x2101) == 0x2101));
}

ROOT::RVec<bool> is_MC_AntiMuon_aFSR(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother_id) {
    return ((pdgId == -13) && ((flags & 0x2101) == 0x2101));
}

bool is_MC_Event(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother_id, const int FSR) {
    
    if (FSR == 1) {
        return 
        ((ROOT::VecOps::Sum(is_MC_Z0(pdgId, flags)) == 1) && 
        (ROOT::VecOps::Sum(is_MC_Muon_bFSR(pdgId, flags, mother_id)) == 1) && 
        (ROOT::VecOps::Sum(is_MC_AntiMuon_bFSR(pdgId, flags, mother_id)) == 1)); 
    } else if (FSR == 2) {
        return 
        ((ROOT::VecOps::Sum(is_MC_Z0(pdgId, flags)) == 1) && 
        (ROOT::VecOps::Sum(is_MC_Muon_aFSR(pdgId, flags, mother_id)) == 1) && 
        (ROOT::VecOps::Sum(is_MC_AntiMuon_aFSR(pdgId, flags, mother_id)) == 1)); 
    }

    return false;
}

