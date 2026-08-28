#include "Checks.h"
#include "Filters.h"
#include "Utils.h"
#include <cstdint>


ROOT::RDF::RNode InvMass(ROOT::RDF::RNode node, const std::string& tag, int FSR) {
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

    return node_fsr;
}



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

        .Define("Muon_inv_mass", mass_leptons<float>, {"good_Muon_pt", "good_Muon_eta", "good_Muon_phi", "good_Muon_mass"});
    return node_rec;
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


    Bitwise mask1: 0, 8, 13 -> 2^13 + 2^8 + 2^0 = 8192 + 256 + 1 = 8449           
    0 x ( 0 0 1 0 )( 0 0 0 1 )( 0 0 0 0 )( 0 0 0 1 ) = 0x2101

    Bitwise mask2: 0, 7, 8 -> 2^8 + 2^7 + 2^0 = 256 + 128 + 1 = 385           
    0 x ( 0 0 0 1 )( 1 0 0 0 )( 0 0 0 1 ) = 0x181
*/
bool is_MC_Event(const ROOT::RVec<int>& pdgId, const ROOT::RVec<int>& flags, const ROOT::RVec<int>& mother_id, const int FSR) {
    const uint16_t num_parts = pdgId.size();
    
    unsigned int count_Z0 = 0;
    unsigned int count_Mu = 0;
    unsigned int count_aMu = 0;

    const uint16_t mask_aFSR = 0x2101;
    const uint16_t mask_bFSR = 0x181;

    if (FSR == 1) {
        int z_idx = -1;

        for (uint8_t i = 0; i < num_parts; i++) {
            if ((pdgId[i] == 23) && ((flags[i] & mask_aFSR) == mask_aFSR)) {
                count_Z0++;
                z_idx = static_cast<int>(i);
            }
        }

        if (count_Z0 != 1 || z_idx < 0) {
            return false;
        }

        for (size_t i = 0; i < num_parts; i++) {
            if (((flags[i] & mask_bFSR) == mask_bFSR) && (mother_id[i] == z_idx)) {
                if (pdgId[i] == 13) {
                    count_Mu++;
                } else if (pdgId[i] == -13) {
                    count_aMu++;
                }
            }
        }

        return ((count_Mu == 1) && (count_aMu == 1));

    } else if (FSR == 2) {
        for (size_t i = 0; i < num_parts; i++) {
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
