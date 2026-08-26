#include "Filters.h"

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

ROOT::RVec<bool> is_MC_Z0(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags) {
    return ((pdgId == 23) && ((flags & 0x2101) == 0x2101));
}

ROOT::RVec<bool> is_MC_Muon(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags) {
    return ((pdgId == 13) && ((flags & 0x2101) == 0x2101));
}

ROOT::RVec<bool> is_MC_AntiMuon(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags) {
    return ((pdgId == -13) && ((flags & 0x2101) == 0x2101));
}

bool is_MC_Event(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags) {

    return ((ROOT::VecOps::Sum(is_MC_Z0(pdgId, flags)) == 1) && 
    (ROOT::VecOps::Sum(is_MC_Muon(pdgId, flags)) == 1) && (ROOT::VecOps::Sum(is_MC_AntiMuon(pdgId, flags)) == 1)); 
}

