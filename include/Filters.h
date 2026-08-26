#ifndef FILTERS_H
#define FILTERS_H

#include <map>
#include <cstdint>
#include <iostream>

#include <Rtypes.h>
#include "ROOT/RVec.hxx"

#include "Config.h"

/**
 * @brief Functor struct: Checks the validated runs from the not ones.
 * @param init_map Validation map.
 * @param run Run of the current event.
 * @param lum_block Luminosity block of the current event.
 * @return Returns true if the event is validated within che CMS standards, false if not.
*/
struct Validation_filter {
    const validation_type& val_map;

    Validation_filter(const validation_type& init_map) : val_map(init_map) {}

    bool operator() (const UInt_t run, const UInt_t lum_block) const;
};

/**
 * @brief Functor struct: Filters the GoodMuons.
 * @param init_structs Struct containing the physical cuts.
 * @param pt ROOT Vector. Stores the transverse momentum.
 * @param eta Pseudorapidity of the particles.
 * @param tight_id Identification muon flag.
 * @param iso Isolation variable.
 * @return Returns the mask for GoodMuons.
*/
struct GoodMuon_filter {
    const config_struct& cfg_struct;

    GoodMuon_filter(const config_struct& init_struct) : cfg_struct(init_struct) {}

    ROOT::RVec<bool> operator() (const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<bool>& tight_id, 
    const ROOT::RVec<float>& iso) const;
};

/**
 * @brief Selects the Z0 events in the MonteCarlo sample.
 * @param pdgId PDG identification index. (23 Z0)
 * @param flags Status flag stored bitwise.
 * @return Returns the mask fot Z0 true events.
*/
ROOT::RVec<bool> is_MC_Z0(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags);

/**
 * @brief Selects muons based on flags requirements.
 * @param pdgId PDG identification index. (13 Muon)
 * @param flags Status flag stored bitwise.
 * @return Returns a mask for true muons.
*/
ROOT::RVec<bool> is_MC_Muon(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags);

/**
 * @brief Accepts antimuons on physics requirements.
 * @param pdgId PDG identification index. (-13 AntiMuon)
 * @param flags Status flag stored bitwise.
 * @return Returns a mask for an antimuon selection.
*/
ROOT::RVec<bool> is_MC_AntiMuon(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags);

/**
 * @brief Generates a mask based on true event selection.
 * @param pdgId PDG identification index.
 * @param flags Status flag stored bitwise.
 * @return Returns true for a Z0 decaying into mu+mu-, else false.
*/
bool is_MC_Event(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags);

#endif