#ifndef FILTERS_H
#define FILTERS_H

#include <map>
#include <cstdint>
#include <iostream>

#include <Rtypes.h>
#include "ROOT/RVec.hxx"

#include "Config.h"

/**
 * @brief Functor struct used to filter the validated runs from the not ones.
 * @param init_map Is the validation map that contains the json file info about the validations runs and the luminosity blocks.
 * @param run Is the run of the current event.
 * @param lum_block Is the luminosity block of the current event.
 * @return Returns true if the event is validated within che CMS standards, false if not.
*/
struct Validation_filter {
    const validation_type& val_map;

    Validation_filter(const validation_type& init_map) : val_map(init_map) {}

    bool operator() (const UInt_t run, const UInt_t lum_block) const;
};

/**
 * @brief Functor struct used to filter the GoodMuons. The filter is based on transeverse momentum, pseudorapidity, thight flag identification
 * and isolation selection.
 * @param init_structs Is the struct containing the physical cut values stored in the config.json file.
 * @param pt Is a RVec storing the transverse momentum of the particles in the event.
 * @param eta Is pseudorapidity of the particles of the event.
 * @param tight_id Is the identification muon flag based on the number of hits in the inner tracker and muon stations, a requirement on the fit quality of the muons tracks and the spatial compatibility with the PV.
 * @param iso Is the isolation variable used to identify isolated leptons tracks.
 * @return Returns an RVec<bool> for all muons in the event. True if the event layed in the fiducial region, false if not.
*/
struct GoodMuon_filter {
    const config_struct& cfg_struct;

    GoodMuon_filter(const config_struct& init_struct) : cfg_struct(init_struct) {}

    ROOT::RVec<bool> operator() (const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<bool>& tight_id, 
    const ROOT::RVec<float>& iso) const;
};

ROOT::RVec<bool> is_MC_Z0(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags);
ROOT::RVec<bool> is_MC_Muon(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags);
ROOT::RVec<bool> is_MC_AntiMuon(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags);

bool is_MC_Event(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags);

#endif