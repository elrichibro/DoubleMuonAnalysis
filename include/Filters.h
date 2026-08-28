#ifndef FILTERS_H
#define FILTERS_H

#include <map>
#include <cstdint>
#include <iostream>

#include <Rtypes.h>
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>


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
 * @return Returns the mask for GoodMuons.
*/
struct GoodMuon_filter {
    const config_struct& cfg_struct;

    GoodMuon_filter(const config_struct& init_struct) : cfg_struct(init_struct) {}

    ROOT::RVec<bool> operator() (const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<bool>& tight_id) const;
};

#endif