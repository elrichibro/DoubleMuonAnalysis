#ifndef FILTERS_H
#define FILTERS_H

#include <map>
#include <cstdint>
#include <iostream>

#include <Rtypes.h>
#include "ROOT/RVec.hxx"

#include "Config.h"

struct Validation_filter {
    const validation_type& val_map;

    Validation_filter(const validation_type& init_map) : val_map(init_map) {}

    bool operator() (const UInt_t run, const UInt_t lum_block) const;
};

struct GoodMuon_filter {
    const config_struct& cfg_struct;

    GoodMuon_filter(const config_struct& init_struct) : cfg_struct(init_struct) {}

    ROOT::RVec<bool> operator() (const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<bool>& tight_id, 
    const ROOT::RVec<float>& iso) const;
};

#endif