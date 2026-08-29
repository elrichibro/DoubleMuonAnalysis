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


ROOT::RDF::RNode node_Kin_cut(const std::string mask_name, ROOT::RDF::RNode node, const std::string pt, const std::string eta,
float pt_cut, float eta_cut);

struct ResultsTagAndProbe {
    std::vector<float> pt_pass;
    std::vector<float> pt_all;
    std::vector<float> eta_pass;
    std::vector<float> eta_all;
};

struct MuonKinematics {
    const ROOT::RVec<float>& pt;
    const ROOT::RVec<float>& eta;
    const ROOT::RVec<float>& phi;
    const ROOT::RVec<float>& mass;
};

struct MuonFlags {
    const ROOT::RVec<bool>& tag;
    const ROOT::RVec<bool>& probe;
    const ROOT::RVec<bool>& global;
    const ROOT::RVec<float>& iso;
};

ResultsTagAndProbe TagAndProbe(const MuonKinematics& kin, const MuonFlags& flags);

#endif