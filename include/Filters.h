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
 * @brief Checks the validated runs from the not ones.
 * @param val_map Validation map.
 * @param run Run of the current event.
 * @param lum_block Luminosity block of the current event.
 * @return Returns the validated node-dataset.
*/
ROOT::RDF::RNode ApplyValidationFilter(ROOT::RDF::RNode node, const validation_type& val_map);


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

/**
 * @brief Defines a new column in the dataset that represents the mask of particles that pass the kinematic cuts.
 * @param mask_name Name mask.
 * @param node RDF node.
 * @param pt Column name for particle transverse momentum.
 * @param eta Column name for particle pseudorapidity.
 * @param pt_cut Minimum transverse momentum threshold.
 * @param eta_cut Pseudorapidity range.
 * @return Returns the node containing the bool mask.
*/
ROOT::RDF::RNode node_Kin_cut(const std::string mask_name, ROOT::RDF::RNode node, const std::string pt, const std::string eta,
float pt_cut, float eta_cut);


/// @brief Transverse momentum and psudorapidity vectors results after T&P
struct ResultsTagAndProbe {
    std::vector<float> pt_pass;// Transverse momentum of probes that pass the selection.
    std::vector<float> pt_all;// Transverse momentum of all probes.
    std::vector<float> eta_pass;// Pseudorapidity of probes that pass the selection.
    std::vector<float> eta_all;// Pseudorapidity of all probes. 
};


/// @brief Muon kinematics variables for selection
struct MuonKinematics {
    const ROOT::RVec<float>& pt;// Transverse momentum
    const ROOT::RVec<float>& eta;// Pseudorapidity
    const ROOT::RVec<float>& phi;// Angular variable
    const ROOT::RVec<float>& mass;// Reconstructed mass
};


/// @brief Muon flags for TagAndProbe selections
struct MuonFlags {
    const ROOT::RVec<bool>& tag;// Muon_tightId flag.
    const ROOT::RVec<bool>& probe;// Muon_isStandalone flag.
    const ROOT::RVec<bool>& global;// Muon_isGlobal flag.
    const ROOT::RVec<float>& iso;// Isolation parameter.
};


/**
 * @brief TagAndProbe function selection
 * @param kin Muon kinematic event values.
 * @param flags Muon event flags.
 * @return Struct containing vectors of passed and total probe kinematics.
*/
ResultsTagAndProbe TagAndProbe(const MuonKinematics& kin, const MuonFlags& flags);

#endif