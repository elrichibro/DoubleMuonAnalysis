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
 * @param run_name Run column name.
 * @param block_name Luminosity block column name.
 * @return Returns the validated node-dataset.
*/
ROOT::RDF::RNode ApplyValidationFilter(ROOT::RDF::RNode node, const validation_type& val_map, const std::string run_name, const std::string block_name);

/**
 * @brief Defines a new column in the dataset that represents the mask of particles that pass the kinematic cuts.
 * @param node RDF node.
 * @param mask_name Name of the new column-mask.
 * @param pt_name Column name for particle transverse momentum.
 * @param eta_name Column name for particle pseudorapidity.
 * @param pt_cut Minimum transverse momentum threshold.
 * @param eta_cut Pseudorapidity range.
 * @return Returns the node containing the bool mask.
*/
ROOT::RDF::RNode ApplyKinMuonFilter(ROOT::RDF::RNode node, const std::string& mask_name, const std::string& pt_name, const std::string& eta_name,
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