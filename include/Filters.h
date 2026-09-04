#ifndef FILTERS_H
#define FILTERS_H

#include <string>
#include <vector>

#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

#include "Config.h"

// ------------------------------------------------------------------------------------------------------------------------------------
// Validation Runs
// ------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Checks the validated runs from the not ones.
 * @param val_map Validation map.
 * @param run_name Run column name.
 * @param block_name Luminosity block column name.
 * @return Returns the validated node-dataset.
*/
ROOT::RDF::RNode ApplyValidationFilter(ROOT::RDF::RNode node, const validation_type& val_map, const std::string& run_name, const std::string& block_name);

// ------------------------------------------------------------------------------------------------------------------------------------
// Kinematical cuts
// ------------------------------------------------------------------------------------------------------------------------------------

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

// ------------------------------------------------------------------------------------------------------------------------------------
// TagAndProbe
// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Transverse momentum, psudorapidity and invariant mass vectors results -> from CalculeteTagAndProbe().
/*
struct ResultsTagAndProbe {
    std::vector<float> pt_pass;// Transverse momentum of probes that pass the selection.
    std::vector<float> pt_all;// Transverse momentum of all probes.
    std::vector<float> eta_pass;// Pseudorapidity of probes that pass the selection.
    std::vector<float> eta_all;// Pseudorapidity of all probes. 
    std::vector<float> mll_pass;// Invariant mass of TagAndProbe pairs that pass the selection.
    std::vector<float> mll_all;// Invariant mass of all TagAndProbe pairs 
};
*/
struct ResultsTagAndProbe {
    ROOT::RVec<float> pt_pass;
    ROOT::RVec<float> pt_fail;

    ROOT::RVec<float> eta_pass;
    ROOT::RVec<float> eta_fail;
    
    ROOT::RVec<float> mll_pass;
    ROOT::RVec<float> mll_fail;
    
    ROOT::RVec<float> tag_pt_pass;
    ROOT::RVec<float> tag_pt_fail;
    
    ROOT::RVec<float> tag_eta_pass;
    ROOT::RVec<float> tag_eta_fail;
};

/// @brief Muon kinematics variables for selection.
struct MuonKinematics_TP {
    const ROOT::RVec<float>& pt;// Transverse momentum
    const ROOT::RVec<float>& eta;// Pseudorapidity
    const ROOT::RVec<float>& phi;// Angular variable
    const ROOT::RVec<float>& mass;// Reconstructed mass
};

/// @brief Muon flags for TagAndProbe selections. Tag muon is a TightId muon. Probe candidate is a StandaloneId muon. Passed probe is a GlobalId muon + Isolation request.
struct MuonFlags_TP {
    const ROOT::RVec<bool>& tight;// Muon_tightId flag.
    const ROOT::RVec<bool>& stand;// Muon_isStandalone flag.
    const ROOT::RVec<bool>& global;// Muon_isGlobal flag.
    const ROOT::RVec<float>& iso;// Isolation parameter.
};

/**
 * @brief TagAndProbe function selection
 * @param kin Muon kinematic event values.
 * @param flags Muon event flags.
 * @param cfg_f Contains the kinematical flags. Passed by value -> small struct -> L1/L2 catche
 * @param cfg_c Stores the cuts values.
 * @return Struct containing vectors of passed and total probe kinematics.
*/
ResultsTagAndProbe CalculateTagAndProbe(const MuonKinematics_TP& kin, const MuonFlags_TP& flags, const flags_config cfg_f, 
const cuts_config cfg_c);

// ------------------------------------------------------------------------------------------------------------------------------------
// ResponseMatrix
// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Results of the Response Matrix calculus.
struct ResultsRespMatrix {
    std::vector<float> pt_gen_RM;// Transverse momentum of the selected generated muons.
    std::vector<float> pt_rec_RM;// Transverse momentum of the selected reconstructed muons.
    std::vector<float> eta_gen_RM;// Pseudorapidity of the selected generated muons.
    std::vector<float> eta_rec_RM;// Pseudorapidity of the selected reconstructed muons.
};

/// @brief Kinematical quantities for Response Matrix calculus -> GENerated and REConstructed muons.
struct MuonKinematics_RM {
    const ROOT::RVec<float>& pt_gen;
    const ROOT::RVec<float>& eta_gen;
    const ROOT::RVec<float>& pt_rec;
    const ROOT::RVec<float>& eta_rec;
};

/// @brief Muon flags for Response Matrix calculus -> REConstructed/GENerated pairing.
struct MuonFlags_RM {
    const ROOT::RVec<int>& gen_flav_rec;// 
    const ROOT::RVec<int>& pair_idx_rec;
    const ROOT::RVec<int>& status_gen;
    const ROOT::RVec<int>& pdg_id_gen;
};

/**
 * @brief TagAndProbe function selection
 * @param kin Muon kinematic event values.
 * @param flags Muon event flags.
 * @param cfg_f
 * @param cfg_c
 * @return Struct containing vectors of passed and total probe kinematics.
*/
ResultsRespMatrix CalculateRespMatrix(const MuonKinematics_RM& kin, const MuonFlags_RM& flags, const flags_config cfg_f, 
const cuts_config cfg_c);


ROOT::RDF::RNode ApplyPt_DataDivision(ROOT::RDF::RNode node, const std::vector<float>& pt_bins, const std::string& column_name);

ROOT::RDF::RNode ApplyEta_DataDivision(ROOT::RDF::RNode node, const std::vector<float>& eta_bins, const std::string& column_name);

#endif