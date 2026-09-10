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

/// @brief 
/// @param node 
/// @param pt_col 
/// @param eta_col 
/// @param mll_col 
/// @param columns_name 
/// @param cfg 
/// @param dataset 
/// @return 
ROOT::RDF::RNode ApplyKinMuonFilter(ROOT::RDF::RNode node, const std::string& pt_col, const std::string& eta_col, const std::string& mll_col, 
std::vector<std::string>& columns_name, const config_struct& cfg, const int dataset, const int succes);

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
const config_struct& cfg);

// -------
// STRUCTS
// -------

/// @brief Transverse momentum, psudorapidity and invariant mass vectors results -> from CalculeteTagAndProbe().
struct ResultsTagAndProbe {
    ROOT::RVec<bool> mask_pass;// Bool mask of probes that pass the selection.
    ROOT::RVec<float> pt;// Transverse momentum of all probes.
    ROOT::RVec<float> eta;// Pseudorapidity of all probes. 
    ROOT::RVec<float> mll;// Invariant mass of all Tag/Probe pairs. 
    //ROOT::RVec<float> tag_pt;
    //ROOT::RVec<float> tag_eta;
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

/// @brief Results of the Response Matrix calculus.
struct ResultsRespMatrix {
    ROOT::RVec<float> pt_gen_RM;// Transverse momentum of the selected generated muons.
    ROOT::RVec<float> pt_rec_RM;// Transverse momentum of the selected reconstructed muons.
    ROOT::RVec<float> eta_gen_RM;// Pseudorapidity of the selected generated muons.
    ROOT::RVec<float> eta_rec_RM;// Pseudorapidity of the selected reconstructed muons.
};

/// @brief Kinematical quantities for Response Matrix calculus -> GENerated and REConstructed muons.
struct MuonKinematics_RM {
    const ROOT::RVec<float>& pt_gen;// Transeverse momentum of generated level muons.
    const ROOT::RVec<float>& eta_gen;// Pseudorapidity of generated level muons.
    const ROOT::RVec<float>& pt_rec;// Transeverse momentum of reconstructed level muons.
    const ROOT::RVec<float>& eta_rec;// Pseudorapidity of reconstructed level muons.
};

/// @brief Muon flags for Response Matrix calculus -> REConstructed/GENerated pairing.
struct MuonFlags_RM {
    const ROOT::RVec<int>& gen_flav_rec;// Flavour of genParticle for MC matching to status==1 muons.
    const ROOT::RVec<int>& pair_idx_rec;// Index into genParticle list for MC matching to status==1 muons
    const ROOT::RVec<int>& status_gen;// Generated particle status -> stable = 1.
    const ROOT::RVec<int>& pdg_id_gen;// PDG id of the particle.
};

// -------
// METHODS
// -------

// ------------------------------------------------------------------------------------------------------------------------------------
// TagAndProbe
// ------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief TagAndProbe function selection for MonteCarlo sample (by adding a DeltaR selection)
 * @param kin Muon kinematic event values.
 * @param flags Muon event flags.
 * @param cfg_f Contains the kinematical flags. Passed by value -> small struct -> L1/L2 catche
 * @param cfg_c Stores the cuts values.
 * @param DeltaR_flags Flags for Generated/Reconstructed muons pairing.
 * @param gen_eta Pseudorapidity of Generated muon used for DeltaR calculus.
 * @param gen_phi Phi angle of Generated muon used for DeltaR calculus.
 * @return Struct containing probe muons kinematical variables.
*/
ResultsTagAndProbe CalculateTagAndProbe_MC(const MuonKinematics_TP& kin, const MuonFlags_TP& flags, const flags_config cfg_f, 
const cuts_config cfg_c, const MuonFlags_RM& DeltaR_flags, const ROOT::RVec<float> gen_eta, const ROOT::RVec<float> gen_phi);

/**
 * @brief TagAndProbe function selection for data sample.
 * @param kin Muon kinematic event values.
 * @param flags Muon event flags.
 * @param cfg_f Contains the kinematical flags. Passed by value -> small struct -> L1/L2 catche
 * @param cfg_c Stores the cuts values.
 * @return Struct containing probe muons kinematical variables.
*/
ResultsTagAndProbe CalculateTagAndProbe_DATA(const MuonKinematics_TP& kin, const MuonFlags_TP& flags, const flags_config cfg_f, 
const cuts_config cfg_c);

// ------------------------------------------------------------------------------------------------------------------------------------
// ResponseMatrix
// ------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Calculates the Response Matrix between Generated and Reconstructed muons in MonteCarlo sample.
 * @param kin Muon kinematical variables.
 * @param flags Muon event flags.
 * @param cfg_f Contains the kinematical flags. Passed by value -> small struct -> L1/L2 catche
 * @param cfg_c Stores the cuts values.
 * @return Struct containing the transverse momentum and pseudorapidity of muons that pass the selection. 
*/
ResultsRespMatrix CalculateRespMatrix(const MuonKinematics_RM& kin, const MuonFlags_RM& flags, const flags_config cfg_f, 
const cuts_config cfg_c);

#endif