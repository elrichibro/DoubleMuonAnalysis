#ifndef FILTERS_H
#define FILTERS_H

#include <string>
#include <vector>

#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

#include "Config.h"

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Apply the validation check to all the events in the dataset.
/// @param node Input RDF node.
/// @param val_map Struct containing the validated run and the relative luminosity blocks.
/// @param run_name Name of the "run" column in the dataset.
/// @param block_name Name of the "luminosity block" column in the dataset.
/// @return A RDF node with all validated events for post analysis.
ROOT::RDF::RNode ApplyValidationFilter(ROOT::RDF::RNode node, const validation_type& val_map, const std::string& run_name, const std::string& block_name);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Applies the kinematical bin division to and unbinned dataset.
/// @param node Input RDF node
/// @param pt_col Name of the transverse momentum column.
/// @param eta_col Name of the pseudorapidity column. 
/// @param mll_col Name of the invariant mass column.
/// @param columns_name Column names for Snapshot operation.
/// @param cfg General configure struct.
/// @param dataset MC or DATA  
/// @param sample PASS or FAIL
/// @return Returns the RDF node with the new defined quantities.
ROOT::RDF::RNode ApplyKinematicalBinDivision(ROOT::RDF::RNode node, const config_struct& cfg, const std::string& pt_col, const std::string& eta_col,
    const std::string& mll_col, std::vector<std::string>& columns_name, const int dataset, const int sample);

// ------------------------------------------------------------------------------------------------------------------------------------

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


// ------------------------------------------------------------------------------------------------------------------------------------

// ---------
// Functions
// ---------

/// @brief TagAndProbe function selection for MonteCarlo sample (by adding a DeltaR selection)
/// @param kin Muon kinematic event values.
/// @param flags Muon event flags.
/// @param cfg_f Contains the kinematical flags. Passed by value -> small struct -> L1/L2 catche
/// @param cfg_c Stores the cuts values.
/// @param DeltaR_flags Flags for Generated/Reconstructed muons pairing.
/// @param gen_eta Pseudorapidity of Generated muon used for DeltaR calculus.
/// @param gen_phi Phi angle of Generated muon used for DeltaR calculus.
/// @return Struct containing probe muons kinematical variables.
ResultsTagAndProbe CalculateTagAndProbe_MC(const MuonKinematics_TP& kin, const MuonFlags_TP& flags, const flags_config cfg_f, 
const cuts_config cfg_c, const MuonFlags_RM& DeltaR_flags, const ROOT::RVec<float> gen_eta, const ROOT::RVec<float> gen_phi);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief TagAndProbe function selection for data sample.
/// @param kin Muon kinematic event values.
/// @param flags Muon event flags.
/// @param cfg_f Contains the kinematical flags. Passed by value -> small struct -> L1/L2 catche
/// @param cfg_c Stores the cuts values.
/// @return Struct containing probe muons kinematical variables.
ResultsTagAndProbe CalculateTagAndProbe_DATA(const MuonKinematics_TP& kin, const MuonFlags_TP& flags, const flags_config cfg_f, 
const cuts_config cfg_c);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Calculates the Response Matrix between Generated and Reconstructed muons in MonteCarlo sample.
/// @param kin Muon kinematical variables.
/// @param flags Muon event flags.
/// @param cfg_f Contains the kinematical flags. Passed by value -> small struct -> L1/L2 catche
/// @param cfg_c Stores the cuts values.
/// @return Struct containing the transverse momentum and pseudorapidity of muons that pass the selection.
ResultsRespMatrix CalculateRespMatrix(const MuonKinematics_RM& kin, const MuonFlags_RM& flags, const flags_config cfg_f, 
const cuts_config cfg_c);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Calculates the acceptance from MonteCarlo sample.
/// @param node Input RDF node.
/// @param tag Tag for Column name.
/// @param FSR FSR Muon states.
/// @return Returns the value of detector acceptance in Z0->mu+mu- event. First element is central value, second value is the stat. sigma.
std::vector<float> CalculateAcceptance(ROOT::RDF::RNode node, const std::string& tag, int FSR);

#endif