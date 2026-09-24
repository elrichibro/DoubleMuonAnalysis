#ifndef CONTROLCHECK_H
#define CONTROLCHECK_H

#include <ROOT/RDataFrame.hxx>

#include "Config.h"

#include <string>
#include <vector>

// ------------------------------------------------------------------------------------------------------------------------------------

// -----------------
// Validation filter
// -----------------

/// @brief Apply the validation check to all the events in the dataset.
/// @param node Input RDF node.
/// @param val_map Struct containing the validated run and the relative luminosity blocks.
/// @param run_name Name of the "run" column in the dataset.
/// @param block_name Name of the "luminosity block" column in the dataset.
/// @return A RDF node with all validated events for post analysis.
ROOT::RDF::RNode ApplyValidationFilter(ROOT::RDF::RNode node, const validation_type& val_map, const std::string& run_name, const std::string& block_name);

// ------------------------------------------------------------------------------------------------------------------------------------

// ----------
// Acceptance
// ----------

/// @brief Calculates the acceptance from MonteCarlo sample. Cuts values are fixed by definition -> can be overwritten manually.
/// @param node Input RDF node.
/// @param tag Tag for Column name.
/// @param FSR FSR Muon states.
/// @return Returns the value of detector acceptance in Z0->mu+mu- event. First element is central value, second value is the stat. sigma.
std::vector<float> CalculateAcceptance(ROOT::RDF::RNode node, const std::string& tag, int FSR);

// ------------------------------------------------------------------------------------------------------------------------------------

// ----------
// Resolution
// ----------

struct ResolutionResults {
    std::vector<double> mean;
    std::vector<double> sigma;
    std::vector<int> events;
};

/// @brief Calculates the resolution of the detector from MC sample
/// @param node RDF input node 
/// @param cfg General configure input
/// @return A struct containing aritmetic mean and standard deviation for each bin.
ResolutionResults CalculateResolution(ROOT::RDF::RNode node, const config_struct& cfg);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Calculates the invariant mass between the 2 muons in the MonteCarlo data.
/// @param node Input node.
/// @param tag FSR status subfix.
/// @param FSR FSR status. 1 Before FSR. 2 After FSR.
/// @return Returns the node containing the new kinematical variables.
ROOT::RDF::RNode CalculateInvMass(ROOT::RDF::RNode node, const std::string& tag, int FSR);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Selects the muon pair and calculates the invariant mass.
/// @param node Input RDF node.
/// @return Updated node with new defined quantities.
ROOT::RDF::RNode node_recMC(ROOT::RDF::RNode node);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Selects the Z0 events in the MonteCarlo sample.
/// @param pdgId PDG identification index. (23 Z0)
/// @param flags Status flag stored bitwise.
/// @return Returns the mask fot Z0 true events.
ROOT::RVec<bool> is_MC_Z0(const ROOT::RVec<int>& pdgId, const ROOT::RVec<int>& flags);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Gets the index in GenPart of the final state of the Z0.
/// @param pdgId Pdf index for particle identification.
/// @param flags Flags of the Generated Particles.
/// @return Index of Z0 in a Z0->mu+mu- event.
int get_MC_Z0_idx(const ROOT::RVec<int>& pdgId, const ROOT::RVec<int>& flags);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Selects muons based on flags requirements.
/// @param pdgId PDG identification index. (13 Muon)
/// @param flags Status flag stored bitwise.
/// @param mother_id Id of the mother particle.
/// @return Returns a mask for true muons.
ROOT::RVec<bool> is_MC_Muon_bFSR(const ROOT::RVec<int>& pdgId, const ROOT::RVec<int>& flags, const ROOT::RVec<int>& mother_id);
ROOT::RVec<bool> is_MC_Muon_aFSR(const ROOT::RVec<int>& pdgId, const ROOT::RVec<int>& flags, const ROOT::RVec<int>&);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Accepts antimuons on physics requirements.
/// @param pdgId PDG identification index. (-13 AntiMuon)
/// @param flags Status flag stored bitwise.
/// @param mother_id Id of the mother particle.
/// @return Returns a mask for an antimuon selection.
ROOT::RVec<bool> is_MC_AntiMuon_bFSR(const ROOT::RVec<int>& pdgId, const ROOT::RVec<int>& flags, const ROOT::RVec<int>& mother_id);
ROOT::RVec<bool> is_MC_AntiMuon_aFSR(const ROOT::RVec<int>& pdgId, const ROOT::RVec<int>& flags, const ROOT::RVec<int>&);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Generates a mask based on true event selection.
/// @param pdgId PDG identification index.
/// @param flags Status flag stored bitwise.
/// @param mother_id 
/// @param FSR FSR status, 1 = Before FSR, 2 = After FSR.
/// @return Returns true for a Z0 decaying into mu+mu-, else false.
bool is_MC_Event(const ROOT::RVec<int>& pdgId, const ROOT::RVec<int>& flags, const ROOT::RVec<int>& mother_id, const int FSR);

#endif