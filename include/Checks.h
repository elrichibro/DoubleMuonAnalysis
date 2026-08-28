#ifndef CHECKS_H
#define CHECKS_H

#include <map>
#include <cstdint>
#include <iostream>
#include <string>

#include <ROOT/RDataFrame.hxx>

ROOT::RDF::RNode InvMass(ROOT::RDF::RNode node, const std::string& tag, int FSR);

ROOT::RDF::RNode node_recMC(ROOT::RDF::RNode node);

/**
 * @brief Selects the Z0 events in the MonteCarlo sample.
 * @param pdgId PDG identification index. (23 Z0)
 * @param flags Status flag stored bitwise.
 * @return Returns the mask fot Z0 true events.
*/
ROOT::RVec<bool> is_MC_Z0(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags);


int get_MC_Z0_idx(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags);

/**
 * @brief Selects muons based on flags requirements.
 * @param pdgId PDG identification index. (13 Muon)
 * @param flags Status flag stored bitwise.
 * @param mother_id Id of the mother particle.
 * @return Returns a mask for true muons.
*/
ROOT::RVec<bool> is_MC_Muon_bFSR(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother_id);
ROOT::RVec<bool> is_MC_Muon_aFSR(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother_id);


/**
 * @brief Accepts antimuons on physics requirements.
 * @param pdgId PDG identification index. (-13 AntiMuon)
 * @param flags Status flag stored bitwise.
 * @param mother_id Id of the mother particle.
 * @return Returns a mask for an antimuon selection.
*/
ROOT::RVec<bool> is_MC_AntiMuon_bFSR(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother_id);
ROOT::RVec<bool> is_MC_AntiMuon_aFSR(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother_id);


/**
 * @brief Generates a mask based on true event selection.
 * @param pdgId PDG identification index.
 * @param flags Status flag stored bitwise.
 * @return Returns true for a Z0 decaying into mu+mu-, else false.
*/
bool is_MC_Event(const ROOT::RVec<Int_t>& pdgId, const ROOT::RVec<Int_t>& flags, const ROOT::RVec<Int_t>& mother_id, const int FSR);


#endif