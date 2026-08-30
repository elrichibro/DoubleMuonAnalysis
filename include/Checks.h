#ifndef CHECKS_H
#define CHECKS_H

#include <string>
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

#include <map>
#include <cstdint>
#include <iostream>

/**
 * @brief Calculates the invariant mass between the 2 muons in the MonteCarlo data.
 * @param node Input node.
 * @param tag FSR status subfix.
 * @param FSR FSR status. 1 Before FSR. 2 After FSR.
 * @return Returns the node containing the new kinematical variables.
 */
ROOT::RDF::RNode CalculateInvMass(ROOT::RDF::RNode node, const std::string& tag, int FSR);




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

ResultsRespMatrix CalculateRespMatrix(const MuonKinematics_RM& kin, const MuonFlags_RM& flags);


#endif