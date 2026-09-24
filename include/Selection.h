#ifndef SELECTION_H
#define SELECTION_H

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>


#include "Config.h"

// ------------------------------------------------------------------------------------------------------------------------------------

// ---------------
// Response Matrix
// ---------------

/// @brief Results from Response Matrix calculus.
struct ResultsRespMatrix {
    bool match = false;// Matched event flag [2 gen : 2 reco] muons
    bool miss = false;// Missed event flag [2 gen : 1-0 reco] muons
    bool fake = false;// Faked event flag [1-0 gen : 2 reco (or 2 gen : 2 reco without match)] muons

    float mll_rec = -1;// Invariant mass from reconstructed pair
    float mll_gen = -1;// Invariant mass from generated pair
    
    float pt_gen = -1;// Transverse momentum of Z0 from reconstructed pair
    float pt_rec = -1;// Transverse momentum of Z0 from generated pair
    
    float y_gen = -99;// Rapidity of Z0 from reconstructed pair
    float y_rec = -99;// Rapidity of Z0 from generated pair
    
    float phis_gen = -1;// Phi* of Z0 from reconstructed pair
    float phis_rec = -1;// Phi* of Z0 from generated pair
};

/// @brief Kinematical quantities for Response Matrix calculus -> GENerated and REConstructed muons.
struct MuonKinematics_REC {
    const ROOT::RVec<float>& pt;
    const ROOT::RVec<float>& eta;
    const ROOT::RVec<float>& phi;
    const ROOT::RVec<float>& mass;
    const ROOT::RVec<int>& charge;
};

/// @brief Kinematical quantities for Response Matrix calculus -> GENerated and REConstructed muons.
struct MuonKinematics_GEN {
    const ROOT::RVec<float>& pt;
    const ROOT::RVec<float>& eta;
    const ROOT::RVec<float>& phi;
    const ROOT::RVec<float>& mass;
};

/// @brief Muon flags for Response Matrix calculus -> REConstructed/GENerated pairing.
struct MuonFlags_RM {
    const ROOT::RVec<bool>& reco_tight;// Muon_tightId flag.
    const ROOT::RVec<int>& reco_idx_gen;// Index into genParticle list for MC matching to status==1 muons
    const ROOT::RVec<int>& gen_pdg_idx;// PDG id of the particle.
    const ROOT::RVec<int>& gen_status;
    const ROOT::RVec<int>& gen_status_flg;
};

/// @brief Calculates the Response Matrix between Generated and Reconstructed muons in MonteCarlo sample.
/// @param kin Muon kinematical variables.
/// @param flags Muon event flags.
/// @param cfg_f Contains the kinematical flags. Passed by value -> small struct -> L1/L2 catche
/// @param cfg_c Stores the cuts values.
/// @return Struct containing the transverse momentum and pseudorapidity of muons that pass the selection.
ResultsRespMatrix CalculateRespMatrix(const MuonKinematics_REC& kin_rec, const MuonKinematics_GEN& kin_gen, const MuonFlags_RM& flags, 
    const flags_config cfg_f, const cuts_config cfg_c);

/// @brief Wrapper of CalculateRespMatrix function -> variables defined within the wrapper.
/// @param node RDF node input.
/// @param flags_RM Response Matrix flags.
/// @param cuts_RM Response Matrix cuts values.
/// @return RDF node with the filtered quantities for Response MAtrix construction.
ROOT::RDF::RNode CalculateRespMatrixWrapper(ROOT::RDF::RNode node, const flags_config& flags_RM, const cuts_config& cuts_RM);

// ------------------------------------------------------------------------------------------------------------------------------------

// ---------------
// Event selection
// ---------------

/// @brief Filters the event selecting Z0/mu+mu- events falling withing the fiducial region.
/// @param node RDF input node: DATA/MC -> Event selection procedure
/// @param cfg General configure struct
/// @return RNode with filtered events and new variables defined.
ROOT::RDF::RNode EventSelection(ROOT::RDF::RNode node, const config_struct& cfg);


#endif