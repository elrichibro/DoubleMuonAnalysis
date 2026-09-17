#ifndef UNFOLD_H
#define UNFOLD_H

#include "Filters.h"
#include "Config.h"

#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

/// @brief Results of the Response Matrix calculus.
struct ResultsRespMatrix {
    bool match = false;
    bool miss = false;
    bool fake = false;

    float mll_rec = -1;
    float mll_gen = -1;
    
    float pt_gen = -1;
    float pt_rec = -1;
    
    float y_gen = -99;
    float y_rec = -99;
    
    float phis_gen = -1;
    float phis_rec = -1;
};

/// @brief Kinematical quantities for Response Matrix calculus -> GENerated and REConstructed muons.
struct MuonKinematics_RM {
    const ROOT::RVec<float>& pt;
    const ROOT::RVec<float>& eta;
    const ROOT::RVec<float>& phi;
    const ROOT::RVec<float>& mass;
    const ROOT::RVec<int>& charge;
};

/// @brief Muon flags for Response Matrix calculus -> REConstructed/GENerated pairing.
struct MuonFlags_RM {
    const ROOT::RVec<bool>& reco_tight;// Muon_tightId flag.
    const ROOT::RVec<int>& reco_idx_gen;// Index into genParticle list for MC matching to status==1 muons
    const ROOT::RVec<int>& gen_pdg_idx;// PDG id of the particle.
    const ROOT::RVec<int>& gen_status_flg;
};

/// @brief Calculates the Response Matrix between Generated and Reconstructed muons in MonteCarlo sample.
/// @param kin Muon kinematical variables.
/// @param flags Muon event flags.
/// @param cfg_f Contains the kinematical flags. Passed by value -> small struct -> L1/L2 catche
/// @param cfg_c Stores the cuts values.
/// @return Struct containing the transverse momentum and pseudorapidity of muons that pass the selection.
ResultsRespMatrix CalculateRespMatrix(const MuonKinematics_RM& kin_rec, const MuonKinematics_RM& kin_gen, const MuonFlags_RM& flags, const flags_config cfg_f, 
const cuts_config cfg_c);

ROOT::RDF::RNode ApplyUnfold(ROOT::RDF::RNode node);


#endif
