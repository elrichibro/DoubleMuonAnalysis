#ifndef UNFOLD_H
#define UNFOLD_H

#include "Filters.h"
#include "Config.h"

#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

#include "TUnfold.h"
#include "TUnfoldDensity.h"

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

struct RespMatrixHisto {
    ROOT::RDF::RResultPtr<TH2D> histo_pt;
    ROOT::RDF::RResultPtr<TH2D> histo_y;
    ROOT::RDF::RResultPtr<TH2D> histo_phis;

    ROOT::RDF::RResultPtr<TH1D> h1_pt_test;
    ROOT::RDF::RResultPtr<TH1D> h1_y_test;
    ROOT::RDF::RResultPtr<TH1D> h1_phis_test;
};

struct UnfoldDensities {
    std::unique_ptr<TUnfoldDensity> pt_unf;
    std::unique_ptr<TUnfoldDensity> y_unf;
    std::unique_ptr<TUnfoldDensity> phis_unf;
};

struct UnfoldResult {
    std::unique_ptr<TUnfoldDensity> unf_density;
    std::unique_ptr<TH1> h1_out_unf;
    std::unique_ptr<TH2> h2_out_cov;
    std::unique_ptr<TGraph> LCurveScan;
    std::unique_ptr<TSpline> logTauX;
    std::unique_ptr<TSpline> logTauY;
    double tau = 0.0;
    int idx_best = -1;
    double chi2A = 0.0;
    double chi2L = 0.0;
    int ndf = 0;
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
ResultsRespMatrix CalculateRespMatrix(const MuonKinematics_REC& kin_rec, const MuonKinematics_GEN& kin_gen, const MuonFlags_RM& flags, const flags_config cfg_f, 
const cuts_config cfg_c);

ROOT::RDF::RNode CalculateRespMatrixWrapper(ROOT::RDF::RNode node, const flags_config& flags_RM, const cuts_config& cuts_RM);

RespMatrixHisto BuildRespMatrixHisto(ROOT::RDF::RNode node, const config_struct& cfg);

UnfoldDensities CreateUnfoldDensity(RespMatrixHisto& histo);

UnfoldResult ApplyUnfold(UnfoldDensities& densities, EventHisto& event_histo, const config_struct& cfg, RespMatrixHisto& resp_histo);

#endif
