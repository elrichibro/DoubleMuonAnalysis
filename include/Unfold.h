#ifndef UNFOLD_H
#define UNFOLD_H

#include "Filters.h"
#include "Config.h"
#include "AnalysisTools.h"

#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

#include "TUnfold.h"
#include "TUnfoldDensity.h"

// ------------------------------------------------------------------------------------------------------------------------------------

// -------
// Structs
// -------

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

/// @brief Response Matrix histograms
struct RespMatrixHisto {
    ROOT::RDF::RResultPtr<TH2D> h2_pt;// Response Matrix P_t_Z0
    ROOT::RDF::RResultPtr<TH2D> h2_y;// Response Matrix Y_Z0
    ROOT::RDF::RResultPtr<TH2D> h2_phis;// Response Matrix Phi*_Z0

    ROOT::RDF::RResultPtr<TH1D> h1_pt_fake;// Distribution of Reco P_t_Z0 in Fake flag case
    ROOT::RDF::RResultPtr<TH1D> h1_y_fake;// Distribution of Reco Y_Z0 in Fake flag case
    ROOT::RDF::RResultPtr<TH1D> h1_phis_fake;// Distribution of Reco Phi*_Z0 in Fake flag case

    ROOT::RDF::RResultPtr<TH1D> h1_pt_test;// Distribution of Reco P_t_Z0 in matched case
    ROOT::RDF::RResultPtr<TH1D> h1_y_test;// Distribution of Reco Y_Z0 in matched case
    ROOT::RDF::RResultPtr<TH1D> h1_phis_test;// Distribution of Reco Phi*_Z0 in matched case
};

/// @brief Control histograms structs
struct ControlHisto {
    std::unique_ptr<TH1D> h1_Eff_pt;
    std::unique_ptr<TH1D> h1_Pur_pt;
    std::unique_ptr<TH1D> h1_Stab_pt;
    
    std::unique_ptr<TH1D> h1_Eff_y;
    std::unique_ptr<TH1D> h1_Pur_y;
    std::unique_ptr<TH1D> h1_Stab_y;

    std::unique_ptr<TH1D> h1_Eff_phis;
    std::unique_ptr<TH1D> h1_Pur_phis;
    std::unique_ptr<TH1D> h1_Stab_phis;
};

// Unfold Density struct
struct UnfoldDensities {
    std::unique_ptr<TUnfoldDensity> pt_unf;
    std::unique_ptr<TUnfoldDensity> y_unf;
    std::unique_ptr<TUnfoldDensity> phis_unf;
};

/// @brief Unfold Result struct
struct UnfoldResult {
    std::unique_ptr<TUnfoldDensity> unf_density;
    
    std::unique_ptr<TH1> h1_out_unf;
    std::unique_ptr<TH2> h2_out_cov;
    std::unique_ptr<TH2> h2_out_corr;

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

// ------------------------------------------------------------------------------------------------------------------------------------

// ---------------
// Response Matrix
// ---------------

/// @brief Calculates the Response Matrix between Generated and Reconstructed muons in MonteCarlo sample.
/// @param kin Muon kinematical variables.
/// @param flags Muon event flags.
/// @param cfg_f Contains the kinematical flags. Passed by value -> small struct -> L1/L2 catche
/// @param cfg_c Stores the cuts values.
/// @return Struct containing the transverse momentum and pseudorapidity of muons that pass the selection.
ResultsRespMatrix CalculateRespMatrix(const MuonKinematics_REC& kin_rec, const MuonKinematics_GEN& kin_gen, const MuonFlags_RM& flags, const flags_config cfg_f, 
const cuts_config cfg_c);

/// @brief Wrapper of CalculateRespMatrix function -> variables defined within the wrapper.
/// @param node RDF node input.
/// @param flags_RM Response Matrix flags.
/// @param cuts_RM Response Matrix cuts values.
/// @return RDF node with the filtered quantities for Response MAtrix construction.
ROOT::RDF::RNode CalculateRespMatrixWrapper(ROOT::RDF::RNode node, const flags_config& flags_RM, const cuts_config& cuts_RM);

/// @brief Books the Matched, Fakes and ResponseMatrix histograms.
/// @param node RDF input node.
/// @param cfg General configure struct.
/// @return Struct containing the result histograms from CalculateRespMatrix function.
RespMatrixHisto BuildRespMatrixHisto(ROOT::RDF::RNode node, const config_struct& cfg);

// ------------------------------------------------------------------------------------------------------------------------------------

// --------------------------------------------
// Control Plots: Efficiency, Purity, Stability
// --------------------------------------------

/// @brief Creates control histograms and stores them into the return struct.
/// @param node RDF input node.
/// @param cfg General configure struct.
/// @return Control histogram struct.
ControlHisto BuildControlHisto(ROOT::RDF::RNode node, const config_struct& cfg);

// ------------------------------------------------------------------------------------------------------------------------------------

// ----------------
// Unfold Procedure
// ----------------

/// @brief Creates the Unfold Density for LScan, ScanTau or ScanSURE application.
/// @param histo Response matrix result container.
/// @param tag Quantity identifier.
/// @return Struct containing the Unfold
UnfoldDensities CreateUnfoldDensity(RespMatrixHisto& histo, const std::string& tag);

/// @brief Apply the Unfold/Scan procedure to a Unfold density.
/// @param density Unfold density input.
/// @param event_histo Signal input -> DATASET -> data.
/// @param resp_histo Signal input for testing purpose.
/// @param fake_histo Histogram for BKG substraction.
/// @param cfg General configure struct.
/// @param tag Unfold quantity identifier.
/// @return Unfold result struct: Unfolded signal histogram.
UnfoldResult ApplyUnfold(std::unique_ptr<TUnfoldDensity> density, TH1D* event_histo, TH1D* resp_histo, TH1D* fake_histo, 
const config_struct& cfg, const std::string& tag);

/// @brief 
/// @param canvas 
/// @param results 
/// @param histo_resp 
/// @param eff_histos 
/// @param tag 
/// @return 
int VisualizeUnfoldResults(std::vector<std::unique_ptr<TCanvas>>& canvas, UnfoldResult& results, RespMatrixHisto& resp_histo,  const std::string& tag);   

/// @brief 
/// @param canvas 
/// @param resp_histo 
/// @param control_histo 
/// @param tag 
/// @return 
int VisualizeControlPlots(std::vector<std::unique_ptr<TCanvas>>& canvas, RespMatrixHisto& resp_histo, ControlHisto& control_histo, const std::string& tag);

/// @brief 
/// @param results 
/// @param cfg 
/// @param tag 
/// @return 
std::unique_ptr<TH1D> BuildFitResultHistogram(const std::vector<EventFitResult>& results, const config_struct& cfg, const std::string& tag);

/// @brief 
/// @param cfg 
/// @return 
std::unique_ptr<TH1D> EventFit_SignalHisto_Wrapper(const config_struct& cfg, EventHisto& event_histo);

#endif
