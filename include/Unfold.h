#ifndef UNFOLD_H
#define UNFOLD_H

#include <ROOT/RDataFrame.hxx>

#include "Config.h"

#include "TUnfoldDensity.h"

// ------------------------------------------------------------------------------------------------------------------------------------

// ------------------------------------
// Selection Response Matrix Histograms
// ------------------------------------

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

/// @brief Books the Matched, Fakes and ResponseMatrix histograms.
/// @param node RDF input node.
/// @param cfg General configure struct.
/// @return Struct containing the result histograms from CalculateRespMatrix function.
RespMatrixHisto BuildRespMatrixHisto(ROOT::RDF::RNode node, const config_struct& cfg);

// ------------------------------------------------------------------------------------------------------------------------------------

// ---------------------------------------------------------
// Control Plots for Unfoding: Efficiency, Purity, Stability
// ---------------------------------------------------------

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

/// @brief Creates control histograms and stores them into the return struct.
/// @param node RDF input node.
/// @param cfg General configure struct.
/// @return Control histogram struct.
ControlHisto BuildControlHisto(ROOT::RDF::RNode node, const config_struct& cfg);

// ------------------------------------------------------------------------------------------------------------------------------------

// --------------------
// Pre-Unfold procedure
// --------------------

// Unfold Density struct
struct UnfoldDensities {
    std::unique_ptr<TUnfoldDensity> pt_unf;
    std::unique_ptr<TUnfoldDensity> y_unf;
    std::unique_ptr<TUnfoldDensity> phis_unf;
};

/// @brief Creates the Unfold Density for LScan, ScanTau or ScanSURE application.
/// @param histo Response matrix result container.
/// @param tag Quantity identifier.
/// @return Struct containing the Unfold
UnfoldDensities CreateUnfoldDensity(RespMatrixHisto& histo, const std::string& tag);

// ------------------------------------------------------------------------------------------------------------------------------------

// ----------------
// Unfold Procedure
// ----------------

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

// ------------------------------------------------------------------------------------------------------------------------------------

// ---------------------
// Visualization options
// ---------------------

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


#endif
