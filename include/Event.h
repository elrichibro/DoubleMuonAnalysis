#ifndef EVENT_H
#define EVENT_H 

#include <ROOT/RDataFrame.hxx>
#include <TTree.h>
#include <TDirectory.h>

#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooAbsData.h>

#include <vector>
#include <string>

#include <TH1D.h>
#include <TH2D.h>

#include "Config.h"

// ------------------------------------------------------------------------------------------------------------------------------------

// -------------------------
// Selection Event histogram
// -------------------------

/// @brief RDF Histograms obtained from DATA sample.
struct EventSelectionHisto {
    ROOT::RDF::RResultPtr<TH1D> h1_mll;// 1D histogram: Z0 invariant mass
    ROOT::RDF::RResultPtr<TH1D> h1_pt;// 1D histogram: Z0 Transverse momentum
    ROOT::RDF::RResultPtr<TH1D> h1_y;// 1D histogram: Z0 rapidity
    ROOT::RDF::RResultPtr<TH1D> h1_phis;// 1D histogram: Z0 phi*

    ROOT::RDF::RResultPtr<TH2D> h2_mll_pt;// 2D histogram: Z0 transverse momentum vs invariant mass
    ROOT::RDF::RResultPtr<TH2D> h2_mll_y;// 2D histogram: Z0 rapidity vs invariant mass
    ROOT::RDF::RResultPtr<TH2D> h2_mll_phis;// 2D histogram: Z0 phi* vs invariant mass
};

/// @brief Creates histograms from RDF input node and save them into EventHisto struct. The selection cuts/flags can be setted by JSON values.
/// @param node RDF input node.
/// @param cfg General config struct.
/// @return The histogram struct filled with event histograms.
EventSelectionHisto BuildEventSelection_Histo(ROOT::RDF::RNode node, const config_struct& cfg);

// ------------------------------------------------------------------------------------------------------------------------------------

// -----------------------------------------
// Building Histograms for fitting procedure
// -----------------------------------------

/// @brief 
/// @param ev_sel_histo 
/// @param tag 
/// @return 
std::vector<std::unique_ptr<TH1D>> BuildEventFit_Histo(EventSelectionHisto& ev_sel_histo, const std::string& tag);

/// @brief 
/// @param tree 
/// @param cfg 
/// @param tag 
/// @return 
std::vector<std::unique_ptr<RooDataSet>> BuildEventFit_SignalModel(TTree* tree, const config_struct& cfg, const std::string& tag);

// ------------------------------------------------------------------------------------------------------------------------------------

// -------------
// Fit porcedure
// -------------

// Event fit results
struct EventFitResult {
    int bin_idx;// Index of the P_t, Y or Phi* quantity.
    
    double n_sig;// Signal yield.
    double n_sig_err;// Signal yield fit error.
    
    double n_bkg;// Backgorund yield.
    double n_bkg_err;// Background yield error.
    
    double lambda;// Exponential lambda result.
    double lambda_err;// Lambda fit error.
    
    int fit_status;// Fit status -> 0 = succes
};

/// @brief Fits binned data with a signal model hardcoded. Saves the fit and residual plot into the output Event file.
/// @param bin_idx Bin index of P_t, Y, Phis star.
/// @param h_mll Histogram -> Binned data for input.
/// @param o_dir Output directory.
/// @param tag Identifies the ortogonal quantity.
/// @return Returns a struct containing the fit results.
EventFitResult EventSingleFit(int bin_idx, TH1D* h_mll, RooDataSet* d_mll_model, TDirectory* o_dir, const std::string& tag, const bool save_plots);

/// @brief Wrapper for EventSingleFit that loops on all P_t, Y or Phi* bins.
/// @param container Container for binned input data.
/// @param o_dir Output directory.
/// @param tag Identifies the ortogonal quantity.
/// @return Returns all the Fit results.
std::vector<EventFitResult> EventSingleFit_Wrapper(std::vector<std::unique_ptr<TH1D>>& container, std::vector<std::unique_ptr<RooDataSet>>& container_model,
const std::string& tag, const bool save_plots, TDirectory* o_dir = nullptr);

// ------------------------------------------------------------------------------------------------------------------------------------

// ---------------------------
// All Event procedure wrapper
// ---------------------------

/// @brief 
/// @param cfg 
/// @return 
std::unique_ptr<TH1D> EventFit_SignalHisto_Wrapper(EventSelectionHisto& ev_sel_histo, const config_struct& cfg);

// ------------------------------------------------------------------------------------------------------------------------------------

// ----------------------------------
// Fit visualization/saving functions
// ----------------------------------

/// @brief 
/// @param results 
/// @param cfg 
/// @param tag 
/// @return 
std::unique_ptr<TH1D> BuildFitResult_Histo(const std::vector<EventFitResult>& results, const config_struct& cfg, const std::string& tag);

/// @brief 
/// @param mll 
/// @param model 
/// @param data 
/// @param bkg_pdf 
/// @param res 
/// @param h_mll 
/// @param o_dir 
/// @param bin_idx 
/// @param tag 
void SaveEventFitCanvas(RooRealVar& mll, RooAbsPdf& model, RooAbsData& data, RooAbsPdf& bkg_pdf, const EventFitResult& res, TH1D* h_mll, 
TDirectory* o_dir, const int bin_idx, const std::string& tag);


#endif