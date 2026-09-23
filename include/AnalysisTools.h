#ifndef ANALYSISTOOLS_H
#define ANALYSISTOOLS_H

#include "Config.h"
#include "Filters.h"


#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include "TTree.h"

#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>

#include <RooDataSet.h>

#include <RooExponential.h>

#include <string>
#include <vector>

// ------------------------------------------------------------------------------------------------------------------------------------

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

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief 
/// @param event_histo 
/// @param tag 
/// @return 
std::vector<std::unique_ptr<TH1D>> PrepareEventFit(EventHisto& event_histo, const std::string& tag);

std::vector<std::unique_ptr<RooDataSet>> PrepareEventFitModel(TTree* tree, const config_struct& cfg, const std::string& tag);

// ------------------------------------------------------------------------------------------------------------------------------------

void SaveEventFitCanvas(RooRealVar& mll, RooAbsPdf& model, RooAbsData& data, RooAbsPdf& bkg_pdf, const EventFitResult& res, TH1D* h_mll, 
    TDirectory* o_dir, const int bin_idx, const std::string& tag);

// ------------------------------------------------------------------------------------------------------------------------------------

// ------
// Fitter
// ------

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
std::vector<EventFitResult> EventFitWrapper(std::vector<std::unique_ptr<TH1D>>& container, std::vector<std::unique_ptr<RooDataSet>>& container_model,
const std::string& tag, const bool save_plots, TDirectory* o_dir = nullptr);


#endif