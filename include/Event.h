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

/// @brief Divides the TH2D into TH1D Invariant mass histograms for each bin of pt, y or phis (selected tag quantity). 
/// @param ev_sel_histo Histogram struct obtained from event selection.
/// @param tag Event quantity: pt, y, phis.
/// @return A vector containing the TH1D histograms of the Invariant Mass for each bin of pt, y, phis.
std::vector<std::unique_ptr<TH1D>> BuildEventFit_Histo(EventSelectionHisto& ev_sel_histo, const std::string& tag);

/// @brief Divides the unbinned data directly from Event selection ttree and fills a std vector of RooDataSet. That vector will be used for model generation with RooKeysPdf.
/// @param tree Input unbinned data ttree.
/// @param cfg General config struct.
/// @param tag Studied quantity.
/// @return A vector of unbinned type containers for model generation in fit procedure.
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
    
    int fit_status;// Fit status -> 0 = success
};

/// @param h_mll Histogram -> Binned data for input.
/// @param o_dir Output directory.
/// @param tag Identifies the ortogonal quantity.
/// @return Returns a struct containing the fit results.

/// @brief Fits binned data with a signal model obtained from MC unbinned data. Saves the fit results and the plot with residuals into the output Event file.
/// @param bin_idx Bin index of P_t, Y, Phis*.
/// @param h_mll Histogram -> Binned data for input.
/// @param d_mll_model Unbinne data for model generation.
/// @param o_dir Output directory.
/// @param tag Studied quantity.
/// @param save_plots Saving plots flag.
/// @return A struct containing the fit results. Quantity of interest is the signal yield.
EventFitResult EventSingleFit(int bin_idx, TH1D* h_mll, RooDataSet* d_mll_model, TDirectory* o_dir, const std::string& tag, const bool save_plots);

/// @brief Wrapper for EventSingleFit that loops on all P_t, Y or Phi* bins.
/// @param container Container for binned input data.
/// @param o_dir Output directory.
/// @param tag Identifies the ortogonal quantity.
/// @return Returns all the Fit results.

/// @brief Wrapper for EventSingleFit that loops on all P_t, Y or Phi* bins.
/// @param container Input histograms container.
/// @param container_model Input unbinned data for model generation container.
/// @param tag Studied quantity tag: P_t, Y, Phi*.
/// @param save_plots Save plots flag.
/// @param o_dir Output directory.
/// @return A vector of Event single fit structs -> all the results for all (P_t, Y or Phi* -> only one quantity under study) bins.
std::vector<EventFitResult> EventSingleFit_Wrapper(std::vector<std::unique_ptr<TH1D>>& container, std::vector<std::unique_ptr<RooDataSet>>& container_model,
const std::string& tag, const bool save_plots, TDirectory* o_dir = nullptr);

// ------------------------------------------------------------------------------------------------------------------------------------

// ---------------------------
// All Event procedure wrapper
// ---------------------------

/// @brief Wrapper of all steps of Event option.
/// @param ev_sel_histo Event selection histograms.
/// @param cfg General configuration struct.
/// @return A sigle TH1D of the signal with subtracted background.
std::unique_ptr<TH1D> EventFit_SignalHisto_Wrapper(EventSelectionHisto& ev_sel_histo, const config_struct& cfg);

// ------------------------------------------------------------------------------------------------------------------------------------

// ----------------------------------
// Fit visualization/saving functions
// ----------------------------------

/// @brief Builds the P_t, Y or Phi* histogram with background subtraction. Obtained from fit procedure.
/// @param results Vector containing all bins fit results.
/// @param cfg General configuration struct.
/// @param tag Studied quantity.
/// @return A TH1D P_t, Y or Phi* histogram of the signal(only) yield.
std::unique_ptr<TH1D> BuildFitResult_Histo(const std::vector<EventFitResult>& results, const config_struct& cfg, const std::string& tag);

/// @param mll RooFit observable.
/// @param model RooFit signal pdf.
/// @param data RooFit data.
/// @param bkg_pdf RooFit background pdf.
/// @param res Fit results struct
/// @param h_mll Input histogram for bin tuning.
/// @param o_dir Output directory.
/// @param bin_idx Bin index for fit identification.
/// @param tag Studied quantity.
/// @return 0 if success, else error code.
int SaveEventFitCanvas(RooRealVar& mll, RooAbsPdf& model, RooAbsData& data, RooAbsPdf& bkg_pdf, const EventFitResult& res, TH1D* h_mll, 
TDirectory* o_dir, const int bin_idx, const std::string& tag);


#endif