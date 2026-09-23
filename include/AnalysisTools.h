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

// ----------------
// Analysis Structs
// ----------------

/// @brief Struct for Template load
struct Template_RooF{
    int eta_bin_idx;// Eta template bin index.
    int pt_bin_idx;// Pt template bin index.
    
    std::unique_ptr<TH1D>h_MC_pass{nullptr};// TH1D container for MonteCarlo (passed sample).
    std::unique_ptr<TH1D>h_MC_fail{nullptr};// TH1D container for MonteCarlo (failed sample).
    std::unique_ptr<TH1D>h_DATA_pass{nullptr};// TH1D container for Data (passed sample).
    std::unique_ptr<TH1D>h_DATA_fail{nullptr};// TH1D container for Data (failed sample).

    std::unique_ptr<RooDataSet>d_MC_pass{nullptr};// RooDataSet container for MonteCarlo (passed sample).
    std::unique_ptr<RooDataSet>d_MC_fail{nullptr};// RooDataSet container for MonteCarlo (failed sample).
    std::unique_ptr<RooDataSet>d_DATA_pass{nullptr};// RooDataSet container for Data (passed sample).
    std::unique_ptr<RooDataSet>d_DATA_fail{nullptr};// RooDataSet container for Data (failed sample).
};

/// @brief Fit results container
struct FitResult {
    int eta_bin_idx;// Eta bin index.
    int pt_bin_idx;// Pt bin index.
    
    double efficiency;// Efficiency value.
    double efficiency_err;// Efficiency fit error.
    
    double n_tot;// Total number of signal events.
    double n_tot_err;// Total number of signal events fit error.
    
    double mu;// Gaussian mean.
    double mu_err;// Gaussian mean fit error.
    
    double sigma;// Gaussian sigma.
    double sigma_err;// Gaussian sigma fit error.

    double lambda_pass;// Lambda PASS of background model.
    double lambda_pass_err;// Lambda PASS fit error.

    double lambda_fail;// Lambda FAIL of background model.
    double lambda_fail_err;// Lambda FAIL fit error.
    
    int fit_status;// Fit status
};

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

// --------------
// Template maker
// --------------

/// @brief Creates an intermediate binned status of data optimizated for analysis process.
/// @param node RDF input node.
/// @param cfg Configure general struct.
/// @param dataset (1)DATA/(2)MC dataset type.
/// @return 0 if succes, else error code.
int BinnedTemplateMaker(ROOT::RDF::RNode node, const config_struct& cfg, const int dataset);

/// @brief Creates a flat tree -> 1 event = 1 muon. For RooFit/RooDataSet input.
/// @param node RDF input node.
/// @param cfg Configuration general struct.
/// @return 0 if succeds, else error code.
int UnbinnedTemplateMaker(ROOT::RDF::RNode node, const config_struct& cfg);

// ------------------------------------------------------------------------------------------------------------------------------------

// ---------------------------
// Template Loader - PreFitter
// ---------------------------

/// @brief Loads the template binned data into the system Template struct.
/// @param cfg Configure general struct.
/// @param container Container for template data (binned/unbinned formats -> lazy option).
/// @return 0 if succes, else error code.
int LoadBinnedTemplate(const config_struct& cfg, std::vector<Template_RooF>& container);

/// @brief Loads the tree unbinned data into the RooFit RooDataSet variable and moves into the template struct.
/// @param tree Input flat tree.
/// @param cfg General configuration struct.
/// @param dataset Dataset identifier: 1 = DATA, 2 = MC.
/// @param container Template system container. 
/// @return 0 if succeds, else error code.
int LoadUnbinnedTemplate(TTree* tree, const config_struct& cfg, const int dataset, std::vector<Template_RooF>& container);

/// @brief 
/// @param event_histo 
/// @param tag 
/// @return 
std::vector<std::unique_ptr<TH1D>> PrepareEventFit(EventHisto& event_histo, const std::string& tag);

std::vector<std::unique_ptr<RooDataSet>> PrepareEventFitModel(TTree* tree, const config_struct& cfg, const std::string& tag);

// ------------------------------------------------------------------------------------------------------------------------------------

// ---------------------
// Visualization methods
// ---------------------

/// @brief Saves the 2D plots of the fitted results in each bin of the bin_settup.
/// @param o_file Input file.
/// @param results Fits results.
/// @param cfg General configuration struct.
/// @param vals Name of the fit values booked.
/// @return 0 if succeds, else error code(1).
int SaveMapFittedValues(TFile* o_file, const std::vector<FitResult>& results, const config_struct& cfg, const std::vector<std::string>& vals);

/// @brief Saves the fit plot for each bin selection as a RooPlot.
/// @param mll Observable.
/// @param sample Pass/Fail category.
/// @param data Data of the sample -> binned/unbinned. 
/// @param simPdf Model.
/// @param res Results struct.
/// @param o_file Output file.
/// @param cfg General configuration struct.
void SaveBinFitCanvas(RooRealVar& mll, RooCategory& sample, RooAbsData& data, RooSimultaneous& simPdf, const FitResult& res, TFile* o_file, 
const config_struct& cfg);

void SaveEventFitCanvas(RooRealVar& mll, RooAbsPdf& model, RooAbsData& data, RooAbsPdf& bkg_pdf, const EventFitResult& res, TH1D* h_mll, 
    TDirectory* o_dir, const int bin_idx, const std::string& tag);

// ------------------------------------------------------------------------------------------------------------------------------------

// ------
// Fitter
// ------

/// @brief Fits the data vs the hardcoded model.
/// @param analysis_struct Data template struct.
/// @param cfg General configuration struct.
/// @param results Struct for fit results.
/// @param o_file Output file for SaveBinFitCanvas function.
/// @return 0 if succeds, else error code.
int EfficiencyFitter(std::vector<Template_RooF>& analysis_struct, const config_struct& cfg, std::vector<FitResult>& results, TFile* o_file);

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