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

#include <string>
#include <vector>

#include "Manager.h"// for snapshot type

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

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Creates an intermediate status of data optimizated for analysis process.
/// @param node RDF input node.
/// @param cfg Configure general struct.
/// @param dataset (1)DATA/(2)MC dataset type.
/// @return 0 if succes, else error code.
int TemplateMaker(ROOT::RDF::RNode node, const config_struct& cfg, const int dataset);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Loads RVectors into RooDataSet (scalar type needed).
/// @param tree Loaded TTree -> template data.
/// @param branch_name Specific bin settup branch.
/// @param mll Invariant mass variable for RooFit usage.
/// @param data_name For RooDataSet variable.
/// @return RooDataSet of a specific bins settup (template struct element).
std::unique_ptr<RooDataSet> LoadRVecIntoDataset(TTree* tree, const std::string& branch_name, RooRealVar& mll, const std::string& data_name);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Loads the template data into the system struct.
/// @param cfg Configure general struct.
/// @param container Container for template data (binned/unbinned formats -> lazy option).
/// @return 0 if succes, else error code.
int LoadTemplate(const config_struct& cfg, std::vector<Template_RooF>& container);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief 
/// @param container 
/// @param cfg 
void CheckPlotsTemplate(const std::vector<Template_RooF>& container, const config_struct& cfg);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief 
/// @param results 
/// @param cfg 
/// @param o_file 
/// @param vals 
/// @return 
int SaveMapFittedValues(const std::vector<FitResult>& results, const config_struct& cfg, TFile* o_file, const std::vector<std::string>& vals);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief 
/// @param mll 
/// @param sample 
/// @param data 
/// @param simPdf 
/// @param res 
/// @param o_file 
/// @param o_dir 
void SaveBinFitCanvas(RooRealVar& mll, RooCategory& sample, RooAbsData& data, RooSimultaneous& simPdf, const FitResult& res, TFile* o_file, 
const config_struct& cfg);

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief 
/// @param analysis_struct 
/// @param cfg 
/// @param results 
/// @param o_file 
/// @return 
int Eff_BinnedFit(std::vector<Template_RooF>& analysis_struct, const config_struct& cfg, std::vector<FitResult>& results, TFile* o_file);

#endif