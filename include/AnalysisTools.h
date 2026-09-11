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

struct Template_RooF{
    int eta_bin_idx;
    int pt_bin_idx;
    
    float eta_min;
    float eta_max;
    
    float pt_min;
    float pt_max;
    
    TH1D* h_MC_pass{nullptr};
    TH1D* h_MC_fail{nullptr};
    TH1D* h_DATA_pass{nullptr};
    TH1D* h_DATA_fail{nullptr};

    RooDataSet* d_MC_pass{nullptr};
    RooDataSet* d_MC_fail{nullptr};
    RooDataSet* d_DATA_pass{nullptr};
    RooDataSet* d_DATA_fail{nullptr};
};

struct FitResult {
    int eta_bin_idx;
    int pt_bin_idx;
    
    double efficiency;
    double efficiency_err;
    
    double n_tot;
    double n_tot_err;
    
    double mu;
    double mu_err;
    
    double sigma;
    double sigma_err;

    double lambda_pass;
    double lambda_pass_err;

    double lambda_fail;
    double lambda_fail_err;
    
    int fit_status;
};

/// @brief 
/// @param node_DATA 
/// @param node_MC 
/// @param cfg 
/// @return 
int TemplateMaker(ROOT::RDF::RNode node, const config_struct& cfg, const int dataset);

/// @brief 
/// @param tree 
/// @param branch_name 
/// @param mll 
/// @param ds_name 
/// @return 
RooDataSet* LoadRVecIntoDataset(TTree* tree, const std::string& branch_name, RooRealVar& mll, const std::string& ds_name);

/// @brief 
/// @param cfg 
/// @param container 
/// @return 
int LoadTemplate(const config_struct& cfg, std::vector<Template_RooF>& container);

/// @brief 
/// @param container 
/// @param cfg 
void CheckPlotsTemplate(const std::vector<Template_RooF>& container, const config_struct& cfg);

/// @brief 
/// @param analysis_struct 
/// @param cfg 
/// @param results 
/// @param o_file 
/// @return 
int Eff_BinnedFit(std::vector<Template_RooF>& analysis_struct, const config_struct& cfg, std::vector<FitResult>& results, TFile* o_file);

/// @brief 
/// @param results 
/// @param cfg 
/// @param o_file 
/// @param vals 
/// @return 
int SaveMapFittedValues(const std::vector<FitResult>& results, const config_struct& cfg, TFile* o_file, const std::vector<std::string>& vals);

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


#endif