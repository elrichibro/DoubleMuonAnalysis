#ifndef ANALYSISTOOLS_H
#define ANALYSISTOOLS_H

#include "Config.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include <ROOT/RDataFrame.hxx>

#include <string>
#include <vector>


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
/// @param node 
/// @param cfg 
/// @param entry_map 
/// @return 
std::vector<ROOT::RDF::RResultPtr<TH3D>> TemplateMaker(ROOT::RDF::RNode node, const config_struct& cfg, std::vector<ROOT::RDF::RResultPtr<TH2D>>& entry_map);

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
/// @param container 
/// @param cfg 
/// @param results 
/// @return 
int Eff_BinnedFit(std::vector<Template_RooF>& container, const analysis_params& params, std::vector<FitResult>& results, const int verb);

/// @brief 
/// @param results 
/// @param cfg 
/// @param vals 
/// @return 
int SaveFitPlots(const std::vector<FitResult>& results, const config_struct& cfg, const std::vector<std::string>& vals);

#endif