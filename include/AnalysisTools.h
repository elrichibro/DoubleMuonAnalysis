#ifndef ANALYSISTOOLS_H
#define ANALYSISTOOLS_H

#include "TH1D.h"
#include "Config.h"
#include <ROOT/RDataFrame.hxx>


struct MC_Template_RooF{
    int pt_bin_idx;
    int eta_bin_idx;
    float eta_min;
    float eta_max;
    float pt_min;
    float pt_max;
    
    TH1D* h_pass{nullptr};
    TH1D* h_fail{nullptr};
};

int LoadMCTemplate(const config_struct& cfg, std::vector<MC_Template_RooF>& container);



// ------------------------------------------------------------------------------------------------------------------------------------
// MC Template Maker
// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RResultPtr<TH3D> TemplateMaker_MC(ROOT::RDF::RNode node, const config_struct& cfg, const bool mask);


#endif