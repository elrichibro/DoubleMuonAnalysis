#ifndef FILTERS_H
#define FILTERS_H

#include <string>
#include <vector>
#include <utility>

#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

#include "Config.h"

// ------------------------------------------------------------------------------------------------------------------------------------

// -----------------
// Validation filter
// -----------------

/// @brief Apply the validation check to all the events in the dataset.
/// @param node Input RDF node.
/// @param val_map Struct containing the validated run and the relative luminosity blocks.
/// @param run_name Name of the "run" column in the dataset.
/// @param block_name Name of the "luminosity block" column in the dataset.
/// @return A RDF node with all validated events for post analysis.
ROOT::RDF::RNode ApplyValidationFilter(ROOT::RDF::RNode node, const validation_type& val_map, const std::string& run_name, const std::string& block_name);

// ------------------------------------------------------------------------------------------------------------------------------------

// ----------
// Acceptance
// ----------

/// @brief Calculates the acceptance from MonteCarlo sample.
/// @param node Input RDF node.
/// @param tag Tag for Column name.
/// @param FSR FSR Muon states.
/// @return Returns the value of detector acceptance in Z0->mu+mu- event. First element is central value, second value is the stat. sigma.
std::vector<float> CalculateAcceptance(ROOT::RDF::RNode node, const std::string& tag, int FSR);

struct ResolutionResults {
    std::vector<double> mean;
    std::vector<double> sigma;
    std::vector<int> events;
};

ResolutionResults CalculateResolution(ROOT::RDF::RNode node, const config_struct& cfg);
// ------------------------------------------------------------------------------------------------------------------------------------

// -------
// Structs
// -------

/// @brief 
struct EventHisto {
    ROOT::RDF::RResultPtr<TH1D> h1_mll;
    ROOT::RDF::RResultPtr<TH1D> h1_pt;
    ROOT::RDF::RResultPtr<TH1D> h1_y;
    ROOT::RDF::RResultPtr<TH1D> h1_phis;

    ROOT::RDF::RResultPtr<TH2D> h2_mll_pt;
    ROOT::RDF::RResultPtr<TH2D> h2_mll_y;
    ROOT::RDF::RResultPtr<TH2D> h2_mll_phis;
};

// ------------------------------------------------------------------------------------------------------------------------------------

// ---------
// Functions
// ---------

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief 
/// @param node 
/// @param cfg 
/// @return 
ROOT::RDF::RNode EventSelection(ROOT::RDF::RNode node, const config_struct& cfg);

/// @brief 
/// @param node 
/// @param cfg 
/// @return 
EventHisto BuildEventHisto(ROOT::RDF::RNode node, const config_struct& cfg);

#endif