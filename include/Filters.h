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

/// @brief Calculates the acceptance from MonteCarlo sample. Cuts values are fixed by definition -> can be overwritten manually.
/// @param node Input RDF node.
/// @param tag Tag for Column name.
/// @param FSR FSR Muon states.
/// @return Returns the value of detector acceptance in Z0->mu+mu- event. First element is central value, second value is the stat. sigma.
std::vector<float> CalculateAcceptance(ROOT::RDF::RNode node, const std::string& tag, int FSR);

// ------------------------------------------------------------------------------------------------------------------------------------

// ----------
// Resolution
// ----------

struct ResolutionResults {
    std::vector<double> mean;
    std::vector<double> sigma;
    std::vector<int> events;
};

/// @brief Calculates the resolution of the detector from MC sample
/// @param node RDF input node 
/// @param cfg General configure input
/// @return A struct containing aritmetic mean and standard deviation for each bin.
ResolutionResults CalculateResolution(ROOT::RDF::RNode node, const config_struct& cfg);

// ------------------------------------------------------------------------------------------------------------------------------------

// ------------
// Event struct
// ------------

/// @brief RDF Histograms obtained from DATA sample.
struct EventHisto {
    ROOT::RDF::RResultPtr<TH1D> h1_mll;// 1D histogram: Z0 invariant mass
    ROOT::RDF::RResultPtr<TH1D> h1_pt;// 1D histogram: Z0 Transverse momentum
    ROOT::RDF::RResultPtr<TH1D> h1_y;// 1D histogram: Z0 rapidity
    ROOT::RDF::RResultPtr<TH1D> h1_phis;// 1D histogram: Z0 phi*

    ROOT::RDF::RResultPtr<TH2D> h2_mll_pt;// 2D histogram: Z0 transverse momentum vs invariant mass
    ROOT::RDF::RResultPtr<TH2D> h2_mll_y;// 2D histogram: Z0 rapidity vs invariant mass
    ROOT::RDF::RResultPtr<TH2D> h2_mll_phis;// 2D histogram: Z0 phi* vs invariant mass
};

// ------------------------------------------------------------------------------------------------------------------------------------

// ---------------
// Event functions
// ---------------

/// @brief Filters the event selecting Z0/mu+mu- events falling withing the fiducial region.
/// @param node RDF input node: DATA/MC -> Event selection procedure
/// @param cfg General configure struct
/// @return RNode with filtered events and new variables defined.
ROOT::RDF::RNode EventSelection(ROOT::RDF::RNode node, const config_struct& cfg);

/// @brief Creates histograms from RDF input node and save them into EventHisto struct. The selection cuts/flags can be setted by JSON values.
/// @param node RDF input node.
/// @param cfg General config struct.
/// @return The histogram struct filled with event histograms.
EventHisto BuildEventHisto(ROOT::RDF::RNode node, const config_struct& cfg);

#endif