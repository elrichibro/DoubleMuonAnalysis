#ifndef UTILS_H
#define UTILS_H

#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

#include <TMath.h>

#include <cmath>
#include <cstdint>

// ------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Calculates the invariant mass of the first two particles in the event.
 * @tparam Template: float, double.
 * @param pt ROOT Vector containing the transverse momentum of the particles.
 * @param eta Pseudorapidity.
 * @param phi Angular variable in cilindrical cordinates.
 * @param mass Mass values of the event particles.
 * @return Returns the invariant mass.
 */
template <typename T>
T CalculateInvariantMAss(const ROOT::RVec<T>& pt, const ROOT::RVec<T>& eta, const ROOT::RVec<T>& phi, const ROOT::RVec<T>& mass) {
    
    return ROOT::VecOps::InvariantMass(ROOT::RVec<T>{pt[0], pt[1]}, ROOT::RVec<T>{eta[0], eta[1]}, ROOT::RVec<T>{phi[0], phi[1]},
        ROOT::RVec<T>{mass[0], mass[1]} );

}

// ------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Calculates the special angular variable between the first two particles in the event.
 * @tparam Template: float, double.
 * @param eta Pseudorapidity of the particles.
 * @param phi Angular variable in cilindrical cordinates.
 * @return Returns special angular variable.
 */
template <typename T>
T CalculatePhiStar(const ROOT::RVec<T>& eta, const ROOT::RVec<T>& phi) {
    T delta_phi = std::abs(ROOT::VecOps::DeltaPhi(phi[0], phi[1]));
    T delta_eta = std::abs(eta[0] - eta[1]);

    T cos = std::tanh(delta_eta / 2.0);
    T sin = std::sqrt(1.0 - (cos * cos));

    return std::tan((TMath::Pi() - delta_phi) / 2.0) * sin;
}

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RNode ApplyPt_DataDivision(ROOT::RDF::RNode node, const std::vector<float>& pt_bins, const std::string& column_name) {
    ROOT::RDF::RNode div_pt_node = node
        .Define("Delta1_Pt", [&pt_bins](const ROOT::RVec<float> pt) {
            return ((pt >= pt_bins[0]) && (pt < pt_bins[1])); 
        }, {column_name})
        .Define("Delta2_Pt", [&pt_bins](const ROOT::RVec<float> pt) {
            return (pt >= pt_bins[1]) && (pt < pt_bins[2]); 
        }, {column_name})
        .Define("Delta3_Pt", [&pt_bins](const ROOT::RVec<float> pt) {
            return (pt >= pt_bins[2]) && (pt < pt_bins[3]); 
        }, {column_name})
        .Define("Delta4_Pt", [&pt_bins](const ROOT::RVec<float> pt) {
            return (pt >= pt_bins[3]) && (pt < pt_bins[4]); 
        }, {column_name})
        .Define("Delta5_Pt", [&pt_bins](const ROOT::RVec<float> pt) {
            return (pt >= pt_bins[4]) && (pt < pt_bins[5]); 
        }, {column_name})
        .Define("Delta6_Pt", [&pt_bins](const ROOT::RVec<float> pt) {
            return (pt >= pt_bins[5]) && (pt < pt_bins[6]); 
        }, {column_name});

    return div_pt_node;
}

ROOT::RDF::RNode ApplyEta_DataDivision(ROOT::RDF::RNode node, const std::vector<float>& eta_bins, const std::string& column_name) {
    ROOT::RDF::RNode div_eta_node = node
        .Define("Delta1_Eta", [&eta_bins](const ROOT::RVec<float> eta) {
            return ((eta >= eta_bins[0]) && (eta < eta_bins[1])); 
        }, {column_name})
        .Define("Delta2_Eta", [&eta_bins](const ROOT::RVec<float> eta) {
            return (eta >= eta_bins[1]) && (eta < eta_bins[2]); 
        }, {column_name})
        .Define("Delta3_Eta", [&eta_bins](const ROOT::RVec<float> eta) {
            return (eta >= eta_bins[2]) && (eta < eta_bins[3]); 
        }, {column_name})
        .Define("Delta4_Eta", [&eta_bins](const ROOT::RVec<float> eta) {
            return (eta >= eta_bins[3]) && (eta < eta_bins[4]); 
        }, {column_name})
        .Define("Delta5_Eta", [&eta_bins](const ROOT::RVec<float> eta) {
            return (eta >= eta_bins[4]) && (eta < eta_bins[5]); 
        }, {column_name})
        .Define("Delta6_Eta", [&eta_bins](const ROOT::RVec<float> eta) {
            return (eta >= eta_bins[5]) && (eta < eta_bins[6]); 
        }, {column_name})
        .Define("Delta7_Eta", [&eta_bins](const ROOT::RVec<float> eta) {
            return (eta >= eta_bins[6]) && (eta < eta_bins[7]); 
        }, {column_name})
        .Define("Delta8_Eta", [&eta_bins](const ROOT::RVec<float> eta) {
            return (eta >= eta_bins[7]) && (eta < eta_bins[8]); 
        }, {column_name});

    return div_eta_node;
}


#endif