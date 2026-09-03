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
    ROOT::RDF::RNode div_pt_node = node;
    
    for(int i = 0; i < pt_bins.size() - 1; i++) {
        std::string name_mask = column_name + "_DeltaMask" + std::to_string(i);
        std::string name_column = column_name + "_Delta" + std::to_string(i);

        float min_pt = pt_bins[i];
        float max_pt = pt_bins[i + 1];
        
        div_pt_node = div_pt_node
            .Define(name_mask, [min_pt, max_pt](const ROOT::RVec<float>& pt) {
                return ((pt >= min_pt) && (pt < max_pt)); 
            }, {column_name});
    
        div_pt_node = div_pt_node
            .Define(name_column, [](const ROOT::RVec<float>& pt, const ROOT::RVec<bool>& mask) {
                return pt[mask];
            }, {column_name, name_mask} );
    }
    return div_pt_node;
}

ROOT::RDF::RNode ApplyEta_DataDivision(ROOT::RDF::RNode node, const std::vector<float>& eta_bins, const std::string& column_name) {
    ROOT::RDF::RNode div_eta_node = node;
    
    for(int i = 0; i < eta_bins.size() - 1; i++) {
        std::string name_mask = column_name + "_DeltaMask" + std::to_string(i);
        std::string name_column = column_name + "_Delta" + std::to_string(i);

        float min_eta = eta_bins[i];
        float max_eta = eta_bins[i + 1];
        
        div_eta_node = div_eta_node
            .Define(name_mask, [min_eta, max_eta](const ROOT::RVec<float>& eta) {
                return ((eta >= min_eta) && (eta < max_eta)); 
            }, {column_name});
    
        div_eta_node = div_eta_node
            .Define(name_column, [](const ROOT::RVec<float>& eta, const ROOT::RVec<bool>& mask) {
                return eta[mask];
            }, {column_name, name_mask} );
    }

    return div_eta_node;
}

ROOT::RDF::RNode ApplyPt_Columns(ROOT::RDF::RNode node, std::string column_name, const int n) {
    ROOT::RDF::RNode pt_columns = node;
    
    for (int i = 0; i < n; i++) {
        std::string name = column_name + "_Delta" + std::to_string(i);
        std::string selection = column_name + "[Delta" + std::to_string(i) + "_Pt]";
        pt_columns.Filter(name, selection);
    }
    
    return pt_columns;
}

#endif