#ifndef UTILS_H
#define UTILS_H

#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

#include <TMath.h>
#include <Math/Vector4D.h>

#include <cmath>
#include <cstdint>

#include <vector>
#include <string>

// ------------------------------------------------------------------------------------------------------------------------------------

/*
This file contains:
    - CalculateInvariantMass()
    - CalculateInvariantMass_Pair()
    - CalculatePhiStar()
*/

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Calculates the invariant mass of the first two particles in the event.
/// @tparam T Template: float, double.
/// @param pt ROOT Vector containing the transverse momentum of the particles.
/// @param eta Pseudorapidity.
/// @param phi Angular variable in cilindrical cordinates.
/// @param mass Mass values of the event particles.
/// @return Returns the invariant mass.
template <typename T>
T CalculateInvariantMass(const ROOT::RVec<T>& pt, const ROOT::RVec<T>& eta, const ROOT::RVec<T>& phi, const ROOT::RVec<T>& mass) {
    
    return ROOT::VecOps::InvariantMass(ROOT::RVec<T>{pt[0], pt[1]}, ROOT::RVec<T>{eta[0], eta[1]}, ROOT::RVec<T>{phi[0], phi[1]},
        ROOT::RVec<T>{mass[0], mass[1]});
}

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Calculates the invariant mass of a pair of particles givin the individual quantities.
/// @tparam T float, double
/// @param pt1 Transverse momentum of first particle.
/// @param pt2 Transverse momentum of second particle.
/// @param eta1 Pseudorapidity of the first particle.
/// @param eta2 Pseudorapidity of the second particle.
/// @param phi1 Angular emission variable of first particle.
/// @param phi2 Angular emission variable of second particle.
/// @param mass1 Mass quantity of first particle.
/// @param mass2 Mass quantity of second particle.
/// @return The invariant mass of the pair.
template <typename T>
T CalculateInvariantMass_Pair(const T pt1, const T pt2, const T eta1, const T eta2, const T phi1, const T phi2, const T mass1, const T mass2) {
    
    const T px1 = pt1 * std::cos(phi1);
    const T py1 = pt1 * std::sin(phi1);
    const T pz1 = pt1 * std::sinh(eta1);
    const T E1 = std::sqrt((px1 * px1) + (py1 * py1) + (pz1 * pz1) + (mass1 * mass1));

    const T px2 = pt2 * std::cos(phi2);
    const T py2 = pt2 * std::sin(phi2);
    const T pz2 = pt2 * std::sinh(eta2);
    const T E2 = std::sqrt((px2 * px2) + (py2 * py2) + (pz2 * pz2) + (mass2 * mass2));

    const T px = px1 + px2;
    const T py = py1 + py2;
    const T pz = pz1 + pz2; 
    const T E = E1 + E2;

    const T mass_sqr = (E * E) - ((px * px) + (py * py) + (pz * pz));
    const T inv_mass = (mass_sqr > static_cast<T>(0)) ? std::sqrt(mass_sqr) : 0; 
    
    return inv_mass;
}

// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Calculates the special angular variable between the first two particles in the event.
/// @tparam T : float, double.
/// @param eta Pseudorapidity of the particles.
/// @param phi Angular variable in cilindrical cordinates.
/// @return Returns special angular variable.
template <typename T>
T CalculatePhiStar(const ROOT::RVec<T>& eta, const ROOT::RVec<T>& phi) {
    T delta_phi = std::abs(ROOT::VecOps::DeltaPhi(phi[0], phi[1]));
    T delta_eta = std::abs(eta[0] - eta[1]);

    T cos = std::tanh(delta_eta / 2.0);
    T sin = std::sqrt(1.0 - (cos * cos));

    return std::tan((TMath::Pi() - delta_phi) / 2.0) * sin;
}

template <typename T>
T CalculatePhiStar_Pair(const T eta1, const T eta2, const T phi1, const T phi2) {
    T delta_phi = std::abs(ROOT::VecOps::DeltaPhi(phi1, phi2));
    T delta_eta = std::abs(eta1 - eta2);

    T cos = std::tanh(delta_eta / 2.0);
    T sin = std::sqrt(1.0 - (cos * cos));

    return std::tan((TMath::Pi() - delta_phi) / 2.0) * sin;
}

// ------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
T CalculatePtZ0(const ROOT::RVec<T>& pt, const ROOT::RVec<T>& eta, const ROOT::RVec<T>& phi, const ROOT::RVec<T>& mass) {
    ROOT::Math::PtEtaPhiMVector mu1(pt[0], eta[0], phi[0], mass[0]);
    ROOT::Math::PtEtaPhiMVector mu2(pt[1], eta[1], phi[1], mass[1]);
        
    auto z0 = mu1 + mu2;   
    T pt_Z0 = static_cast<T>(z0.Pt()); 
    
    return pt_Z0;
}

template <typename T>
T CalculatePtZ0_Raw(const ROOT::RVec<T>& pt, const ROOT::RVec<T>& phi) {
    T px_Z0 = (pt[0] * std::cos(phi[0])) + (pt[1] * std::cos(phi[1]));
    T py_Z0 = (pt[0] * std::sin(phi[0])) + (pt[1] * std::sin(phi[1]));

    T pt_Z0 = std::sqrt((px_Z0 * px_Z0) + (py_Z0 * py_Z0));

    return pt_Z0;
}

template <typename T>
T CalculatePtZ0_Raw_Pair(const T pt1, const T pt2, const T phi1, const T phi2) {
    T px_Z0 = (pt1 * std::cos(phi1)) + (pt2 * std::cos(phi2));
    T py_Z0 = (pt1 * std::sin(phi1)) + (pt2 * std::sin(phi2));

    T pt_Z0 = std::sqrt((px_Z0 * px_Z0) + (py_Z0 * py_Z0));

    return pt_Z0;
}

// ------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
T CalculateRapidityZ0(const ROOT::RVec<T>& pt, const ROOT::RVec<T>& eta, const ROOT::RVec<T>& phi, const ROOT::RVec<T>& mass) {
    ROOT::Math::PtEtaPhiMVector mu1(pt[0], eta[0], phi[0], mass[0]);
    ROOT::Math::PtEtaPhiMVector mu2(pt[1], eta[1], phi[1], mass[1]);
        
    auto z0 = mu1 + mu2;

    T y_Z0 = static_cast<T>(z0.Rapidity());
    return y_Z0; 
}

template <typename T>
T CalculateRapidityZ0_Raw(const ROOT::RVec<T>& pt, const ROOT::RVec<T>& eta, const ROOT::RVec<T>& phi, const ROOT::RVec<T>& mass) {
    ROOT::RVec<T> pz_mu = pt * std::sinh(eta);

    ROOT::RVec<T> E = std::sqrt((pt * pt) + (pz_mu * pz_mu) + (mass * mass));
    
    T E_tot = E[0] + E[1];
    T pz_tot = pz_mu[0] + pz_mu[1];

    T y_Z0 = static_cast<T>(0.5 * std::log((E_tot + pz_tot) / (E_tot - pz_tot)));
    
    return y_Z0; 
}

template <typename T>
T CalculateRapidityZ0_Raw_Pair(const T pt1, const T pt2, const T eta1, const T eta2, const T phi1, const T phi2, const T mass1, const T mass2) {
    T pz_mu1 = pt1 * std::sinh(eta1);
    T pz_mu2 = pt2 * std::sinh(eta2);

    T E1 = std::sqrt((pt1 * pt1) + (pz_mu1 * pz_mu1) + (mass1 * mass1));
    T E2 = std::sqrt((pt2 * pt2) + (pz_mu2 * pz_mu2) + (mass2 * mass2));
    
    T E_tot = E1 + E2;
    T pz_tot = pz_mu1 + pz_mu2;

    T y_Z0 = static_cast<T>(0.5 * std::log((E_tot + pz_tot) / (E_tot - pz_tot)));
    
    return y_Z0; 
}

std::vector<double> CreateBins(int nbins, double min, double max, const std::string& distribution, double split = 0.0);

#endif