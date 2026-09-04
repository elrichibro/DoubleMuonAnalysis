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
T CalculateInvariantMass(const ROOT::RVec<T>& pt, const ROOT::RVec<T>& eta, const ROOT::RVec<T>& phi, const ROOT::RVec<T>& mass) {
    
    return ROOT::VecOps::InvariantMass(ROOT::RVec<T>{pt[0], pt[1]}, ROOT::RVec<T>{eta[0], eta[1]}, ROOT::RVec<T>{phi[0], phi[1]},
        ROOT::RVec<T>{mass[0], mass[1]} );

}

// ------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
T CalculateInvariantMass_Pair(const float pt1, const float pt2, const float eta1, const float eta2, const float phi1, const float phi2, 
    const float mass1, const float mass2) {
    
    float pt = pt1 - pt2;
    float phi = phi1 - phi2;
    float eta = eta1 - eta2;
    
    float px = pt * std::cos(phi);
    float py = pt * std::sin(phi);
    float pz = pt * std::sinh(eta);
    
    return std::sqrt((px * px) + (py * py) + (pz * pz) + (mass1 - mass1));
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

#endif