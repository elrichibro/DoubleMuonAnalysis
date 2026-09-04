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