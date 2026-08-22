#ifndef UTILS_H
#define UTILS_H

#include <ROOT/RVec.hxx>
#include <cmath>

/**
 * @brief Calculate the invariant mass of the first two particles in the event.
 * @tparam Template: float, double.
 * @param pt Is the RVec containing the transverse momentum of the particles.
 * @param eta Is the pseudorapidity of the particles.
 * @param phi Is the angular variable in cilindrical cordinates.
 * @param mass Mass values of the particles.
 * @return Returns the invariant mass depending on the input type.
 */
template <typename T>
T mass_leptons(const ROOT::RVec<T>& pt, const ROOT::RVec<T>& eta, const ROOT::RVec<T>& phi, const ROOT::RVec<T>& mass) {
    
    return ROOT::VecOps::InvariantMass(ROOT::RVec<T>{pt[0], pt[1]}, ROOT::RVec<T>{eta[0], eta[1]}, ROOT::RVec<T>{phi[0], phi[1]},
        ROOT::RVec<T>{mass[0], mass[1]} );

}

/**
 * @brief Calculate the special angular variable between the first two particles in the event.
 * @tparam Template: float, double.
 * @param eta Is the pseudorapidity of the particles.
 * @param phi Is the angular variable in cilindrical cordinates.
 * @return Returns special angular variable depending on the input type.
 */
template <typename T>
T phi_star(const ROOT::RVec<T>& eta, const ROOT::RVec<T>& phi) {
    T delta_phi = std::abs(ROOT::VecOps::DeltaPhi(phi[0], phi[1]));
    T delta_eta = std::abs(eta[0] - eta[1]);

    T cos = std::tanh(delta_eta / 2.0);
    T sin = std::sqrt(1.0 - (cos * cos));

    return std::tan((M_PI - delta_phi) / 2.0) * sin;
}


#endif
