#ifndef UTILS_H
#define UTILS_H

#include <ROOT/RVec.hxx>
#include <cmath>
#include <cstdint>
#include <fstream>

#include <nlohmann/json.hpp>


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

using validation_type = std::unordered_map<std::uint32_t, std::vector<std::pair<std::uint16_t, std::uint16_t>>>;

inline validation_type Validation_load(const std::string& json_path) {
    validation_type validation_map;
    
    std::ifstream file(json_path);
    if (!file.is_open()) {
        throw std::runtime_error("Error: cannot open " + json_path);
    }

    nlohmann::json json_data = nlohmann::json::parse(file);// JSON CORE !!!

    for (const auto& [run_number, lum_block] : json_data.items()) {
        std::uint32_t run = static_cast<std::uint32_t>(std::stoul(run_number));

        std::vector<std::pair<std::uint16_t, std::uint16_t>> blocks;
        blocks.reserve(lum_block.size());

        for (const auto& iter : lum_block) {
            std::uint16_t start = static_cast<std::uint16_t>(iter[0].get<unsigned int>());
            std::uint16_t end   = static_cast<std::uint16_t>(iter[1].get<unsigned int>());
            blocks.emplace_back(start, end);
        }
        validation_map[run] = std::move(blocks);
    }
    return validation_map;
}

#endif
