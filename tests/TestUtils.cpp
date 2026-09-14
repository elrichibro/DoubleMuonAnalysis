#include <gtest/gtest.h>
#include <cmath>
#include <TMath.h>
#include <ROOT/RVec.hxx>

#include "Utils.h"

TEST(KinematicTest, InvariantMassPair_particles_at_rest) {
    constexpr double m0 = 0.105;
    double m_inv = CalculateInvariantMass_Pair<double>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, m0, m0);
    
    EXPECT_NEAR(m_inv, 2.0 * m0, 1e-6);
}

TEST(KinematicTest, InvariantMassPair_back_to_back) {
    constexpr float pt = 50.0;
    constexpr float mass = 0.105;
    
    float m_inv = CalculateInvariantMass_Pair<float>(pt, pt, 0.0, 0.0, 0.0, 0.0 + static_cast<float>(M_PI), mass, mass);
    float m_inv_exp = 2.0 * std::sqrt((pt * pt) + (mass * mass));

    EXPECT_NEAR(m_inv, m_inv_exp, 1e-4);
}

template <typename T>
void TestInvMass_both_functions() {
    ROOT::RVec<T> pt{static_cast<T>(45.0), static_cast<T>(35.0)};
    ROOT::RVec<T> eta{static_cast<T>(0.1), static_cast<T>(-0.5)};
    ROOT::RVec<T> phi{static_cast<T>(0.2), static_cast<T>(3.0)};
    ROOT::RVec<T> mass{static_cast<T>(0.105), static_cast<T>(0.105)};


    T m_inv_vec = CalculateInvariantMass<T>(pt, eta, phi, mass);
    T m_inv_pair = CalculateInvariantMass_Pair<T>(pt[0], pt[1], eta[0], eta[1], phi[0], phi[1], mass[0], mass[1]);

    EXPECT_NEAR(m_inv_vec, m_inv_pair, static_cast<T>(1e-4));
}

TEST(KinematicTest, InvariantMass_both_functions) {
    TestInvMass_both_functions<float>();
    TestInvMass_both_functions<double>();
}

TEST(KinematicTest, PhiStar_back_to_back) {
    ROOT::RVec<double> eta{0.5, -0.5};
    ROOT::RVec<double> phi{0.0, M_PI};

    double phi_star = CalculatePhiStar<double>(eta, phi);
    
    EXPECT_NEAR(phi_star, 0.0, 1e-6);
}

TEST(KinematicTest, PhiStar_positive) {
    ROOT::RVec<double> eta{0.2, 0.8};
    ROOT::RVec<double> phi{0.1, 2.5};

    double phi_star = CalculatePhiStar<double>(eta, phi);
    
    EXPECT_GT(phi_star, 0.0);

    EXPECT_FALSE(std::isnan(phi_star));
}