#include <gtest/gtest.h>
#include <ROOT/RVec.hxx>

#include "ControlCheck.h"

TEST(MonteCarloChecksTest, is_MC_Z0_valid_input) {
    ROOT::RVec<Int_t> pdgId{11, 23, -11, 23, -23, 21, 20};
    ROOT::RVec<Int_t> flags{0x2100, 0x2101, 0x0, 0x2020, 0x1902, 0x202, 0x1234};

    auto mask = is_MC_Z0(pdgId, flags);

    EXPECT_FALSE(mask[0]);
    EXPECT_TRUE(mask[1]);
    EXPECT_FALSE(mask[2]);
    EXPECT_FALSE(mask[3]);
    EXPECT_FALSE(mask[4]);
    EXPECT_FALSE(mask[5]);
    EXPECT_FALSE(mask[6]);

    EXPECT_EQ(get_MC_Z0_idx(pdgId, flags), 1);
}

TEST(MonteCarloChecksTest, get_MC_Z0_idx_not_found) {
    ROOT::RVec<Int_t> pdgId{13, -13, 23, 11, 10, 1, 23, -23};
    ROOT::RVec<Int_t> flags{0x2101, 0x2101, 0x1000, 0x2101, 0x0, 0x2020, 0x1902, 0x000};

    EXPECT_EQ(get_MC_Z0_idx(pdgId, flags), -1);
}

TEST(MonteCarloChecksTest, is_MC_mother_idx) {
    ROOT::RVec<Int_t> pdgId{1, 2, 3, 23, 13, -13};
    ROOT::RVec<Int_t> flags{0x0, 0x0, 0x0, 0x2101, 0x181, 0x181};
    ROOT::RVec<Int_t> mother_id{0, 0, 0, 3, 3, 3};

    auto res_mu = is_MC_Muon_bFSR(pdgId, flags, mother_id);
    auto res_anti = is_MC_AntiMuon_bFSR(pdgId, flags, mother_id);

    EXPECT_FALSE(res_mu[0]);
    EXPECT_FALSE(res_mu[1]);
    EXPECT_FALSE(res_mu[2]);
    EXPECT_FALSE(res_mu[3]);
    EXPECT_TRUE(res_mu[4]);
    EXPECT_FALSE(res_mu[5]);

    EXPECT_FALSE(res_anti[0]);
    EXPECT_FALSE(res_anti[1]);
    EXPECT_FALSE(res_anti[2]);
    EXPECT_FALSE(res_anti[3]);
    EXPECT_FALSE(res_anti[4]);
    EXPECT_TRUE(res_anti[5]);
}

TEST(MonteCarloChecksTest, is_MC_Z0_repeated) {
    ROOT::RVec<Int_t> pdgId{23, 23, 23, 13, -13};
    ROOT::RVec<Int_t> flags{0x2101, 0x2101, 0x2101, 0x181, 0x181};
    ROOT::RVec<Int_t> mother_id{0, 1, 2, 0, 3};

    auto res_mu = is_MC_Muon_bFSR(pdgId, flags, mother_id);
    auto res_anti = is_MC_AntiMuon_bFSR(pdgId, flags, mother_id);

    EXPECT_FALSE(res_mu[0]);
    EXPECT_FALSE(res_mu[1]);
    EXPECT_FALSE(res_mu[2]);
    EXPECT_TRUE(res_mu[3]);
    EXPECT_FALSE(res_mu[4]);
        
    EXPECT_FALSE(res_anti[0]);
    EXPECT_FALSE(res_anti[1]);
    EXPECT_FALSE(res_anti[2]);
    EXPECT_FALSE(res_anti[3]);
    EXPECT_FALSE(res_anti[4]);
}

TEST(MCChecksTest, is_MC_wrong_mother_idx) {
    ROOT::RVec<Int_t> pdgId{13, -13, 13};
    ROOT::RVec<Int_t> flags{0x2101, 0x2101, 0x0000};
    ROOT::RVec<Int_t> mother_id{-1, -1, -1};

    auto res_mu = is_MC_Muon_aFSR(pdgId, flags, mother_id);
    auto res_anti = is_MC_AntiMuon_aFSR(pdgId, flags, mother_id);

    EXPECT_TRUE(res_mu[0]);
    EXPECT_FALSE(res_mu[1]);
    EXPECT_FALSE(res_mu[2]);

    EXPECT_FALSE(res_anti[0]);
    EXPECT_TRUE(res_anti[1]);
    EXPECT_FALSE(res_anti[2]);
}