#include "Unfold.h"
#include "Utils.h"
#include <vector>

ResultsRespMatrix CalculateRespMatrix(const MuonKinematics_RM& kin_rec, const MuonKinematics_RM& kin_gen, const MuonFlags_RM& flags, const flags_config cfg_f, 
const cuts_config cfg_c) {
    
    ResultsRespMatrix results;

    int n_tight = ROOT::VecOps::Sum(flags.reco_tight);

    // Fast exit
    if (n_tight > 2) {
        return results;
    }
    
    std::vector<int> reco_idx;
    std::vector<int> gen_idx;

    // Number reconstructed muons
    const int n_muons_rec = kin_rec.pt.size();
    const int n_muons_gen = kin_gen.pt.size();

    float m_rec = 0;
    bool pass_reco = false;
    int rec1 = -1;
    int rec2 = -1;

    float m_gen = 0;
    bool pass_gen = false;
    int gen1 = -1;
    int gen2 = -1;

    for (int i = 0; i < n_muons_rec; i++) {
        if (flags.reco_tight[i] != 1) {
            continue;
        }
        bool pass_kin = (((kin_rec.pt[i] > cfg_c.pt_cut) && (std::abs(kin_rec.eta[i]) < cfg_c.eta_cut)) || (!cfg_f.en_kinematics));
        
        if (pass_kin) {
            reco_idx.push_back(i);
        }
    }

    if (reco_idx.size() == 2) {
        rec1 = reco_idx[0];
        rec2 = reco_idx[1];
        
        if (kin_rec.charge[rec1] != kin_rec.charge[rec2]) {
            m_rec = CalculateInvariantMass_Pair<float>(kin_rec.pt[rec1], kin_rec.pt[rec2], kin_rec.eta[rec1], kin_rec.eta[rec2], kin_rec.phi[rec1],
            kin_rec.phi[rec2], kin_rec.mass[rec1], kin_rec.mass[rec2]);
            
            if ((m_rec > 60) && (m_rec < 120)) {
                pass_reco = true;
            }
        }
    }


    for (int i = 0; i < n_muons_gen; i++) {
        if (abs(flags.gen_pdg_idx[i]) != 13) {
            continue;
        }
    
        bool mask = ((flags.gen_status_flg[i] & (1 << 0)) > 0) && (((flags.gen_status_flg[i] & (1 << 7)) > 0) || ((flags.gen_status_flg[i] & (1 << 8)) > 0)) && 
        (((flags.gen_status_flg[i] & (1 << 13)) > 0) || ((flags.gen_status_flg[i] & (1 << 14)) > 0));

        if (!mask) {
            continue;
        }
        
        bool pass_kin = (((kin_gen.pt[i] > cfg_c.pt_cut) && (std::abs(kin_gen.eta[i]) < cfg_c.eta_cut)) || (!cfg_f.en_kinematics));
        
        if (pass_kin) {
            gen_idx.push_back(i);
        }
    }


    if (gen_idx.size() == 2) {
        gen1 = gen_idx[0];
        gen2 = gen_idx[1];
        
        if (flags.gen_pdg_idx[gen1] != flags.gen_pdg_idx[gen2]) {
            m_gen = CalculateInvariantMass_Pair<float>(kin_gen.pt[gen1], kin_gen.pt[gen2], kin_gen.eta[gen1], kin_gen.eta[gen2],
             kin_gen.phi[gen1], kin_gen.phi[gen2], kin_gen.mass[gen1], kin_gen.mass[gen2]);
            
            if ((m_gen > 60) && (m_gen < 120)) {
                pass_gen = true;
            }
        }
    }

    bool matched = false;

    if (pass_reco && pass_gen) {
        bool match1 = ((flags.reco_idx_gen[rec1] == gen1) && (flags.reco_idx_gen[rec2] == gen2));
        bool match2 = ((flags.reco_idx_gen[rec1] == gen2) && (flags.reco_idx_gen[rec2] == gen1));
        
        matched = (match1 || match2);
    }

    if (matched) {
        results.match = true;
        
        results.mll_rec = m_rec;
        results.mll_gen = m_gen;
        
        results.pt_gen = CalculatePtZ0_Raw_Pair<float>(kin_gen.pt[gen1], kin_gen.pt[gen2], kin_gen.phi[gen1], kin_gen.phi[gen2]);
        results.pt_rec = CalculatePtZ0_Raw_Pair<float>(kin_rec.pt[rec1], kin_rec.pt[rec2], kin_rec.phi[rec1], kin_rec.phi[rec2]);
        
        results.y_gen = CalculateRapidityZ0_Raw_Pair<float>(kin_rec.pt[gen1], kin_rec.pt[gen2], kin_rec.eta[gen1], kin_rec.eta[gen2],
            kin_rec.phi[gen1], kin_rec.phi[gen2],  kin_rec.mass[gen1], kin_rec.mass[gen2]);
        
        results.y_rec = CalculateRapidityZ0_Raw_Pair<float>(kin_rec.pt[rec1], kin_rec.pt[rec2], kin_rec.eta[rec1], kin_rec.eta[rec2],
            kin_rec.phi[rec1], kin_rec.phi[rec2],  kin_rec.mass[rec1], kin_rec.mass[rec2]);

        results.phis_gen = CalculatePhiStar_Pair<float>(kin_rec.eta[gen1], kin_rec.eta[gen2], kin_rec.phi[gen1], kin_rec.phi[gen2]);
        results.phis_rec = CalculatePhiStar_Pair<float>(kin_rec.eta[rec1], kin_rec.eta[rec1], kin_rec.phi[rec2], kin_rec.phi[rec2]);
    
    } else if ((pass_gen) && (!pass_reco)) {
        results.miss = true;
        results.mll_gen = m_gen;

        results.pt_gen = CalculatePtZ0_Raw_Pair<float>(kin_gen.pt[gen1], kin_gen.pt[gen2], kin_gen.phi[gen1], kin_gen.phi[gen2]);

        results.y_gen = CalculateRapidityZ0_Raw_Pair<float>(kin_rec.pt[gen1], kin_rec.pt[gen2], kin_rec.eta[gen1], kin_rec.eta[gen2],
            kin_rec.phi[gen1], kin_rec.phi[gen2],  kin_rec.mass[gen1], kin_rec.mass[gen2]);

        results.phis_gen = CalculatePhiStar_Pair<float>(kin_rec.eta[gen1], kin_rec.eta[gen2], kin_rec.phi[gen1], kin_rec.phi[gen2]);

    } else if (pass_reco && !matched) {
        results.fake = true;
        results.mll_rec = m_rec;

        results.pt_rec = CalculatePtZ0_Raw_Pair<float>(kin_rec.pt[rec1], kin_rec.pt[rec2], kin_rec.phi[rec1], kin_rec.phi[rec2]);

        results.y_rec = CalculateRapidityZ0_Raw_Pair<float>(kin_rec.pt[rec1], kin_rec.pt[rec2], kin_rec.eta[rec1], kin_rec.eta[rec2],
            kin_rec.phi[rec1], kin_rec.phi[rec2],  kin_rec.mass[rec1], kin_rec.mass[rec2]);

        results.phis_rec = CalculatePhiStar_Pair<float>(kin_rec.eta[rec1], kin_rec.eta[rec1], kin_rec.phi[rec2], kin_rec.phi[rec2]);
    }

    if (gen_idx.size() > 2) {
        std::cout << "Rare situation: " << gen_idx.size() << std::endl;
    }

    return results;
}

ROOT::RDF::RNode ApplyUnfold(ROOT::RDF::RNode node) {
    
}

