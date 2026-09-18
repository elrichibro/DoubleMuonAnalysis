#include "Unfold.h"
#include "Utils.h"

#include <vector>

/*
Reconstructed Muon:
    - Tight flag
    - pt > 25 GeV
    - |eta| < 2.4 

Generated Muon:
    - Pdg index = 13/-13
    - Generator status flags: 0 + (7 || 8) + (13 || 14)
    - pt > 25 GeV
    - |eta| < 2.4

Reconstructed pair:
    - 2 reconstructed muons
    - Opposite charge
    - Invariant mass range -> 60 < m_ll < 120 GeV

Generated pair:
    - 2 Generated muons
    - Opposite pdg index0
    - Invariant mass range -> 60 < m_ll < 120 GeV

Match :
    - 2 posibilities -> match between Generator index and mother index for reconstructed particle (that is the muon in the GenPart collection)
*/


/// @brief 
/// @param kin_rec 
/// @param kin_gen 
/// @param flags 
/// @param cfg_f 
/// @param cfg_c 
/// @return 
ResultsRespMatrix CalculateRespMatrix(const MuonKinematics_REC& kin_rec, const MuonKinematics_GEN& kin_gen, const MuonFlags_RM& flags, const flags_config cfg_f, 
const cuts_config cfg_c) {
    
    ResultsRespMatrix results;

    int n_tight = ROOT::VecOps::Sum(flags.reco_tight);

    // Fast exit -> Experiment cut
    if (n_tight > 2) {
        return results;
    }
    
    std::vector<int> reco_idx;
    std::vector<int> gen_idx;

    // Number reconstructed muons
    const int n_muons_rec = kin_rec.pt.size();
    const int n_muons_gen = kin_gen.pt.size();

    float m_rec = 0;
    int rec1 = -1;
    int rec2 = -1;

    float m_gen = 0;
    int gen1 = -1;
    int gen2 = -1;

    bool matched = false;
    bool pass_gen = false;
    bool pass_reco = false;

    // Loop on reconstructed muons -> searching for a muon in fiducial region.
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
        if (abs(flags.gen_pdg_idx[i]) != 13 || (flags.gen_status[i] != 1)) {
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
        
        results.y_gen = CalculateRapidityZ0_Raw_Pair<float>(kin_gen.pt[gen1], kin_gen.pt[gen2], kin_gen.eta[gen1], kin_gen.eta[gen2],
            kin_gen.phi[gen1], kin_gen.phi[gen2],  kin_gen.mass[gen1], kin_gen.mass[gen2]);
        
        results.y_rec = CalculateRapidityZ0_Raw_Pair<float>(kin_rec.pt[rec1], kin_rec.pt[rec2], kin_rec.eta[rec1], kin_rec.eta[rec2],
            kin_rec.phi[rec1], kin_rec.phi[rec2],  kin_rec.mass[rec1], kin_rec.mass[rec2]);

        results.phis_gen = CalculatePhiStar_Pair<float>(kin_gen.eta[gen1], kin_gen.eta[gen2], kin_gen.phi[gen1], kin_gen.phi[gen2]);
        results.phis_rec = CalculatePhiStar_Pair<float>(kin_rec.eta[rec1], kin_rec.eta[rec2], kin_rec.phi[rec1], kin_rec.phi[rec2]);
    
    } else if ((pass_gen) && (!pass_reco)) {
        results.miss = true;
        results.mll_gen = m_gen;

        results.pt_gen = CalculatePtZ0_Raw_Pair<float>(kin_gen.pt[gen1], kin_gen.pt[gen2], kin_gen.phi[gen1], kin_gen.phi[gen2]);

        results.y_gen = CalculateRapidityZ0_Raw_Pair<float>(kin_gen.pt[gen1], kin_gen.pt[gen2], kin_gen.eta[gen1], kin_gen.eta[gen2],
            kin_gen.phi[gen1], kin_gen.phi[gen2],  kin_gen.mass[gen1], kin_gen.mass[gen2]);

        results.phis_gen = CalculatePhiStar_Pair<float>(kin_gen.eta[gen1], kin_gen.eta[gen2], kin_gen.phi[gen1], kin_gen.phi[gen2]);

    } else if (pass_reco && !matched) {
        results.fake = true;
        results.mll_rec = m_rec;

        results.pt_rec = CalculatePtZ0_Raw_Pair<float>(kin_rec.pt[rec1], kin_rec.pt[rec2], kin_rec.phi[rec1], kin_rec.phi[rec2]);

        results.y_rec = CalculateRapidityZ0_Raw_Pair<float>(kin_rec.pt[rec1], kin_rec.pt[rec2], kin_rec.eta[rec1], kin_rec.eta[rec2],
            kin_rec.phi[rec1], kin_rec.phi[rec2],  kin_rec.mass[rec1], kin_rec.mass[rec2]);

        results.phis_rec = CalculatePhiStar_Pair<float>(kin_rec.eta[rec1], kin_rec.eta[rec2], kin_rec.phi[rec1], kin_rec.phi[rec2]);
    }

    return results;
}

ROOT::RDF::RNode CalculateRespMatrixWrapper(ROOT::RDF::RNode node, const flags_config& flags_RM, const cuts_config& cuts_RM) {
    ROOT::RDF::RNode node_RM = node;

    node_RM = node_RM
        .Define("RespMatrix_mask", [flags_RM, cuts_RM](const ROOT::RVec<float>& pt_rec, const ROOT::RVec<float>& eta_rec, const ROOT::RVec<float>& phi_rec, const ROOT::RVec<float>& mass_rec, const ROOT::RVec<int>& charge_rec,
        const ROOT::RVec<float>& pt_gen, const ROOT::RVec<float>& eta_gen, const ROOT::RVec<float>& phi_gen, const ROOT::RVec<float>& mass_gen,
        const ROOT::RVec<bool>& reco_tight, const ROOT::RVec<int>& reco_idx_gen, const ROOT::RVec<int>& gen_pdg_idx, const ROOT::RVec<int>& gen_status, const ROOT::RVec<int>& gen_status_flg) {
            
            MuonKinematics_REC kin_rec{pt_rec, eta_rec, phi_rec, mass_rec, charge_rec};
            MuonKinematics_GEN kin_gen{pt_gen, eta_gen, phi_gen, mass_gen};
            MuonFlags_RM val{reco_tight, reco_idx_gen, gen_pdg_idx, gen_status, gen_status_flg};
            
            return CalculateRespMatrix(kin_rec, kin_gen, val, flags_RM, cuts_RM);
        
        }, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge", "GenPart_pt", "GenPart_eta", "GenPart_phi", 
            "GenPart_mass", "Muon_tightId", "Muon_genPartIdx", "GenPart_pdgId", "GenPart_status", "GenPart_statusFlags"});
    
    node_RM = node_RM
        .Define("Matched", [](const ResultsRespMatrix& resp){ return resp.match; }, {"RespMatrix_mask"})
        .Define("Missed", [](const ResultsRespMatrix& resp){ return resp.miss; }, {"RespMatrix_mask"})
        .Define("Faked", [](const ResultsRespMatrix& resp){ return resp.fake; }, {"RespMatrix_mask"})
        .Filter("Matched || Missed || Faked"); 

    node_RM = node_RM
        .Define("Rec_InvMass", [](const ResultsRespMatrix& res) { return res.mll_rec; }, {"RespMatrix_mask"})
        .Define("Gen_InvMass", [](const ResultsRespMatrix& res) { return res.mll_gen; }, {"RespMatrix_mask"})
            
        .Define("Rec_Pt", [](const ResultsRespMatrix& res) { return res.pt_rec; }, {"RespMatrix_mask"})
        .Define("Gen_Pt", [](const ResultsRespMatrix& res) { return res.pt_gen; }, {"RespMatrix_mask"})
        
        .Define("Rec_Y", [](const ResultsRespMatrix& res) { return res.y_rec; }, {"RespMatrix_mask"})
        .Define("Gen_Y", [](const ResultsRespMatrix& res) { return res.y_gen; }, {"RespMatrix_mask"})

        .Define("Rec_Phis", [](const ResultsRespMatrix& res) { return res.phis_rec; }, {"RespMatrix_mask"})
        .Define("Gen_Phis", [](const ResultsRespMatrix& res) { return res.phis_gen; }, {"RespMatrix_mask"});

    node_RM = node_RM
        .Define("weight", [](){ return 1.0; }, {})
        .Define("Rec_Pt_Unf", [](bool miss, float pt_rec){ return miss ? -1.0 : pt_rec; }, {"Missed","Rec_Pt"})
        .Define("Gen_Pt_Unf", [](bool fake, float pt_gen){ return fake ? -1.0 : pt_gen; }, {"Faked","Gen_Pt"})

        .Define("Rec_Y_Unf",  [](bool miss, float y_rec){ return miss ? -100.0 : y_rec; }, {"Missed","Rec_Y"})
        .Define("Gen_Y_Unf",  [](bool fake, float y_gen){ return fake ? -100.0 : y_gen; }, {"Faked","Gen_Y"})

        .Define("Rec_Phis_Unf",[](bool miss, float phis_rec){ return miss ? -1.0 : phis_rec; }, {"Missed","Rec_Phis"})
        .Define("Gen_Phis_Unf",[](bool fake, float phis_gen){ return fake ? -1.0 : phis_gen; }, {"Faked","Gen_Phis"});
    
    return node_RM;
}


RespMatrixHisto BuildRespMatrixHisto(ROOT::RDF::RNode node) {
    RespMatrixHisto histo;
    
    ROOT::RDF::RNode node_unf = node;

    int reco_bins_pt = 30;
    int gen_bins_pt = 30;
    
    int reco_bins_y = 30;
    int gen_bins_y = 30;
    
    int reco_bins_phis = 30;
    int gen_bins_phis = 30;

    auto hResp_pt = node_unf.Histo2D({"hResp_pt","", reco_bins_pt, 0.0, 100.0, gen_bins_pt, 0.0, 100.0}, "Rec_Pt_Unf", "Gen_Pt_Unf", "weight");

    auto hResp_y = node_unf.Histo2D({"hResp_y","", reco_bins_y, -2.5, 2.5, gen_bins_y, -2.5, 2.5}, "Rec_Y_Unf", "Gen_Y_Unf", "weight");

    auto hResp_phis = node_unf.Histo2D({"hResp_phis","", reco_bins_phis, 0.0, 3.0, gen_bins_phis, 0.0, 3.0}, "Rec_Phis_Unf", 
    "Gen_Phis_Unf", "weight");

    histo.histo_pt = hResp_pt;
    histo.histo_y = hResp_y;
    histo.histo_phis = hResp_phis;

    return histo;
}

UnfoldDensities CreateUnfoldDensity(RespMatrixHisto& histo) {
    UnfoldDensities densities;

    densities.pt_unf = std::make_unique<TUnfoldDensity>(histo.histo_pt.GetPtr(), TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature, TUnfold::kEConstraintNone,
        TUnfoldDensity::kDensityModeBinWidth);

    densities.y_unf = std::make_unique<TUnfoldDensity>(histo.histo_y.GetPtr(), TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature, TUnfold::kEConstraintNone,
        TUnfoldDensity::kDensityModeBinWidth);

    densities.phis_unf = std::make_unique<TUnfoldDensity>(histo.histo_phis.GetPtr(), TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature, TUnfold::kEConstraintNone,
        TUnfoldDensity::kDensityModeBinWidth);

    return densities;

}
