#include "Unfold.h"
#include "Utils.h"

#include <vector>

#include <TGraph.h>
#include <TSpline.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMarker.h>

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

        .Define("Rec_Y_Unf", [](bool miss, float y_rec){ return miss ? -100.0 : y_rec; }, {"Missed","Rec_Y"})
        .Define("Gen_Y_Unf", [](bool fake, float y_gen){ return fake ? -100.0 : y_gen; }, {"Faked","Gen_Y"})

        .Define("Rec_Phis_Unf",[](bool miss, float phis_rec){ return miss ? -1.0 : phis_rec; }, {"Missed","Rec_Phis"})
        .Define("Gen_Phis_Unf",[](bool fake, float phis_gen){ return fake ? -1.0 : phis_gen; }, {"Faked","Gen_Phis"});
    
    return node_RM;
}


RespMatrixHisto BuildRespMatrixHisto(ROOT::RDF::RNode node, const config_struct& cfg) {
    
    ROOT::RDF::RNode node_unf = node;
    RespMatrixHisto resp_histo;

    std::vector<double> reco_bins_pt, gen_bins_pt, reco_bins_y, gen_bins_y, reco_bins_phis, gen_bins_phis;

    if (cfg.unfold.use_bins == true) {

        const auto& pt = cfg.unfold.pt_bins;
        const auto& y = cfg.unfold.y_bins;
        const auto& phis = cfg.unfold.phis_bins;

        reco_bins_pt = CreateBins(pt.reco_bins, pt.min, pt.max, pt.distribution, pt.split);
        gen_bins_pt = CreateBins(pt.gen_bins, pt.min, pt.max, pt.distribution, pt.split);
            
        reco_bins_y = CreateBins(y.reco_bins, y.min, y.max, y.distribution, y.split);
        gen_bins_y = CreateBins(y.gen_bins, y.min, y.max, y.distribution, y.split);
            
        reco_bins_phis = CreateBins(phis.reco_bins, phis.min, phis.max, phis.distribution, phis.split);
        gen_bins_phis = CreateBins(phis.gen_bins, phis.min, phis.max, phis.distribution, phis.split);
        
        // ------------------
        // Matched Histograms
        // ------------------

        auto node_matched = node_unf.Filter("Matched");
        resp_histo.h1_pt_test = node_matched.Histo1D({"h_rec_pt_matched", "", pt.reco_bins, reco_bins_pt.data()}, "Rec_Pt");
        resp_histo.h1_y_test = node_matched.Histo1D({"h_rec_y_matched", "", y.reco_bins, reco_bins_y.data()}, "Rec_Y");
        resp_histo.h1_phis_test = node_matched.Histo1D({"h_rec_phis_matched", "", phis.reco_bins, reco_bins_phis.data()}, "Rec_Phis");

        // --------------------
        // Faked BKG Histograms
        // --------------------

        auto node_faked = node_unf.Filter("Faked");
        resp_histo.h1_pt_fake = node_faked.Histo1D({"h_rec_pt_faked", "", pt.reco_bins, reco_bins_pt.data()}, "Rec_Pt");
        resp_histo.h1_y_fake = node_faked.Histo1D({"h_rec_y_faked", "", y.reco_bins, reco_bins_y.data()}, "Rec_Y");
        resp_histo.h1_phis_fake = node_faked.Histo1D({"h_rec_phis_faked", "", phis.reco_bins, reco_bins_phis.data()}, "Rec_Phis");

        // ------------------
        // Response Histogram
        // ------------------

        resp_histo.h2_pt = node_unf.Histo2D({"hResp_pt","", pt.reco_bins, reco_bins_pt.data(), pt.gen_bins , gen_bins_pt.data()}, "Rec_Pt_Unf",
        "Gen_Pt_Unf", "weight");

        resp_histo.h2_y = node_unf.Histo2D({"hResp_y","", y.reco_bins, reco_bins_y.data(), y.gen_bins, gen_bins_y.data()}, "Rec_Y_Unf", 
        "Gen_Y_Unf", "weight");

        resp_histo.h2_phis = node_unf.Histo2D({"hResp_phis","", phis.reco_bins, reco_bins_phis.data(), phis.gen_bins, gen_bins_phis.data()}, 
        "Rec_Phis_Unf", "Gen_Phis_Unf", "weight");

    } else {
        reco_bins_pt = cfg.unfold.pt_bins.reco_vec;
        gen_bins_pt = cfg.unfold.pt_bins.gen_vec;
        
        reco_bins_y = cfg.unfold.y_bins.reco_vec;
        gen_bins_y = cfg.unfold.y_bins.gen_vec;
        
        reco_bins_phis = cfg.unfold.phis_bins.reco_vec;
        gen_bins_phis = cfg.unfold.phis_bins.gen_vec;

        int n_reco_pt = static_cast<int>(reco_bins_pt.size()) - 1;
        int n_gen_pt  = static_cast<int>(gen_bins_pt.size()) - 1;

        int n_reco_y = static_cast<int>(reco_bins_y.size()) - 1;
        int n_gen_y  = static_cast<int>(gen_bins_y.size()) - 1;

        int n_reco_phis = static_cast<int>(reco_bins_phis.size()) - 1;
        int n_gen_phis  = static_cast<int>(gen_bins_phis.size()) - 1;

        // ------------------
        // Matched Histograms
        // ------------------

        auto node_matched = node_unf.Filter("Matched");
        resp_histo.h1_pt_test = node_matched.Histo1D({"h_rec_pt_matched", "", n_reco_pt, reco_bins_pt.data()}, "Rec_Pt");
        resp_histo.h1_y_test = node_matched.Histo1D({"h_rec_y_matched", "", n_reco_y, reco_bins_y.data()}, "Rec_Y");
        resp_histo.h1_phis_test = node_matched.Histo1D({"h_rec_phis_matched", "", n_reco_phis, reco_bins_phis.data()}, "Rec_Phis");

        // --------------------
        // Faked BKG Histograms
        // --------------------

        auto node_faked = node_unf.Filter("Faked");
        resp_histo.h1_pt_fake = node_faked.Histo1D({"h_rec_pt_faked", "", n_reco_pt, reco_bins_pt.data()}, "Rec_Pt");
        resp_histo.h1_y_fake = node_faked.Histo1D({"h_rec_y_faked", "", n_reco_y, reco_bins_y.data()}, "Rec_Y");
        resp_histo.h1_phis_fake = node_faked.Histo1D({"h_rec_phis_faked", "", n_reco_phis, reco_bins_phis.data()}, "Rec_Phis");

        // ------------------
        // Response Histogram
        // ------------------

        resp_histo.h2_pt = node_unf.Histo2D({"hResp_pt","", n_reco_pt, reco_bins_pt.data(), n_gen_pt, gen_bins_pt.data()}, "Rec_Pt_Unf",
        "Gen_Pt_Unf", "weight");

        resp_histo.h2_y = node_unf.Histo2D({"hResp_y","", n_reco_y, reco_bins_y.data(), n_gen_y, gen_bins_y.data()}, "Rec_Y_Unf", "Gen_Y_Unf",
        "weight");

        resp_histo.h2_phis = node_unf.Histo2D({"hResp_phis","", n_reco_phis, reco_bins_phis.data(), n_gen_phis, gen_bins_phis.data()}, 
        "Rec_Phis_Unf", "Gen_Phis_Unf", "weight");
    }

    return resp_histo;
}

EffPurHisto BuildEffPurHisto(ROOT::RDF::RNode node, const config_struct& cfg) {
    ROOT::RDF::RNode node_histo = node;
    EffPurHisto eff_pur_histo;

    std::vector<double> reco_bins_pt, gen_bins_pt, reco_bins_y, gen_bins_y, reco_bins_phis, gen_bins_phis;

    const auto& pt = cfg.unfold.pt_bins;
    const auto& y = cfg.unfold.y_bins;
    const auto& phis = cfg.unfold.phis_bins;

    if (cfg.unfold.use_bins) {

        reco_bins_pt = CreateBins(pt.reco_bins, pt.min, pt.max, pt.distribution, pt.split);
        gen_bins_pt = CreateBins(pt.gen_bins, pt.min, pt.max, pt.distribution, pt.split);
            
        reco_bins_y = CreateBins(y.reco_bins, y.min, y.max, y.distribution, y.split);
        gen_bins_y = CreateBins(y.gen_bins, y.min, y.max, y.distribution, y.split);
            
        reco_bins_phis = CreateBins(phis.reco_bins, phis.min, phis.max, phis.distribution, phis.split);
        gen_bins_phis = CreateBins(phis.gen_bins, phis.min, phis.max, phis.distribution, phis.split);
    
    } else {
        reco_bins_pt = pt.reco_vec;
        gen_bins_pt = pt.gen_vec;
        
        reco_bins_y = y.reco_vec;
        gen_bins_y = y.gen_vec;
        
        reco_bins_phis = phis.reco_vec;
        gen_bins_phis = phis.gen_vec;

    }

    int n_reco_pt = static_cast<int>(reco_bins_pt.size()) - 1;
    int n_gen_pt = static_cast<int>(gen_bins_pt.size()) - 1;
    
    int n_reco_y = static_cast<int>(reco_bins_y.size()) - 1;
    int n_gen_y = static_cast<int>(gen_bins_y.size()) - 1;
    
    int n_reco_phis = static_cast<int>(reco_bins_phis.size()) - 1;
    int n_gen_phis = static_cast<int>(gen_bins_phis.size()) - 1;

    auto FindBinIndex = [](float value, const std::vector<double>& bins) -> int {
        for (int i = 0; i < bins.size() - 1; ++i) {
            if ((value > bins[i]) && (value <= bins[i + 1])) {
                
                return static_cast<int>(i);
            }
        }
        
        return -1;
    };

    node_histo = node_histo
        .Define("Purity_Pt_Pass", [FindBinIndex, reco_bins_pt](bool matched, float gen, float reco) {
            if (!matched) {
                return false;
            }

            int bin_gen = FindBinIndex(gen, reco_bins_pt);
            int bin_reco = FindBinIndex(reco, reco_bins_pt);
            
            return (bin_reco != -1) && (bin_reco == bin_gen);
        }, {"Matched", "Gen_Pt", "Rec_Pt"});

    node_histo = node_histo
        .Define("Purity_Y_Pass", [FindBinIndex, reco_bins_y](bool matched, float gen, float reco) {
            if (!matched) {
                return false;
            }

            int bin_gen = FindBinIndex(gen, reco_bins_y);
            int bin_reco = FindBinIndex(reco, reco_bins_y);
            
            return (bin_reco != -1) && (bin_reco == bin_gen);
        }, {"Matched", "Gen_Y", "Rec_Y"});

    node_histo = node_histo
        .Define("Purity_Phis_Pass", [FindBinIndex, reco_bins_phis](bool matched, float gen, float reco) {
            if (!matched) {
                return false;
            }

            int bin_gen = FindBinIndex(gen, reco_bins_phis);
            int bin_reco = FindBinIndex(reco, reco_bins_phis);
            
            return (bin_reco != -1) && (bin_reco == bin_gen);
        }, {"Matched", "Gen_Phis", "Rec_Phis"});

    // -----
    // Pt_Z0
    // -----

    // Efficiency
    auto h_eff_den_pt = node_histo
        .Filter("Matched || Missed")
        .Histo1D({"h_eff_den_pt", "", n_gen_pt, gen_bins_pt.data()}, "Gen_Pt", "weight");
    
    auto h_eff_num_pt = node_histo
        .Filter("Matched")
        .Histo1D({"h_eff_num_pt", "", n_gen_pt, gen_bins_pt.data()}, "Gen_Pt", "weight");

    // Purity
    auto h_pur_den_pt = node_histo
        .Filter("Matched || Faked")
        .Histo1D({"h_pur_den_pt", "", n_reco_pt, reco_bins_pt.data()}, "Rec_Pt", "weight");
    auto h_pur_num_pt = node_histo
        .Filter("Matched && Purity_Pt_Pass")
        .Histo1D({"h_pur_num_pt", "", n_reco_pt, reco_bins_pt.data()}, "Rec_Pt", "weight");

    // ----
    // Y_Z0
    // ----

    // Efficiency
    auto h_eff_den_y = node_histo
        .Filter("Matched || Missed")
        .Histo1D({"h_eff_den_y", "", n_gen_y, gen_bins_y.data()}, "Gen_Y", "weight");
    auto h_eff_num_y = node_histo
        .Filter("Matched")
        .Histo1D({"h_eff_num_y", "", n_gen_y, gen_bins_y.data()}, "Gen_Y", "weight");
    
    // Purity
    auto h_pur_den_y = node_histo
        .Filter("Matched || Faked")
        .Histo1D({"h_pur_den_y", "", n_reco_y, reco_bins_y.data()}, "Rec_Y", "weight");
    auto h_pur_num_y = node_histo
        .Filter("Matched && Purity_Y_Pass")
        .Histo1D({"h_pur_num_y", "", n_reco_y, reco_bins_y.data()}, "Rec_Y", "weight");

    // -------
    // Phis_Z0
    // -------

    // Efficiency
    auto h_eff_den_phis = node_histo
        .Filter("Matched || Missed")
        .Histo1D({"h_eff_den_phis", "", n_gen_phis, gen_bins_phis.data()}, "Gen_Phis", "weight");
    auto h_eff_num_phis = node_histo
        .Filter("Matched")
        .Histo1D({"h_eff_num_phis", "", n_gen_phis, gen_bins_phis.data()}, "Gen_Phis", "weight");
    
    // Purity
    auto h_pur_den_phis = node_histo
        .Filter("Matched || Faked")
        .Histo1D({"h_pur_den_phis", "", n_reco_phis, reco_bins_phis.data()}, "Rec_Phis", "weight");
    auto h_pur_num_phis = node_histo
        .Filter("Matched && Purity_Phis_Pass")
        .Histo1D({"h_pur_num_phis", "", n_reco_phis, reco_bins_phis.data()}, "Rec_Phis", "weight");


    auto CreateHistogrm = [](ROOT::RDF::RResultPtr<TH1D>& num, ROOT::RDF::RResultPtr<TH1D>& den, const std::string& name, const std::string& title) {
        auto h_result = std::unique_ptr<TH1D>(static_cast<TH1D*>(num->Clone(name.c_str())));
        
        h_result->SetDirectory(nullptr);
        h_result->SetTitle(title.c_str());
        
        // Event Loop
        h_result->Divide(num.GetPtr(), den.GetPtr(), 1.0, 1.0, "B");
        
        return h_result;
    };

    eff_pur_histo.h1_Eff_pt = CreateHistogrm(h_eff_num_pt, h_eff_den_pt, "h1_Eff_pt", "Efficiency; Gen_Pt_Z0 [GeV]; Efficiency");
    eff_pur_histo.h1_Pur_pt = CreateHistogrm(h_pur_num_pt,h_pur_den_pt, "h1_Pur_pt", "Purity; Rec_Pt_Z0 [GeV]; Purity");

    eff_pur_histo.h1_Eff_y = CreateHistogrm(h_eff_num_y, h_eff_den_y, "h1_Eff_y", "Efficiency; Gen_Y_Z0; Efficiency");
    eff_pur_histo.h1_Pur_y = CreateHistogrm(h_pur_num_y, h_pur_den_y, "h1_Pur_y", "Purity; Rec_Y_Z0; Purity");

    eff_pur_histo.h1_Eff_phis = CreateHistogrm(h_eff_num_phis, h_eff_den_phis, "h1_Eff_phis", "Efficiency; Gen_Phis_Z0; Efficiency");
    eff_pur_histo.h1_Pur_phis = CreateHistogrm(h_pur_num_phis, h_pur_den_phis, "h1_Pur_phis", "Purity; Rec_Phis_Z0; Purity");


    return eff_pur_histo;
}


UnfoldDensities CreateUnfoldDensity(RespMatrixHisto& resp_histo) {
    UnfoldDensities densities;

    densities.pt_unf = std::make_unique<TUnfoldDensity>(resp_histo.h2_pt.GetPtr(), TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature, TUnfold::kEConstraintNone,
        TUnfoldDensity::kDensityModeBinWidth);

    densities.y_unf = std::make_unique<TUnfoldDensity>(resp_histo.h2_y.GetPtr(), TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature, TUnfold::kEConstraintNone,
        TUnfoldDensity::kDensityModeBinWidth);

    densities.phis_unf = std::make_unique<TUnfoldDensity>(resp_histo.h2_phis.GetPtr(), TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature, TUnfold::kEConstraintNone,
        TUnfoldDensity::kDensityModeBinWidth);

    return densities;
}

UnfoldResult ApplyUnfold(std::unique_ptr<TUnfoldDensity> density, TH1D* event_histo, TH1D* resp_histo, TH1D* fake_histo, const config_struct& cfg, const std::string& tag) {
    UnfoldResult results;

    results.unf_density = std::move(density); 

    // Starting the Second Event Loop on DATA !
    results.unf_density->SetInput(resp_histo);
    
    // Subtracting fake background
    //results.unf_density->SubtractBackground(fake_histo, "Fake signal", 1.0, 0.05);
    
    TGraph *lc = nullptr;
    TSpline *sx = nullptr;
    TSpline *sy = nullptr;
    //TSpline *scanResult = nullptr;


    // Options to Scan are: ScanSURE(), ScanLCurve(), ScanTau()
    results.idx_best = results.unf_density->ScanLcurve(cfg.unfold.scan.n_iter, cfg.unfold.scan.tau_min, cfg.unfold.scan.tau_max, &lc, &sx, &sy);// CORE

    results.LCurveScan.reset(lc);
    
    results.logTauX.reset(sx);
    results.logTauY.reset(sy);

    results.tau = results.unf_density->GetTau();
    
    std::string name_h1_unf = "hUnfolded_" + tag;
    std::string name_h2_cov = "hCov_" + tag;
    std::string name_h2_corr = "hCorr_" + tag;

    results.h1_out_unf.reset(results.unf_density->GetOutput(name_h1_unf.c_str()));
    results.h2_out_cov.reset(results.unf_density->GetEmatrixTotal(name_h2_cov.c_str()));
    results.h2_out_corr.reset(results.unf_density->GetRhoIJtotal(name_h2_corr.c_str()));

    results.chi2A = results.unf_density->GetChi2A();
    results.chi2L = results.unf_density->GetChi2L();
    results.ndf = results.unf_density->GetNdf();

    return results;
}


int VisualizeUnfoldResults(std::vector<std::unique_ptr<TCanvas>>& canvas, UnfoldResult& results, RespMatrixHisto& resp_histo, const std::string& tag) {
    std::cout << "Initializing visualization " << tag << " sample." << std::endl;
    
    std::cout << "tau scelto = " << results.tau << " (indice " << results.idx_best << ")\n";
    std::cout << "chi2A = " << results.chi2A << "  chi2L = " << results.chi2L << "  ndf = "   << results.ndf << "\n";

    // --------------------------
    // Canvas 1 - Response Matrix
    // --------------------------

    std::string name_c1 = "Response matrix " + tag;
    auto c1 = std::make_unique<TCanvas>(name_c1.c_str(), name_c1.c_str(), 800, 600);
    
    TH2D* histo_resp;
    if (tag == "pt") {
        histo_resp = resp_histo.h2_pt.GetPtr();
    } else if (tag == "y") {
        histo_resp = resp_histo.h2_y.GetPtr();
    } else if (tag == "phis") {
        histo_resp = resp_histo.h2_phis.GetPtr();
    }

    if (histo_resp == nullptr) {
        std::cout << "ERROR: invalid Response Matrix histogram, exiting..." << std::endl;
        return 1;
    }

    //c1->SetLogz();
    histo_resp->Draw("COLZ");

    std::string title_reco = tag + " reco [GeV]";
    std::string title_gen = tag + " gen [GeV]";
    histo_resp->GetXaxis()->SetTitle(title_reco.c_str());
    histo_resp->GetYaxis()->SetTitle(title_gen.c_str());
    
    // ----------------
    // Canvas 2 - LScan
    // ----------------
    
    std::string name_c2 = "L-curve " + tag;
    auto c2 = std::make_unique<TCanvas>(name_c2.c_str(), name_c2.c_str(), 800, 600);
    
    results.LCurveScan->SetTitle("L-curve;log_{10}(#chi^{2}_{data});log_{10}(curvature)");
    results.LCurveScan->Draw("ALP");

    double xBest;
    double yBest;
    
    results.LCurveScan->GetPoint(results.idx_best, xBest, yBest);
    
    TMarker *mBest = new TMarker(xBest, yBest, 20);
    mBest->SetMarkerColor(kRed);
    mBest->SetMarkerSize(1.5);
    mBest->Draw("SAME");

    // ------------------------
    // Canvas 3 - Unfold Output
    // ------------------------

    TH1D* h_truth_GEN = static_cast<TH1D*>(histo_resp->ProjectionY(("h_truth_GEN_" + tag).c_str()));
    h_truth_GEN->SetDirectory(nullptr);

    std::cout << "Binning MC:" << h_truth_GEN->GetNbinsX() << std::endl;
    std::cout << "Binning Unfolded Truth:" << results.h1_out_unf->GetNbinsX() << std::endl;

    std::string name_c3 = "Unfolded vs Truth MC " + tag;
    auto c3 = std::make_unique<TCanvas>(name_c3.c_str(), name_c3.c_str(), 800, 600);
    results.h1_out_unf->SetLineColor(kRed);
    results.h1_out_unf->SetMarkerColor(kRed);
    results.h1_out_unf->SetMarkerStyle(20);
    results.h1_out_unf->Scale(1,"width");
    results.h1_out_unf->Draw("E");
    results.h1_out_unf->SetStats(0);

    h_truth_GEN->SetLineColor(kBlue);
    h_truth_GEN->SetLineWidth(2);
    h_truth_GEN->Scale(1,"width");
    h_truth_GEN->Draw("HIST SAME");
    h_truth_GEN->SetStats(0);

    TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(results.h1_out_unf.get(), "Unfolded", "lep");
    leg->AddEntry(h_truth_GEN, "MC Truth (gen)", "l");
    leg->Draw();
    
    // ----------------------------
    // Canvas 4 - Covariance Matrix
    // ----------------------------

    std::string name_c4 = "Covariance matrix " + tag;
    auto c4 = std::make_unique<TCanvas>(name_c4.c_str(), name_c4.c_str(), 800, 600);
    results.h2_out_cov->Draw("COLZ");

    // ----------------------------
    // Canvas 6 - Correlation Matrix
    // ----------------------------

    std::string name_c6 = "Correlation matrix " + tag;
    auto c6 = std::make_unique<TCanvas>(name_c6.c_str(), name_c6.c_str(), 800, 600);
    results.h2_out_corr->Draw("COLZ");

    canvas.push_back(std::move(c1));
    canvas.push_back(std::move(c2));
    canvas.push_back(std::move(c3));
    canvas.push_back(std::move(c4));
    canvas.push_back(std::move(c6));

    return 0;
}


int VisualizeControlPlots(std::vector<std::unique_ptr<TCanvas>>& canvas, RespMatrixHisto& resp_histo, EffPurHisto& eff_histo, const std::string& tag) {
    // --------------------------
    // Canvas 1 - Response Matrix
    // --------------------------

    std::string name_c1 = "Response matrix " + tag;
    auto c1 = std::make_unique<TCanvas>(name_c1.c_str(), name_c1.c_str(), 800, 600);
    
    TH2D* histo_resp;
    if (tag == "pt") {
        histo_resp = resp_histo.h2_pt.GetPtr();
    } else if (tag == "y") {
        histo_resp = resp_histo.h2_y.GetPtr();
    } else if (tag == "phis") {
        histo_resp = resp_histo.h2_phis.GetPtr();
    }

    if (histo_resp == nullptr) {
        std::cout << "ERROR: invalid Response Matrix histogram, exiting..." << std::endl;
        return 1;
    }

    //c1->SetLogz();
    histo_resp->Draw("COLZ");

    std::string title_reco = tag + " reco [GeV]";
    std::string title_gen = tag + " gen [GeV]";
    histo_resp->GetXaxis()->SetTitle(title_reco.c_str());
    histo_resp->GetYaxis()->SetTitle(title_gen.c_str());
    
    std::string name_c5 = "Efficiency/Purity_" + tag; 
    auto c5 = std::make_unique<TCanvas>(name_c5.c_str(), name_c5.c_str(), 800, 600);
    TH1D* h_eff;
    TH1D* h_pur;

    if (tag == "pt") {
        h_eff = eff_histo.h1_Eff_pt.get();
        h_pur = eff_histo.h1_Pur_pt.get();
    } else if (tag == "y") {
        h_eff = eff_histo.h1_Eff_y.get();
        h_pur = eff_histo.h1_Pur_y.get();
    } else if (tag == "phis") {
        h_eff = eff_histo.h1_Eff_phis.get();
        h_pur = eff_histo.h1_Pur_phis.get();
    }

    if (h_eff && h_pur) {
        c5->cd();
        h_eff->SetLineColor(kRed);
        h_eff->SetLineWidth(2);
        h_eff->SetMinimum(0.0);
        h_eff->SetMaximum(1.15);
        h_eff->SetStats(0);

        h_pur->SetLineColor(kBlue);
        h_pur->SetLineWidth(2);
        h_pur->SetStats(0);

        h_eff->Draw("HIST");
        h_pur->Draw("HIST SAME");

        auto leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_eff, "Efficiency", "l");
        leg->AddEntry(h_pur, "Purity", "l");
        leg->Draw();

        c5->Update();
        canvas.push_back(std::move(c5));
    }

    canvas.push_back(std::move(c1));

    return 0;

}