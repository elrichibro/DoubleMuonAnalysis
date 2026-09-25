#include "Selection.h"

#include "Utils.h"

#include <vector>

// ------------------------------------------------------------------------------------------------------------------------------------

ResultsRespMatrix CalculateRespMatrix(const MuonKinematics_REC& kin_rec, const MuonKinematics_GEN& kin_gen, const MuonFlags_RM& flags, const flags_config cfg_f, 
const cuts_config cfg_c) {
    
    ResultsRespMatrix results;

    // Summing collection objects -> reconstructions flag -> fast exit if sum(tight) > 2
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

    // Defualt values
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

        // Fiducial region -> [ pt > JSON_cut && |eta| < JSON_cut ]
        bool pass_kin = (((kin_rec.pt[i] > cfg_c.pt_cut) && (std::abs(kin_rec.eta[i]) < cfg_c.eta_cut)) || (!cfg_f.en_kinematics));
        
        if (pass_kin) {
            reco_idx.push_back(i);
        }
    }

    // Two reconstructed muons
    if (reco_idx.size() == 2) {
        rec1 = reco_idx[0];
        rec2 = reco_idx[1];
        
        // Opposite charge condition
        if (kin_rec.charge[rec1] != kin_rec.charge[rec2]) {

            // Invariant mass of reconstructed pair
            m_rec = CalculateInvariantMass_Pair<float>(kin_rec.pt[rec1], kin_rec.pt[rec2], kin_rec.eta[rec1], kin_rec.eta[rec2], 
                kin_rec.phi[rec1], kin_rec.phi[rec2], kin_rec.mass[rec1], kin_rec.mass[rec2]);
            
            // Invariant mass fiducial region
            if ((m_rec > 60) && (m_rec < 120)) {
                pass_reco = true;
            }
        }
    }

    // Loop on generated muons
    for (int i = 0; i < n_muons_gen; i++) {
        if (abs(flags.gen_pdg_idx[i]) != 13 || (flags.gen_status[i] != 1)) {
            continue;
        }
    
        // Mask identification
        bool mask = ((flags.gen_status_flg[i] & (1 << 0)) > 0) && (((flags.gen_status_flg[i] & (1 << 7)) > 0) || ((flags.gen_status_flg[i] & (1 << 8)) > 0)) && 
        (((flags.gen_status_flg[i] & (1 << 13)) > 0) || ((flags.gen_status_flg[i] & (1 << 14)) > 0));

        if (!mask) {
            continue;
        }
        
        // Kinematical region (same as reconstructed)
        bool pass_kin = (((kin_gen.pt[i] > cfg_c.pt_cut) && (std::abs(kin_gen.eta[i]) < cfg_c.eta_cut)) || (!cfg_f.en_kinematics));
        
        if (pass_kin) {
            gen_idx.push_back(i);
        }
    }

    // Two generated muons
    if (gen_idx.size() == 2) {
        gen1 = gen_idx[0];
        gen2 = gen_idx[1];
        
        // PDG Id condition
        if (flags.gen_pdg_idx[gen1] != flags.gen_pdg_idx[gen2]) {
            
            // Invariant mass
            m_gen = CalculateInvariantMass_Pair<float>(kin_gen.pt[gen1], kin_gen.pt[gen2], kin_gen.eta[gen1], kin_gen.eta[gen2],
             kin_gen.phi[gen1], kin_gen.phi[gen2], kin_gen.mass[gen1], kin_gen.mass[gen2]);
            
            // Selection in fiducial region
            if ((m_gen > 60) && (m_gen < 120)) {
                pass_gen = true;
            }
        }
    }

    // Matched muons condition
    if (pass_reco && pass_gen) {
        bool match1 = ((flags.reco_idx_gen[rec1] == gen1) && (flags.reco_idx_gen[rec2] == gen2));
        bool match2 = ((flags.reco_idx_gen[rec1] == gen2) && (flags.reco_idx_gen[rec2] == gen1));
        
        matched = (match1 || match2);
    }


    if (matched) {
        // Flag
        results.match = true;
        
        // Invariant masses
        results.mll_rec = m_rec;
        results.mll_gen = m_gen;
        
        // Interesting quantities
        results.pt_gen = CalculatePtZ0_Raw_Pair<float>(kin_gen.pt[gen1], kin_gen.pt[gen2], kin_gen.phi[gen1], kin_gen.phi[gen2]);
        results.pt_rec = CalculatePtZ0_Raw_Pair<float>(kin_rec.pt[rec1], kin_rec.pt[rec2], kin_rec.phi[rec1], kin_rec.phi[rec2]);
        
        results.y_gen = CalculateRapidityZ0_Raw_Pair<float>(kin_gen.pt[gen1], kin_gen.pt[gen2], kin_gen.eta[gen1], kin_gen.eta[gen2],
            kin_gen.phi[gen1], kin_gen.phi[gen2],  kin_gen.mass[gen1], kin_gen.mass[gen2]);
        
        results.y_rec = CalculateRapidityZ0_Raw_Pair<float>(kin_rec.pt[rec1], kin_rec.pt[rec2], kin_rec.eta[rec1], kin_rec.eta[rec2],
            kin_rec.phi[rec1], kin_rec.phi[rec2],  kin_rec.mass[rec1], kin_rec.mass[rec2]);

        results.phis_gen = CalculatePhiStar_Pair<float>(kin_gen.eta[gen1], kin_gen.eta[gen2], kin_gen.phi[gen1], kin_gen.phi[gen2]);
        results.phis_rec = CalculatePhiStar_Pair<float>(kin_rec.eta[rec1], kin_rec.eta[rec2], kin_rec.phi[rec1], kin_rec.phi[rec2]);
    
    // Missed muons
    } else if ((pass_gen) && (!pass_reco)) {
        // Flag
        results.miss = true;

        // Invariant mass
        results.mll_gen = m_gen;

        // Only generated muons quantities
        results.pt_gen = CalculatePtZ0_Raw_Pair<float>(kin_gen.pt[gen1], kin_gen.pt[gen2], kin_gen.phi[gen1], kin_gen.phi[gen2]);

        results.y_gen = CalculateRapidityZ0_Raw_Pair<float>(kin_gen.pt[gen1], kin_gen.pt[gen2], kin_gen.eta[gen1], kin_gen.eta[gen2],
            kin_gen.phi[gen1], kin_gen.phi[gen2],  kin_gen.mass[gen1], kin_gen.mass[gen2]);

        results.phis_gen = CalculatePhiStar_Pair<float>(kin_gen.eta[gen1], kin_gen.eta[gen2], kin_gen.phi[gen1], kin_gen.phi[gen2]);

    // Faked muons
    } else if (pass_reco && !matched) {
        // Flag
        results.fake = true;
        
        // Reconstructed invariant mass
        results.mll_rec = m_rec;

        // Only reconstructed quantities
        results.pt_rec = CalculatePtZ0_Raw_Pair<float>(kin_rec.pt[rec1], kin_rec.pt[rec2], kin_rec.phi[rec1], kin_rec.phi[rec2]);

        results.y_rec = CalculateRapidityZ0_Raw_Pair<float>(kin_rec.pt[rec1], kin_rec.pt[rec2], kin_rec.eta[rec1], kin_rec.eta[rec2],
            kin_rec.phi[rec1], kin_rec.phi[rec2],  kin_rec.mass[rec1], kin_rec.mass[rec2]);

        results.phis_rec = CalculatePhiStar_Pair<float>(kin_rec.eta[rec1], kin_rec.eta[rec2], kin_rec.phi[rec1], kin_rec.phi[rec2]);
    }
    return results;
}

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RNode CalculateRespMatrixWrapper(ROOT::RDF::RNode node, const flags_config& flags_RM, const cuts_config& cuts_RM) {
    ROOT::RDF::RNode node_RM = node;

    node_RM = node_RM
        // Filling result struct
        .Define("RespMatrix_mask", [flags_RM, cuts_RM](const ROOT::RVec<float>& pt_rec, const ROOT::RVec<float>& eta_rec, const ROOT::RVec<float>& phi_rec, const ROOT::RVec<float>& mass_rec, const ROOT::RVec<int>& charge_rec,
        const ROOT::RVec<float>& pt_gen, const ROOT::RVec<float>& eta_gen, const ROOT::RVec<float>& phi_gen, const ROOT::RVec<float>& mass_gen,
        const ROOT::RVec<bool>& reco_tight, const ROOT::RVec<int>& reco_idx_gen, const ROOT::RVec<int>& gen_pdg_idx, const ROOT::RVec<int>& gen_status, const ROOT::RVec<int>& gen_status_flg) {
            
            MuonKinematics_REC kin_rec{pt_rec, eta_rec, phi_rec, mass_rec, charge_rec};
            MuonKinematics_GEN kin_gen{pt_gen, eta_gen, phi_gen, mass_gen};
            MuonFlags_RM val{reco_tight, reco_idx_gen, gen_pdg_idx, gen_status, gen_status_flg};
            
            return CalculateRespMatrix(kin_rec, kin_gen, val, flags_RM, cuts_RM);
        
        }, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge", "GenPart_pt", "GenPart_eta", "GenPart_phi", 
            "GenPart_mass", "Muon_tightId", "Muon_genPartIdx", "GenPart_pdgId", "GenPart_status", "GenPart_statusFlags"})
    
        // Defining flags
        .Define("Matched", [](const ResultsRespMatrix& resp){ return resp.match; }, {"RespMatrix_mask"})
        .Define("Missed", [](const ResultsRespMatrix& resp){ return resp.miss; }, {"RespMatrix_mask"})
        .Define("Faked", [](const ResultsRespMatrix& resp){ return resp.fake; }, {"RespMatrix_mask"})
        
        // Filtering events
        .Filter("Matched || Missed || Faked")

        // Defining response matrix quantities
        .Define("Rec_InvMass", [](const ResultsRespMatrix& res) { return res.mll_rec; }, {"RespMatrix_mask"})
        .Define("Gen_InvMass", [](const ResultsRespMatrix& res) { return res.mll_gen; }, {"RespMatrix_mask"})
            
        .Define("Rec_Pt", [](const ResultsRespMatrix& res) { return res.pt_rec; }, {"RespMatrix_mask"})
        .Define("Gen_Pt", [](const ResultsRespMatrix& res) { return res.pt_gen; }, {"RespMatrix_mask"})
        
        .Define("Rec_Y", [](const ResultsRespMatrix& res) { return res.y_rec; }, {"RespMatrix_mask"})
        .Define("Gen_Y", [](const ResultsRespMatrix& res) { return res.y_gen; }, {"RespMatrix_mask"})

        .Define("Rec_Phis", [](const ResultsRespMatrix& res) { return res.phis_rec; }, {"RespMatrix_mask"})
        .Define("Gen_Phis", [](const ResultsRespMatrix& res) { return res.phis_gen; }, {"RespMatrix_mask"});
    
    return node_RM;
}

// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RNode EventSelection(ROOT::RDF::RNode node, const config_struct& cfg) {
    ROOT::RDF::RNode node_event = node;

    // Defining Event cuts from JSON info
    std::string good_muon = "true";
    if (cfg.flag_ES.en_kinematics) {
        good_muon += " && (Muon_pt > " + std::to_string(cfg.cut_ES.pt_cut) + " && abs(Muon_eta) < " + std::to_string(cfg.cut_ES.eta_cut)
         + ")";
    }

    if (cfg.flag_ES.en_tight_muon) {
        good_muon += " && (Muon_tightId == true)";
    }

    std::string event_cut = "GoodEvent";
    if (cfg.flag_ES.en_mass_window) {
        event_cut += " && (InvariantMass > " + std::to_string(cfg.cut_ES.mass_min) + " && InvariantMass < " 
        + std::to_string(cfg.cut_ES.mass_max) + ")";
    }

    node_event = node_event
        // Defiing mask with kinematical cuts
        .Define("GoodMuon", good_muon)

        // Selecting variables in collection
        .Define("GM_Pt", "Muon_pt[GoodMuon]")
        .Define("GM_Eta", "Muon_eta[GoodMuon]")
        .Define("GM_Phi", "Muon_phi[GoodMuon]")
        .Define("GM_Mass", "Muon_mass[GoodMuon]")
        .Define("GM_Charge", "Muon_charge[GoodMuon]")

        .Define("EventPair", "Sum(GoodMuon) == 2")
        .Define("GoodEvent", "EventPair && (GM_Charge[0] != GM_Charge[1])")

        // Invariant mass booking
        .Define("InvariantMass", [] (const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi,
         const ROOT::RVec<float>& mass) -> float {
            
            return static_cast<float>(ROOT::VecOps::InvariantMass(pt, eta, phi, mass));
        
        }, {"GM_Pt", "GM_Eta", "GM_Phi", "GM_Mass"})

        // Invariant mass selection
        .Filter(event_cut, "InvMass selection -> Good Event")

        // Defining Z0 quantities of interest
        .Define("Pt_Z", [] (const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi,
         const ROOT::RVec<float>& mass) {
            
            float pt_Z = CalculatePtZ0<float>(pt, eta, phi, mass);
            return pt_Z;
        
        }, {"GM_Pt", "GM_Eta", "GM_Phi", "GM_Mass"})
        
        .Define("Y_Z", [] (const ROOT::RVec<float>& pt, const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi,
         const ROOT::RVec<float>& mass) {
            
            float y_Z = CalculateRapidityZ0<float>(pt, eta, phi, mass);
            return y_Z;
        
        }, {"GM_Pt", "GM_Eta", "GM_Phi", "GM_Mass"})
        
        .Define("Phis_Z", [] (const ROOT::RVec<float>& eta, const ROOT::RVec<float>& phi) {
            
            float phis_Z = CalculatePhiStar<float>(eta, phi);
            return phis_Z;
        
        }, {"GM_Eta", "GM_Phi"});

    return node_event;
}