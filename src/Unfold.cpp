#include "Unfold.h"

#include "Utils.h"

#include <TCanvas.h>
#include <TMarker.h>
#include <TLegend.h>

// ------------------------------------------------------------------------------------------------------------------------------------

RespMatrixHisto BuildRespMatrixHisto(ROOT::RDF::RNode node, const config_struct& cfg) {
    
    ROOT::RDF::RNode node_unf = node;
    
    // Initializing container
    RespMatrixHisto resp_histo;

    // Defining Unfold quantities with over/underflow bins
    node_unf = node_unf
        .Define("weight", [](){ return 1.0; }, {})// Standard weight -> for TagAndProbe implementation
        .Define("Rec_Pt_Unf", [](bool miss, float pt_rec){ return miss ? -1.0 : pt_rec; }, {"Missed","Rec_Pt"})
        .Define("Gen_Pt_Unf", [](bool fake, float pt_gen){ return fake ? -1.0 : pt_gen; }, {"Faked","Gen_Pt"})

        .Define("Rec_Y_Unf", [](bool miss, float y_rec){ return miss ? -100.0 : y_rec; }, {"Missed","Rec_Y"})
        .Define("Gen_Y_Unf", [](bool fake, float y_gen){ return fake ? -100.0 : y_gen; }, {"Faked","Gen_Y"})

        .Define("Rec_Phis_Unf",[](bool miss, float phis_rec){ return miss ? -1.0 : phis_rec; }, {"Missed","Rec_Phis"})
        .Define("Gen_Phis_Unf",[](bool fake, float phis_gen){ return fake ? -1.0 : phis_gen; }, {"Faked","Gen_Phis"});

    // Bin vectors initialization
    std::vector<double> reco_bins_pt, gen_bins_pt, reco_bins_y, gen_bins_y, reco_bins_phis, gen_bins_phis;

    const auto& pt = cfg.unfold.pt_bins;
    const auto& y = cfg.unfold.y_bins;
    const auto& phis = cfg.unfold.phis_bins;

    // Selection between custom/CreateBins bins -> JSON settup

    // -----
    // Pt_Z0
    // -----

    reco_bins_pt = (cfg.unfold.use_custom_bins) ? pt.reco_vec : CreateBins(pt.reco_bins, pt.min, pt.max, pt.distribution);
    gen_bins_pt = (cfg.unfold.use_custom_bins) ? pt.gen_vec : CreateBins(pt.gen_bins, pt.min, pt.max, pt.distribution);
    
    // ----
    // Y_Z0
    // ----

    reco_bins_y = (cfg.unfold.use_custom_bins) ? y.reco_vec : CreateBins(y.reco_bins, y.min, y.max, y.distribution);
    gen_bins_y = (cfg.unfold.use_custom_bins) ? y.gen_vec : CreateBins(y.gen_bins, y.min, y.max, y.distribution);

    // -------
    // Phi*_Z0
    // -------

    reco_bins_phis = (cfg.unfold.use_custom_bins) ? phis.reco_vec : CreateBins(phis.reco_bins, phis.min, phis.max, phis.distribution);
    gen_bins_phis = (cfg.unfold.use_custom_bins) ? phis.gen_vec : CreateBins(phis.gen_bins, phis.min, phis.max, phis.distribution);

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

    resp_histo.h2_y = node_unf.Histo2D({"hResp_y","", n_reco_y, reco_bins_y.data(), n_gen_y, gen_bins_y.data()}, "Rec_Y_Unf", 
    "Gen_Y_Unf", "weight");

    resp_histo.h2_phis = node_unf.Histo2D({"hResp_phis","", n_reco_phis, reco_bins_phis.data(), n_gen_phis, gen_bins_phis.data()}, 
    "Rec_Phis_Unf", "Gen_Phis_Unf", "weight");

    return resp_histo;
}

// ------------------------------------------------------------------------------------------------------------------------------------

ControlHisto BuildControlHisto(ROOT::RDF::RNode node, const config_struct& cfg) {
    ROOT::RDF::RNode node_histo = node;
    
    // Initializating container
    ControlHisto control_histo;

    node_histo = node_histo.Define("weight", [](){ return 1.0; }, {});


    // Bins initialization
    const auto& pt = cfg.unfold.pt_bins;
    const auto& y = cfg.unfold.y_bins;
    const auto& phis = cfg.unfold.phis_bins;

    std::vector<double> reco_bins_pt = (cfg.unfold.use_custom_bins) ? pt.reco_vec : CreateBins(pt.reco_bins, pt.min, pt.max, pt.distribution);
    std::vector<double> gen_bins_pt = (cfg.unfold.use_custom_bins) ? pt.gen_vec : CreateBins(pt.gen_bins, pt.min, pt.max, pt.distribution);
    
    std::vector<double> reco_bins_y = (cfg.unfold.use_custom_bins) ? y.reco_vec : CreateBins(y.reco_bins, y.min, y.max, y.distribution);
    std::vector<double> gen_bins_y = (cfg.unfold.use_custom_bins) ? y.gen_vec : CreateBins(y.gen_bins, y.min, y.max, y.distribution);
    
    std::vector<double> reco_bins_phis = (cfg.unfold.use_custom_bins) ? phis.reco_vec : CreateBins(phis.reco_bins, phis.min, phis.max, phis.distribution);
    std::vector<double> gen_bins_phis = (cfg.unfold.use_custom_bins) ? phis.gen_vec : CreateBins(phis.gen_bins, phis.min, phis.max, phis.distribution);

    int n_reco_pt = static_cast<int>(reco_bins_pt.size()) - 1;
    int n_gen_pt = static_cast<int>(gen_bins_pt.size()) - 1;
    
    int n_reco_y = static_cast<int>(reco_bins_y.size()) - 1;
    int n_gen_y = static_cast<int>(gen_bins_y.size()) - 1;
    
    int n_reco_phis = static_cast<int>(reco_bins_phis.size()) - 1;
    int n_gen_phis = static_cast<int>(gen_bins_phis.size()) - 1;


    // Find bin index in a bin vector -> used for purity/stability
    auto FindBinIndex = [](float value, const std::vector<double>& bins) -> int {
        for (int i = 0; i < bins.size() - 1; ++i) {
            if ((value > bins[i]) && (value <= bins[i + 1])) {
                
                return static_cast<int>(i);
            }
        }
        
        return -1;
    };

    // Lambda for Histogram creation
    auto CreateHistogrm = [](ROOT::RDF::RResultPtr<TH1D>& num, ROOT::RDF::RResultPtr<TH1D>& den, const std::string& name, const std::string& title) {
        auto h_result = std::unique_ptr<TH1D>(static_cast<TH1D*>(num->Clone(name.c_str())));
        
        h_result->SetDirectory(nullptr);
        h_result->SetTitle(title.c_str());
        
        // Event Loop
        h_result->Divide(num.GetPtr(), den.GetPtr(), 1.0, 1.0, "B");
        
        return h_result;
    };

    node_histo = node_histo
        .Define("Purity_Pt_Pass", [FindBinIndex, reco_bins_pt](bool matched, float gen, float reco) {
            if (!matched) {
                return false;
            }

            int bin_gen = FindBinIndex(gen, reco_bins_pt);
            int bin_reco = FindBinIndex(reco, reco_bins_pt);
            
            return (bin_reco != -1) && (bin_reco == bin_gen);
        }, {"Matched", "Gen_Pt", "Rec_Pt"})

        .Define("Purity_Y_Pass", [FindBinIndex, reco_bins_y](bool matched, float gen, float reco) {
            if (!matched) {
                return false;
            }

            int bin_gen = FindBinIndex(gen, reco_bins_y);
            int bin_reco = FindBinIndex(reco, reco_bins_y);
            
            return (bin_reco != -1) && (bin_reco == bin_gen);
        }, {"Matched", "Gen_Y", "Rec_Y"})

        .Define("Purity_Phis_Pass", [FindBinIndex, reco_bins_phis](bool matched, float gen, float reco) {
            if (!matched) {
                return false;
            }

            int bin_gen = FindBinIndex(gen, reco_bins_phis);
            int bin_reco = FindBinIndex(reco, reco_bins_phis);
            
            return (bin_reco != -1) && (bin_reco == bin_gen);
        }, {"Matched", "Gen_Phis", "Rec_Phis"})

        .Define("Stability_Pt_Pass", [FindBinIndex, gen_bins_pt](bool matched, float gen, float reco) {
            if (!matched) {
                return false;
            }

            int bin_gen = FindBinIndex(gen, gen_bins_pt);
            int bin_reco = FindBinIndex(reco, gen_bins_pt);
            
            return (bin_gen != -1) && (bin_reco == bin_gen);
        }, {"Matched", "Gen_Pt", "Rec_Pt"})

        .Define("Stability_Y_Pass", [FindBinIndex, gen_bins_y](bool matched, float gen, float reco) {
            if (!matched) {
                return false;
            }

            int bin_gen = FindBinIndex(gen, gen_bins_y);
            int bin_reco = FindBinIndex(reco, gen_bins_y);
            
            return (bin_gen != -1) && (bin_reco == bin_gen);
        }, {"Matched", "Gen_Y", "Rec_Y"})

        .Define("Stability_Phis_Pass", [FindBinIndex, gen_bins_phis](bool matched, float gen, float reco) {
            if (!matched) {
                return false;
            }

            int bin_gen = FindBinIndex(gen, gen_bins_phis);
            int bin_reco = FindBinIndex(reco, gen_bins_phis);
            
            return (bin_gen != -1) && (bin_reco == bin_gen);
        }, {"Matched", "Gen_Phis", "Rec_Phis"});


    // -------------------------------
    // Efficiency / Purity / Stability
    // -------------------------------
    
    // node_eff_num / node_stab_den
    ROOT::RDF::RNode node_matched = node_histo.Filter("Matched");
    auto h_matched_pt = node_matched.Histo1D({"h_matched_pt", "", n_gen_pt, gen_bins_pt.data()}, "Gen_Pt", "weight");
    auto h_matched_y = node_matched.Histo1D({"h_matched_y", "", n_gen_y, gen_bins_y.data()}, "Gen_Y", "weight");
    auto h_matched_phis = node_matched.Histo1D({"h_matched_phis", "", n_gen_phis, gen_bins_phis.data()}, "Gen_Phis", "weight");
    
    ROOT::RDF::RNode node_eff_den = node_histo.Filter("Matched || Missed");
    auto h_eff_den_pt = node_eff_den.Histo1D({"h_eff_den_pt", "", n_gen_pt, gen_bins_pt.data()}, "Gen_Pt", "weight");
    auto h_eff_den_y = node_eff_den.Histo1D({"h_eff_den_y", "", n_gen_y, gen_bins_y.data()}, "Gen_Y", "weight");
    auto h_eff_den_phis = node_eff_den.Histo1D({"h_eff_den_phis", "", n_gen_phis, gen_bins_phis.data()}, "Gen_Phis", "weight");

    ROOT::RDF::RNode node_pur_den = node_histo.Filter("Matched || Faked");
    auto h_pur_den_pt = node_pur_den.Histo1D({"h_pur_den_pt", "", n_reco_pt, reco_bins_pt.data()}, "Rec_Pt", "weight");
    auto h_pur_den_y = node_pur_den.Histo1D({"h_pur_den_y", "", n_reco_y, reco_bins_y.data()}, "Rec_Y", "weight");
    auto h_pur_den_phis = node_pur_den.Histo1D({"h_pur_den_phis", "", n_reco_phis, reco_bins_phis.data()}, "Rec_Phis", "weight");

    auto h_pur_num_pt = node_histo.Filter("Matched && Purity_Pt_Pass").Histo1D({"h_pur_num_pt", "", n_reco_pt, reco_bins_pt.data()}, "Rec_Pt", "weight");
    auto h_pur_num_y = node_histo.Filter("Matched && Purity_Y_Pass").Histo1D({"h_pur_num_y", "", n_reco_y, reco_bins_y.data()}, "Rec_Y", "weight");
    auto h_pur_num_phis = node_histo.Filter("Matched && Purity_Phis_Pass").Histo1D({"h_pur_num_phis", "", n_reco_phis, reco_bins_phis.data()}, "Rec_Phis", "weight");
    
    auto h_stab_num_pt = node_histo.Filter("Matched && Stability_Pt_Pass").Histo1D({"h_stab_num_pt", "", n_gen_pt, gen_bins_pt.data()}, "Gen_Pt", "weight");
    auto h_stab_num_y = node_histo.Filter("Matched && Stability_Y_Pass").Histo1D({"h_stab_num_y", "", n_gen_y, gen_bins_y.data()}, "Gen_Y", "weight");
    auto h_stab_num_phis = node_histo.Filter("Matched && Stability_Phis_Pass").Histo1D({"h_stab_num_phis", "", n_gen_phis, gen_bins_phis.data()}, "Gen_Phis", "weight");

    // Histogram creation
    control_histo.h1_Eff_pt = CreateHistogrm(h_matched_pt, h_eff_den_pt, "h1_Eff_pt", "Efficiency; Gen_Pt_Z0 [GeV]; Efficiency");
    control_histo.h1_Pur_pt = CreateHistogrm(h_pur_num_pt,h_pur_den_pt, "h1_Pur_pt", "Purity; Rec_Pt_Z0 [GeV]; Purity");
    control_histo.h1_Stab_pt = CreateHistogrm(h_stab_num_pt, h_matched_pt, "h1_Stab_pt", "Stability; Gen_Pt_Z0 [GeV]; Stability");

    control_histo.h1_Eff_y = CreateHistogrm(h_matched_y, h_eff_den_y, "h1_Eff_y", "Efficiency; Gen_Y_Z0; Efficiency");
    control_histo.h1_Pur_y = CreateHistogrm(h_pur_num_y, h_pur_den_y, "h1_Pur_y", "Purity; Rec_Y_Z0; Purity");
    control_histo.h1_Stab_y = CreateHistogrm(h_stab_num_y, h_matched_y, "h1_Stab_y", "Stability; Gen_Y_Z0; Stability");
    
    control_histo.h1_Eff_phis = CreateHistogrm(h_matched_phis, h_eff_den_phis, "h1_Eff_phis", "Efficiency; Gen_Phis_Z0; Efficiency");
    control_histo.h1_Pur_phis = CreateHistogrm(h_pur_num_phis, h_pur_den_phis, "h1_Pur_phis", "Purity; Rec_Phis_Z0; Purity");
    control_histo.h1_Stab_phis = CreateHistogrm(h_stab_num_phis, h_matched_phis, "h1_Stab_phis", "Stability; Gen_Phis_Z0; Stability");

    return control_histo;
}

// ------------------------------------------------------------------------------------------------------------------------------------

UnfoldDensities CreateUnfoldDensity(RespMatrixHisto& resp_histo, const std::string& tag) {
    UnfoldDensities densities;

    if (tag == "pt") {
        densities.pt_unf = std::make_unique<TUnfoldDensity>(resp_histo.h2_pt.GetPtr(), TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature, TUnfold::kEConstraintNone,
        TUnfoldDensity::kDensityModeBinWidth);
    } else if (tag == "y") {
        densities.y_unf = std::make_unique<TUnfoldDensity>(resp_histo.h2_y.GetPtr(), TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature, TUnfold::kEConstraintNone,
        TUnfoldDensity::kDensityModeBinWidth);
    } else if (tag == "phis") {
        densities.phis_unf = std::make_unique<TUnfoldDensity>(resp_histo.h2_phis.GetPtr(), TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature, TUnfold::kEConstraintNone,
        TUnfoldDensity::kDensityModeBinWidth);
    } else {
        std::cout << "ERROR: invalid tag input, exiting..." << std::endl;
        return densities;
    }

    return densities;
}

// ------------------------------------------------------------------------------------------------------------------------------------

UnfoldResult ApplyUnfold(std::unique_ptr<TUnfoldDensity> density, TH1D* event_histo, TH1D* resp_histo, TH1D* fake_histo, const config_struct& cfg, const std::string& tag) {
    UnfoldResult results;

    results.unf_density = std::move(density); 

    if (cfg.unfold.closure_test == true) {
        results.unf_density->SetInput(resp_histo);
    } else {
        results.unf_density->SetInput(event_histo);
    }
    // Starting the Second Event Loop on DATA !
    
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

// ------------------------------------------------------------------------------------------------------------------------------------

int VisualizeUnfoldResults(std::vector<std::unique_ptr<TCanvas>>& canvas, UnfoldResult& results, RespMatrixHisto& resp_histo, const std::string& tag) {
    std::cout << "Initializing visualization " << tag << " sample." << std::endl;
    
    std::cout << "tau scelto = " << results.tau << " (indice " << results.idx_best << ")\n";
    std::cout << "chi2A = " << results.chi2A << "  chi2L = " << results.chi2L << "  ndf = "   << results.ndf << "\n";

    // --------------------------
    // Canvas 1 - Response Matrix
    // --------------------------

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
    /*
    
    std::string name_c1 = "Response matrix " + tag;
    auto c1 = std::make_unique<TCanvas>(name_c1.c_str(), name_c1.c_str(), 800, 600);
    
    histo_resp->Draw("COLZ");

    std::string title_reco = tag + " reco [GeV]";
    std::string title_gen = tag + " gen [GeV]";
    histo_resp->GetXaxis()->SetTitle(title_reco.c_str());
    histo_resp->GetYaxis()->SetTitle(title_gen.c_str());
    */
    
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
/*
    std::string name_c4 = "Covariance matrix " + tag;
    auto c4 = std::make_unique<TCanvas>(name_c4.c_str(), name_c4.c_str(), 800, 600);
    results.h2_out_cov->Draw("COLZ");

    // ----------------------------
    // Canvas 6 - Correlation Matrix
    // ----------------------------

    std::string name_c6 = "Correlation matrix " + tag;
    auto c6 = std::make_unique<TCanvas>(name_c6.c_str(), name_c6.c_str(), 800, 600);
    results.h2_out_corr->Draw("COLZ");
*/
    //canvas.push_back(std::move(c1));
    canvas.push_back(std::move(c2));
    canvas.push_back(std::move(c3));
    //canvas.push_back(std::move(c4));
    //canvas.push_back(std::move(c6));

    return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------

int VisualizeControlPlots(std::vector<std::unique_ptr<TCanvas>>& canvas, RespMatrixHisto& resp_histo, ControlHisto& control_histo, const std::string& tag) {
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
    histo_resp->Draw("TEXTS COLZ");

    std::string title_reco = tag + " reco [GeV]";
    std::string title_gen = tag + " gen [GeV]";
    histo_resp->GetXaxis()->SetTitle(title_reco.c_str());
    histo_resp->GetYaxis()->SetTitle(title_gen.c_str());
    
    std::string name_c5 = "Efficiency/Purity/Stability_" + tag; 
    auto c5 = std::make_unique<TCanvas>(name_c5.c_str(), name_c5.c_str(), 800, 600);
    TH1D* h_eff;
    TH1D* h_pur;
    TH1D* h_stab;

    if (tag == "pt") {
        h_eff = control_histo.h1_Eff_pt.get();
        h_pur = control_histo.h1_Pur_pt.get();
        h_stab = control_histo.h1_Stab_pt.get();

    } else if (tag == "y") {
        h_eff = control_histo.h1_Eff_y.get();
        h_pur = control_histo.h1_Pur_y.get();
        h_stab = control_histo.h1_Stab_y.get();

    } else if (tag == "phis") {
        h_eff = control_histo.h1_Eff_phis.get();
        h_pur = control_histo.h1_Pur_phis.get();
        h_stab = control_histo.h1_Stab_phis.get();
    
    }

    if (h_eff && h_pur && h_stab) {
        c5->cd();
        h_eff->SetLineColor(kRed);
        h_eff->SetLineWidth(2);
        h_eff->SetMinimum(0.0);
        h_eff->SetMaximum(1.15);
        h_eff->SetStats(0);

        h_pur->SetLineColor(kBlue);
        h_pur->SetLineWidth(2);
        h_pur->SetStats(0);

        h_stab->SetLineColor(6);
        h_stab->SetLineWidth(2);
        h_stab->SetStats(0);

        h_eff->Draw("HIST");
        h_pur->Draw("HIST SAME");
        h_stab->Draw("HIST SAME");

        auto leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_eff, "Efficiency", "l");
        leg->AddEntry(h_pur, "Purity", "l");
        leg->AddEntry(h_stab, "Stability", "l");
        leg->Draw();

        c5->Update();
        canvas.push_back(std::move(c5));
    }

    canvas.push_back(std::move(c1));

    return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------
