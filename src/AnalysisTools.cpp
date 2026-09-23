#include "AnalysisTools.h"
#include "Utils.h"

#include "TH1D.h"

#include <map>
#include <memory>

#include <RooHistPdf.h>
#include <RooFormulaVar.h>

#include <RooGaussian.h>
#include <RooAddPdf.h>

#include <RooFFTConvPdf.h>
#include <RooKeysPdf.h>

#include <RooFitResult.h>
#include <RooMinimizer.h>

#include <RooPlot.h>
#include <RooHist.h>
#include <TPaveText.h>

#include <TCanvas.h>
#include <TLine.h>

// ------------------------------------------------------------------------------------------------------------------------------------

void SaveEventFitCanvas(RooRealVar& mll, RooAbsPdf& model, RooAbsData& data, RooAbsPdf& bkg_pdf, const EventFitResult& res, TH1D* h_mll,
TDirectory* o_dir, const int bin_idx, const std::string& tag) {
    
    if (!o_dir) {
        return;
    }

    std::cout << "Saving histogram for bin: " << bin_idx << std::endl;
    o_dir->cd();

    int n_bins_frame = h_mll->GetNbinsX();

    TCanvas canvas(("canvas_" + tag + "_bin_" + std::to_string(bin_idx)).c_str(), "", 800, 800);
    
    TPad* pad_main  = new TPad("pad_main", "", 0.0, 0.30, 1.0, 1.0);
    TPad* pad_resid = new TPad("pad_resid", "", 0.0, 0.0, 1.0, 0.30);
    
    pad_main->SetBottomMargin(0.02);
    pad_resid->SetTopMargin(0.02);
    pad_resid->SetBottomMargin(0.35);
    pad_main->Draw();
    pad_resid->Draw();

    // ----------------
    // MAIN PAD (Fit)
    // ----------------
    pad_main->cd();

    RooPlot* frame = mll.frame(
        RooFit::Title(("mll_" + tag + "_bin_" + std::to_string(bin_idx)).c_str()), 
        RooFit::Bins(n_bins_frame)
    );

    data.plotOn(frame, RooFit::Name("data_hist"));
    model.plotOn(frame, RooFit::Name("model"), RooFit::LineColor(kBlue));
    model.plotOn(frame, RooFit::Components(bkg_pdf), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));

    frame->GetXaxis()->SetLabelSize(0);
    frame->GetXaxis()->SetTitle("");
    frame->Draw();

    // ---------------------
    // RESIDUALS PAD (Pulls)
    // ---------------------
    pad_resid->cd();

    RooHist* hpull = frame->pullHist("data_hist", "model");
    
    RooPlot* frame_pull = mll.frame(RooFit::Title(""));
    
    if (hpull) {
        frame_pull->addPlotable(hpull, "P");
    }

    frame_pull->GetYaxis()->SetTitle("Pull");
    frame_pull->GetYaxis()->SetTitleSize(0.11);
    frame_pull->GetYaxis()->SetLabelSize(0.09);
    frame_pull->GetYaxis()->SetTitleOffset(0.45);
    frame_pull->GetYaxis()->SetNdivisions(504);
    frame_pull->GetXaxis()->SetTitleSize(0.12);
    frame_pull->GetXaxis()->SetLabelSize(0.10);
    frame_pull->GetYaxis()->SetRangeUser(-5.0, 5.0);
    frame_pull->Draw();

    TLine* zero_line = new TLine(mll.getMin(), 0.0, mll.getMax(), 0.0);
    zero_line->SetLineColor(kGray+2);
    zero_line->SetLineStyle(2);
    zero_line->Draw("same");

    canvas.Write(("canvas_" + tag + "_bin_" + std::to_string(bin_idx)).c_str(), TObject::kOverwrite);

    delete frame;
    delete frame_pull;
}

// ------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::unique_ptr<TH1D>> PrepareEventFit(EventHisto& event_histo, const std::string& tag) {

    // Binned data container -> input
    std::vector<std::unique_ptr<TH1D>> container;
    
    // Central dataset [ X(P_t_Z0/Y_Z0/Phi*_Z0) , Y(Invariant Mass) ]
    TH2D* histo;

    // Kinematical quantity selection.
    if (tag == "pt") {
        histo = event_histo.h2_mll_pt.GetPtr();
    } else if (tag == "y") {
        histo = event_histo.h2_mll_y.GetPtr();
    } else if (tag == "phis") {
        histo = event_histo.h2_mll_phis.GetPtr();
    } else {
        std::cout << "ERROR: invalid tag input, exiting..." << std::endl;
    }

    int n_bins = histo->GetNbinsX();
    
    container.reserve(n_bins);

    // Saving Proyections of Invariant Mass (mll)
    for (int i = 0; i < n_bins; i++) {
        int bin_root = i + 1;

        std::string name = "h_mll_" + tag + "_bin" + std::to_string(i);
        TH1D* histo_mll = histo->ProjectionY(name.c_str(), bin_root, bin_root);

        histo_mll->SetDirectory(nullptr);

        container.emplace_back(std::unique_ptr<TH1D>(histo_mll));
    }

    return container;
}

// ------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::unique_ptr<RooDataSet>> PrepareEventFitModel(TTree* tree, const config_struct& cfg, const std::string& tag) {
    std::vector<std::unique_ptr<RooDataSet>> container;
    std::vector<double> vector_bins;

    if (cfg.unfold.use_custom_bins == true) {
        if (tag == "pt") {
            vector_bins = cfg.unfold.pt_bins.reco_vec;
        } else if (tag == "y") {
            vector_bins = cfg.unfold.y_bins.reco_vec;
        } else if (tag == "phis") {
            vector_bins = cfg.unfold.phis_bins.reco_vec;
        }
    } else {
        if (tag == "pt") {
            const auto& pt = cfg.unfold.pt_bins;
            vector_bins = CreateBins(pt.reco_bins, pt.min, pt.max, pt.distribution);

        } else if (tag == "y") {
            const auto& y = cfg.unfold.y_bins;
            vector_bins = CreateBins(y.reco_bins, y.min, y.max, y.distribution);            

        } else if (tag == "phis") {
            const auto& phis = cfg.unfold.phis_bins;
            vector_bins = CreateBins(phis.reco_bins, phis.min, phis.max, phis.distribution);
        }
    }

    float min_pt = vector_bins.front();
    float max_pt = vector_bins.back();

    int n_bins = vector_bins.size() - 1;
    
    container.reserve(n_bins);

    std::vector<std::vector<double>> mll_buffers(n_bins);
    
    
    for (int i = 0; i < n_bins; i++) {
        mll_buffers[i].reserve(20000);
    }

    float t_val, t_mll;

    tree->SetBranchAddress("InvariantMass", &t_mll);
    
    if (tag == "pt") {
        tree->SetBranchAddress("Pt_Z", &t_val);
    } else if (tag == "y") {
        tree->SetBranchAddress("Y_Z", &t_val);
    } else if (tag == "phis") {
        tree->SetBranchAddress("Phis_Z", &t_val);
    }

    Long64_t nentries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nentries; entry++) {
        tree->GetEntry(entry);
        
        int i_vector = std::distance(vector_bins.begin(),  std::upper_bound(vector_bins.begin(), vector_bins.end(), t_val)) - 1;
        if ((i_vector < 0) || (i_vector >= n_bins)) {
           continue;
        }

        mll_buffers[i_vector].push_back(t_mll);
    }

    RooRealVar mll("mll", "Invariant Mass", 60.0, 120.0);
    RooArgSet vars(mll);

    for (int i = 0; i < n_bins; i++) {
        if (mll_buffers[i].empty()) {
            continue;
        }

        std::string sample_tag = "MC";
        std::string data_name = "d_" + sample_tag + "_idx" + std::to_string(i);

        auto data = std::make_unique<RooDataSet>(data_name.c_str(), data_name.c_str(), vars);

        for (float mass : mll_buffers[i]) {
            mll.setVal(mass);
            data->add(vars);
        }

        container.push_back(std::move(data));
    }

    return container;
}

// ------------------------------------------------------------------------------------------------------------------------------------

EventFitResult EventSingleFit(int bin_idx, TH1D* h_mll, RooDataSet* d_mll_model, TDirectory* o_dir, const std::string& tag, const bool save_plots) {
    EventFitResult result;

    result.bin_idx = bin_idx;
    
    // Unpacking histogram
    if (!h_mll || !d_mll_model) {
        std::cout << "ERROR: invalid input histogram, exiting..." << std::endl;
        return result;
    }

    // Defining variable
    RooRealVar mll("mll", "m_{#mu#mu}", 60.0, 120.0, "GeV");

    // Defining data
    RooDataHist data_hist(("h_mll_" + tag + "_bin_" + std::to_string(bin_idx)).c_str(), "", mll, h_mll);

    // Signal distribution
    RooKeysPdf sig_pdf("sig_pdf", "Signal PDF from MC (KDE)", mll, *d_mll_model, RooKeysPdf::MirrorBoth);
    
    //RooHistPdf sig_pdf("sig_pdf", "Signal PDF from MC", mll, hist_mc_model);
    //sig_pdf.setInterpolationOrder(3);    
    
    // Bkg distribution -> exp
    RooRealVar lambda("lambda", "Decay constant", -0.02, -0.5, 0.0);
    RooExponential bkg_pdf("bkg_pdf", "Exponential Bkg", mll, lambda);

    double total_entries = h_mll->Integral();
    // Yields
    RooRealVar n_sig("n_sig", "Signal Yield", total_entries * 0.9, 0.0, total_entries * 1.5);
    RooRealVar n_bkg("n_bkg", "Bkg Yield", total_entries * 0.1, 0.0, total_entries * 1.5);

    // MODEL
    RooAddPdf model("model", "Signal + Bkg", RooArgList(sig_pdf, bkg_pdf), RooArgList(n_sig, n_bkg));

    // Fit
    std::unique_ptr<RooFitResult> fit_res(model.fitTo(
        data_hist, 
        RooFit::Extended(true),
        RooFit::NumCPU(8), 
        RooFit::Save(true)
    ));

    // Saving values
    result.n_sig = n_sig.getVal();
    result.n_sig_err = n_sig.getError();
    
    result.n_bkg = n_bkg.getVal();
    result.n_bkg_err = n_bkg.getError();
    
    result.lambda = lambda.getVal();
    result.lambda_err = lambda.getError();
    
    result.fit_status = fit_res ? fit_res->status() : -1;

    if (save_plots) {
        SaveEventFitCanvas(mll, model, data_hist, bkg_pdf, result, h_mll, o_dir, bin_idx, tag);
    }
    
    return result;
}

// ------------------------------------------------------------------------------------------------------------------------------------

std::vector<EventFitResult> EventFitWrapper(std::vector<std::unique_ptr<TH1D>>& container, std::vector<std::unique_ptr<RooDataSet>>& container_model, const std::string& tag, const bool save_plots, TDirectory* o_dir) {
    // Fit global container
    std::vector<EventFitResult> results;

    results.reserve(container.size());

    int succed_fits = 0;
    int failed_fits = 0;

    // Loop on Input histograms -> check needed(?)
    for (int i = 0; i < container.size(); i++) {
        
        // Getting smart pointer
        TH1D* histo = container[i].get();
        RooDataSet* data_model = container_model[i].get();

        // Fit
        EventFitResult res = EventSingleFit(i, histo, data_model, o_dir, tag, save_plots);

        // Fit status + counting
        if (res.fit_status == 0) {
            succed_fits++;
        } else {
            failed_fits++;
        }

        results.push_back(res);
    }

    std::cout << "Fits succeded: " << succed_fits << ", failed: " << failed_fits << std::endl;
    
    return results;
}