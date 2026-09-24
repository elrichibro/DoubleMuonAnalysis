#include "Event.h"// main include

#include "Utils.h"

#include <RooAddPdf.h>
#include <RooKeysPdf.h>
#include <RooExponential.h>

#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooHist.h>

#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>

#include <iostream>

// ------------------------------------------------------------------------------------------------------------------------------------

EventSelectionHisto BuildEventSelection_Histo(ROOT::RDF::RNode node, const config_struct& cfg) {
    ROOT::RDF::RNode node_event = node;
    EventSelectionHisto histo;
      
    std::vector<double> bins_pt, bins_y, bins_phis;
    int n_pt, n_y, n_phis;

    // Choosing binning options
    if (cfg.unfold.use_custom_bins == true) {
        bins_pt = cfg.unfold.pt_bins.reco_vec;
        bins_y = cfg.unfold.y_bins.reco_vec;
        bins_phis = cfg.unfold.phis_bins.reco_vec;

    } else {
        const auto& pt = cfg.unfold.pt_bins;
        const auto& y = cfg.unfold.y_bins;
        const auto& phis = cfg.unfold.phis_bins;

        bins_pt = CreateBins(pt.reco_bins, pt.min, pt.max, pt.distribution);
        bins_y = CreateBins(y.reco_bins, y.min, y.max, y.distribution);            
        bins_phis = CreateBins(phis.reco_bins, phis.min, phis.max, phis.distribution);
    }

    // Mll bins - HARDCODED
    int n_mll = 70;
    
    std::vector<double> mll_bins(n_mll + 1);
    double step = (120.0 - 60.0) / n_mll;

    for (int i = 0; i < n_mll + 1; i++) {
        mll_bins[i] = 60.0 + (i * step);
    }

    n_pt = static_cast<int>(bins_pt.size()) - 1;
    n_y = static_cast<int>(bins_y.size()) - 1;
    n_phis = static_cast<int>(bins_phis.size()) - 1;

    ROOT::RDF::TH2DModel model_mll_pt("h2_mll_pt", "Pt vs Mll", n_pt, bins_pt.data(), n_mll, mll_bins.data());
    ROOT::RDF::TH2DModel model_mll_y("h2_mll_y", "Y vs Mll", n_y, bins_y.data(), n_mll, mll_bins.data());
    ROOT::RDF::TH2DModel model_mll_phis("h2_mll_phis", "Phis vs Mll", n_phis, bins_phis.data(), n_mll, mll_bins.data());

    histo.h1_mll = node_event.Histo1D({"hInvMass_fit","", 70, 60.0, 120.0}, "InvariantMass");
        
    histo.h1_pt = node_event.Histo1D({"hPt_event", "", n_pt, bins_pt.data()}, "Pt_Z");
    histo.h1_y = node_event.Histo1D({"hRapidity_event", "", n_y, bins_y.data()}, "Y_Z");
    histo.h1_phis = node_event.Histo1D({"hPhis_event", "", n_phis, bins_phis.data()}, "Phis_Z");

    histo.h2_mll_pt = node_event.Histo2D(model_mll_pt, "Pt_Z", "InvariantMass");
    histo.h2_mll_y = node_event.Histo2D(model_mll_y, "Y_Z", "InvariantMass");
    histo.h2_mll_phis = node_event.Histo2D(model_mll_phis, "Phis_Z", "InvariantMass");

    return histo;
}

// ------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::unique_ptr<TH1D>> BuildEventFit_Histo(EventSelectionHisto& ev_sel_histo, const std::string& tag) {

    // Binned data container -> input
    std::vector<std::unique_ptr<TH1D>> container;
    
    // Central dataset [ X(P_t_Z0/Y_Z0/Phi*_Z0) , Y(Invariant Mass) ]
    TH2D* histo;

    // Kinematical quantity selection.
    if (tag == "pt") {
        histo = ev_sel_histo.h2_mll_pt.GetPtr();
    } else if (tag == "y") {
        histo = ev_sel_histo.h2_mll_y.GetPtr();
    } else if (tag == "phis") {
        histo = ev_sel_histo.h2_mll_phis.GetPtr();
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



std::vector<std::unique_ptr<RooDataSet>> BuildEventFit_SignalModel(TTree* tree, const config_struct& cfg, const std::string& tag) {
    
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



std::vector<EventFitResult> EventSingleFit_Wrapper(std::vector<std::unique_ptr<TH1D>>& container, std::vector<std::unique_ptr<RooDataSet>>& container_model, const std::string& tag, const bool save_plots, TDirectory* o_dir) {
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

// ------------------------------------------------------------------------------------------------------------------------------------

std::unique_ptr<TH1D> EventFit_SignalHisto_Wrapper(EventSelectionHisto& ev_sel_histo, const config_struct& cfg) {
    std::string tag = cfg.event.event_quantity;
    
    // MC DATASET -> Model
    ROOT::RDataFrame mc_frame("MC_Event_Tree", cfg.selection.o_sel_file_data);
    ROOT::RDF::RNode node_event_model = mc_frame;

    EventSelectionHisto event_model_histo = BuildEventSelection_Histo(node_event_model, cfg);

    std::vector<std::unique_ptr<TH1D>> event_container = BuildEventFit_Histo(ev_sel_histo, tag);

    // Unpacking model
    std::unique_ptr<TFile> model_file(TFile::Open(cfg.selection.o_sel_file_data.c_str(), "READ"));
    TTree* model_tree = model_file->Get<TTree>("MC_Event_Tree");

    if (!model_file || !model_file) {
        std::cout << "ERROR: invalid read operation" << std::endl;
        return nullptr;
    }

    std::vector<std::unique_ptr<RooDataSet>> event_model_container = BuildEventFit_SignalModel(model_tree, cfg, tag);
    
    std::unique_ptr<TFile> o_fit_file = nullptr;
    TDirectory* fit_dir;

    if (cfg.event.save_fit_plots) {
        // Output file
        o_fit_file = std::make_unique<TFile>(cfg.event.o_event_file.c_str(), "UPDATE");        
        
        if (o_fit_file->IsZombie()) {
            std::cout << "ERROR: invalid output file, exiting..." << cfg.event.o_event_file << std::endl;
            return nullptr;
        }

        o_fit_file->cd();

        std::string dir_name = "Event_InvMass_DATA_fits_" + tag;
        fit_dir = o_fit_file->GetDirectory(dir_name.c_str());

        if (!fit_dir) {
            fit_dir = o_fit_file->mkdir(dir_name.c_str());
        }
    }

    std::vector<EventFitResult> event_results = EventSingleFit_Wrapper(event_container, event_model_container, tag, cfg.event.save_fit_plots, fit_dir);

    if (cfg.event.save_fit_plots) {
        o_fit_file->cd();
        o_fit_file->Write();
        o_fit_file->Close();
    }

    std::unique_ptr<TH1D> h_sig = BuildFitResult_Histo(event_results, cfg, tag);

    return h_sig;
}

// ------------------------------------------------------------------------------------------------------------------------------------

std::unique_ptr<TH1D> BuildFitResult_Histo(const std::vector<EventFitResult>& results, const config_struct& cfg, const std::string& tag) {
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
    
    int n_bins = vector_bins.size() - 1;

    std::string name = "h_signal_" + tag;
    auto histo = std::make_unique<TH1D>(name.c_str(), name.c_str(), n_bins, vector_bins.data());

    histo->SetDirectory(nullptr);

    for (int i = 0; i < n_bins; ++i) {
        int bin_root = i + 1; 

        histo->SetBinContent(bin_root, results[i].n_sig);
        histo->SetBinError(bin_root, results[i].n_sig_err);
    }
    return histo;
}

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
