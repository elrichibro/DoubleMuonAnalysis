#ifndef MANAGER_H
#define MANAGER_H

#include <vector>
#include <memory>
#include <string>
#include <algorithm>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"

#include "Config.h"

#include <any>

#include <ROOT/RDataFrame.hxx>

// ------------------------------------------------------------------------------------------------------------------------------------
// PipelineObj class
// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Is the Pipeline object: TH1D, TH2D, TEfficiency... used by the OutputManager class to book histograms, save or print them.
class PipelineObj {
    public:
        PipelineObj(){};// Constructor
        virtual ~PipelineObj(){};// Destructor
        
        virtual void Write(TFile& file) = 0;// Writes the object into a .root file.
        virtual void Draw(TCanvas& canvas) = 0;// Prints the object on screen within a canvas.
        
        virtual void Process(){};// Used now only by ObjectTEff
        
        virtual std::string GetName() const = 0;
};

// ------------------------------------------------------------------------------------------------------------------------------------
// ObjectTH1 class
// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Child class of PipelineObj -> Histogram 1 dimensional.
class ObjectTH1 : public PipelineObj {
    private:
        std::string name_th1;
        ROOT::RDF::RResultPtr<TH1D> th1;
    public:
        /// @brief ObjectTH1 class constructor.  
        /// @param name Name of the booked object.
        /// @param ptr Histogram 1D pointer.
        ObjectTH1(std::string name, ROOT::RDF::RResultPtr<TH1D> ptr) : name_th1(name), th1(ptr) {}
        
        std::string GetName() const override {return name_th1;}// Name getter
        
        void Write(TFile& file) override; 
        void Draw(TCanvas& canvas) override;
};

// ------------------------------------------------------------------------------------------------------------------------------------
// ObjectTH2 class
// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Child class of PipelineObj -> Histogram 2D. 
class ObjectTH2 : public PipelineObj {
    private:
        std::string name_th2;
        ROOT::RDF::RResultPtr<TH2D> th2;

    public:
        /// @brief ObjectTH2 constructor.
        /// @param name Name of the object to be booked.
        /// @param ptr Pointer to the 2D histogram.
        ObjectTH2(std::string name, ROOT::RDF::RResultPtr<TH2D> ptr) : name_th2(name), th2(ptr) {}
        
        std::string GetName() const override {return name_th2;};// Name getter
        
        void Write(TFile& file) override; 
        void Draw(TCanvas& canvas) override;
};

// ------------------------------------------------------------------------------------------------------------------------------------
// ObjectTEff class
// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief Child template class of PipelineObj. TEfficiency can accepts bouth: TH1D or TH2D. Used for complex analysis plots. 
template <typename T>
class ObjectTEff : public PipelineObj {
    private:
        std::string name_eff;
        ROOT::RDF::RResultPtr<T> hist_num;// Numerator histogram.
        ROOT::RDF::RResultPtr<T> hist_den;// Denominator histogram.

        std::unique_ptr<TEfficiency> eff_obj;// Efficiency pointer -> needed for write and draw operation. Is not a smar pointer -> starts the event loop !!!.

    public:
        /// @brief Constructor of ObjectTEff class
        /// @param name Name of the pipeline object to book.
        /// @param ptr1 Numerator histogram pointer.
        /// @param ptr2 Denominator histogram pointer.
        ObjectTEff(std::string name, ROOT::RDF::RResultPtr<T> ptr1, ROOT::RDF::RResultPtr<T> ptr2) : name_eff(name), 
        hist_num(ptr1), hist_den(ptr2) {}
        
        std::string GetName() const override {return name_eff;};
        
        /// @brief Write method for ObjectTEff
        /// @param file Output file path
        void Write(TFile& file) override {
            if (eff_obj) {
                file.cd();
                eff_obj->Write();
            }
        }

        /// @brief Prints the TEfficiency object.
        /// @param canvas Canvas pointer used.
        void Draw(TCanvas& canvas) override {
            // Check
            if (!eff_obj) {
                return;
            }
            canvas.cd();
            
            // Checked at runtime -> template usage to determinate the histogram dimension.
            if constexpr (std::is_same<T, TH2D>::value) {
                canvas.SetRightMargin(0.15);
                eff_obj->Draw("COLZ");
            } else {
                eff_obj->Draw("AP");
            }
        }
        
        /// @brief Makes a TEfficiency pointer needed for Write and Draw operations. Starts the event loop. Only used in OutputManager->Run() method. 
        void Process() override {
            if (TEfficiency::CheckConsistency(*hist_num, *hist_den)) {
                eff_obj = std::make_unique<TEfficiency>(*hist_num, *hist_den);
                eff_obj->SetName(name_eff.c_str());
            } else {
                std::cerr << "ERROR: problem with histograms." << std::endl;
            }
        }
};

// ------------------------------------------------------------------------------------------------------------------------------------
// OutputManager class
// ------------------------------------------------------------------------------------------------------------------------------------

/// @brief OutputManager class -> The instance controls the output of the program -> saves pipeline objects and print plots.
class OutputManager {
    private: 
        std::vector<std::unique_ptr<PipelineObj>> pipeline;
        
        using snapshot_type = decltype(std::declval<ROOT::RDF::RNode>().Snapshot("", "", std::vector<std::string>{}, ROOT::RDF::RSnapshotOptions{}));
        
        std::vector<snapshot_type> snapshot_vec;

        bool visualize = false;
        
        bool save_sel_plots = false;
        bool save_sel_data = false;

        std::string o_file_plots = "";
        std::string o_file_data = "";

    public:
        /// @brief OutputManager class constructor.
        /// @param output Output file path
        /// @param vis Flag for visualization option.
        /// @param sav Flag for saving the .root file containing all objects booked under request.
        OutputManager(const config_struct& cfg) : o_file_plots(cfg.io.o_file_plots), o_file_data(cfg.io.o_file_data), visualize(cfg.general.visualize), 
        save_sel_plots(cfg.general.save_sel_plots), save_sel_data(cfg.general.save_sel_data) {};
        ~OutputManager(){};
    
        // Overload method: used to add PipelineObjs to the pipe.
        void AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH1D> hist);
        void AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH2D> hist);
        
        template <typename T>
        void AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<T> hist1, ROOT::RDF::RResultPtr<T> hist2) {
            pipeline.push_back(std::make_unique<ObjectTEff<T>>(name, hist1, hist2));
        }
        
        void Run();
    
        void BookAnalysis(ROOT::RDF::RNode node, const config_struct& cfg);

        void Clear() {pipeline.clear();}
};

// ----------------
// ANALYSIS WRAPPER
// ----------------




#endif