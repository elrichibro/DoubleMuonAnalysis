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

#include <ROOT/RDataFrame.hxx>

class PipelineObj {
    public:
        PipelineObj(){};
        virtual ~PipelineObj(){};
        
        virtual void Write(TFile& file) = 0;// Pure virtuals
        virtual void Draw(TCanvas& canvas) = 0;
        
        virtual void Process(){};
        
        virtual std::string GetName() const = 0;
};


class ObjectTH1 : public PipelineObj {
    private:
        std::string name_th1;
        ROOT::RDF::RResultPtr<TH1D> th1;
    public:
        ObjectTH1(std::string name, ROOT::RDF::RResultPtr<TH1D> ptr) : name_th1(name), th1(ptr) {}
        std::string GetName() const override {return name_th1;}
        
        void Write(TFile& file) override; 
        void Draw(TCanvas& canvas) override;
};

class ObjectTH2 : public PipelineObj {
    private:
        std::string name_th2;
        ROOT::RDF::RResultPtr<TH2D> th2;

    public:
        ObjectTH2(std::string name, ROOT::RDF::RResultPtr<TH2D> ptr) : name_th2(name), th2(ptr) {}
        std::string GetName() const override {return name_th2;};
        
        void Write(TFile& file) override; 
        void Draw(TCanvas& canvas) override;
};


template <typename T>
class ObjectTEff : public PipelineObj {
    private:
        std::string name_eff;
        ROOT::RDF::RResultPtr<T> hist_num;
        ROOT::RDF::RResultPtr<T> hist_den;

        std::unique_ptr<TEfficiency> eff_obj;

    public:
        ObjectTEff(std::string name, ROOT::RDF::RResultPtr<T> ptr1, ROOT::RDF::RResultPtr<T> ptr2) : name_eff(name), 
        hist_num(ptr1), hist_den(ptr2) {}
        
        std::string GetName() const override {return name_eff;};


        void Write(TFile& file) override {
            if (eff_obj) {
                file.cd();
                eff_obj->Write();
            }
        }

        void Draw(TCanvas& canvas) override {
            if (!eff_obj) {
                return;
            }
            canvas.cd();
            if constexpr (std::is_same<T, TH2D>::value) {
                //canvas.SetRightMargin(0.15);
                eff_obj->Draw("COLZ");
            } else {
                eff_obj->Draw("AP");
            }
        }
        
        void Process() override {
            if (TEfficiency::CheckConsistency(*hist_num, *hist_den)) {
                eff_obj = std::make_unique<TEfficiency>(*hist_num, *hist_den);
                eff_obj->SetName(name_eff.c_str());
            } else {
                std::cerr << "ERROR: problem with histograms." << std::endl;
            }
        }
};

class OutputManager {
    private: 
        std::vector<std::unique_ptr<PipelineObj>> pipeline;
        
        bool visualize = false;
        bool save = false;
    
        std::string output_file = "";

    public:
        OutputManager(std::string output, bool vis, bool sav) : output_file(output), visualize(vis), save(sav) {};
        ~OutputManager(){};
    
        void AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH1D> hist);
        void AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH2D> hist);

        void GetVisualize(bool value) {visualize = value;}
        bool GetVisualize() {return visualize;}

        void GetSave(bool value) {save = value;}
        bool GetSave() {return save;}
        
        template <typename T>
        void AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<T> hist1, ROOT::RDF::RResultPtr<T> hist2) {
            pipeline.push_back(std::make_unique<ObjectTEff<T>>(name, hist1, hist2));
        }
        
        void Run();
    
        void Clear() {pipeline.clear();}
};




#endif