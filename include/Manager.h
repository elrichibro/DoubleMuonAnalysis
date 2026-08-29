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
            
            virtual void Process(TCanvas& canvas, TFile& file) = 0;
            virtual std::string GetName() const = 0;
    };


    class ObjectTH1 : public PipelineObj {
        private:
            std::string name_th1;
            ROOT::RDF::RResultPtr<TH1D> th1;
        public:
            ObjectTH1(std::string name, ROOT::RDF::RResultPtr<TH1D> ptr) : name_th1(name), th1(ptr) {}
            std::string GetName() const override {return name_th1;}
            void Process(TCanvas& canvas, TFile& file) override;
    };

    class ObjectTH2 : public PipelineObj {
        private:
            std::string name_th2;
            ROOT::RDF::RResultPtr<TH2D> th2;

        public:
            ObjectTH2(std::string name, ROOT::RDF::RResultPtr<TH2D> ptr) : name_th2(name), th2(ptr) {}
            std::string GetName() const override {return name_th2;};
            void Process(TCanvas& canvas, TFile& file) override;
    };

    class ObjectTEff : public PipelineObj {
        private:
            std::string name_eff;
            ROOT::RDF::RResultPtr<TH1D> hist_num;
            ROOT::RDF::RResultPtr<TH1D> hist_den;

        public:
            ObjectTEff(std::string name, ROOT::RDF::RResultPtr<TH1D> ptr1, ROOT::RDF::RResultPtr<TH1D> ptr2) : name_eff(name), 
            hist_num(ptr1), hist_den(ptr2) {}
            
            std::string GetName() const override {return name_eff;};
            
            void Process(TCanvas& canvas, TFile& file) override;
    };

    class OutputManager {
        private: 
            std::vector<std::unique_ptr<PipelineObj>> pipeline;
        
        public:
            OutputManager(){};
            ~OutputManager(){};
        
            void AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH1D> hist);
            void AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH2D> hist);
            void AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH1D> hist1, ROOT::RDF::RResultPtr<TH1D> hist2);
            
            void Run(const std::string& filename);
        
            void Clear() {pipeline.clear();}
    };




    #endif