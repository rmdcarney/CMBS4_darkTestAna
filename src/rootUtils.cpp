#include "rootUtils.h"

namespace rootUtils{


    TF1* fit(std::vector<double> &x, std::vector<double> &y, unsigned lo_i, unsigned hi_i, std::string plotTitle){

        TCanvas *c = new TCanvas("c","Fitting DAC calib curve");
        c->SetGrid();

        std::vector<double> yErr(x.size());
        for( unsigned i=0; i<yErr.size(); i++){
            yErr[i] = 1e-2;     //TODO: what should the uncertainty on the SQUID measuremnt be?
        }

        TGraphErrors *g = new TGraphErrors(x.size(), &(x[0]), &(y[0]), 0, &(yErr[0]) );
        g->Draw("a*");

        //Fit the graph with the predefined "pol3" function
        double lo = x[lo_i];
        double hi = x[hi_i];
        g->Fit("pol1","q","",lo,hi);

        //Access the fit resuts
        TF1 *f = g->GetFunction("pol1");
        f->SetLineWidth(3);
        f->SetLineColor(kRed);
        gStyle->SetOptFit(0);
        f->Draw("SAME");
        c->SaveAs(plotTitle.c_str());

        return f;
    }

    void plot(std::vector<double> &x, std::vector<double> &y, std::string filename, std::string axesTitles, TLine* l){

        TCanvas *c = new TCanvas("c","PSat");
        c->SetGrid();

        TGraph *g = new TGraph(x.size(), &(x[0]), &(y[0]));
        g->Draw("a*");
        g->SetTitle(axesTitles.c_str());
        gStyle->SetOptFit(0);
        if(l){
            l->SetLineColor(kRed);
            l->SetLineWidth(3);    
            l->Draw();
        }

        c->SaveAs(filename.c_str());

        return;
    }
    
    void plot(std::vector<double> &x, std::vector<double> &y, std::string filename, std::string axesTitles, std::vector<TLine*> l){

        TCanvas *c = new TCanvas("c","PSat");
        c->SetGrid();

        TGraph *g = new TGraph(x.size(), &(x[0]), &(y[0]));
        g->Draw("a*");
        g->SetTitle(axesTitles.c_str());
        gStyle->SetOptFit(0);
        
        for( unsigned i=0; i<l.size(); i++){
                        
            l.at(i)->SetLineColor(i+2);
            l.at(i)->SetLineWidth(3);    
            l.at(i)->Draw();
        }

        c->SaveAs(filename.c_str());

        return;
    }


    TLine* createLine( std::vector<double>& data, double m ){

        double max = *(std::max_element(data.begin(), data.end()));
        TLine* l = new TLine(0,0,max,max*m);
        return l;
    }
    
    TLine* createVertLine( std::vector<double>& y, std::vector<double>& x, double point ){

        double min = *(std::min_element(y.begin(), y.end()));
        double max = *(std::max_element(y.begin(), y.end()));
        TLine* l = new TLine(x[point],min,x[point],max); //Constructor: x1, y1, x2, y2
        return l;
    }

}
