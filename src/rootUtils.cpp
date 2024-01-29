#include "rootUtils.h"

namespace rootUtils{


    TF1* fit(std::vector<double> &x, std::vector<double> &y, double lo, double hi, std::string plotTitle){

        TCanvas *c = new TCanvas("c","Fitting DAC calib curve");
        c->SetGrid();

        std::vector<double> yErr(x.size());
        for( unsigned i=0; i<yErr.size(); i++){
            yErr[i] = 1e-2;     //TODO: what should the uncertainty on the SQUID measuremnt be?
        }

        TGraphErrors *g = new TGraphErrors(x.size(), &(x[0]), &(y[0]), 0, &(yErr[0]) );
        g->Draw("a*");

        //Fit the graph with the predefined "pol3" function
        g->Fit("pol1","","",lo,hi);

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

    TLine* createLine( std::vector<double>& data, double m ){

        double max = *(std::max_element(data.begin(), data.end()));
        TLine* l = new TLine(0,0,max,max*m);
        return l;
    }

}
