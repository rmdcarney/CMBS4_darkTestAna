#include "rootUtils.h"

namespace rootUtils{


    TFitResultPtr sfit(std::vector<double> &x, std::vector<double> &y, unsigned lo_i, unsigned hi_i, std::string plotTitle){

        TCanvas *c = new TCanvas("c","Fitting DAC calib curve");
        c->SetGrid();

        std::vector<double> yErr(x.size());
        for( unsigned i=0; i<yErr.size(); i++){
            yErr[i] = 1e-2;     //TODO: what should the uncertainty on the SQUID measuremnt be?
        }

        TGraphErrors *g = new TGraphErrors(x.size(), &(x[0]), &(y[0]), 0, &(yErr[0]) );
        g->GetYaxis()->SetRangeUser(0.,50.);  //   Y
        g->Draw("a*");

        //Fit the graph with the predefined "pol3" function
        double lo = x[lo_i];
        double hi = x[hi_i];
        TFitResultPtr r = g->Fit("pol1","qS","",lo,hi); //S fills the FitResultPtr, q is "quiet mode"

        //std::cout<<"Fit lo, hi: "<<lo<<", "<<hi<<" which had indices: "<<lo_i<<", "<<hi_i<<std::endl;

        //Access the fit resuts
        TF1 *f = g->GetFunction("pol1");
        f->SetLineWidth(3);
        f->SetLineColor(kRed);
        gStyle->SetOptFit(0);
        f->Draw("SAME");
        c->SaveAs(plotTitle.c_str());

        return r;
    }
    
    TF1* fit(std::vector<double> &x, std::vector<double> &y, unsigned lo_i, unsigned hi_i, std::string plotTitle){

        TCanvas *c = new TCanvas("c","Fitting DAC calib curve");
        c->SetGrid();

        std::vector<double> yErr(x.size());
        for( unsigned i=0; i<yErr.size(); i++){
            yErr[i] = 1e-2;     //TODO: what should the uncertainty on the SQUID measuremnt be?
        }

        TGraphErrors *g = new TGraphErrors(x.size(), &(x[0]), &(y[0]), 0, &(yErr[0]) );
        g->GetYaxis()->SetRangeUser(0.,50.);  //   Y
        g->Draw("a*");

        //Fit the graph with the predefined "pol3" function
        double lo = x[lo_i];
        double hi = x[hi_i];
        g->Fit("pol1","q","",lo,hi);
        std::cout<<"Fit lo, hi: "<<lo<<", "<<hi<<" which had indices: "<<lo_i<<", "<<hi_i<<std::endl;

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
    
//    unsigned extractFitParsForTiming( std::filename inputData, double &R0_init, double &L_I_init ){
//
//        //Load IV measurement result, calculate R0 and L_I at this bias point!
//        double R0 = 0, L_I_IV=0;
//        TGraphErrors *gIV = (TGraphErrors*)gIVCurves->GetListOfGraphs()->FindObject(
//                Form("ch%d_%lfnF_%lfmK",
//                    channelID,
//                    configuration->run->cFilter,
//                    configuration->run->TBath));
//        if(gIV==nullptr){
//            cout<<"IV curve for this channel configuration is not found..."<<endl;
//            delete gIT;
//            delete gVT;
//            return;
//        }else{
//            if(VBias>gIV->GetX()[0]){
//                double IBias = gIV->Eval(VBias); //This is a linear interpolation - extrctuing the bias current for this voltage
//                R0 = VBias/(IBias+0.04) - configuration->Rp[channelID];//What is the 0.04, it's just an offset to push into the lienar region, doesn't apply?
//                gIV->Fit("pol1","","",VBias-0.001, VBias+0.001);//Fits a 1D line around the original IV curve
//                double dIdV = ((TF1*)gIV->GetListOfFunctions()->Last())->GetParameter(1);//Extracts gadient from linear fit
//                dIdV = 1./(1./dIdV-configuration->Rp[channelID]);//Subtracts parasitic resistance
//                L_I_IV = (1.-R0*dIdV)/(1.+R0*dIdV);//Not sure where thos is coming from - is this the loop gain? Yes
//                L_I_IV = L_I_IV<0? 100 : L_I_IV; //And this is setting a minimum.
//            }else{
//                cout<<"below IV curve..."<<endl;
//                double IBias = gIV->GetY()[0] * gIV->GetX()[0] / VBias; 
//                R0 = VBias/(IBias+0.04) - configuration->Rp[channelID];
//                R0 = R0<0?0:R0;
//                L_I_IV = R0==0?0:100;
//            }
//        }
//    }

}
