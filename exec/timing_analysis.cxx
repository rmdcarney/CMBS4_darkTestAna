//Local, project
#include "num.h"
#include "rootUtils.h"
#include "util.h"

//Local, rootlib
#include <TError.h>
#include <TLine.h>
#include <TF1.h>
#include <TVirtualFFT.h>
#include <TCanvas.h>
#include <TGraph.h>

//std lib
#include <algorithm>
#include <complex>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#define DEBUG 1
#define genWaveFreq 1.5 //Hz
/*
 * Compile like: g++ -std=c++1y -Wall $(root-config --cflags --libs) linearFit.cxx -o linFit
 * Run like: ./linFit [path_to_file]
 */

//Read in data file, calculate mean values. 
int main(int argc, char *argv[]){

    //Supress root output
    gErrorIgnoreLevel = kFatal;
    double Rp = std::atof(argv[4]);

    //====================================
    // F I L E   I O 
    //====================================
    //Get data from file
    std::cout<<".\n.\n.\nReading data file: "<<argv[1]<<"\t"<<std::flush;
    std::vector<double> time, V_gen, V_SQUID;
    std::vector< std::vector<double> > data;
    data.push_back(time); //push_back is a copy
    data.push_back(time);
    data.push_back(time);
    std::string datafile = argv[1];
    if( util::importData( datafile, data, 't' ) ) return 0;
    time    = data[0];
    V_gen   = data[1];
    V_SQUID = data[2];
    std::cout<<"[done]"<<std::endl;

    int Npoints = time.size();

    //=========================
    // Container setup
    //========================
    TGraph* g_funcGenV  = new TGraph( Npoints, &time[0], &V_gen[0] );
    TGraph* g_squidV    = new TGraph( Npoints, &time[0], &V_SQUID[0] );
//    if(DEBUG){ 
//        TCanvas* c_tmp1     = new TCanvas( "c1", "c1", 600, 400 );
//        std::cout<<"Size: "<<time.size()<<std::endl;
//        c_tmp1->cd();
//        g_funcGenV->Draw("AL");
//        c_tmp1->SaveAs("./timingRaw1.pdf");
//        c_tmp1->SaveAs("./timingRaw1.root");
//        g_squidV->Draw("AL");
//        c_tmp1->SaveAs("./timingRaw2.pdf");
//        c_tmp1->SaveAs("./timingRaw2.root");
//    }

    double *fft_cmplxRe     = new double[Npoints/2+1];      // Real waveforms will have symmetric FFTs, so just look at half (hence n/2)
    double *fft_cmplxIm     = new double[Npoints/2+1];
    double *fft_cmplxAbs    = new double[Npoints/2+1];
    double *fft_cmplxArg    = new double[Npoints/2+1];
    double *freq = new double[Npoints/2+1];
    double dF = 1./(20e-6*Npoints);                  // Why 20e-6? TODO: ask Xinran. It's for the units and given the type of FFT he's doing the scale factor is differet (I vaguely rememebr this from my DFT lectures).  
                                                    
    //=========================
    // FFT - Real to complex 
    //========================
    //Set up FFT class from Root
    TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &Npoints, "R2C");  // "R2C" - a real-input/complex-output discrete Fourier transform (DFT) in one or more dimensions,
    
    //------------------------------
    // Process stimulus waveform
    //------------------------------
    //TODO should be a function
    fftr2c->SetPoints(g_funcGenV->GetY());                      // Perform FFT on generated waveform. This function creates complex number objects from the real input values
    fftr2c->Transform();                                        // Executes the FFT
    fftr2c->GetPointsComplex(fft_cmplxRe,fft_cmplxIm);          // Returns real and imaginary parts of complex numbers

    std::cout<<"FFT computed. Ndim: "<<fftr2c->GetNdim()<<std::endl;

    //Extract the complex magnitude (z) and angle (theta) from the real and imaginary coefficients of the complex output of the fft
    for(int iFreq=0; iFreq<Npoints/2+1; iFreq++){               
//    for(int iFreq=0; iFreq<2000; iFreq++){               
        double Abs = TMath::Sqrt(fft_cmplxRe[iFreq]*fft_cmplxRe[iFreq]+fft_cmplxIm[iFreq]*fft_cmplxIm[iFreq]); // |z|    = sqrt( a^2 + b^2 ), where z = a + bi 
        double Arg = TMath::ATan(fft_cmplxIm[iFreq]/fft_cmplxRe[iFreq])                                        // arg(z) = atan( b/a ), where z = a + bi  
            + (fft_cmplxRe[iFreq]<=0 ? fft_cmplxIm[iFreq]<=0 ? - TMath::Pi() : TMath::Pi() : 0); //Ternary to figure out what quadrant we're in. 
        //cout<<">>"<<fft_cmplxRe[iFreq]<<" + i "<<fft_cmplxIm[iFreq]<<" = "<<Abs<<" e^( i "<<Arg<<")"<<endl;
        fft_cmplxAbs[iFreq] = Abs;
        fft_cmplxArg[iFreq] = Arg;
        freq[iFreq] = iFreq * dF;

    }

    //Copy data into TGraphs for later
    TGraph *g_genFFT_cmplxAbs = new TGraph(Npoints/2+1, freq, fft_cmplxAbs);
    TGraph *g_genFFT_cmplxArg = new TGraph(Npoints/2+1, freq, fft_cmplxArg);
    g_genFFT_cmplxAbs->SetNameTitle(Form("genFFT_cmplxAbs "));
    g_genFFT_cmplxArg->SetNameTitle(Form("genFFT_cmplxArg "));
//    if(DEBUG){ 
//        TCanvas* c_tmp2     = new TCanvas( "c2", "c2", 600, 400 );
//        c_tmp2->cd();
//        c_tmp2->SetLogy();
//        g_genFFT_cmplxAbs->Draw("AP");
//        gPad->SetLogy();
//        c_tmp2->SaveAs("./genFFT_cmplxAbs.pdf");
//        c_tmp2->SaveAs("./genFFT_cmplxAbs.root");
//        TCanvas* c_tmp3     = new TCanvas( "c3", "c3", 600, 400 );
//        c_tmp3->cd();
//        g_genFFT_cmplxArg->Draw("AP");
//        c_tmp3->SaveAs("./genFFT_cmplxArg.pdf");
//        c_tmp3->SaveAs("./genFFT_cmplxArg.root");
//    }
    //------------------------------
    // Process response waveform
    //------------------------------
    //TODO should be a function
    fftr2c->SetPoints(g_squidV->GetY());                      // Perform FFT on generated waveform. This function creates complex number objects from the real input values
    fftr2c->Transform();                                        // Executes the FFT
    fftr2c->GetPointsComplex(fft_cmplxRe,fft_cmplxIm);          // Returns real and imaginary parts of complex numbers

    std::cout<<"FFT computed. Ndim: "<<fftr2c->GetNdim()<<std::endl;

    //Extract the complex magnitude (z) and angle (theta) from the real and imaginary coefficients of the complex output of the fft
    for(int iFreq=0; iFreq<Npoints/2+1; iFreq++){               
//    for(int iFreq=0; iFreq<2000; iFreq++){               
        double Abs = TMath::Sqrt(fft_cmplxRe[iFreq]*fft_cmplxRe[iFreq]+fft_cmplxIm[iFreq]*fft_cmplxIm[iFreq]); // |z|    = sqrt( a^2 + b^2 ), where z = a + bi 
        double Arg = TMath::ATan(fft_cmplxIm[iFreq]/fft_cmplxRe[iFreq])                                        // arg(z) = atan( b/a ), where z = a + bi  
            + (fft_cmplxRe[iFreq]<=0 ? fft_cmplxIm[iFreq]<=0 ? - TMath::Pi() : TMath::Pi() : 0); //Ternary to figure out what quadrant we're in. 
        //cout<<">>"<<fft_cmplxRe[iFreq]<<" + i "<<fft_cmplxIm[iFreq]<<" = "<<Abs<<" e^( i "<<Arg<<")"<<endl;
        fft_cmplxAbs[iFreq] = Abs;
        fft_cmplxArg[iFreq] = Arg;
        freq[iFreq] = iFreq * dF;

    }

    //Copy data into TGraphs for later
    TGraph *g_detFFT_cmplxAbs = new TGraph(Npoints/2+1, freq, fft_cmplxAbs);
    TGraph *g_detFFT_cmplxArg = new TGraph(Npoints/2+1, freq, fft_cmplxArg);
    g_detFFT_cmplxAbs->SetNameTitle(Form("detFFT_cmplxAbs "));
    g_detFFT_cmplxArg->SetNameTitle(Form("detFFT_cmplxArg "));
//    if(DEBUG){ 
//        TCanvas* c_tmp2     = new TCanvas( "c2", "c2", 600, 400 );
//        c_tmp2->cd();
//        c_tmp2->SetLogy();
//        g_detFFT_cmplxAbs->Draw("AP");
//        gPad->SetLogy();
//        c_tmp2->SaveAs("./detFFT_cmplxAbs.pdf");
//        c_tmp2->SaveAs("./detFFT_cmplxAbs.root");
//        TCanvas* c_tmp3     = new TCanvas( "c3", "c3", 600, 400 );
//        c_tmp3->cd();
//        g_detFFT_cmplxArg->Draw("AP");
//        c_tmp3->SaveAs("./detFFT_cmplxArg.pdf");
//        c_tmp3->SaveAs("./detFFT_cmplxArg.root");
//    }


    //---------------------------------------------------------
    //Extract complex response at multiples of base frequency
    //---------------------------------------------------------
    double Freq = genWaveFreq; 
    std::cout<<"square wave frequency is "<<Freq<<"Hz"<<std::endl;
    int NFitPoints = (int)((Npoints/2+1)*(dF/Freq/2))-2;
    TGraph *gZAbs    = new TGraph( NFitPoints );
    TGraph *gZArg    = new TGraph( NFitPoints );
    double *ErrorAbs = new double[ NFitPoints];
    
    //Divide response by stimulus
    for(int iFreq=0; iFreq< NFitPoints ; iFreq++){
        int iFreqAll = (iFreq*2+1)*(Freq/dF);
        gZAbs->SetPoint(iFreq, g_detFFT_cmplxAbs->GetX()[iFreqAll], g_genFFT_cmplxAbs->GetY()[iFreqAll] / g_detFFT_cmplxAbs->GetY()[iFreqAll] );
        gZArg->SetPoint(iFreq, g_detFFT_cmplxArg->GetX()[iFreqAll], g_genFFT_cmplxArg->GetY()[iFreqAll] - g_detFFT_cmplxArg->GetY()[iFreqAll] );
        ErrorAbs[iFreq] = 0;
        for(int idf=-(int)(Freq/dF)+2; idf<(int)(Freq/dF)-2; idf++){
            if(idf==0){}
            else{
                ErrorAbs[iFreq] += TMath::Power(g_genFFT_cmplxAbs->GetY()[iFreqAll+idf] / g_detFFT_cmplxAbs->GetY()[iFreqAll],2);
                ErrorAbs[iFreq] += TMath::Power(g_detFFT_cmplxAbs->GetY()[iFreqAll+idf] * gZAbs->GetY()[iFreq] / g_detFFT_cmplxAbs->GetY()[iFreqAll],2);
            }
        }
        ErrorAbs[iFreq]/=2*(int)(Freq/dF)-4.;
        ErrorAbs[iFreq] = TMath::Sqrt(ErrorAbs[iFreq]);
    }

    std::cout<<"Some examples of errors: "<<ErrorAbs[0]<<", "<<ErrorAbs[200]<<", "<<ErrorAbs[500]<<std::endl;
//    if(DEBUG){ 
//        TCanvas* c_tmp4     = new TCanvas( "c4", "c4", 600, 400 );
//        c_tmp4->cd();
//        c_tmp4->SetLogy();
//        gZAbs->Draw("AP");
//        gPad->SetLogy();
//        c_tmp4->SaveAs("./norm_cmplxAbs.pdf");
//        c_tmp4->SaveAs("./norm_cmplxAbs.root");
//        TCanvas* c_tmp5     = new TCanvas( "c5", "c5", 600, 400 );
//        c_tmp5->cd();
//        gZArg->Draw("AP");
//        c_tmp5->SaveAs("./norm_cmplxArg.pdf");
//        c_tmp5->SaveAs("./norm_cmplxArg.root");
//    }

    //========
    // Fit 
    //========
    //-------------------------
    // Parameter setup
    //-------------------------
    //Open IV output file and extract relevant data
    std::cout<<".\n.\n.\nReading data file: "<<argv[2]<<"\t"<<std::flush;
    std::vector<double> IV_Vps, IV_Vtes, IV_Ites;
    std::vector< std::vector<double> > IV_data;
    IV_data.push_back(IV_Vps); //push_back is a copy
    IV_data.push_back(IV_Vps);
    IV_data.push_back(IV_Vps);
    std::string IVfile = argv[2];
    if( util::importData( IVfile, IV_data, 'i' ) ) return 0;
    IV_Vps      = IV_data[0]; // Voltage supplied by power supply [V]
    IV_Vtes     = IV_data[1]; // Corresponding voltage across TES [nV]
    IV_Ites     = IV_data[2]; // Corresponding measured current flowing through TES at that voltage [uA]
    std::cout<<"[done]"<<std::endl;

    //Interpolate between the IV curve points and look-up the corresponding TES voltage and current
    double V0(1.), I0(1.), R0(1.); //Small signal TES voltage and current and Vbias 
    double V_functionGeneratorOffset = std::atof(argv[3]);
    std::pair<unsigned, unsigned> loHiLU_indices;

    std::cout<<"V_fgo   = "<<V_functionGeneratorOffset<<std::endl;
    std::cout<<"IV_Vps  = "<<IV_Vps[0]<<std::endl;

    if( util::lookUpIV( std::atof(argv[3]), V0, I0, IV_Vps, IV_Vtes, IV_Ites, loHiLU_indices )) return 0;
    R0 = V0/I0 - Rp;// Make sure to remove parasitic resistance too. [nV]/[uA] = [mOhm]

    //Fit tangent at IV_Vtes (+/- 1nV) to obtain dI/dV
    TFitResultPtr rFitRslt  = rootUtils::sfit( IV_Vtes, IV_Ites, loHiLU_indices.first+5, loHiLU_indices.second-1, "TES_IV_fitatR0.pdf"); 
    if( !rFitRslt->IsValid() ){
        TFile* fitPtrFile = new TFile("timingAna_fitResultPtrFAIL_info.root", "RECREATE");
        rFitRslt->Write();
        fitPtrFile->Close();
        std::cout<<"ERROR: Fitting the tangent to the IV curve during paramter estimation failed. Please check the fit output curve (pdf) and fit data (root file) and try again.\nExiting."<<std::endl;
        return 0;
    }    

    //Calculate imaginary loop gain
    double dIdV             = rFitRslt->Parameter(1); //Get gradient from fit
    double L_im_IV          = (1-R0*dIdV) / (1+R0*dIdV); //See write-up in readMe for derivation
   
    if( DEBUG ){
        std::cout<<"------------ PARAMETER SETUP ------------"<<std::endl;
        std::cout<<"V0      : "<<V0<<std::endl;
        std::cout<<"I0      : "<<I0<<std::endl;
        std::cout<<"R0      : "<<R0<<std::endl;
        std::cout<<"dIdV    : "<<dIdV<<std::endl;
        std::cout<<"L_im_IV : "<<L_im_IV<<std::endl;
        std::cout<<"------------------------------------------"<<std::endl;
    }

    //fitting function DOI: 10.1063/1.4759111
    ////[0] -> Rp+R0(1+beta) [1]->L            
    ////[2] -> tau/(L_I*R0*(2+beta))           
    ////[3] -> tau1                            
    ////[4] -> (1-L_I)/ (L_I*R0*(2+beta))      
    ////[5] ->g1T1/(g1T1+g1b)/(L_I*R0*(2+beta))
    double R0_IV = R0;
    double par[6] = {Rp+R0*(1+0.1), 2e-7, 0.001/(L_im_IV*R0*(2+0.1)), 
        0.001, (1.-L_im_IV)/(L_im_IV*R0*(2+0.1)), 0.3/(L_im_IV*R0*(2+0.1))};
    //complexZ *impediance = new complexZ(par, Rp, R0);
    util::complexChi2 Chi2(gZAbs->GetX(), gZAbs->GetY(), gZArg->GetY(), ErrorAbs, gZAbs->GetN(), par, Rp, R0);

    //-------------------------
    // Fit FFT 
    //-------------------------
    ROOT::Fit::Fitter fitter;
    fitter.Config().SetParamsSettings(6,par);
    fitter.Config().ParSettings(0).SetLimits(0,50);
    //fitter.Config().ParSettings(0).Fix();
    fitter.Config().ParSettings(0).SetName("Rp+R0(1+beta)");
    fitter.Config().ParSettings(1).SetLimits(0,0.01);
    fitter.Config().ParSettings(1).SetName("L");
    fitter.Config().ParSettings(2).SetLimits(0,1000);
    fitter.Config().ParSettings(2).SetName("tau/(L_I*R0*(2+beta))");
    fitter.Config().ParSettings(3).SetLimits(0.0001,1);
    fitter.Config().ParSettings(3).SetName("tau_a");
    fitter.Config().ParSettings(4).SetLimits(-20, gZArg->GetY()[0]<-1? 0 : 1000);
    fitter.Config().ParSettings(4).SetName("(1-L_I)/ (L_I*R0*(2+beta))");
    fitter.Config().ParSettings(5).SetLimits(0,10000);
    fitter.Config().ParSettings(5).SetName("GTa/(GTa+Gab)/(L_I*R0*(2+beta))");
    fitter.Config().SetMinimizer("Minuit2","Migrad");
    fitter.Config().MinimizerOptions().SetPrintLevel(2);
    fitter.Config().MinimizerOptions().SetMaxIterations(1000);
    fitter.Config().MinimizerOptions().SetMaxFunctionCalls(1e4);
    std::cout<<"fit start..."<<std::endl;
    fitter.FitFCN(6,Chi2,0, gZAbs->GetN()/2, false);
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);
     

    //-------------------------
    // Parameter extraction
    //-------------------------
    //Force R0 to the IV curve value:
    double beta = R0==0? 0 : (result.GetParams()[0] - Rp)/R0 - 1.;
    double L_I  = R0==0? 0 : 1./(R0*(2.+beta)*result.GetParams()[4] + 1.);
    double tau  = result.GetParams()[2]*(L_I*R0*(2+beta));
    double L    = result.GetParams()[1];
    double tauA = result.GetParams()[3];
    double gRatio =R0==0? 0 : 1./(result.GetParams()[5]*L_I*R0*(2+beta))-1;//g1b/g1T1
    std::cout<<"beta	="<<beta<<std::endl;
    std::cout<<"L_I	="<<L_I<<std::endl;
    std::cout<<"tau	="<<tau<<std::endl;
    std::cout<<"Gab/GTa	="<<gRatio<<std::endl;
    std::cout<<"R0	="<<R0<<std::endl;
    std::cout<<"Rp	="<<Rp<<std::endl;
    std::cout<<"L_im_IV	="<<L_im_IV<<std::endl;
    Chi2.fimpediance->setPar(result.GetParams());


    return 0;

}

