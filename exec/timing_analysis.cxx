//Local, project
#include "num.h"
#include "rootUtils.h"
#include "util.h"

//Local, rootlib
#include <TComplex.h>
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

#define DEBUG 0
#define genWaveFreq 1.5 //Hz
#define sampleRate 150e3 //Hz 
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

    //int Npoints     = time.size();
    int Npoints     = 800e3;
    const int fSize = (Npoints/2)+1;

    std::cout<<"Number of data points: "<<Npoints<<std::endl;
    std::cout<<"Number of frequency poitns for fft: "<<fSize<<std::endl;

    //=========================
    // Container setup
    //========================
    TGraph g_funcGenV( Npoints, &time[0], &V_gen[0] );
    TGraph g_squidV( Npoints, &time[0], &V_SQUID[0] );
    if(DEBUG){ 
        TCanvas c_tmp1( "c1", "c1", 600, 400 );
        std::cout<<"Size: "<<time.size()<<std::endl;
        c_tmp1.cd();
        g_funcGenV.Draw("AL");
        c_tmp1.SaveAs("./timing_VFuncGen.pdf");
        c_tmp1.SaveAs("./timing_VFuncGen.root");
        g_squidV.Draw("AL");
        c_tmp1.SaveAs("./timing_VSQUID.pdf");
        c_tmp1.SaveAs("./timing_VSQUID.root");
    }

    std::vector<double> fft_cmplxRe(fSize,0.); 
    std::vector<double> fft_cmplxIm(fSize,0.);      // Real waveforms will have symmetric FFTs, so just look at half (hence n/2)
    std::vector<double> fft_cmplxAbs(fSize,0.);    
    std::vector<double> fft_cmplxArg(fSize,0.);    
    std::vector<double> freq(fSize,0.);            
    
    double dF = sampleRate/Npoints;                         // dF in Hz 
                                                    
    //=========================
    // FFT - Real to complex 
    //========================
    //Set up FFT class from Root
    TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &Npoints, "R2C");  // "R2C" - a real-input/complex-output discrete Fourier transform (DFT) in one or more dimensions,
    
    //------------------------------
    // Process stimulus waveform
    //------------------------------
    //TODO should be a function
    fftr2c->SetPoints(g_funcGenV.GetY());                      // Perform FFT on generated waveform. This function creates complex number objects from the real input values
    fftr2c->Transform();                                        // Executes the FFT
    fftr2c->GetPointsComplex(&fft_cmplxRe[0],&fft_cmplxIm[0]);          // Returns real and imaginary parts of complex numbers

    std::cout<<"FFT computed. Ndim: "<<fftr2c->GetNdim()<<std::endl;

    //Extract the complex magnitude (z) and angle (theta) from the real and imaginary coefficients of the complex output of the fft
    for(int iFreq=0; iFreq<fSize; iFreq++){               
        double Abs = TMath::Sqrt(fft_cmplxRe[iFreq]*fft_cmplxRe[iFreq]+fft_cmplxIm[iFreq]*fft_cmplxIm[iFreq]); // |z|    = sqrt( a^2 + b^2 ), where z = a + bi 
        double Arg = TMath::ATan(fft_cmplxIm[iFreq]/fft_cmplxRe[iFreq])                                        // arg(z) = atan( b/a ), where z = a + bi  
            + (fft_cmplxRe[iFreq]<=0 ? fft_cmplxIm[iFreq]<=0 ? - TMath::Pi() : TMath::Pi() : 0); //Ternary to figure out what quadrant we're in. 
        //cout<<">>"<<fft_cmplxRe[iFreq]<<" + i "<<fft_cmplxIm[iFreq]<<" = "<<Abs<<" e^( i "<<Arg<<")"<<endl;
        fft_cmplxAbs[iFreq] = Abs;
        fft_cmplxArg[iFreq] = Arg;
        freq[iFreq] = iFreq * dF;

    }

    std::cout<<"DEBUG: Freq values: "<<freq[0]<<", "<<freq[1]<<", "<<freq[2]<<std::endl;

    //Copy data into TGraphs for later
    TGraph g_genFFT_cmplxAbs(fSize, &freq[0], &fft_cmplxAbs[0]);
    TGraph g_genFFT_cmplxArg(fSize, &freq[0], &fft_cmplxArg[0]);
    g_genFFT_cmplxAbs.SetNameTitle(Form("genFFT_cmplxAbs "));
    g_genFFT_cmplxArg.SetNameTitle(Form("genFFT_cmplxArg "));
    if(DEBUG){ 
        TCanvas c_tmp2( "c2", "c2", 600, 400 );
        c_tmp2.cd();
        c_tmp2.SetLogy();
        g_genFFT_cmplxAbs.SetMarkerStyle(20);
        g_genFFT_cmplxAbs.Draw("ALP");
        gPad->SetLogy();
        gPad->SetLogx();
        c_tmp2.SaveAs("./genFFT_cmplxAbs.pdf");
        c_tmp2.SaveAs("./genFFT_cmplxAbs.root");
        TCanvas c_tmp3( "c3", "c3", 600, 400 );
        c_tmp3.cd();
        g_genFFT_cmplxArg.Draw("AP");
        c_tmp3.SaveAs("./genFFT_cmplxArg.pdf");
        c_tmp3.SaveAs("./genFFT_cmplxArg.root");
    }
    //------------------------------
    // Process response waveform
    //------------------------------
    //TODO should be a function
    fftr2c->SetPoints(g_squidV.GetY());                      // Perform FFT on generated waveform. This function creates complex number objects from the real input values
    fftr2c->Transform();                                        // Executes the FFT
    fftr2c->GetPointsComplex(&fft_cmplxRe[0],&fft_cmplxIm[0]);          // Returns real and imaginary parts of complex numbers

    std::cout<<"FFT computed. Ndim: "<<fftr2c->GetNdim()<<std::endl;

    //Extract the complex magnitude (z) and angle (theta) from the real and imaginary coefficients of the complex output of the fft
    for(int iFreq=0; iFreq<fSize; iFreq++){               
        double Abs = TMath::Sqrt(fft_cmplxRe[iFreq]*fft_cmplxRe[iFreq]+fft_cmplxIm[iFreq]*fft_cmplxIm[iFreq]); // |z|    = sqrt( a^2 + b^2 ), where z = a + bi 
        double Arg = TMath::ATan(fft_cmplxIm[iFreq]/fft_cmplxRe[iFreq])                                        // arg(z) = atan( b/a ), where z = a + bi  
            + (fft_cmplxRe[iFreq]<=0 ? fft_cmplxIm[iFreq]<=0 ? - TMath::Pi() : TMath::Pi() : 0); //Ternary to figure out what quadrant we're in. 
        //cout<<">>"<<fft_cmplxRe[iFreq]<<" + i "<<fft_cmplxIm[iFreq]<<" = "<<Abs<<" e^( i "<<Arg<<")"<<endl;
        fft_cmplxAbs[iFreq] = Abs;
        fft_cmplxArg[iFreq] = Arg;
        freq[iFreq] = iFreq * dF;

    }

    //Copy data into TGraphs for later
    TGraph g_detFFT_cmplxAbs(fSize, &freq[0], &fft_cmplxAbs[0]);
    TGraph g_detFFT_cmplxArg(fSize, &freq[0], &fft_cmplxArg[0]);
    g_detFFT_cmplxAbs.SetNameTitle(Form("detFFT_cmplxAbs "));
    g_detFFT_cmplxArg.SetNameTitle(Form("detFFT_cmplxArg "));
    if(DEBUG){ 
        TCanvas c_tmp2( "c2", "c2", 600, 400 );
        c_tmp2.cd();
        c_tmp2.SetLogy();
        gPad->SetLogx();
        g_detFFT_cmplxAbs.SetMarkerStyle(20);
        g_detFFT_cmplxAbs.Draw("ALP");
        gPad->SetLogy();
        c_tmp2.SaveAs("./detFFT_cmplxAbs.pdf");
        c_tmp2.SaveAs("./detFFT_cmplxAbs.root");
        TCanvas c_tmp3( "c3", "c3", 600, 400 );
        c_tmp3.cd();
        g_detFFT_cmplxArg.Draw("AP");
        c_tmp3.SaveAs("./detFFT_cmplxArg.pdf");
        c_tmp3.SaveAs("./detFFT_cmplxArg.root");
    }


    //---------------------------------------------------------
    //Extract complex response at multiples of base frequency
    //---------------------------------------------------------
    std::cout<<"square wave frequency is "<<genWaveFreq<<"Hz"<<std::endl;
    int NFitPoints = (int)((fSize)*(dF/genWaveFreq/2))-2; 

    std::cout<<"======= FIT DEBUG ======\ndF = "<<dF<<"\ndF/genWaveFreq/2 = "<<dF/genWaveFreq/2<<"\nNFitPoint = "<<NFitPoints<<std::endl;

    TGraph gZAbs( NFitPoints );
    TGraph gZArg ( NFitPoints );
    std::vector<double> ErrorAbs(NFitPoints);
    
    //Divide response by stimulus
    for(int iFreq=0; iFreq< NFitPoints ; iFreq++){
        int iFreqAll = (iFreq*2+1)*(genWaveFreq/dF);
        if( iFreq<20 ){
            std::cout<<"iFreqAll: "<<iFreqAll<<", genFFT: "<<g_genFFT_cmplxAbs.GetY()[iFreqAll]<<", detFFT: "<<g_detFFT_cmplxAbs.GetY()[iFreqAll]<<std::endl;
        }
        gZAbs.SetPoint(iFreq, g_detFFT_cmplxAbs.GetX()[iFreqAll], g_genFFT_cmplxAbs.GetY()[iFreqAll] / g_detFFT_cmplxAbs.GetY()[iFreqAll] );
        gZArg.SetPoint(iFreq, g_detFFT_cmplxArg.GetX()[iFreqAll], g_genFFT_cmplxArg.GetY()[iFreqAll] - g_detFFT_cmplxArg.GetY()[iFreqAll] );
        ErrorAbs[iFreq] = 0;
        for(int idf=-(int)(genWaveFreq/dF)+2; idf<(int)(genWaveFreq/dF)-2; idf++){
            if(idf==0){}
            else{
                ErrorAbs[iFreq] += TMath::Power(g_genFFT_cmplxAbs.GetY()[iFreqAll+idf] / g_detFFT_cmplxAbs.GetY()[iFreqAll],2);
                ErrorAbs[iFreq] += TMath::Power(g_detFFT_cmplxAbs.GetY()[iFreqAll+idf] * gZAbs.GetY()[iFreq] / g_detFFT_cmplxAbs.GetY()[iFreqAll],2);
            }
        }
        ErrorAbs[iFreq]/=2*(int)(genWaveFreq/dF)-4.;
        ErrorAbs[iFreq] = TMath::Sqrt(ErrorAbs[iFreq]);
    }

    std::cout<<"Some examples of errors: "<<ErrorAbs[0]<<", "<<ErrorAbs[200]<<", "<<ErrorAbs[500]<<std::endl;
    if(DEBUG){ 
        TCanvas c_tmp4( "c4", "c4", 600, 400 );
        c_tmp4.cd();
        c_tmp4.SetLogy();
        gPad->SetLogx();
        gZAbs.SetMarkerStyle(20);
        gZAbs.Draw("ALP");
        gPad->SetLogy();
        c_tmp4.SaveAs("./norm_cmplxAbs.pdf");
        c_tmp4.SaveAs("./norm_cmplxAbs.root");
        TCanvas c_tmp5( "c5", "c5", 600, 400 );
        c_tmp5.cd();
        gZArg.Draw("AP");
        c_tmp5.SaveAs("./norm_cmplxArg.pdf");
        c_tmp5.SaveAs("./norm_cmplxArg.root");
    }

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

    std::cout<<"V_functionGeneratorOffset   = "<<V_functionGeneratorOffset<<std::endl;
    std::cout<<"IV_Vps  = "<<IV_Vps[0]<<std::endl;

    if( util::lookUpIV( std::atof(argv[3]), V0, I0, IV_Vps, IV_Vtes, IV_Ites, loHiLU_indices )) return 0;
    R0 = V0/I0 - Rp;// Make sure to remove parasitic resistance too. [nV]/[uA] = [mOhm]

    //Fit tangent at IV_Vtes (+/- 1nV) to obtain dI/dV
    TFitResultPtr rFitRslt  = rootUtils::sfit( IV_Vtes, IV_Ites, loHiLU_indices.first+5, loHiLU_indices.second-1, "TES_IV_fitatR0.pdf"); 
    if( !rFitRslt->IsValid() ){
        TFile fitPtrFile("timingAna_fitResultPtrFAIL_info.root", "RECREATE");
        rFitRslt->Write();
        fitPtrFile.Close();
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
    util::complexChi2 Chi2(gZAbs.GetX(), gZAbs.GetY(), gZArg.GetY(), &ErrorAbs[0], gZAbs.GetN(), par, Rp, R0);

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
    fitter.Config().ParSettings(4).SetLimits(-20, gZArg.GetY()[0]<-1? 0 : 1000);
    fitter.Config().ParSettings(4).SetName("(1-L_I)/ (L_I*R0*(2+beta))");
    fitter.Config().ParSettings(5).SetLimits(0,10000);
    fitter.Config().ParSettings(5).SetName("GTa/(GTa+Gab)/(L_I*R0*(2+beta))");
    fitter.Config().SetMinimizer("Minuit2","Migrad");
    fitter.Config().MinimizerOptions().SetPrintLevel(2);
    fitter.Config().MinimizerOptions().SetMaxIterations(1000);
    fitter.Config().MinimizerOptions().SetMaxFunctionCalls(1e4);
    std::cout<<"fit start..."<<std::endl;
    fitter.FitFCN(6,Chi2,0, gZAbs.GetN()/2, false);
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);
     

    //-------------------------
    // Parameter extraction
    //-------------------------
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

    //-------------------------------------------
    // Draw results in time and frequency domain
    //-------------------------------------------
    TGraph gZFitAbs;
    TGraph gZFitArg;
    for(int ip=0; ip<gZAbs.GetN(); ip++){
        double freq_i =  gZAbs.GetX()[ip];
        gZFitAbs.SetPoint(ip, freq_i, std::abs(Chi2.fimpediance->Z(freq_i)));//fZAbs.Eval(freq));
        gZFitArg.SetPoint(ip, freq_i, std::arg(Chi2.fimpediance->Z(freq_i)));//fZArg.Eval(freq));
    }

    TCanvas cZFitAbs("cZFitAbs","cFitAbs",600,400);
    cZFitAbs.cd();
        gPad->SetLogx();
    cZFitAbs.SetLogy();
    gZFitAbs.SetMarkerStyle(20);
    gZFitAbs.Draw("ALP");
    gPad->SetLogy();
    cZFitAbs.Update();
    cZFitAbs.SaveAs("./gZFitAbs.pdf");
    cZFitAbs.SaveAs("./gZFitAbs.root");
    TCanvas cZFitArg("cZFitArg","cFitArg",600,400);
    cZFitArg.cd();
    cZFitArg.SetLogy();
    gZFitArg.Draw("ALP");
        gPad->SetLogx();
    gPad->SetLogy();
    cZFitArg.Update();
    cZFitArg.SaveAs("./gZFitArg.pdf");
    cZFitArg.SaveAs("./gZFitArg.root");

    //111111111111111111111111111111
    //draw fitted time domain I-T
    //111111111111111111111111111111
    TVirtualFFT *fftc2r = TVirtualFFT::FFT(1, &Npoints, "C2R");
    int ih=0;
    unsigned counter(0), counterAll(0);
    for(int ip=0; ip<g_genFFT_cmplxArg.GetN(); ip++){
        double freq_i = g_genFFT_cmplxArg.GetX()[ip];
        std::complex<double> Zfit = Chi2.fimpediance->Z(freq_i);
        //complex<double> Vsquare = TMath::Abs(freq_i/Freq.-(ih*2+1))<0.001 ? g_genFFT_cmplxAbs.GetY()[ip] : 0;
        std::complex<double> Vsquare = g_genFFT_cmplxAbs.GetY()[ip]*TMath::Cos(g_genFFT_cmplxArg.GetY()[ip]) 
            + std::complex<double>{0.0,1.0} * (g_genFFT_cmplxAbs.GetY()[ip]*TMath::Sin(g_genFFT_cmplxArg.GetY()[ip]));
        if(TMath::Abs(freq_i/genWaveFreq-(ih*2+1))>0.001){
            Vsquare=0;
            counter++;
        }else{
            ih++;
        }
        counterAll++;
        if( ip < 20 ){
            std::cout<<freq_i<<" "<<Vsquare.real()<<std::endl;
        }
        std::complex<double> Ifit = Vsquare/((double)Npoints);
        Ifit/=Zfit;
        TComplex Ifit_root(Ifit.real(),  Ifit.imag()); //TODO what is this for?
        //TComplex Ifit_con (Ifit.real(), -Ifit.imag());
        fftc2r->SetPointComplex(ip, Ifit_root);
        //fftc2r->SetPointComplex(Npoint-ip-1, Ifit_con);
    }
        std::cout<<"Number of zero instances: "<<counter<<" / "<<counterAll<<std::endl;
    fftc2r->Transform();

    //Instantiate a fit output graph
    TGraph gVTFit(Npoints, g_squidV.GetX(), g_squidV.GetY());
    //GetPoints copies the fft content into the supplied array
    fftc2r->GetPoints(gVTFit.GetY());

    double offsetFit = TMath::MinElement(gVTFit.GetN(),gVTFit.GetY());
    for( unsigned i=0; i<gVTFit.GetN(); i++ ) gVTFit.GetY()[i] -= offsetFit;

    double maxFit = TMath::MaxElement(gVTFit.GetN(),gVTFit.GetY());
    std::cout<<"Max for gIT: "<<maxFit<<std::endl;
    std::cout<<"Point 10 before norm is: "<<gVTFit.GetY()[9]<<std::endl;
    for( unsigned i=0; i<gVTFit.GetN(); i++ ) gVTFit.GetY()[i] *= (1./maxFit);
    std::cout<<"Point 10 after norm is: "<<gVTFit.GetY()[9]<<std::endl;
    
    offsetFit = TMath::MinElement(g_funcGenV.GetN(),g_funcGenV.GetY());
    for( unsigned i=0; i<g_funcGenV.GetN(); i++ ) g_funcGenV.GetY()[i] -= offsetFit;
    maxFit = TMath::MaxElement(g_funcGenV.GetN(),g_funcGenV.GetY());
    std::cout<<"Max for gfunGen: "<<maxFit<<std::endl;
    std::cout<<"Point 10 before norm is: "<<g_funcGenV.GetY()[9]<<std::endl;
    for( unsigned i=0; i<g_funcGenV.GetN(); i++ ) g_funcGenV.GetY()[i] *= (1./maxFit);
    std::cout<<"Point 10 after norm is: "<<g_funcGenV.GetY()[9]<<std::endl;
//    if( DEBUG ){

        TCanvas cTimeFit( "cTimeFit", "cTimeFit", 600, 400 );
        cTimeFit.cd();
        g_funcGenV.SetLineColor(kRed);
        g_funcGenV.Draw("AL");
        gVTFit.SetMarkerStyle(7);
        gVTFit.Draw("Psame");
        //gVTFit.Print();
        cTimeFit.Update();
        cTimeFit.SaveAs("./fittedTime_functionGen.pdf");
        cTimeFit.SaveAs("./fittedTime_functionGen.root");

  //  }
    //222222222222222222222222222222
    //draw fitted time domain I-T
    //222222222222222222222222222222
     ih=0;
    counter= 0, counterAll = 0;
    for(int ip=0; ip<g_detFFT_cmplxArg.GetN(); ip++){
        double freq_i = g_detFFT_cmplxArg.GetX()[ip];
        std::complex<double> Zfit = Chi2.fimpediance->Z(freq_i);
        //complex<double> Vsquare = TMath::Abs(freq_i/Freq.-(ih*2+1))<0.001 ? g_genFFT_cmplxAbs.GetY()[ip] : 0;
        std::complex<double> Vsquare = g_detFFT_cmplxAbs.GetY()[ip]*TMath::Cos(g_detFFT_cmplxArg.GetY()[ip]) 
            + std::complex<double>{0.0,1.0} * (g_detFFT_cmplxAbs.GetY()[ip]*TMath::Sin(g_detFFT_cmplxArg.GetY()[ip]));
        if(TMath::Abs(freq_i/genWaveFreq-(ih*2+1))>0.001){
            Vsquare=0;
            counter++;
        }else{
            ih++;
        }
        counterAll++;
        if( ip < 20 ){
            std::cout<<freq_i<<" "<<Vsquare.real()<<std::endl;
        }
        std::complex<double> Ifit = Vsquare/((double)Npoints);
        Ifit/=Zfit;
        TComplex Ifit_root(Ifit.real(),  Ifit.imag()); //TODO what is this for?
        //TComplex Ifit_con (Ifit.real(), -Ifit.imag());
        fftc2r->SetPointComplex(ip, Ifit_root);
        //fftc2r->SetPointComplex(Npoint-ip-1, Ifit_con);
    }
    fftc2r->Transform();

    //Instantiate a fit output graph
    TGraph gITFit(Npoints, g_funcGenV.GetX(), g_squidV.GetY());
    //GetPoints copies the fft content into the supplied array
    fftc2r->GetPoints(gITFit.GetY());

    offsetFit = TMath::MinElement(gITFit.GetN(),gITFit.GetY());
    for( unsigned i=0; i<gITFit.GetN(); i++ ) gITFit.GetY()[i] -= offsetFit;

    maxFit = TMath::MaxElement(gITFit.GetN(),gITFit.GetY());
    std::cout<<"Max for gIT: "<<maxFit<<std::endl;
    std::cout<<"Point 10 before norm is: "<<gITFit.GetY()[9]<<std::endl;
    for( unsigned i=0; i<gITFit.GetN(); i++ ) gITFit.GetY()[i] *= (1./maxFit);
    std::cout<<"Point 10 after norm is: "<<gITFit.GetY()[9]<<std::endl;
    
    offsetFit = TMath::MinElement(g_squidV.GetN(),g_squidV.GetY());
    for( unsigned i=0; i<g_squidV.GetN(); i++ ) g_squidV.GetY()[i] -= offsetFit;
    maxFit = TMath::MaxElement(g_squidV.GetN(),g_squidV.GetY());
    std::cout<<"Max for gfunGen: "<<maxFit<<std::endl;
    std::cout<<"Point 10 before norm is: "<<g_squidV.GetY()[9]<<std::endl;
    for( unsigned i=0; i<g_squidV.GetN(); i++ ) g_squidV.GetY()[i] *= (1./maxFit);
    std::cout<<"Point 10 after norm is: "<<g_squidV.GetY()[9]<<std::endl;
//    if( DEBUG ){

        TCanvas cTimeFit_TES( "cTimeFit_TES", "cTimeFit_TES", 600, 400 );
        cTimeFit_TES.cd();
        g_squidV.SetLineColor(kRed);
        g_squidV.Draw("AL");
        gITFit.SetMarkerStyle(7);
        //gITFit.Print();
        cTimeFit_TES.Update();
        cTimeFit_TES.SaveAs("./fittedTime_TES.pdf");
        cTimeFit_TES.SaveAs("./fittedTime_TES.root");

  //  }
        

        return 0;

}

