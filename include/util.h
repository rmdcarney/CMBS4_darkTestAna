#ifndef UTIL_h
#define UTIL_h

#include <num.h>

#include <algorithm>
#include <complex>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

#define _USE_MATH_DEFINES //Used for pi        

namespace util{

    struct complexZ{
        
        double fpar[6];
        double fRp;
        double fR0;
        
        //Constructor
        complexZ(double (&par)[6], double Rp, double R0){
            for(int i=0; i<5; i++){
                fpar[i] = par[i];
            }
            fRp = Rp;
            fR0 = R0;
        }

        //Copy constructor
        complexZ(complexZ &a){
            complexZ(a.fpar, a.fRp, a.fR0);
        }

        //Setters
        void setPar(const double* par){
            for(int i=0; i<5; i++){
                fpar[i] = par[i];
            }
        }
        
        //Getters
        double absZ(double freq){
            return std::abs(Z(freq));
        }
        double argZ(double freq){
            return std::arg(Z(freq));
        }
        double getAbsZTF1(double *x, double *p){
            //cout<<"@"<<x[0]<<" Z("<<Z(x[0]).real()<<","<<Z(x[0]).imag()<<")"<<endl;
            //cout<<"pars {"<<fpar[0]<<","<<fpar[1]<<","<<fpar[2]<<","<<fpar[3]<<","<<fpar[4]<<","<<fpar[5]<<","<<endl;
            return absZ(x[0]);
        }
        double getArgZTF1(double *x, double *p){
            return argZ(x[0]);
        }

        //This is the fit function used //TODO check definition is correct
        std::complex<double> Z(double freq, const double*par = nullptr) const {

            const double *parInUse = (par==nullptr) ? fpar : par;
            freq*=2*M_PI;
            std::complex<double> A = parInUse[0] + std::complex<double>{0.0,1.0} * (double)(parInUse[1]*freq); //[0] -> Rp+R0(1+beta) [1]->L
            std::complex<double> B = parInUse[4] + std::complex<double>{0.0,1.0} * (double)(parInUse[2]*freq); //[2] -> tau/(L_I*R0*(2+beta))
            std::complex<double> C = 1. + std::complex<double>{0.0,1.0} * (double)(parInUse[3]*freq);          //[3] -> tau1
            std::complex<double> Z = A + 1./(B - parInUse[5]/C);			        //[4] -> (1-L_I)/ (L_I*R0*(2+beta))
                                                                                    //[5] ->g1T1/(g1T1+g1b)/(L_I*R0*(2+beta))
            return Z;
        }
    };
    
    //==================================
    // Chi2 - minimization for fitting
    //==================================
    struct complexChi2{
    
        //Constructor
        complexChi2(double* freq, double* abs, double* arg, double* errAbs, int Npoints, double (&par)[6], double Rp, double R0) :
            ffreq(freq), fabs(abs), farg(arg), ferr(errAbs), fN(Npoints){ 
               //TODO: decide if you want to do this with new or not. 
                fimpediance = new util::complexZ(par, Rp, R0);
            }

        //Destructor
        ~complexChi2(){ delete fimpediance;}

        //overload ()
        double operator() (const double *par) const {

            //double dF=20;
            double chi2 = 0; // no error considered!

            for(int i=0; i<fN/2; i++){
                double freq = ffreq[i];
                std::complex<double> Zmeasure = fabs[i]*std::cos(farg[i]) + std::complex<double>{0.0,1.0}*fabs[i]*std::sin(farg[i]);
                std::complex<double> dZ = Zmeasure - fimpediance->Z(freq,par);
                //cout<<i<<endl;
                chi2 += std::norm(dZ)/ferr[i]/ferr[i];
            }	
            return chi2;
        }

        //Bunch of stuff on the heap here, TODO: figure out memory management/appropriate containers.
        complexZ *fimpediance;
        const double* ffreq;
        const double* fabs;
        const double* farg;
        const double* ferr;
        int fN;
        double *fpar;
    };

    unsigned importTempCurve( std::string& filename, std::vector<double>& t, std::vector<double>& r);
    unsigned importData( std::string& filename, std::vector< std::vector<double> >& v, char filetype );
    unsigned sortByX( std::vector<double>& x, std::vector<double>& y );
    unsigned lookUpThermTemp( const double& measuredResistance, double& measuredTemp, std::vector<double>& r, std::vector<double>& t);
    unsigned lookUpIV( const double& V_set, double& V_TES, double& I_TES, std::vector<double>& vps, std::vector<double>& vtes, std::vector<double>& ites, std::pair<unsigned,unsigned>& i );
    unsigned getParasiticRegion( std::vector<double>& d2, unsigned& lo, unsigned& hi, double threshold = 1. );
    unsigned getNormalRegion( std::vector<double>& d2, unsigned& lo, unsigned& hi, double threshold_hi = 0.01, double threshold_lo = 0.01 );

}
#endif
 


