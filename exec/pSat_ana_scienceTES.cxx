//Local, project
#include "num.h"
#include "rootUtils.h"
#include "util.h"

//Local, rootlib
#include <TError.h>
#include <TLine.h>
#include <TF1.h>

//std lib
#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

/*
 * Compile like: g++ -std=c++1y -Wall $(root-config --cflags --libs) linearFit.cxx -o linFit
 * Run like: ./linFit [path_to_file]
 */

//Read in data file, calculate mean values. 
int main(int argc, char *argv[]){

    //Supress root output
    gErrorIgnoreLevel = kFatal;
    std::string prefix = argv[2];

    //====================================
    // F I L E   I O 
    //====================================
    //Get data from file
    std::cout<<".\n.\n.\nReading data file: "<<argv[1]<<std::flush;
    std::vector<double> I_keithley, SQUIDCtrl_V;
    std::string datafile = argv[1];
    if( util::importData( datafile, I_keithley, SQUIDCtrl_V, 'p' ) ) return 0;
    std::cout<<"[done]"<<std::endl;

    //Sort by x
    std::cout<<"Sorting data in x, ascending.."<<std::flush;
    if( util::sortByX( I_keithley, SQUIDCtrl_V ) ) return 0;
    std::cout<<"[done]"<<std::endl;

    //====================================
    // Transform to correct units
    //====================================
    //Reflect curve in x-axis amd convert from SQUID readout output: volts to current measured by SQUID: uA
    //TODO import conversion factor from lut
    std::cout<<"Converting SQUID readout from voltage to current...";
    const double conversion = -1./0.03617; //Conversion factor for SQUID 2 taken from here: https://commons.lbl.gov/display/CMB/BlueFors+Dilution+Refrigerator#BlueForsDilutionRefrigerator-UsingQuantumDesignSQUIDsfor10mOhmTES
    std::vector<double> Imeas;
    num::transform( SQUIDCtrl_V, conversion, Imeas );
    std::cout<<"[done]"<<std::endl;
   
    //=======================================
    // Find fit ranges based on shape of IV
    //=======================================
    //Fit linear region
    std::cout<<"Using moments to set fit ranges..."<<std::flush;

    //Calculate approximations of the 1st and 2nd moments
    std::vector<double> d;
    if( num::deriv( I_keithley, Imeas, d ) ){
        std::cout<<"ERROR: could not differentiate. Exiting."<<std::endl;
        return 0;
    }
    std::string plotTitle = prefix + "00_diff.root";
    rootUtils::plot( I_keithley, d, plotTitle, "");
    std::vector<double> d2;
    if( num::deriv( I_keithley, d, d2 ) ){
        std::cout<<"ERROR: could not differentiate. Exiting."<<std::endl;
        return 0;
    }
    plotTitle = prefix + "0_diff2.root";
    rootUtils::plot( I_keithley, d2, plotTitle, "");

    //Use derivatives to select fit regions
    unsigned par_i_lo, par_i_hi;
    if( util::getParasiticRegion( d2, par_i_lo, par_i_hi ) ) return 0;
    unsigned norm_i_lo, norm_i_hi;
    if( util::getNormalRegion( d2, norm_i_lo, norm_i_hi, 0.02 ) ) return 0;
    std::cout<<"[done]"<<std::endl;

    //=======================================
    // Tether SQUID output
    //=======================================
    std::cout<<"Fitting normal region to tether SQUID..."<<std::flush;
    TF1* f = rootUtils::fit( I_keithley, Imeas, norm_i_lo, norm_i_hi, "normalFit.pdf" ); 
    double c = f->GetParameter(0);
    double m = f->GetParameter(1);

    //Subtract offset - this sets an absolute range for the SQUID (which was untethered before)
    num::removeOffset( Imeas, c );
    TLine* l_sq = rootUtils::createLine( I_keithley, m );
    TLine* l_sq_hi = rootUtils::createVertLine( Imeas, I_keithley,  norm_i_lo );
    TLine* l_sq_lo = rootUtils::createVertLine( Imeas, I_keithley, norm_i_hi );
    std::vector<TLine *> lv;
    lv.push_back( l_sq );
    lv.push_back( l_sq_lo );
    lv.push_back( l_sq_hi );

    plotTitle = prefix + "1_SQUIDtethered.pdf";
    std::string title_plot1 = "IV measurement;Current measured by sourcemeter [mA];Current measured in TES branch by SQUID [uA];";
    rootUtils::plot(I_keithley , Imeas, plotTitle, title_plot1, lv); 
    std::cout<<"[done]"<<std::endl;
   
    //================================================================
    // Convert recorded current on sourcemeter to voltage across TES 
    //================================================================
    std::cout<<"Calculating voltage across TES..."<<std::flush;
    const double shuntResistance = 0.5; //0.5 mOhm - TODO: this actual resistance still needs to be measured. 
    std::vector<double> V1;
    num::transform( I_keithley, 1e3 ); //Convert mA to uA
    V1 = num::minus(I_keithley,Imeas); //I_Shunt = I_keithley - I_TES [uA]
    num::transform( V1, shuntResistance ); //V1 = I_shunt * shunt_resistance. For now will use this as votlage across TES.  Result is in nV: V = I.R, nV[1e-9] = uA[1e-6] . mOhm [1e-3]
    //for(unsigned i=0;i<V1.size();i++) std::cout<<V1[i]<<std::endl;

    plotTitle = prefix + "2_unitsConverted.pdf";
    std::string title_plots34 = "IV measurement;Voltage across TES [nV];Current measured in TES branch by SQUID [uA];";
    rootUtils::plot(V1 , Imeas, plotTitle, title_plots34); 
    std::cout<<"[done]"<<std::endl;

    //====================================
    // Measure Rp and correct offset
    //====================================
    //Fit parasitic region and remove offset
    std::cout<<"Fitting parasitic region to measure resistance and remove offset..."<<std::flush;
    TF1* fp = rootUtils::fit( V1, Imeas, par_i_lo, par_i_hi, "parasiticFit.pdf" ); 
    double cp = fp->GetParameter(0);
    double mp = fp->GetParameter(1);
    double Rp = 1./mp;
    num::removeOffset( Imeas, cp, par_i_lo, par_i_hi ); //Remove offset for the SQUID part of the data
    TLine* lp = rootUtils::createLine( V1, mp );
    TLine* lp_hi = rootUtils::createVertLine( Imeas, V1,  par_i_lo );
    TLine* lp_lo = rootUtils::createVertLine( Imeas, V1, par_i_hi );
    std::vector<TLine *> lp_v;
    lp_v.push_back( lp );
    lp_v.push_back( lp_lo );
    lp_v.push_back( lp_hi );
    plotTitle = prefix + "3_Rp_corrected.pdf";
    rootUtils::plot(V1 , Imeas, plotTitle, title_plots34, lp_v); 
    std::cout<<"[done]"<<std::endl;
    
    //====================================
    //Fit normal region and calculate normal re
    //====================================
    std::cout<<"Fitting normal region to measure normal resistance..."<<std::flush;
    TF1* fn = rootUtils::fit( V1, Imeas, norm_i_lo, norm_i_hi, "normalFit.pdf" ); 
    double cn = fn->GetParameter(0);
    double mn = fn->GetParameter(1);
    TLine* ln = rootUtils::createLine( V1, mn );
    TLine* ln_hi = rootUtils::createVertLine( Imeas, V1,  norm_i_lo );
    TLine* ln_lo = rootUtils::createVertLine( Imeas, V1, norm_i_hi );
    std::vector<TLine *> ln_v;
    ln_v.push_back( ln );
    ln_v.push_back( ln_lo );
    ln_v.push_back( ln_hi );
    double Rn = (1./mn)-Rp;
    plotTitle = prefix + "4_normalFit.pdf";
    plotTitle = prefix + "4i_normalFit.root";
    title_plots34 = "SAT MF2 R135a-1: Pixel 155 I-V @ 90 mK;Voltage across TES [nV];Current measured in TES branch by SQUID [uA];";
    rootUtils::plot(V1 , Imeas, plotTitle, title_plots34, ln_v); 
    rootUtils::plot(V1 , Imeas, plotTitle, title_plots34); 
    std::cout<<"[done]"<<std::endl;

    //====================================
    // Fit and find Psat
    //====================================
    // Make a P=I^2.R_tes vs R_tes plot
    std::cout<<"Plotting TES power vs resistance..."<<std::flush;
   
    //Recover the TES resistance by subtracting parasitic
    std::vector<double> TES_resistance = num::divide(V1,Imeas);
    num::removeOffset( TES_resistance, Rp ); 

    //P = Imeas * I meas * TES_resistance
    std::vector<double> currentSq = num::multiply(Imeas, Imeas); //I^2
    num::transform( currentSq, 1e-3); //fW to pW
    std::vector<double> power = num::multiply( currentSq, TES_resistance ); //P = I^2 . R

    //Plot
    plotTitle = prefix + "5_PR.pdf";
    std::string title_plot5 = "Power vs resistance;TES resistance [mOhm];TES power [pW];";
    rootUtils::plot( TES_resistance, power, plotTitle, title_plot5 );
    plotTitle = prefix + "5_PR.root";
    rootUtils::plot( TES_resistance, power, plotTitle, title_plot5 );
    std::cout<<"[done]"<<std::endl;

    std::cout<<"\n=================== R E S U L T S ==================="<<std::endl;
    std::cout<<"Parasitic resistance calculated as: "<<Rp<<" [mOhm] "<<std::endl;
    std::cout<<"Normal resistance (parasitic subtracted) calculated as: "<<Rn<<" [mOhm] "<<std::endl;
    std::cout<<"[Program complete]"<<std::endl;

    return 0;

}
