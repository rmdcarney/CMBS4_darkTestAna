//Local, project
#include "num.h"
#include "rootUtils.h"
#include "util.h"

//Local, rootlib
#include <TLine.h>
#include <TF1.h>

//std lib
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

/*
 * Compile like: g++ -std=c++1y -Wall $(root-config --cflags --libs) linearFit.cxx -o linFit
 * Run like: ./linFit [path_to_file]
 */

//Read in data file, calculate mean values. 
int main(int argc, char *argv[]){

    //====================================
    // F I L E   I O 
    //====================================
    //Get data from file
    std::ifstream infile(argv[1]);
    std::string line;
    std::vector<double> I_keithley, SQUIDCtrl_V;

    std::cout<<"Reading data file: "<<argv[1]<<std::endl;
    while (std::getline(infile, line)){
        std::istringstream iss(line);
        double sqV, sqV_stdDev;
        double na1, na2;
        double SMU_I_i, SMU_I_RMS;
        if( !(iss >> sqV >> sqV_stdDev >> na1 >> na2 >> SMU_I_i >> SMU_I_RMS ) ){ 
            std::cout<<"ERROR: could not read file. Formatting?"<<std::endl;
            return 1; 
        } //ss failed

        I_keithley.push_back( SMU_I_RMS );
        SQUIDCtrl_V.push_back( sqV );
    }

    std::cout<<"Number of data points: "<<std::endl;
    std::cout<<I_keithley.size()<<", "<<SQUIDCtrl_V.size()<<std::endl;

    //Sort by x
    if( util::sortByX( I_keithley, SQUIDCtrl_V ) ) return 0;

    //====================================
    // Transform to correct units
    //====================================
    //Reflect curve in x-axis amd convert from SQUID readout output: volts to current measured by SQUID: uA
    //TODO import conversion factor from lut
    const double conversion = -1./0.03617; //Conversion factor for SQUID 2 taken from here: https://commons.lbl.gov/display/CMB/BlueFors+Dilution+Refrigerator#BlueForsDilutionRefrigerator-UsingQuantumDesignSQUIDsfor10mOhmTES
    std::vector<double> Imeas;
    num::transform( SQUIDCtrl_V, conversion, Imeas );
   
    //====================================
    // Fit and recover PSat
    //====================================
    //Fit linear region
    // TODO filename
    std::cout<<"Fitting normal region to recover offset..."<<std::endl;

    //Calculate approximations of the 1st and 2nd moments
    std::vector<double> d;
    if( num::deriv( I_keithley, Imeas, d ) ){
        std::cout<<"ERROR: could not differentiate. Exiting."<<std::endl;
        return 0;
    }
    rootUtils::plot( I_keithley, d, "0_diff.root", "");
    std::vector<double> d2;
    if( num::deriv( I_keithley, d, d2 ) ){
        std::cout<<"ERROR: could not differentiate. Exiting."<<std::endl;
        return 0;
    }
    rootUtils::plot( I_keithley, d2, "0_diff2.root", "");

    //Use derivatives to select fit regions


    TF1* f = rootUtils::fit( I_keithley, Imeas, 1.05, 1.15, "normalFit.pdf" ); 
    double c = f->GetParameter(0);
    double m = f->GetParameter(1);

    //Subtract offset - this sets an absolute range for the SQUID (which was untethered before)
    //TODO: filename
    num::removeOffset( Imeas, c );
    TLine* l = rootUtils::createLine( I_keithley, m );
    std::cout<<"Trying to draw line at ("<<l->GetX1()<<","<<l->GetY1()<<") to ("<<l->GetX2()<<","<<l->GetY2()<<")"<<std::endl;
    std::string title_plot1 = "IV measurement;Current measured by sourcemeter [mA];Current measured in TES branch by SQUID [uA];";
    rootUtils::plot(I_keithley , Imeas, "1_offset.pdf", title_plot1, l); 
   
    //================================================================
    // Convert recorded current on sourcemeter to voltage across TES 
    //================================================================
    const double shuntResistance = 0.5; //0.5 mOhm - TODO: this actual resistance still needs to be measured. 
    std::vector<double> V1;
    num::transform( I_keithley, 1e3 ); //Convert mA to uA
    V1 = num::minus(I_keithley,Imeas); //I_Shunt = I_keithley - I_TES [uA]
    num::transform( V1, shuntResistance ); //V1 = I_shunt * shunt_resistance. For now will use this as votlage across TES.  Result is in nV: V = I.R, nV[1e-9] = uA[1e-6] . mOhm [1e-3]
    //for(unsigned i=0;i<V1.size();i++) std::cout<<V1[i]<<std::endl;

    std::cout<<"V1 size: "<<V1.size()<<std::endl;
    std::cout<<"Imeas size: "<<Imeas.size()<<std::endl;

    //TODO filename
    std::string title_plots34 = "IV measurement;Voltage across TES [nV];Current measured in TES branch by SQUID [uA];";
    rootUtils::plot(V1 , Imeas, "2_unitConversion.pdf", title_plots34); 
    
    //====================================
    // Remove Rp
    //====================================
    //Fit parasitic region and remove offset
    //TODO: find fit region
    //TODO: filename
    std::cout<<"Fitting parasitic region to measure resistance and remove offset..."<<std::endl;
    TF1* fp = rootUtils::fit( V1, Imeas, 10., 140., "parasiticFit.pdf" ); 
    double cp = fp->GetParameter(0);
    double mp = fp->GetParameter(1);
    TLine* lp = rootUtils::createLine( V1, mp );
    double Rp = 1./mp;
    std::cout<<"Parasitic resistance calculated as: "<<Rp<<" [mOhm] "<<std::endl;
    num::removeOffset( Imeas, cp, 324, 465 ); //Remove offset for the SQUID part of the data
    
    //TODO: fit region
    //TODO: filename
    //Fit normal region and calculate normal re
    std::cout<<"Fitting normal region to measure normal resistance..."<<std::endl;
    TF1* fn = rootUtils::fit( V1, Imeas, 500., 560., "normalFit.pdf" ); 
    double cn = fn->GetParameter(0);
    double mn = fn->GetParameter(1);
    TLine* ln = rootUtils::createLine( V1, mn );
    double Rn = (1./mn)-Rp;
    std::cout<<"Trying to draw line at ("<<ln->GetX1()<<","<<ln->GetY1()<<") to ("<<ln->GetX2()<<","<<ln->GetY2()<<")"<<std::endl;
    std::cout<<"Normal resistance calculated as: "<<Rn<<" [mOhm] "<<std::endl;
    
    rootUtils::plot(V1 , Imeas, "3_Rp_corrected.pdf", title_plots34, lp); 
    rootUtils::plot(V1 , Imeas, "4_normalFit.pdf", title_plots34, ln); 

    //====================================
    // Fit and find Psat
    //====================================
    // Make a P=I^2.R_tes vs R_tes plot
   
    //TODO fit region
    //TODO filename
    //Recover the TES resistance by subtracting parasitic
    std::vector<double> TES_resistance = num::divide(V1,Imeas);
    num::removeOffset( TES_resistance, Rp ); 

    //P = Imeas * I meas * TES_resistance
    std::vector<double> currentSq = num::multiply(Imeas, Imeas); //I^2
    num::transform( currentSq, 1e-3); //fW to pW
    std::vector<double> power = num::multiply( currentSq, TES_resistance ); //P = I^2 . R

    //Plot
    std::string title_plot5 = "Power vs resistance;TES resistance [mOhm];TES power [pW];";
    rootUtils::plot( TES_resistance, power, "5_PR.pdf", title_plot5 );


    return 0;

}
