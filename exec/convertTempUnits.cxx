//Local, project
#include "num.h"
#include "util.h"

//std lib
#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

//TODO: maybe make the temperature calibration curve a compile-time choice?

//Read in data file, calculate mean values. 
int main(int argc, char *argv[]){

//    std::string prefix = argv[2];

    //====================================
    // F I L E   I O 
    //====================================
    //Get data from file
    std::string thermCurve = "/Users/seacow/Dropbox/seacow_13Inch2014_LBNL/workspace/CMB_S4/thermometerCurves/R30748.txt"; 
    std::cout<<".\n.\n.\nReading thermometer curve: "<<thermCurve<<std::flush;
    std::vector<double> thermR, thermT;
    if( util::importTempCurve( thermCurve, thermT, thermR ) ) return 0;
    std::cout<<"[done]"<<std::endl;
    std::cout<<"\n\n Temp curve: "<<std::endl;
    for(unsigned i=0; i<thermT.size(); i++ ){
        std::cout<<thermR[i]<<"\t"<<thermT[i]<<std::endl;
    }

    std::cout<<".\n.\n.\nReading Tc: "<<argv[1]<<std::flush;
    std::vector<double> therm_r, SQUIDCtrl_V;
    std::string datafile = argv[1];
    if( util::importData( datafile, therm_r, SQUIDCtrl_V, 't') ) return 0;
    std::cout<<"[done]"<<std::endl;

    //Sort by x
    std::cout<<"Sorting data in x, ascending.."<<std::flush;
    if( util::sortByX( therm_r, SQUIDCtrl_V ) ) return 0;
    std::cout<<"[done]"<<std::endl;

    //====================================
    // Transform to correct units
    //===================================
    std::vector<double> therm_t; 
    for( unsigned i=0; i<therm_r.size(); i++ ){

        double t = 0.0;
        util::lookUpThermTemp( therm_r[i], t, thermR, thermT );
        therm_t.push_back(t);

    }
    for(unsigned i=0; i<therm_r.size(); i++ ){

        std::cout<<therm_t[i]<<"\t"<<SQUIDCtrl_V[i]<<std::endl;
//        std::cout<<i<<"\t"<<therm_r[i]<<"\t"<<therm_t[i]<<std::endl;

    }
    return 0;
}
