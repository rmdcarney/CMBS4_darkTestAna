#ifndef UTIL_h
#define UTIL_h

#include <num.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

namespace util{

    unsigned importTempCurve( std::string& filename, std::vector<double>& t, std::vector<double>& r);
    unsigned importData( std::string& filename, std::vector<double>& x, std::vector<double>& y, char filetype );
    unsigned sortByX( std::vector<double>& x, std::vector<double>& y );
    unsigned lookUpThermTemp( const double& measuredResistance, double& measuredTemp, std::vector<double>& r, std::vector<double>& t);
    unsigned getParasiticRegion( std::vector<double>& d2, unsigned& lo, unsigned& hi, double threshold = 1. );
    unsigned getNormalRegion( std::vector<double>& d2, unsigned& lo, unsigned& hi, double threshold_hi = 0.01, double threshold_lo = 0.01 );

}
#endif
 


