#ifndef UTIL_h
#define UTIL_h

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

namespace util{

    unsigned importData( std::string& filename, std::vector<double>& x, std::vector<double>& y );
    unsigned sortByX( std::vector<double>& x, std::vector<double>& y );
    unsigned getParasiticRegion( std::vector<double>& d2, unsigned& lo, unsigned& hi, double threshold = 2. );
    unsigned getNormalRegion( std::vector<double>& d2, unsigned& lo, unsigned& hi, double threshold_hi = 0.02, double threshold_lo = 0.05 );

}
#endif
 


