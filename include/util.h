#ifndef UTIL_h
#define UTIL_h

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

namespace util{

    unsigned sortByX( std::vector<double>& x, std::vector<double>& y );
    unsigned getParasiticRegion( std::vector<double>& d2, unsigned& lo, unsigned& hi );
    unsigned getNormalRegion( std::vector<double>& d2, unsigned& lo, unsigned& hi );

}
#endif
 


