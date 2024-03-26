#ifndef NUM_h
#define NUM_h

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace num{

    double getSum( std::vector<double>& data );
    double getMean( std::vector<double>& data );
    double getFrac( const double& x, const double& hi, const double& lo );
    double lerp(const double& f, const double& hi, const double& lo);
    unsigned getClosestPair(const double& x, std::vector<double>& v, std::pair<unsigned, unsigned> &i, std::pair<double,double> &p );
    void transform( std::vector<double>& data, double scalar );
    //Overloaded copy version
    void transform( std::vector<double>& data, double scalar, std::vector<double>& results);
    
    std::vector<double> add( std::vector<double>& a, std::vector<double>& b );
    std::vector<double> minus( std::vector<double>& a, std::vector<double>& b );
    std::vector<double> divide( std::vector<double>& a, std::vector<double>& b );
    std::vector<double> multiply( std::vector<double>& a, std::vector<double>& b );
    
    void removeOffset( std::vector<double>& data, double offset );
    //Overloaded with bounds
    void removeOffset( std::vector<double>& data, double offset, unsigned start, unsigned end );
    unsigned deriv( std::vector<double>& x, std::vector<double>& y, std::vector<double>& result );
}
#endif
 
