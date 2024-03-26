#include "num.h"

namespace num{

    double getFrac( const double& x, const double& hi, const double& lo ){
        return ( (x-lo)/(hi-lo) );
    }

    double lerp(const double& f, const double& hi, const double& lo){
        return (lo * (1.0 - f)) + (hi * f);
    }

    unsigned getClosestPair(const double& x, std::vector<double>& v, std::pair<unsigned, unsigned> &i, std::pair<double,double> &p ){
        
        if (v.empty()){
            std::cout<<"ERROR: LUT in getClosestPair() empty. Exiting. "<<std::endl;
            return 1; // Can't help you if we're empty.
        }

        //Find first value (assuming vector already sorted) that is closest
        std::vector<double>::const_iterator it = std::lower_bound(v.begin(), v.end(), x);

        //Make a pair with the closest values. 
        if( it == v.end() ){
            i = std::make_pair(v.size()-1, v.size()-1);
            p = std::make_pair(v.back(), v.back());
        } else if (it == v.begin()){
            i =  std::make_pair(0, 0);
            p = std::make_pair(v.front(), v.front());
        } else {
            unsigned index = it - v.begin();
            i = std::make_pair(index-1, index);
            p = std::make_pair(*(it-1), *(it));
        }

//        std::cout<<"\n------\nPair indices are: "<<i.first<<", "<<i.second<<std::endl;
//        std::cout<<"Value is: "<<x<<std::endl;
//        std::cout<<"Pair content are: "<<p.first<<", "<<p.second<<std::endl;
        return 0;
    }

    double getSum( std::vector<double>& data ){

        double sum = 0.;
        std::for_each( data.begin(), data.end(), [&](const double x){ sum += x; });
        return sum;
    }

    double getMean( std::vector<double>& data ){

        double mean = getSum( data )/(data.size());
        return mean;
    }

    void transform( std::vector<double>& data, double scalar ){

        std::transform(data.begin(), data.end(), data.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
        return;
    }

    //Overloaded copy version
    void transform( std::vector<double>& data, double scalar, std::vector<double>& results){

        results = data;
        std::transform(results.begin(), results.end(), results.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
        return;
    }

    std::vector<double> add( std::vector<double>& a, std::vector<double>& b ){

        std::vector<double> result;
        for(unsigned i=0; i<a.size(); i++){
            result.push_back(a[i] + b[i]);
        }

        return result;
    }
    std::vector<double>  minus( std::vector<double>& a, std::vector<double>& b ){

        std::vector<double> result;
        for(unsigned i=0; i<a.size(); i++){
            result.push_back(a[i] - b[i]);
        }

        return result;
    }

    std::vector<double> divide( std::vector<double>& a, std::vector<double>& b ){

        std::vector<double> result;
        for(unsigned i=0; i<a.size(); i++){
            result.push_back(a[i]/b[i]);
        }
        return result;
    }

    std::vector<double> multiply( std::vector<double>& a, std::vector<double>& b ){

        std::vector<double> result;
        for(unsigned i=0; i<a.size(); i++){
            result.push_back(a[i]*b[i]);
        }
        return result;
    }

    void removeOffset( std::vector<double>& data, double offset ){

        for(auto& element : data) element -= offset;
        return;
    }

    //Overloaded with bounds
    void removeOffset( std::vector<double>& data, double offset, unsigned start, unsigned end ){

        for(unsigned i=0; i<data.size();i++){
            if( i<start ) continue;
            if( i>end   ) continue;
            data[i] -= offset;
        }
        return;
    }

    unsigned deriv( std::vector<double>& x, std::vector<double>& y, std::vector<double>& result ){

        if( x.size() != y.size() ){
            std::cout<<"ERROR: x and y need to be the same length."<<std::endl;
            return 1;
        }

        //Have the first dervative result be a constant
        result.push_back(0);

        //Use the first data point derviatve (which is in teh oparasitic resitantce region) to normalize the rest of the derivative plot
        //This gives us a standard way to set thresholds for the fit. 
        double dx_n = x[1]-x[0];
        double dy_n = y[1]-y[0];
        double norm = dy_n/dx_n;

        //Normalized dx
        for(unsigned i=1; i<x.size(); i++ ){
            double dy = y[i] - y[i-1];
            double dx = x[i] - x[i-1];
            double d = dy/dx;
            result.push_back( d/norm );
        }
        return 0;
    }

}
