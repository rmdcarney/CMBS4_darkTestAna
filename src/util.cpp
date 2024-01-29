#include "util.h"

namespace util{

    unsigned sortByX( std::vector<double>& x, std::vector<double>& y ){

        if( x.size() != y.size() ) return 1;

        //Make a vector of pairs to sort together
        std::vector< std::pair<double, double> > tmp;
        for(size_t i=0; i<x.size(); ++i){
            tmp.push_back( std::make_pair( y[i], x[i] ));
        }
        std::sort( std::begin(tmp), std::end(tmp),[](const std::pair<double,double> &a, const std::pair<double,double> &b){ return a.second < b.second; });

        //Put back
        for(size_t i=0; i<x.size(); ++i){
            x[i] = tmp[i].second;
            y[i] = tmp[i].first;
        }
        return 0;
    }


}
