#include "util.h"

namespace util{

    unsigned importTempCurve( std::string& filename, std::vector<double>& t, std::vector<double>& r){

        std::ifstream infile(filename);
        std::string line;

        //Import data from text file into two vectors
        while (std::getline(infile, line)){
            std::istringstream iss(line);
            double temp, resistance;
            if( !(iss >> temp >> resistance) ){ 
                std::cout<<"ERROR: could not read file. Formatting?"<<std::endl;
                return 1; 
            } //ss failed

            r.push_back( resistance );
            t.push_back( temp );
        }

        //Sort
        if( util::sortByX( r, t ) ){
            std::cout<<"ERROR: during temp curve import could not sort vectors. Exiting."<<std::endl;
            return 1;
        }

        return 0;
    }


    unsigned lookUpThermTemp( const double& measuredResistance, double& measuredTemp, std::vector<double>& r, std::vector<double>& t){

        if( r.empty() || t.empty() ){
            std::cout<<"ERROR: called LUT before filling. Try running importTempCurve first. Exiting."<<std::endl;
            return 1;
        }
        
        if( r.size() != t.size() ){
            std::cout<<"ERROR: temperature and resistance containers are different sizes. Exiting."<<std::endl;
            return 1;
        }

        //Get the closest 
        std::pair<unsigned,unsigned> indices;
        std::pair<double,double> luValues;
        num::getClosestPair( measuredResistance, r, indices, luValues );
       
        //Look up equivalent temp - NB: this assumes the vectors were sorted together on import, see import function and util::sortByX()
        double upperT = t[indices.first];
        double lowerT = t[indices.second];
        double upperR = r[indices.first];
        double lowerR = r[indices.second];

        //Linear interpolation
        measuredTemp = num::lerp( num::getFrac( measuredResistance, upperR, lowerR ), upperT, lowerT );
        return 0;
    }


    unsigned importData( std::string& filename, std::vector<double>& x, std::vector<double>& y, char filetype ){
        
        std::ifstream infile(filename);
        std::string line;

        switch( filetype ){
            case 'p': {
                          while (std::getline(infile, line)){
                              std::istringstream iss(line);
                              double sqV, sqV_stdDev;
                              double na1, na2;
                              double SMU_I_i, SMU_I_RMS;
                              if( !(iss >> sqV >> sqV_stdDev >> na1 >> na2 >> SMU_I_i >> SMU_I_RMS ) ){ 
                                  std::cout<<"ERROR: could not read file. Formatting?"<<std::endl;
                                  return 1; 
                              } //ss failed

                              x.push_back( SMU_I_RMS );
                              y.push_back( sqV );
                          } break; }
            case 't': {
                          while (std::getline(infile, line)){
                              std::istringstream iss(line);
                              double r, t;
                              double na1, na2;
                              if( !(iss >> t >> na1 >> r >> na2  ) ){ 
                                  std::cout<<"ERROR: could not read file. Formatting?"<<std::endl;
                                  return 1; 
                              } //ss failed

                              x.push_back( t );
                              y.push_back( r );
                          } break; }
        }
        return 0;
    }

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

    //Using the supplied second derivative vector, choose the upper bound of the parasitic fit region
    unsigned getParasiticRegion( std::vector<double>& d2, unsigned& lo, unsigned& hi, double threshold ){

        //Starting from the lowest x-value, find when d2 exceeds 2. Tis assumes d2 has been sorted and normalized
        unsigned tmp = 2e6;
        for(unsigned i=0; i<d2.size(); i++){
            
            if( std::fabs(d2[i]) > threshold && i>0 ){
                tmp = i-1;
                break;
            }
        }
        if( tmp > 1e6 ){
            std::cout<<"ERROR: something went wrong trying to find parasitic fit region. Second derivative did not exceed "<<threshold<<". Exiting."<<std::endl;
            return 1;
        }

        //Select the middle 90% of that range, just in case the border is weird.
        //TODO does it make sense to do a 90% here?
        lo = 0;
//        unsigned range = tmp - lo;
//        double range90 = 0.9*range;
//        unsigned adjRange = std::rint(range90);
//        unsigned diff = range-adjRange;
//        lo += unsigned( diff/2. );
//        hi = lo + adjRange;
        hi = tmp;
        
        return 0;
    }

    //Using the supplied second derivative vector, choose the upper bound of the normal fit region
    unsigned getNormalRegion( std::vector<double>& d2, unsigned& lo, unsigned& hi, double threshold_hi, double threshold_lo ){

        unsigned tmpHi(2e6), tmpLo(2e6);
        for(unsigned i=d2.size()-1; i>=0 && i<d2.size(); i--){
            
            if( std::fabs(d2[i]) < threshold_hi ){
                tmpHi = i;
                break;
            }
        }
        for(unsigned i=tmpHi-1; i>=0 && i<tmpHi; i--){

            if( std::fabs(d2[i]) > threshold_lo ){
                tmpLo = i-1;
                break;
            }
        }

        if( tmpLo > 1e6 || tmpHi > 1e6 ){
            std::cout<<"ERROR: something went wrong trying to find normal fit region. Second derivative did not exceed "<<threshold_hi<<". Exiting."<<std::endl;
            return 1;
        }

        //Select the middle 90% of that range, just in case the border is weird.
        //TODO does it make sense to do a 90% here?
        unsigned range = tmpHi - tmpLo;
        double range90 = 0.9*range;
        unsigned adjRange = std::rint(range90);
        unsigned diff = unsigned((range-adjRange)/2.);
        lo = tmpLo + diff;
        hi = tmpHi - diff;
        
        return 0;
    }
}
