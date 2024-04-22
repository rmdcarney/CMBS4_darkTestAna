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

    //Look up TES voltage and current given fixed bias voltage. 
    unsigned lookUpIV( const double& V_set, double& V_TES, double& I_TES, std::vector<double>& vps, std::vector<double>& vtes, std::vector<double>& ites, std::pair<unsigned,unsigned>& i ){


        if( vps.empty() || vtes.empty() || ites.empty() ){
            std::cout<<"ERROR: called LUT before filling. Exiting."<<std::endl;
            return 1;
        }
        if( vps.size() != vtes.size() || vps.size() != vtes.size() ){
            std::cout<<"ERROR: containers are different sizes. Exiting."<<std::endl;
            return 1;
        }

        //Get the closest pair of values in Vps to V_set
        std::pair<unsigned,unsigned> indices;
        std::pair<double,double> luValues;
        num::getClosestPair( V_set, vps, indices, luValues );
       
        //Look up equivalent TES voltage/current
        //NB: this assumes the vectors were sorted together. If not,  see import function and util::sortByX()
        double lowerV = vtes[indices.first];
        double upperV = vtes[indices.second];
        double lowerI = ites[indices.first];
        double upperI = ites[indices.second];
        i.first  = indices.first;
        i.second = indices.second; 

        //Linear interpolation
        double frac = num::getFrac( V_set, luValues.first, luValues.second );

        // Apply same fraction to retrive TES values
        V_TES = num::lerp( frac, upperV, lowerV );
        I_TES = num::lerp( frac, upperI, lowerI );

        std::cout<<" ilo    = "<<indices.first<<std::endl;
        std::cout<<" ihi    = "<<indices.second<<std::endl;

        std::cout<<" V_lo   = "<<lowerV<<std::endl;
        std::cout<<" V_hi   = "<<upperV<<std::endl;
        std::cout<<" I_lo   = "<<lowerI<<std::endl;
        std::cout<<" I_hi   = "<<upperI<<std::endl;
        std::cout<<" Frac   = "<<frac<<std::endl;
        std::cout<<" VTES   = "<<V_TES<<std::endl;
        std::cout<<" ITES   = "<<I_TES<<std::endl;


        return 0;
    }
    
    //LU Thermometer temperautree based on resistance recorded
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
        double lowerT = t[indices.first];
        double upperT = t[indices.second];
        double lowerR = r[indices.first];
        double upperR = r[indices.second];

        //Linear interpolation
        measuredTemp = num::lerp( num::getFrac( measuredResistance, upperR, lowerR ), upperT, lowerT );
        return 0;
    }


    unsigned importData( std::string& filename, std::vector< std::vector<double> >& v, char filetype ){
        
        std::ifstream infile(filename);
        std::string line;

        switch( filetype ){
            case 'p': { //IV or Power
                          if( v.size() < 2 ){
                              std::cout<<"ERROR: need at least two vectors to store IV data. Exiting."<<std::endl;
                              return 1;
                          }
                          while (std::getline(infile, line)){
                              std::istringstream iss(line);
                              double sqV, sqV_stdDev;
                              double na1, na2;
                              double SMU_I_i, SMU_I_RMS;
                              if( !(iss >> sqV >> sqV_stdDev >> na1 >> na2 >> SMU_I_i >> SMU_I_RMS ) ){ 
                                  std::cout<<"ERROR: could not read file. Formatting?"<<std::endl;
                                  return 1; 
                              } //ss failed

                              v[0].push_back( SMU_I_RMS );
                              v[1].push_back( sqV );
                          } break; }
            case 'c': { //Tc
                          if( v.size() < 2 ){
                              std::cout<<"ERROR: need at least two vectors to store Tc data. Exiting."<<std::endl;
                              return 1;
                          }
                          while (std::getline(infile, line)){
                              std::istringstream iss(line);
                              double r, t;
                              double na1, na2;
                              if( !(iss >> t >> na1 >> r >> na2  ) ){ 
                                  std::cout<<"ERROR: could not read file. Formatting?"<<std::endl;
                                  return 1; 
                              } //ss failed

                              v[0].push_back( t );
                              v[1].push_back( r );
                          } break; }
            case 't': { //Timing
                          if( v.size() < 3 ){
                              std::cout<<"ERROR: need at least three vectors to store timing data. Exiting."<<std::endl;
                              return 1;
                          }
                          while (std::getline(infile, line)){
                              std::istringstream iss(line);
                              double t, vGen, vSQUID;
                              double na1, na2;
                              if( !(iss >> t >> vGen >> vSQUID  ) ){ 
                                  std::cout<<"ERROR: could not read file. Formatting?"<<std::endl;
                                  return 1; 
                              } //ss failed

                              v[0].push_back( t );
                              v[1].push_back( vGen );
                              v[2].push_back( vSQUID );
                          } 
                          break; }
            case 'i': { //IV curve (saved after running psat analysis)
                          if( v.size() < 3 ){
                              std::cout<<"ERROR: need at least three vectors to store IV data. Exiting."<<std::endl;
                              return 1;
                          }
                          //Ignore the first line (comment line) 
                          std::getline(infile, line);
                          //Read in the rest of the file
                          while (std::getline(infile, line)){
                              std::istringstream iss(line);
                              double V_ps, V_TES, I_TES;
                              double na;
                              if( !(iss >> na >> V_ps >> V_TES >> I_TES  ) ){ 
                                  std::cout<<"ERROR: could not read file. Formatting?"<<std::endl;
                                  return 1; 
                              } //ss failed

                              v[0].push_back( V_ps );
                              v[1].push_back( V_TES );
                              v[2].push_back( I_TES );
                          }
                          
                          //This data is read in descending order, so needs to be reversed.
                          std::reverse( v[0].begin(), v[0].end() );
                          std::reverse( v[1].begin(), v[1].end() );
                          std::reverse( v[2].begin(), v[2].end() );
                          break; }
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
