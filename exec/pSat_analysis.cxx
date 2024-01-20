#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFunction.h>
#include <TF1.h>
#include <TLine.h>

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

/*
 * Compile like: g++ -std=c++1y -Wall $(root-config --cflags --libs) linearFit.cxx -o linFit
 * Run like: ./linFit [path_to_file]
 */
TF1* fit(std::vector<double> &x, std::vector<double> &y, double lo, double hi, std::string plotTitle){
    
    TCanvas *c = new TCanvas("c","Fitting DAC calib curve");
    c->SetGrid();

    std::vector<double> yErr(x.size());
    for( unsigned i=0; i<yErr.size(); i++){
        yErr[i] = 1e-2;     //TODO: what should the uncertainty on the SQUID measuremnt be?
    }

    TGraphErrors *g = new TGraphErrors(x.size(), &(x[0]), &(y[0]), 0, &(yErr[0]) );
    g->Draw("a*");
    
    //Fit the graph with the predefined "pol3" function
    g->Fit("pol1","","",lo,hi);
    
    //Access the fit resuts
    TF1 *f = g->GetFunction("pol1");
    f->SetLineWidth(3);
    f->SetLineColor(kRed);
    gStyle->SetOptFit(0);
    f->Draw("SAME");
    c->SaveAs(plotTitle.c_str());

    return f;
}

void plot(std::vector<double> &x, std::vector<double> &y, std::string filename, std::string axesTitles, TLine* l = 0){

    TCanvas *c = new TCanvas("c","PSat");
    c->SetGrid();
    
    TGraph *g = new TGraph(x.size(), &(x[0]), &(y[0]));
    g->Draw("a*");
    g->SetTitle(axesTitles.c_str());
    gStyle->SetOptFit(0);
   if(l){
       l->SetLineColor(kRed);
       l->SetLineWidth(3);    
       l->Draw();
    }

    c->SaveAs(filename.c_str());
    
    return;
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

TLine* createLine( std::vector<double>& data, double m ){

    double max = *(std::max_element(data.begin(), data.end()));
    TLine* l = new TLine(0,0,max,max*m);
    return l;
}


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
    
    //====================================
    // Transform to correct units
    //====================================
    //Reflect curve in x-axis amd convert from SQUID readout output: volts to current measured by SQUID: uA
    const double conversion = -1./0.03617; //Conversion factor for SQUID 2 taken from here: https://commons.lbl.gov/display/CMB/BlueFors+Dilution+Refrigerator#BlueForsDilutionRefrigerator-UsingQuantumDesignSQUIDsfor10mOhmTES
    std::vector<double> Imeas;
    transform( SQUIDCtrl_V, conversion, Imeas );
   
    //====================================
    // Fit and recover PSat
    //====================================
    //Fit linear region
    std::cout<<"Fitting normal region to recover offset..."<<std::endl;
    TF1* f = fit( I_keithley, Imeas, 1.05, 1.15, "normalFit.pdf" ); 
    double c = f->GetParameter(0);
    double m = f->GetParameter(1);

    //Subtract offset - this sets an absolute range for the SQUID (which was untethered before)
    removeOffset( Imeas, c );
    TLine* l = createLine( I_keithley, m );
    std::cout<<"Trying to draw line at ("<<l->GetX1()<<","<<l->GetY1()<<") to ("<<l->GetX2()<<","<<l->GetY2()<<")"<<std::endl;
    std::string title_plot1 = "IV measurement;Current measured by sourcemeter [mA];Current measured in TES branch by SQUID [uA];";
    plot(I_keithley , Imeas, "1_offset.pdf", title_plot1, l); 
   
    //================================================================
    // Convert recorded current on sourcemeter to voltage across TES 
    //================================================================
    const double shuntResistance = 0.5; //0.5 mOhm - TODO: this actual resistance still needs to be measured. 
    std::vector<double> V1;
    transform( I_keithley, 1e3 ); //Convert mA to uA
    V1 = minus(I_keithley,Imeas); //I_Shunt = I_keithley - I_TES [uA]
    transform( V1, shuntResistance ); //V1 = I_shunt * shunt_resistance. For now will use this as votlage across TES.  Result is in nV: V = I.R, nV[1e-9] = uA[1e-6] . mOhm [1e-3]
    //for(unsigned i=0;i<V1.size();i++) std::cout<<V1[i]<<std::endl;

    std::cout<<"V1 size: "<<V1.size()<<std::endl;
    std::cout<<"Imeas size: "<<Imeas.size()<<std::endl;

    std::string title_plots34 = "IV measurement;Voltage across TES [nV];Current measured in TES branch by SQUID [uA];";
    plot(V1 , Imeas, "2_unitConversion.pdf", title_plots34); 
    
    //====================================
    // Remove Rp
    //====================================
    //Fit parasitic region and remove offset
    std::cout<<"Fitting parasitic region to measure resistance and remove offset..."<<std::endl;
    TF1* fp = fit( V1, Imeas, 10., 140., "parasiticFit.pdf" ); 
    double cp = fp->GetParameter(0);
    double mp = fp->GetParameter(1);
    TLine* lp = createLine( V1, mp );
    double Rp = 1./mp;
    std::cout<<"Parasitic resistance calculated as: "<<Rp<<" [mOhm] "<<std::endl;
    removeOffset( Imeas, cp, 324, 465 ); //Remove offset for the SQUID part of the data
    
    //Fit normal region and calculate normal re
    std::cout<<"Fitting normal region to measure normal resistance..."<<std::endl;
    TF1* fn = fit( V1, Imeas, 500., 560., "normalFit.pdf" ); 
    double cn = fn->GetParameter(0);
    double mn = fn->GetParameter(1);
    TLine* ln = createLine( V1, mn );
    double Rn = (1./mn)-Rp;
    std::cout<<"Trying to draw line at ("<<ln->GetX1()<<","<<ln->GetY1()<<") to ("<<ln->GetX2()<<","<<ln->GetY2()<<")"<<std::endl;
    std::cout<<"Normal resistance calculated as: "<<Rn<<" [mOhm] "<<std::endl;
    
    plot(V1 , Imeas, "3_Rp_corrected.pdf", title_plots34, lp); 
    plot(V1 , Imeas, "4_normalFit.pdf", title_plots34, ln); 

    //====================================
    // Fit and find Psat
    //====================================
    // Make a P=I^2.R_tes vs R_tes plot
   
    //Recover the TES resistance by subtracting parasitic
    std::vector<double> TES_resistance = divide(V1,Imeas);
    removeOffset( TES_resistance, Rp ); 
    std::vector<double> currentSq = multiply(Imeas, Imeas);
    transform( currentSq, 1e-3);
    std::vector<double> power = multiply( currentSq, TES_resistance );

    std::string title_plot5 = "Power vs resistance;TES resistance [mOhm];TES power [pW];";
    plot( TES_resistance, power, "5_PR.pdf", title_plot5 );


    return 0;

}
