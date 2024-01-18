#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFunction.h>
#include <TF1.h>

#include <math.h>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>



/*
 * Compile like: g++ -std=c++1y -Wall $(root-config --cflags --libs) linearFit.cxx -o linFit
 * Run like: ./linFit [path_to_file]
 */


// Errorfunction
// par[0] = Mean
// par[1] = Sigma
// par[2] = Normlization
// par[3] = Offset
#define SQRT2 1.414213562
double scurveFct(double *x, double *par) {
    return par[3] + 0.5*( 2-std::erfc( (x[0]-par[0])/(par[1]*SQRT2) ) )*par[2];
}

double reverseScurveFct(double *x, double *par) {
    return par[3] + 0.5*( std::erfc( (x[0]-par[0])/(par[1]*SQRT2) ) )*par[2];
}

void fit(std::vector<double> &x, std::vector<double> &y){
    
    TCanvas *c = new TCanvas("c","Fitting DAC calib curve");
    c->SetGrid();

    std::vector<double> yErr(x.size());
    for( unsigned i=0; i<yErr.size(); i++){
        yErr[i] = 5e-9;
    }

    TGraphErrors *g = new TGraphErrors(x.size(), &(x[0]), &(y[0]), 0, &(yErr[0]) );
    g->Draw("a*");
    
    //Define the fit function
    //TF1 *fitFunc = new TF1("scurveFct",scurveFct,0.1,0.2,4);
    TF1 *fitFunc = new TF1("reverseScurveFct",reverseScurveFct,0.1,0.2,4);
    

    //Fit the graph with the predefined "pol3" function
    //g->Fit("pol1");
    g->Fit("reverseScurveFct");
    
    //Access the fit resuts
    //TF1 *f = g->GetFunction("pol1");
    //f->SetLineWidth(1);
    gStyle->SetOptFit(0);
    c->SaveAs("out.pdf");

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

//Read in data file, calculate mean values. 
int main(int argc, char *argv[]){

    double p0 = 3.234;
    double p1 = -0.3711;

    std::ifstream infile(argv[1]);
    std::string line;

    std::vector<double> Vset, Vmeas;

    std::cout<<"Reading data file: "<<argv[1]<<std::endl;
    while (std::getline(infile, line)){
        std::istringstream iss(line);
        double temp, sqV, x0, x1;
        if( !(iss >> temp >> x0 >> sqV >> x1) ){ 
            std::cout<<"ERROR: could not read file. Formatting?"<<std::endl;
            return 1; 
        } //ss failed

        double tmp = (std::log10(temp)-p0)/p1;
        double t = std::pow(10, tmp);

        Vset.push_back( t );
        Vmeas.push_back( sqV );
    }

    std::cout<<Vset.size()<<", "<<Vmeas.size()<<std::endl;
    std::cout<<"Fitting... "<<std::endl;
    fit( Vset, Vmeas );
    
    return 0;

}
