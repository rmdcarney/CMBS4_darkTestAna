#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFunction.h>
#include <TF1.h>

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

/*
 * Compile like: g++ -std=c++1y -Wall $(root-config --cflags --libs) linearFit.cxx -o linFit
 * Run like: ./linFit [path_to_file]
 */
void fit(std::vector<double> &x, std::vector<double> &y){
    
    TCanvas *c = new TCanvas("c","Fitting DAC calib curve");
    c->SetGrid();

    std::vector<double> yErr(x.size());
    for( unsigned i=0; i<yErr.size(); i++){
        yErr[i] = 5e-4;
    }

    TGraphErrors *g = new TGraphErrors(x.size(), &(x[0]), &(y[0]), 0, &(yErr[0]) );
    g->Draw("a*");
    
    //Fit the graph with the predefined "pol3" function
    g->Fit("pol1");
    
    //Access the fit resuts
    TF1 *f = g->GetFunction("pol1");
    f->SetLineWidth(1);
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

    std::ifstream infile(argv[1]);
    std::string line;

    std::vector<double> Vset, Vmeas;

    std::cout<<"Reading data file: "<<argv[1]<<std::endl;
    while (std::getline(infile, line)){
        std::istringstream iss(line);
        double Vset_i, Vmeas_i;
        if( !(iss >> Vset_i >> Vmeas_i) ){ 
            std::cout<<"ERROR: could not read file. Formatting?"<<std::endl;
            return 1; 
        } //ss failed

        Vset.push_back( Vset_i );
        Vmeas.push_back( Vmeas_i );
    }

    std::cout<<Vset.size()<<", "<<Vmeas.size()<<std::endl;
    std::cout<<"Fitting... "<<std::endl;
    fit( Vset, Vmeas );
    
    return 0;

}
