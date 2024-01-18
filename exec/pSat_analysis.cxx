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
TF1* fit(std::vector<double> &x, std::vector<double> &y){
    
    TCanvas *c = new TCanvas("c","Fitting DAC calib curve");
    c->SetGrid();

    std::vector<double> yErr(x.size());
    for( unsigned i=0; i<yErr.size(); i++){
        yErr[i] = 5e-4;
    }

    TGraphErrors *g = new TGraphErrors(x.size(), &(x[0]), &(y[0]), 0, &(yErr[0]) );
    g->Draw("a*");
    
    //Fit the graph with the predefined "pol3" function
    g->Fit("pol1","","",0.7,1.0);
    
    //Access the fit resuts
    TF1 *f = g->GetFunction("pol1");
    f->SetLineWidth(1);
    gStyle->SetOptFit(0);
    f->Draw("SAME");
    c->SaveAs("fit.pdf");

    return f;
}

void plot(std::vector<double> &x, std::vector<double> &y, std::string filename, TLine* l = 0){

    TCanvas *c = new TCanvas("c","PSat");
    c->SetGrid();
    
    TGraph *g = new TGraph(x.size(), &(x[0]), &(y[0]));
    g->Draw("a*");
    g->SetTitle("IV measurement;Voltage across TES [V];Current measured in TES branch by SQUID [uA];");
    gStyle->SetOptFit(0);
   l->SetLineColor(kRed);
   l->SetLineWidth(3);    
   if(l){
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

void removeOffset( std::vector<double>& data, double offset ){
    
    for(auto& element : data) element -= offset;
    return;
}

TLine* createLine( std::vector<double>& data, double m ){

    double max = *(std::max_element(data.begin(), data.end()));
    TLine* l = new TLine(0,0,max,max*m);
    return l;
}


//Read in data file, calculate mean values. 
int main(int argc, char *argv[]){

    //Get data from file
    std::ifstream infile(argv[1]);
    std::string line;
    std::vector<double> Vset, Imeas;

    std::cout<<"Reading data file: "<<argv[1]<<std::endl;
    while (std::getline(infile, line)){
        std::istringstream iss(line);
        double sqV, sqV_stdDev;
        double na1, na2;
        double SMU_I_i, SMU_I_RMS;
        double Vset_i;
        if( !(iss >> sqV >> sqV_stdDev >> na1 >> na2 >> SMU_I_i >> SMU_I_RMS ) ){ 
            std::cout<<"ERROR: could not read file. Formatting?"<<std::endl;
            return 1; 
        } //ss failed

        Vset.push_back( SMU_I_RMS );
        Imeas.push_back( sqV );
    }

    std::cout<<"Number of data points: "<<std::endl;
    std::cout<<Vset.size()<<", "<<Imeas.size()<<std::endl;
    
    //Reflect curve in x-axis amd convert from volts to uA
    const double conversion = -1./0.03617; //Conversion factor for SQUID 2 taken from here: https://commons.lbl.gov/display/CMB/BlueFors+Dilution+Refrigerator#BlueForsDilutionRefrigerator-UsingQuantumDesignSQUIDsfor10mOhmTES
    transform( Imeas, conversion );
    
    //Fit linear region
    std::cout<<"Fitting linear portion..."<<std::endl;
    TF1* f = fit( Vset, Imeas ); 
    std::cout<<"Plotted in fit.pdf"<<std::endl;
    double c = f->GetParameter(0);
    double m = f->GetParameter(1);

    //Subtract offset
    removeOffset( Imeas, c );
    TLine* l = createLine( Vset, m );
    std::cout<<"Trying to draw line at ("<<l->GetX1()<<","<<l->GetY1()<<") to ("<<l->GetX2()<<","<<l->GetY2()<<")"<<std::endl;
    plot( Vset, Imeas, "offset.pdf", l ); 
    
   

    //plot(Vset , Imeas);

    return 0;

}
