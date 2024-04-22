#ifndef ROOTUTILS_h
#define ROOTUTILS_h

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFunction.h>
#include <TF1.h>
#include <Fit/Fitter.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TFile.h>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace rootUtils{

    TF1* fit(std::vector<double> &x, std::vector<double> &y, unsigned lo_i, unsigned hi_i, std::string plotTitle);
    TFitResultPtr sfit(std::vector<double> &x, std::vector<double> &y, unsigned lo_i, unsigned hi_i, std::string plotTitle);
    void plot(std::vector<double> &x, std::vector<double> &y, std::string filename, std::string axesTitles, TLine* l = 0);
    void plot(std::vector<double> &x, std::vector<double> &y, std::string filename, std::string axesTitles, std::vector<TLine*> l);
    TLine* createLine( std::vector<double>& data, double m );
    TLine* createVertLine( std::vector<double>& y, std::vector<double>& x, double point );

}
#endif
 

