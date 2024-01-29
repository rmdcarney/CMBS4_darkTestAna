#ifndef ROOTUTILS_h
#define ROOTUTILS_h

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFunction.h>
#include <TF1.h>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace rootUtils{

    TF1* fit(std::vector<double> &x, std::vector<double> &y, unsigned lo_i, unsigned hi_i, std::string plotTitle);
    void plot(std::vector<double> &x, std::vector<double> &y, std::string filename, std::string axesTitles, TLine* l = 0);
    TLine* createLine( std::vector<double>& data, double m );

}
#endif
 

