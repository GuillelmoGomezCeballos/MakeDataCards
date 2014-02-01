#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TFile.h>                  // File handle
#include <TH1D.h>                   // 1D histogram class
#endif

void extrapolate(const Double_t mass, TString xsfname = "enhancement_sm4_brww.txt") 
{ 
  vector<Double_t> massv;
  vector<Double_t> value;
 
  //
  // parse xsec file
  //
  ifstream ifs;
  ifs.open(xsfname.Data());
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    Double_t m, val;
    stringstream ss(line);
    ss >> m >> val;
    massv.push_back(m);
    value.push_back(val);
  }
  ifs.close();

  //
  // Make extrapolation
  //
  for(UInt_t im=1; im<massv.size(); im++) {    
    if(mass > massv[im-1] && mass <= massv[im]) {
      double rel = (massv[im]-mass)/(massv[im]-massv[im-1]);
      double res = (1.0 - rel)*value[im] + rel*value[im-1];
      printf("res %6.1f %10.6f\n",mass,res);
    }
  }
}
