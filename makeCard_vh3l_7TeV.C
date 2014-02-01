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
#include "/home/ceballos/releases/CMSSW_5_2_8/src/Smurf/Analysis/HWWlvlv/factors.h"
#endif

void makeCard_vh3l_7TeV(const Double_t mass) 
{
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool isFermioPhobic = false;
  bool isSM4          = false;
  double scaleFactor[2] = {1.0, 1.0};

  // cross section data
  const TString xsfname("xs_br_h_ecm7tev.txt");
  
  // reference mass points
  const Int_t nmass=13;
  const Double_t massPoints[nmass] = {110,115,120,125,130,135,140,150,160,170,180,190,200};

  // directory of input cards for reference mass points
  TString inputDir("/data/smurf/ceballos/inputLimits/limits_wh3l");
  if     (isFermioPhobic == true) {
    inputDir = "/home/ceballos/releases/CMSSW_5_2_3_patch3/src/Ana/nt_scripts/limits_wh3l/fermiophobic";
  }
  else if(isSM4 == true) {
    inputDir = "/home/ceballos/releases/CMSSW_5_2_3_patch3/src/Ana/nt_scripts/limits_wh3l/sm4";
  }

  // directory of output cards
  const TString outputDir("test");

  // card names
  vector<TString> cardnamesv;
  if     (isFermioPhobic == true) {
    cardnamesv.push_back("FF_vh3l_cut.txt");
  }
  else if(isSM4 == true) {
    cardnamesv.push_back("SM4_vh3l_cut.txt");
  }
  else {
    cardnamesv.push_back("vh3l1_shape_7TeV.txt");
    cardnamesv.push_back("vh3l2_shape_7TeV.txt");
    cardnamesv.push_back("vh3l1_cut_7TeV.txt");
    cardnamesv.push_back("vh3l2_cut_7TeV.txt");
  }

  // shape files
  vector<TString> shapenamesv;
  shapenamesv.push_back("vh3l1_input_7TeV.root");
  shapenamesv.push_back("vh3l2_input_7TeV.root");

  // type
  vector<TString> type;
  type.push_back("sssf");
  type.push_back("ossf");

  //--------------------------------------------------------------------------------------------------------------
  // Main code 
  //==============================================================================================================  
 
  TString newDir = outputDir + TString("/");
  newDir += mass;
  gSystem->mkdir(newDir,kTRUE);
  
  vector<Double_t> massv;
  vector<Double_t> ggHv;
  vector<Double_t> qqHv;
  vector<Double_t> WHv;
  vector<Double_t> ZHv;
  vector<Double_t> ttHv;
  vector<Double_t> BRHWWv;
  vector<Double_t> BRHZZv;
  vector<Double_t> BRHBBv;
  vector<Double_t> BRHTTv;
  
  //
  // parse xsec file
  //
  ifstream ifs;
  ifs.open(xsfname.Data());
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    Double_t m, ggh, qqh, wh, zh, tth, brhww, brhzz, brhbb, brhtt, brhgg,brhzg;
    stringstream ss(line);
    ss >> m >> ggh >> qqh >> wh >> zh >> tth >> brhww >> brhzz >> brhbb >> brhtt >> brhgg >> brhzg;
    massv.push_back(m);
    ggHv.push_back(ggh);
    qqHv.push_back(qqh);
    WHv.push_back(wh);
    ZHv.push_back(zh);
    ttHv.push_back(tth);
    BRHWWv.push_back(brhww);
    BRHZZv.push_back(brhzz);
    BRHBBv.push_back(brhbb);
    BRHTTv.push_back(brhtt);
  }
  ifs.close();

  //
  // Find reference mass point
  //
  Double_t refmass=0;
  if(mass <= massPoints[0]) { 
    refmass=massPoints[0];
  } else if(mass >= massPoints[nmass-1]) {
    refmass=massPoints[nmass-1];
  } else { 
    for(Int_t im=1; im<nmass; im++) {    
      if(mass > massPoints[im-1] && mass <= massPoints[im]) {
        if(mass < 0.5*(massPoints[im-1]+massPoints[im]) && (mass <= 300 || mass >= 325))
	  refmass=massPoints[im-1];
	else
	  refmass=massPoints[im];
      }
    }
  }
  assert(refmass>0);
  
  //
  // Get cross sections
  //
  Int_t imass=-1, iref=-1;
  for(UInt_t j=0; j<massv.size(); j++) {
    if(mass==massv[j])    imass=j;
    if(refmass==massv[j]) iref=j;
  }
  assert(imass>-1);
  assert(iref>-1);

  if     (isFermioPhobic == true) {
    scaleFactor[0] = 0.0;
    scaleFactor[1] = enhancementFactor(mass,2)/enhancementFactor(refmass,2); // FF BR(H->WW) enhancement factor
  }
  else if(isSM4 == true) {
    scaleFactor[0] = enhancementFactor(mass,3)/enhancementFactor(refmass,3); // SM4 BR(H->tautau) enhancement factor
    scaleFactor[1] = enhancementFactor(mass,1)/enhancementFactor(refmass,1); // SM4 BR(H->WW) enhancement factor
  }

  //
  // Read in reference card and produce new card
  //
  for(UInt_t icard=0; icard<cardnamesv.size(); icard++) {
    TString refname  = inputDir + TString("/");
    refname += refmass;
    refname += TString("/") + cardnamesv[icard];
    TString outname = newDir + TString("/") + cardnamesv[icard];
    ifstream refcard; refcard.open(refname.Data()); assert(refcard.is_open());
    ofstream outcard; outcard.open(outname.Data()); assert(outcard.is_open());
    while(getline(refcard,line)) {
      string label = line.substr(0,4);
      if(label==string("rate")) {
        string str;
	Double_t VHtt,VHww,WZ,ZZ,Wjets,Wgamma,VVV,VHtt_SM,VHww_SM;
	stringstream ss(line);
	ss >> str >> VHtt >> VHww >> WZ >> ZZ >> Wjets >> Wgamma >> VVV >> VHtt_SM >> VHww_SM;
	outcard << str << "   ";
	if((ZHv[imass]+WHv[imass]+ttHv[imass])>0) outcard << VHtt*(ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHTTv[imass]/BRHTTv[iref]*scaleFactor[0] << "   ";
	else             outcard << 0 << "   ";
	
	if((ZHv[imass]+WHv[imass]+ttHv[imass])>0) outcard << VHww*(ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHWWv[imass]/BRHWWv[iref]*scaleFactor[1] << "   ";
	else             outcard << 0 << "   ";

	outcard << WZ << "   " << ZZ << "   " << Wjets << "   " << Wgamma << "   " << VVV << "   " << VHtt_SM << "   " << VHww_SM << endl;
      } else {
        outcard << line << endl;
      }
    }
    refcard.close();
    outcard.close();
    cout << outname << " created!" << endl;
    if(shapenamesv.size()>icard) {
      TString refshapename = inputDir + TString("/");
      refshapename += refmass;
      refshapename += TString("/") + shapenamesv[icard];
      TFile refshapefile(refshapename);
      TH1D *histo_Data	     = (TH1D*)refshapefile.Get("histo_Data");        assert(histo_Data);  
      TH1D *histo_WH_htt     = (TH1D*)refshapefile.Get("histo_WH_htt");      assert(histo_WH_htt);
      TH1D *histo_WH_hww     = (TH1D*)refshapefile.Get("histo_WH_hww");      assert(histo_WH_hww);
      TH1D *histo_WZ	     = (TH1D*)refshapefile.Get("histo_WZ");          assert(histo_WZ);    
      TH1D *histo_ZZ	     = (TH1D*)refshapefile.Get("histo_ZZ");          assert(histo_ZZ);    
      TH1D *histo_Wjets      = (TH1D*)refshapefile.Get("histo_Wjets");       assert(histo_Wjets); 
      TH1D *histo_Wgamma     = (TH1D*)refshapefile.Get("histo_Wgamma");      assert(histo_Wgamma);
      TH1D *histo_VVV	     = (TH1D*)refshapefile.Get("histo_VVV");         assert(histo_VVV);
      TH1D *histo_WH_htt_SM  = (TH1D*)refshapefile.Get("histo_WH_htt_SM");   assert(histo_WH_htt_SM);
      TH1D *histo_WH_hww_SM  = (TH1D*)refshapefile.Get("histo_WH_hww_SM");   assert(histo_WH_hww_SM);
      
      TH1D *out_histo_WH_htt = (TH1D*)histo_WH_htt->Clone("histo_WH_htt");   
      out_histo_WH_htt->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHTTv[imass]/BRHTTv[iref]*scaleFactor[0]);
      
      TH1D *out_histo_WH_hww = (TH1D*)histo_WH_hww->Clone("histo_WH_hww");
      out_histo_WH_hww->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHWWv[imass]/BRHWWv[iref]*scaleFactor[1]);
      
      TH1D *out_histo_Data      = (TH1D*)histo_Data	->Clone("histo_Data");     
      TH1D *out_histo_WZ        = (TH1D*)histo_WZ	->Clone("histo_WZ");	   
      TH1D *out_histo_ZZ        = (TH1D*)histo_ZZ	->Clone("histo_ZZ");	   
      TH1D *out_histo_Wjets     = (TH1D*)histo_Wjets	->Clone("histo_Wjets");    
      TH1D *out_histo_Wgamma    = (TH1D*)histo_Wgamma	->Clone("histo_Wgamma");   
      TH1D *out_histo_VVV       = (TH1D*)histo_VVV	->Clone("histo_VVV");	   
      TH1D *out_histo_WH_htt_SM = (TH1D*)histo_WH_htt_SM->Clone("histo_WH_htt_SM");
      TH1D *out_histo_WH_hww_SM = (TH1D*)histo_WH_hww_SM->Clone("histo_WH_hww_SM");

      // sssf
      TH1D *out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Up          ;TH1D *histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Up;    
      TH1D *out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Down 	     ;TH1D *histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Down; 
      TH1D *out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Up   	     ;TH1D *histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Up;   
      TH1D *out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Down 	     ;TH1D *histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Down; 
      TH1D *out_histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Up   	   	     ;TH1D *histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Up;	    
      TH1D *out_histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Down 	   	     ;TH1D *histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Down;	    
      TH1D *out_histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Up   	   	     ;TH1D *histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Up;	    
      TH1D *out_histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Down 	   	     ;TH1D *histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Down;	    
      TH1D *out_histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Up     	     ;TH1D *histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Up;     
      TH1D *out_histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Down   	     ;TH1D *histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Down;   
      TH1D *out_histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Up   	     ;TH1D *histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Up;   
      TH1D *out_histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Down 	     ;TH1D *histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Down; 
      TH1D *out_histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Up  		     ;TH1D *histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Up;  
      TH1D *out_histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Down       	     ;TH1D *histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Down;       
      TH1D *out_histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Up    ;TH1D *histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Up;   
      TH1D *out_histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Down  ;TH1D *histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Down; 
      TH1D *out_histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Up    ;TH1D *histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Up;   
      TH1D *out_histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Down  ;TH1D *histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Down; 
      // ossf
      TH1D *out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Up          ;TH1D *histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Up;    
      TH1D *out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Down 	     ;TH1D *histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Down; 
      TH1D *out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Up   	     ;TH1D *histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Up;   
      TH1D *out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Down 	     ;TH1D *histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Down; 
      TH1D *out_histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Up   	   	     ;TH1D *histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Up;	    
      TH1D *out_histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Down 	   	     ;TH1D *histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Down;	    
      TH1D *out_histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Up   	   	     ;TH1D *histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Up;	    
      TH1D *out_histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Down 	   	     ;TH1D *histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Down;	    
      TH1D *out_histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Up     	     ;TH1D *histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Up;     
      TH1D *out_histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Down   	     ;TH1D *histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Down;   
      TH1D *out_histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Up   	     ;TH1D *histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Up;   
      TH1D *out_histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Down 	     ;TH1D *histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Down; 
      TH1D *out_histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Up  		     ;TH1D *histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Up;  
      TH1D *out_histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Down       	     ;TH1D *histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Down;       
      TH1D *out_histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Up    ;TH1D *histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Up;   
      TH1D *out_histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Down  ;TH1D *histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Down; 
      TH1D *out_histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Up    ;TH1D *histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Up;   
      TH1D *out_histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Down  ;TH1D *histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Down; 
      if       (type[icard] == "sssf"){
       histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Up = (TH1D*)refshapefile.Get("histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Up");assert(histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Up);
       histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Down = (TH1D*)refshapefile.Get("histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Down");assert(histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Down);
       histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Up = (TH1D*)refshapefile.Get("histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Up");assert(histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Up);
       histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Down = (TH1D*)refshapefile.Get("histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Down");assert(histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Down);
       histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Up = (TH1D*)refshapefile.Get("histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Up");assert(histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Up);
       histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Down = (TH1D*)refshapefile.Get("histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Down");assert(histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Down);
       histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Up = (TH1D*)refshapefile.Get("histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Up");assert(histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Up);
       histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Down = (TH1D*)refshapefile.Get("histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Down");assert(histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Down);
       histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Up = (TH1D*)refshapefile.Get("histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Up");assert(histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Up);
       histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Down = (TH1D*)refshapefile.Get("histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Down");assert(histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Down);
       histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Up = (TH1D*)refshapefile.Get("histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Up");assert(histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Up);
       histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Down = (TH1D*)refshapefile.Get("histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Down");assert(histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Down);
       histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Up = (TH1D*)refshapefile.Get("histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Up");assert(histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Up);
       histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Down = (TH1D*)refshapefile.Get("histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Down");assert(histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Down);
       histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Up = (TH1D*)refshapefile.Get("histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Up");assert(histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Up);
       histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Down = (TH1D*)refshapefile.Get("histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Down");assert(histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Down);
       histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Up = (TH1D*)refshapefile.Get("histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Up");assert(histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Up);
       histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Down = (TH1D*)refshapefile.Get("histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Down");assert(histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Down);
      } else if(type[icard] == "ossf"){
       histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Up = (TH1D*)refshapefile.Get("histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Up");assert(histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Up);
       histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Down = (TH1D*)refshapefile.Get("histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Down");assert(histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Down);
       histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Up = (TH1D*)refshapefile.Get("histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Up");assert(histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Up);
       histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Down = (TH1D*)refshapefile.Get("histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Down");assert(histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Down);
       histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Up = (TH1D*)refshapefile.Get("histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Up");assert(histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Up);
       histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Down = (TH1D*)refshapefile.Get("histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Down");assert(histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Down);
       histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Up = (TH1D*)refshapefile.Get("histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Up");assert(histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Up);
       histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Down = (TH1D*)refshapefile.Get("histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Down");assert(histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Down);
       histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Up = (TH1D*)refshapefile.Get("histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Up");assert(histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Up);
       histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Down = (TH1D*)refshapefile.Get("histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Down");assert(histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Down);
       histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Up = (TH1D*)refshapefile.Get("histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Up");assert(histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Up);
       histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Down = (TH1D*)refshapefile.Get("histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Down");assert(histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Down);
       histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Up = (TH1D*)refshapefile.Get("histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Up");assert(histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Up);
       histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Down = (TH1D*)refshapefile.Get("histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Down");assert(histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Down);
       histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Up = (TH1D*)refshapefile.Get("histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Up");assert(histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Up);
       histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Down = (TH1D*)refshapefile.Get("histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Down");assert(histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Down);
       histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Up = (TH1D*)refshapefile.Get("histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Up");assert(histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Up);
       histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Down = (TH1D*)refshapefile.Get("histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Down");assert(histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Down);
      }
      TH1D *histo_WH_htt_CMS_MVALepEffBoundingUp = (TH1D*)refshapefile.Get("histo_WH_htt_CMS_MVALepEffBoundingUp");assert(histo_WH_htt_CMS_MVALepEffBoundingUp);
      TH1D *histo_WH_htt_CMS_MVALepEffBoundingDown = (TH1D*)refshapefile.Get("histo_WH_htt_CMS_MVALepEffBoundingDown");assert(histo_WH_htt_CMS_MVALepEffBoundingDown);
      TH1D *histo_WH_hww_CMS_MVALepEffBoundingUp = (TH1D*)refshapefile.Get("histo_WH_hww_CMS_MVALepEffBoundingUp");assert(histo_WH_hww_CMS_MVALepEffBoundingUp);
      TH1D *histo_WH_hww_CMS_MVALepEffBoundingDown = (TH1D*)refshapefile.Get("histo_WH_hww_CMS_MVALepEffBoundingDown");assert(histo_WH_hww_CMS_MVALepEffBoundingDown);
      TH1D *histo_WZ_CMS_MVALepEffBoundingUp = (TH1D*)refshapefile.Get("histo_WZ_CMS_MVALepEffBoundingUp");assert(histo_WZ_CMS_MVALepEffBoundingUp);
      TH1D *histo_WZ_CMS_MVALepEffBoundingDown = (TH1D*)refshapefile.Get("histo_WZ_CMS_MVALepEffBoundingDown");assert(histo_WZ_CMS_MVALepEffBoundingDown);
      TH1D *histo_ZZ_CMS_MVALepEffBoundingUp = (TH1D*)refshapefile.Get("histo_ZZ_CMS_MVALepEffBoundingUp");assert(histo_ZZ_CMS_MVALepEffBoundingUp);
      TH1D *histo_ZZ_CMS_MVALepEffBoundingDown = (TH1D*)refshapefile.Get("histo_ZZ_CMS_MVALepEffBoundingDown");assert(histo_ZZ_CMS_MVALepEffBoundingDown);
      TH1D *histo_Wgamma_CMS_MVALepEffBoundingUp = (TH1D*)refshapefile.Get("histo_Wgamma_CMS_MVALepEffBoundingUp");assert(histo_Wgamma_CMS_MVALepEffBoundingUp);
      TH1D *histo_Wgamma_CMS_MVALepEffBoundingDown = (TH1D*)refshapefile.Get("histo_Wgamma_CMS_MVALepEffBoundingDown");assert(histo_Wgamma_CMS_MVALepEffBoundingDown);
      TH1D *histo_VVV_CMS_MVALepEffBoundingUp = (TH1D*)refshapefile.Get("histo_VVV_CMS_MVALepEffBoundingUp");assert(histo_VVV_CMS_MVALepEffBoundingUp);
      TH1D *histo_VVV_CMS_MVALepEffBoundingDown = (TH1D*)refshapefile.Get("histo_VVV_CMS_MVALepEffBoundingDown");assert(histo_VVV_CMS_MVALepEffBoundingDown);
      TH1D *histo_WH_htt_SM_CMS_MVALepEffBoundingUp = (TH1D*)refshapefile.Get("histo_WH_htt_SM_CMS_MVALepEffBoundingUp");assert(histo_WH_htt_SM_CMS_MVALepEffBoundingUp);
      TH1D *histo_WH_htt_SM_CMS_MVALepEffBoundingDown = (TH1D*)refshapefile.Get("histo_WH_htt_SM_CMS_MVALepEffBoundingDown");assert(histo_WH_htt_SM_CMS_MVALepEffBoundingDown);
      TH1D *histo_WH_hww_SM_CMS_MVALepEffBoundingUp = (TH1D*)refshapefile.Get("histo_WH_hww_SM_CMS_MVALepEffBoundingUp");assert(histo_WH_hww_SM_CMS_MVALepEffBoundingUp);
      TH1D *histo_WH_hww_SM_CMS_MVALepEffBoundingDown = (TH1D*)refshapefile.Get("histo_WH_hww_SM_CMS_MVALepEffBoundingDown");assert(histo_WH_hww_SM_CMS_MVALepEffBoundingDown);
      TH1D *histo_WH_htt_CMS_MVALepResBoundingUp = (TH1D*)refshapefile.Get("histo_WH_htt_CMS_MVALepResBoundingUp");assert(histo_WH_htt_CMS_MVALepResBoundingUp);
      TH1D *histo_WH_htt_CMS_MVALepResBoundingDown = (TH1D*)refshapefile.Get("histo_WH_htt_CMS_MVALepResBoundingDown");assert(histo_WH_htt_CMS_MVALepResBoundingDown);
      TH1D *histo_WH_hww_CMS_MVALepResBoundingUp = (TH1D*)refshapefile.Get("histo_WH_hww_CMS_MVALepResBoundingUp");assert(histo_WH_hww_CMS_MVALepResBoundingUp);
      TH1D *histo_WH_hww_CMS_MVALepResBoundingDown = (TH1D*)refshapefile.Get("histo_WH_hww_CMS_MVALepResBoundingDown");assert(histo_WH_hww_CMS_MVALepResBoundingDown);
      TH1D *histo_WZ_CMS_MVALepResBoundingUp = (TH1D*)refshapefile.Get("histo_WZ_CMS_MVALepResBoundingUp");assert(histo_WZ_CMS_MVALepResBoundingUp);
      TH1D *histo_WZ_CMS_MVALepResBoundingDown = (TH1D*)refshapefile.Get("histo_WZ_CMS_MVALepResBoundingDown");assert(histo_WZ_CMS_MVALepResBoundingDown);
      TH1D *histo_ZZ_CMS_MVALepResBoundingUp = (TH1D*)refshapefile.Get("histo_ZZ_CMS_MVALepResBoundingUp");assert(histo_ZZ_CMS_MVALepResBoundingUp);
      TH1D *histo_ZZ_CMS_MVALepResBoundingDown = (TH1D*)refshapefile.Get("histo_ZZ_CMS_MVALepResBoundingDown");assert(histo_ZZ_CMS_MVALepResBoundingDown);
      TH1D *histo_Wgamma_CMS_MVALepResBoundingUp = (TH1D*)refshapefile.Get("histo_Wgamma_CMS_MVALepResBoundingUp");assert(histo_Wgamma_CMS_MVALepResBoundingUp);
      TH1D *histo_Wgamma_CMS_MVALepResBoundingDown = (TH1D*)refshapefile.Get("histo_Wgamma_CMS_MVALepResBoundingDown");assert(histo_Wgamma_CMS_MVALepResBoundingDown);
      TH1D *histo_VVV_CMS_MVALepResBoundingUp = (TH1D*)refshapefile.Get("histo_VVV_CMS_MVALepResBoundingUp");assert(histo_VVV_CMS_MVALepResBoundingUp);
      TH1D *histo_VVV_CMS_MVALepResBoundingDown = (TH1D*)refshapefile.Get("histo_VVV_CMS_MVALepResBoundingDown");assert(histo_VVV_CMS_MVALepResBoundingDown);
      TH1D *histo_WH_htt_SM_CMS_MVALepResBoundingUp = (TH1D*)refshapefile.Get("histo_WH_htt_SM_CMS_MVALepResBoundingUp");assert(histo_WH_htt_SM_CMS_MVALepResBoundingUp);
      TH1D *histo_WH_htt_SM_CMS_MVALepResBoundingDown = (TH1D*)refshapefile.Get("histo_WH_htt_SM_CMS_MVALepResBoundingDown");assert(histo_WH_htt_SM_CMS_MVALepResBoundingDown);
      TH1D *histo_WH_hww_SM_CMS_MVALepResBoundingUp = (TH1D*)refshapefile.Get("histo_WH_hww_SM_CMS_MVALepResBoundingUp");assert(histo_WH_hww_SM_CMS_MVALepResBoundingUp);
      TH1D *histo_WH_hww_SM_CMS_MVALepResBoundingDown = (TH1D*)refshapefile.Get("histo_WH_hww_SM_CMS_MVALepResBoundingDown");assert(histo_WH_hww_SM_CMS_MVALepResBoundingDown);
      TH1D *histo_WH_htt_CMS_MVAMETResBoundingUp = (TH1D*)refshapefile.Get("histo_WH_htt_CMS_MVAMETResBoundingUp");assert(histo_WH_htt_CMS_MVAMETResBoundingUp);
      TH1D *histo_WH_htt_CMS_MVAMETResBoundingDown = (TH1D*)refshapefile.Get("histo_WH_htt_CMS_MVAMETResBoundingDown");assert(histo_WH_htt_CMS_MVAMETResBoundingDown);
      TH1D *histo_WH_hww_CMS_MVAMETResBoundingUp = (TH1D*)refshapefile.Get("histo_WH_hww_CMS_MVAMETResBoundingUp");assert(histo_WH_hww_CMS_MVAMETResBoundingUp);
      TH1D *histo_WH_hww_CMS_MVAMETResBoundingDown = (TH1D*)refshapefile.Get("histo_WH_hww_CMS_MVAMETResBoundingDown");assert(histo_WH_hww_CMS_MVAMETResBoundingDown);
      TH1D *histo_WZ_CMS_MVAMETResBoundingUp = (TH1D*)refshapefile.Get("histo_WZ_CMS_MVAMETResBoundingUp");assert(histo_WZ_CMS_MVAMETResBoundingUp);
      TH1D *histo_WZ_CMS_MVAMETResBoundingDown = (TH1D*)refshapefile.Get("histo_WZ_CMS_MVAMETResBoundingDown");assert(histo_WZ_CMS_MVAMETResBoundingDown);
      TH1D *histo_ZZ_CMS_MVAMETResBoundingUp = (TH1D*)refshapefile.Get("histo_ZZ_CMS_MVAMETResBoundingUp");assert(histo_ZZ_CMS_MVAMETResBoundingUp);
      TH1D *histo_ZZ_CMS_MVAMETResBoundingDown = (TH1D*)refshapefile.Get("histo_ZZ_CMS_MVAMETResBoundingDown");assert(histo_ZZ_CMS_MVAMETResBoundingDown);
      TH1D *histo_Wgamma_CMS_MVAMETResBoundingUp = (TH1D*)refshapefile.Get("histo_Wgamma_CMS_MVAMETResBoundingUp");assert(histo_Wgamma_CMS_MVAMETResBoundingUp);
      TH1D *histo_Wgamma_CMS_MVAMETResBoundingDown = (TH1D*)refshapefile.Get("histo_Wgamma_CMS_MVAMETResBoundingDown");assert(histo_Wgamma_CMS_MVAMETResBoundingDown);
      TH1D *histo_VVV_CMS_MVAMETResBoundingUp = (TH1D*)refshapefile.Get("histo_VVV_CMS_MVAMETResBoundingUp");assert(histo_VVV_CMS_MVAMETResBoundingUp);
      TH1D *histo_VVV_CMS_MVAMETResBoundingDown = (TH1D*)refshapefile.Get("histo_VVV_CMS_MVAMETResBoundingDown");assert(histo_VVV_CMS_MVAMETResBoundingDown);
      TH1D *histo_WH_htt_SM_CMS_MVAMETResBoundingUp = (TH1D*)refshapefile.Get("histo_WH_htt_SM_CMS_MVAMETResBoundingUp");assert(histo_WH_htt_SM_CMS_MVAMETResBoundingUp);
      TH1D *histo_WH_htt_SM_CMS_MVAMETResBoundingDown = (TH1D*)refshapefile.Get("histo_WH_htt_SM_CMS_MVAMETResBoundingDown");assert(histo_WH_htt_SM_CMS_MVAMETResBoundingDown);
      TH1D *histo_WH_hww_SM_CMS_MVAMETResBoundingUp = (TH1D*)refshapefile.Get("histo_WH_hww_SM_CMS_MVAMETResBoundingUp");assert(histo_WH_hww_SM_CMS_MVAMETResBoundingUp);
      TH1D *histo_WH_hww_SM_CMS_MVAMETResBoundingDown = (TH1D*)refshapefile.Get("histo_WH_hww_SM_CMS_MVAMETResBoundingDown");assert(histo_WH_hww_SM_CMS_MVAMETResBoundingDown);
      TH1D *histo_WH_htt_CMS_MVAJESBoundingUp = (TH1D*)refshapefile.Get("histo_WH_htt_CMS_MVAJESBoundingUp");assert(histo_WH_htt_CMS_MVAJESBoundingUp);
      TH1D *histo_WH_htt_CMS_MVAJESBoundingDown = (TH1D*)refshapefile.Get("histo_WH_htt_CMS_MVAJESBoundingDown");assert(histo_WH_htt_CMS_MVAJESBoundingDown);
      TH1D *histo_WH_hww_CMS_MVAJESBoundingUp = (TH1D*)refshapefile.Get("histo_WH_hww_CMS_MVAJESBoundingUp");assert(histo_WH_hww_CMS_MVAJESBoundingUp);
      TH1D *histo_WH_hww_CMS_MVAJESBoundingDown = (TH1D*)refshapefile.Get("histo_WH_hww_CMS_MVAJESBoundingDown");assert(histo_WH_hww_CMS_MVAJESBoundingDown);
      TH1D *histo_WZ_CMS_MVAJESBoundingUp = (TH1D*)refshapefile.Get("histo_WZ_CMS_MVAJESBoundingUp");assert(histo_WZ_CMS_MVAJESBoundingUp);
      TH1D *histo_WZ_CMS_MVAJESBoundingDown = (TH1D*)refshapefile.Get("histo_WZ_CMS_MVAJESBoundingDown");assert(histo_WZ_CMS_MVAJESBoundingDown);
      TH1D *histo_ZZ_CMS_MVAJESBoundingUp = (TH1D*)refshapefile.Get("histo_ZZ_CMS_MVAJESBoundingUp");assert(histo_ZZ_CMS_MVAJESBoundingUp);
      TH1D *histo_ZZ_CMS_MVAJESBoundingDown = (TH1D*)refshapefile.Get("histo_ZZ_CMS_MVAJESBoundingDown");assert(histo_ZZ_CMS_MVAJESBoundingDown);
      TH1D *histo_Wgamma_CMS_MVAJESBoundingUp = (TH1D*)refshapefile.Get("histo_Wgamma_CMS_MVAJESBoundingUp");assert(histo_Wgamma_CMS_MVAJESBoundingUp);
      TH1D *histo_Wgamma_CMS_MVAJESBoundingDown = (TH1D*)refshapefile.Get("histo_Wgamma_CMS_MVAJESBoundingDown");assert(histo_Wgamma_CMS_MVAJESBoundingDown);
      TH1D *histo_VVV_CMS_MVAJESBoundingUp = (TH1D*)refshapefile.Get("histo_VVV_CMS_MVAJESBoundingUp");assert(histo_VVV_CMS_MVAJESBoundingUp);
      TH1D *histo_VVV_CMS_MVAJESBoundingDown = (TH1D*)refshapefile.Get("histo_VVV_CMS_MVAJESBoundingDown");assert(histo_VVV_CMS_MVAJESBoundingDown);
      TH1D *histo_WH_htt_SM_CMS_MVAJESBoundingUp = (TH1D*)refshapefile.Get("histo_WH_htt_SM_CMS_MVAJESBoundingUp");assert(histo_WH_htt_SM_CMS_MVAJESBoundingUp);
      TH1D *histo_WH_htt_SM_CMS_MVAJESBoundingDown = (TH1D*)refshapefile.Get("histo_WH_htt_SM_CMS_MVAJESBoundingDown");assert(histo_WH_htt_SM_CMS_MVAJESBoundingDown);
      TH1D *histo_WH_hww_SM_CMS_MVAJESBoundingUp = (TH1D*)refshapefile.Get("histo_WH_hww_SM_CMS_MVAJESBoundingUp");assert(histo_WH_hww_SM_CMS_MVAJESBoundingUp);
      TH1D *histo_WH_hww_SM_CMS_MVAJESBoundingDown = (TH1D*)refshapefile.Get("histo_WH_hww_SM_CMS_MVAJESBoundingDown");assert(histo_WH_hww_SM_CMS_MVAJESBoundingDown);
      TH1D *histo_Wjets_CMS_MVAWBoundingUp = (TH1D*)refshapefile.Get("histo_Wjets_CMS_MVAWBoundingUp");assert(histo_Wjets_CMS_MVAWBoundingUp);
      TH1D *histo_Wjets_CMS_MVAWBoundingDown = (TH1D*)refshapefile.Get("histo_Wjets_CMS_MVAWBoundingDown");assert(histo_Wjets_CMS_MVAWBoundingDown);
      TH1D *histo_WZ_CMS_MVAWZNLOBoundingUp = (TH1D*)refshapefile.Get("histo_WZ_CMS_MVAWZNLOBoundingUp");assert(histo_WZ_CMS_MVAWZNLOBoundingUp);
      TH1D *histo_WZ_CMS_MVAWZNLOBoundingDown = (TH1D*)refshapefile.Get("histo_WZ_CMS_MVAWZNLOBoundingDown");assert(histo_WZ_CMS_MVAWZNLOBoundingDown);
      TH1D *histo_ZZ_CMS_MVAZZNLOBoundingUp = (TH1D*)refshapefile.Get("histo_ZZ_CMS_MVAZZNLOBoundingUp");assert(histo_ZZ_CMS_MVAZZNLOBoundingUp);
      TH1D *histo_ZZ_CMS_MVAZZNLOBoundingDown = (TH1D*)refshapefile.Get("histo_ZZ_CMS_MVAZZNLOBoundingDown");assert(histo_ZZ_CMS_MVAZZNLOBoundingDown);

      if       (type[icard] == "sssf"){
       out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Up = (TH1D*)histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Up->Clone("histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Up");
       out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Down = (TH1D*)histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Down->Clone("histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Down");
       out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Up = (TH1D*)histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Up->Clone("histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Up");
       out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Down = (TH1D*)histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Down->Clone("histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Down");
       out_histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Up = (TH1D*)histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Up->Clone("histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Up");
       out_histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Down = (TH1D*)histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Down->Clone("histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Down");
       out_histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Up = (TH1D*)histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Up->Clone("histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Up");
       out_histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Down = (TH1D*)histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Down->Clone("histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Down");
       out_histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Up = (TH1D*)histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Up->Clone("histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Up");
       out_histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Down = (TH1D*)histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Down->Clone("histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Down");
       out_histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Up = (TH1D*)histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Up->Clone("histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Up");
       out_histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Down = (TH1D*)histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Down->Clone("histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Down");
       out_histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Up = (TH1D*)histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Up->Clone("histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Up");
       out_histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Down = (TH1D*)histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Down->Clone("histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Down");
       out_histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Up = (TH1D*)histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Up->Clone("histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Up");
       out_histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Down = (TH1D*)histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Down->Clone("histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Down");
       out_histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Up = (TH1D*)histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Up->Clone("histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Up");
       out_histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Down = (TH1D*)histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Down->Clone("histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Down");
      } else if(type[icard] == "ossf"){
       out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Up = (TH1D*)histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Up->Clone("histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Up");
       out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Down = (TH1D*)histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Down->Clone("histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Down");
       out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Up = (TH1D*)histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Up->Clone("histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Up");
       out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Down = (TH1D*)histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Down->Clone("histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Down");
       out_histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Up = (TH1D*)histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Up->Clone("histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Up");
       out_histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Down = (TH1D*)histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Down->Clone("histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Down");
       out_histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Up = (TH1D*)histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Up->Clone("histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Up");
       out_histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Down = (TH1D*)histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Down->Clone("histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Down");
       out_histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Up = (TH1D*)histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Up->Clone("histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Up");
       out_histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Down = (TH1D*)histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Down->Clone("histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Down");
       out_histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Up = (TH1D*)histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Up->Clone("histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Up");
       out_histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Down = (TH1D*)histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Down->Clone("histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Down");
       out_histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Up = (TH1D*)histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Up->Clone("histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Up");
       out_histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Down = (TH1D*)histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Down->Clone("histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Down");
       out_histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Up = (TH1D*)histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Up->Clone("histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Up");
       out_histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Down = (TH1D*)histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Down->Clone("histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Down");
       out_histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Up = (TH1D*)histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Up->Clone("histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Up");
       out_histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Down = (TH1D*)histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Down->Clone("histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Down");
      }
      TH1D *out_histo_WH_htt_CMS_MVALepEffBoundingUp = (TH1D*)histo_WH_htt_CMS_MVALepEffBoundingUp->Clone("histo_WH_htt_CMS_MVALepEffBoundingUp");
      TH1D *out_histo_WH_htt_CMS_MVALepEffBoundingDown = (TH1D*)histo_WH_htt_CMS_MVALepEffBoundingDown->Clone("histo_WH_htt_CMS_MVALepEffBoundingDown");
      TH1D *out_histo_WH_hww_CMS_MVALepEffBoundingUp = (TH1D*)histo_WH_hww_CMS_MVALepEffBoundingUp->Clone("histo_WH_hww_CMS_MVALepEffBoundingUp");
      TH1D *out_histo_WH_hww_CMS_MVALepEffBoundingDown = (TH1D*)histo_WH_hww_CMS_MVALepEffBoundingDown->Clone("histo_WH_hww_CMS_MVALepEffBoundingDown");
      TH1D *out_histo_WZ_CMS_MVALepEffBoundingUp = (TH1D*)histo_WZ_CMS_MVALepEffBoundingUp->Clone("histo_WZ_CMS_MVALepEffBoundingUp");
      TH1D *out_histo_WZ_CMS_MVALepEffBoundingDown = (TH1D*)histo_WZ_CMS_MVALepEffBoundingDown->Clone("histo_WZ_CMS_MVALepEffBoundingDown");
      TH1D *out_histo_ZZ_CMS_MVALepEffBoundingUp = (TH1D*)histo_ZZ_CMS_MVALepEffBoundingUp->Clone("histo_ZZ_CMS_MVALepEffBoundingUp");
      TH1D *out_histo_ZZ_CMS_MVALepEffBoundingDown = (TH1D*)histo_ZZ_CMS_MVALepEffBoundingDown->Clone("histo_ZZ_CMS_MVALepEffBoundingDown");
      TH1D *out_histo_Wgamma_CMS_MVALepEffBoundingUp = (TH1D*)histo_Wgamma_CMS_MVALepEffBoundingUp->Clone("histo_Wgamma_CMS_MVALepEffBoundingUp");
      TH1D *out_histo_Wgamma_CMS_MVALepEffBoundingDown = (TH1D*)histo_Wgamma_CMS_MVALepEffBoundingDown->Clone("histo_Wgamma_CMS_MVALepEffBoundingDown");
      TH1D *out_histo_VVV_CMS_MVALepEffBoundingUp = (TH1D*)histo_VVV_CMS_MVALepEffBoundingUp->Clone("histo_VVV_CMS_MVALepEffBoundingUp");
      TH1D *out_histo_VVV_CMS_MVALepEffBoundingDown = (TH1D*)histo_VVV_CMS_MVALepEffBoundingDown->Clone("histo_VVV_CMS_MVALepEffBoundingDown");
      TH1D *out_histo_WH_htt_SM_CMS_MVALepEffBoundingUp = (TH1D*)histo_WH_htt_SM_CMS_MVALepEffBoundingUp->Clone("histo_WH_htt_SM_CMS_MVALepEffBoundingUp");
      TH1D *out_histo_WH_htt_SM_CMS_MVALepEffBoundingDown = (TH1D*)histo_WH_htt_SM_CMS_MVALepEffBoundingDown->Clone("histo_WH_htt_SM_CMS_MVALepEffBoundingDown");
      TH1D *out_histo_WH_hww_SM_CMS_MVALepEffBoundingUp = (TH1D*)histo_WH_hww_SM_CMS_MVALepEffBoundingUp->Clone("histo_WH_hww_SM_CMS_MVALepEffBoundingUp");
      TH1D *out_histo_WH_hww_SM_CMS_MVALepEffBoundingDown = (TH1D*)histo_WH_hww_SM_CMS_MVALepEffBoundingDown->Clone("histo_WH_hww_SM_CMS_MVALepEffBoundingDown");
      TH1D *out_histo_WH_htt_CMS_MVALepResBoundingUp = (TH1D*)histo_WH_htt_CMS_MVALepResBoundingUp->Clone("histo_WH_htt_CMS_MVALepResBoundingUp");
      TH1D *out_histo_WH_htt_CMS_MVALepResBoundingDown = (TH1D*)histo_WH_htt_CMS_MVALepResBoundingDown->Clone("histo_WH_htt_CMS_MVALepResBoundingDown");
      TH1D *out_histo_WH_hww_CMS_MVALepResBoundingUp = (TH1D*)histo_WH_hww_CMS_MVALepResBoundingUp->Clone("histo_WH_hww_CMS_MVALepResBoundingUp");
      TH1D *out_histo_WH_hww_CMS_MVALepResBoundingDown = (TH1D*)histo_WH_hww_CMS_MVALepResBoundingDown->Clone("histo_WH_hww_CMS_MVALepResBoundingDown");
      TH1D *out_histo_WZ_CMS_MVALepResBoundingUp = (TH1D*)histo_WZ_CMS_MVALepResBoundingUp->Clone("histo_WZ_CMS_MVALepResBoundingUp");
      TH1D *out_histo_WZ_CMS_MVALepResBoundingDown = (TH1D*)histo_WZ_CMS_MVALepResBoundingDown->Clone("histo_WZ_CMS_MVALepResBoundingDown");
      TH1D *out_histo_ZZ_CMS_MVALepResBoundingUp = (TH1D*)histo_ZZ_CMS_MVALepResBoundingUp->Clone("histo_ZZ_CMS_MVALepResBoundingUp");
      TH1D *out_histo_ZZ_CMS_MVALepResBoundingDown = (TH1D*)histo_ZZ_CMS_MVALepResBoundingDown->Clone("histo_ZZ_CMS_MVALepResBoundingDown");
      TH1D *out_histo_Wgamma_CMS_MVALepResBoundingUp = (TH1D*)histo_Wgamma_CMS_MVALepResBoundingUp->Clone("histo_Wgamma_CMS_MVALepResBoundingUp");
      TH1D *out_histo_Wgamma_CMS_MVALepResBoundingDown = (TH1D*)histo_Wgamma_CMS_MVALepResBoundingDown->Clone("histo_Wgamma_CMS_MVALepResBoundingDown");
      TH1D *out_histo_VVV_CMS_MVALepResBoundingUp = (TH1D*)histo_VVV_CMS_MVALepResBoundingUp->Clone("histo_VVV_CMS_MVALepResBoundingUp");
      TH1D *out_histo_VVV_CMS_MVALepResBoundingDown = (TH1D*)histo_VVV_CMS_MVALepResBoundingDown->Clone("histo_VVV_CMS_MVALepResBoundingDown");
      TH1D *out_histo_WH_htt_SM_CMS_MVALepResBoundingUp = (TH1D*)histo_WH_htt_SM_CMS_MVALepResBoundingUp->Clone("histo_WH_htt_SM_CMS_MVALepResBoundingUp");
      TH1D *out_histo_WH_htt_SM_CMS_MVALepResBoundingDown = (TH1D*)histo_WH_htt_SM_CMS_MVALepResBoundingDown->Clone("histo_WH_htt_SM_CMS_MVALepResBoundingDown");
      TH1D *out_histo_WH_hww_SM_CMS_MVALepResBoundingUp = (TH1D*)histo_WH_hww_SM_CMS_MVALepResBoundingUp->Clone("histo_WH_hww_SM_CMS_MVALepResBoundingUp");
      TH1D *out_histo_WH_hww_SM_CMS_MVALepResBoundingDown = (TH1D*)histo_WH_hww_SM_CMS_MVALepResBoundingDown->Clone("histo_WH_hww_SM_CMS_MVALepResBoundingDown");
      TH1D *out_histo_WH_htt_CMS_MVAMETResBoundingUp = (TH1D*)histo_WH_htt_CMS_MVAMETResBoundingUp->Clone("histo_WH_htt_CMS_MVAMETResBoundingUp");
      TH1D *out_histo_WH_htt_CMS_MVAMETResBoundingDown = (TH1D*)histo_WH_htt_CMS_MVAMETResBoundingDown->Clone("histo_WH_htt_CMS_MVAMETResBoundingDown");
      TH1D *out_histo_WH_hww_CMS_MVAMETResBoundingUp = (TH1D*)histo_WH_hww_CMS_MVAMETResBoundingUp->Clone("histo_WH_hww_CMS_MVAMETResBoundingUp");
      TH1D *out_histo_WH_hww_CMS_MVAMETResBoundingDown = (TH1D*)histo_WH_hww_CMS_MVAMETResBoundingDown->Clone("histo_WH_hww_CMS_MVAMETResBoundingDown");
      TH1D *out_histo_WZ_CMS_MVAMETResBoundingUp = (TH1D*)histo_WZ_CMS_MVAMETResBoundingUp->Clone("histo_WZ_CMS_MVAMETResBoundingUp");
      TH1D *out_histo_WZ_CMS_MVAMETResBoundingDown = (TH1D*)histo_WZ_CMS_MVAMETResBoundingDown->Clone("histo_WZ_CMS_MVAMETResBoundingDown");
      TH1D *out_histo_ZZ_CMS_MVAMETResBoundingUp = (TH1D*)histo_ZZ_CMS_MVAMETResBoundingUp->Clone("histo_ZZ_CMS_MVAMETResBoundingUp");
      TH1D *out_histo_ZZ_CMS_MVAMETResBoundingDown = (TH1D*)histo_ZZ_CMS_MVAMETResBoundingDown->Clone("histo_ZZ_CMS_MVAMETResBoundingDown");
      TH1D *out_histo_Wgamma_CMS_MVAMETResBoundingUp = (TH1D*)histo_Wgamma_CMS_MVAMETResBoundingUp->Clone("histo_Wgamma_CMS_MVAMETResBoundingUp");
      TH1D *out_histo_Wgamma_CMS_MVAMETResBoundingDown = (TH1D*)histo_Wgamma_CMS_MVAMETResBoundingDown->Clone("histo_Wgamma_CMS_MVAMETResBoundingDown");
      TH1D *out_histo_VVV_CMS_MVAMETResBoundingUp = (TH1D*)histo_VVV_CMS_MVAMETResBoundingUp->Clone("histo_VVV_CMS_MVAMETResBoundingUp");
      TH1D *out_histo_VVV_CMS_MVAMETResBoundingDown = (TH1D*)histo_VVV_CMS_MVAMETResBoundingDown->Clone("histo_VVV_CMS_MVAMETResBoundingDown");
      TH1D *out_histo_WH_htt_SM_CMS_MVAMETResBoundingUp = (TH1D*)histo_WH_htt_SM_CMS_MVAMETResBoundingUp->Clone("histo_WH_htt_SM_CMS_MVAMETResBoundingUp");
      TH1D *out_histo_WH_htt_SM_CMS_MVAMETResBoundingDown = (TH1D*)histo_WH_htt_SM_CMS_MVAMETResBoundingDown->Clone("histo_WH_htt_SM_CMS_MVAMETResBoundingDown");
      TH1D *out_histo_WH_hww_SM_CMS_MVAMETResBoundingUp = (TH1D*)histo_WH_hww_SM_CMS_MVAMETResBoundingUp->Clone("histo_WH_hww_SM_CMS_MVAMETResBoundingUp");
      TH1D *out_histo_WH_hww_SM_CMS_MVAMETResBoundingDown = (TH1D*)histo_WH_hww_SM_CMS_MVAMETResBoundingDown->Clone("histo_WH_hww_SM_CMS_MVAMETResBoundingDown");
      TH1D *out_histo_WH_htt_CMS_MVAJESBoundingUp = (TH1D*)histo_WH_htt_CMS_MVAJESBoundingUp->Clone("histo_WH_htt_CMS_MVAJESBoundingUp");
      TH1D *out_histo_WH_htt_CMS_MVAJESBoundingDown = (TH1D*)histo_WH_htt_CMS_MVAJESBoundingDown->Clone("histo_WH_htt_CMS_MVAJESBoundingDown");
      TH1D *out_histo_WH_hww_CMS_MVAJESBoundingUp = (TH1D*)histo_WH_hww_CMS_MVAJESBoundingUp->Clone("histo_WH_hww_CMS_MVAJESBoundingUp");
      TH1D *out_histo_WH_hww_CMS_MVAJESBoundingDown = (TH1D*)histo_WH_hww_CMS_MVAJESBoundingDown->Clone("histo_WH_hww_CMS_MVAJESBoundingDown");
      TH1D *out_histo_WZ_CMS_MVAJESBoundingUp = (TH1D*)histo_WZ_CMS_MVAJESBoundingUp->Clone("histo_WZ_CMS_MVAJESBoundingUp");
      TH1D *out_histo_WZ_CMS_MVAJESBoundingDown = (TH1D*)histo_WZ_CMS_MVAJESBoundingDown->Clone("histo_WZ_CMS_MVAJESBoundingDown");
      TH1D *out_histo_ZZ_CMS_MVAJESBoundingUp = (TH1D*)histo_ZZ_CMS_MVAJESBoundingUp->Clone("histo_ZZ_CMS_MVAJESBoundingUp");
      TH1D *out_histo_ZZ_CMS_MVAJESBoundingDown = (TH1D*)histo_ZZ_CMS_MVAJESBoundingDown->Clone("histo_ZZ_CMS_MVAJESBoundingDown");
      TH1D *out_histo_Wgamma_CMS_MVAJESBoundingUp = (TH1D*)histo_Wgamma_CMS_MVAJESBoundingUp->Clone("histo_Wgamma_CMS_MVAJESBoundingUp");
      TH1D *out_histo_Wgamma_CMS_MVAJESBoundingDown = (TH1D*)histo_Wgamma_CMS_MVAJESBoundingDown->Clone("histo_Wgamma_CMS_MVAJESBoundingDown");
      TH1D *out_histo_VVV_CMS_MVAJESBoundingUp = (TH1D*)histo_VVV_CMS_MVAJESBoundingUp->Clone("histo_VVV_CMS_MVAJESBoundingUp");
      TH1D *out_histo_VVV_CMS_MVAJESBoundingDown = (TH1D*)histo_VVV_CMS_MVAJESBoundingDown->Clone("histo_VVV_CMS_MVAJESBoundingDown");
      TH1D *out_histo_WH_htt_SM_CMS_MVAJESBoundingUp = (TH1D*)histo_WH_htt_SM_CMS_MVAJESBoundingUp->Clone("histo_WH_htt_SM_CMS_MVAJESBoundingUp");
      TH1D *out_histo_WH_htt_SM_CMS_MVAJESBoundingDown = (TH1D*)histo_WH_htt_SM_CMS_MVAJESBoundingDown->Clone("histo_WH_htt_SM_CMS_MVAJESBoundingDown");
      TH1D *out_histo_WH_hww_SM_CMS_MVAJESBoundingUp = (TH1D*)histo_WH_hww_SM_CMS_MVAJESBoundingUp->Clone("histo_WH_hww_SM_CMS_MVAJESBoundingUp");
      TH1D *out_histo_WH_hww_SM_CMS_MVAJESBoundingDown = (TH1D*)histo_WH_hww_SM_CMS_MVAJESBoundingDown->Clone("histo_WH_hww_SM_CMS_MVAJESBoundingDown");
      TH1D *out_histo_Wjets_CMS_MVAWBoundingUp = (TH1D*)histo_Wjets_CMS_MVAWBoundingUp->Clone("histo_Wjets_CMS_MVAWBoundingUp");
      TH1D *out_histo_Wjets_CMS_MVAWBoundingDown = (TH1D*)histo_Wjets_CMS_MVAWBoundingDown->Clone("histo_Wjets_CMS_MVAWBoundingDown");
      TH1D *out_histo_WZ_CMS_MVAWZNLOBoundingUp = (TH1D*)histo_WZ_CMS_MVAWZNLOBoundingUp->Clone("histo_WZ_CMS_MVAWZNLOBoundingUp");
      TH1D *out_histo_WZ_CMS_MVAWZNLOBoundingDown = (TH1D*)histo_WZ_CMS_MVAWZNLOBoundingDown->Clone("histo_WZ_CMS_MVAWZNLOBoundingDown");
      TH1D *out_histo_ZZ_CMS_MVAZZNLOBoundingUp = (TH1D*)histo_ZZ_CMS_MVAZZNLOBoundingUp->Clone("histo_ZZ_CMS_MVAZZNLOBoundingUp");
      TH1D *out_histo_ZZ_CMS_MVAZZNLOBoundingDown = (TH1D*)histo_ZZ_CMS_MVAZZNLOBoundingDown->Clone("histo_ZZ_CMS_MVAZZNLOBoundingDown");

      if       (type[icard] == "sssf"){
        out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Up->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHTTv[imass]/BRHTTv[iref]*scaleFactor[0]);   
        out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Down->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHTTv[imass]/BRHTTv[iref]*scaleFactor[0]); 
        out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Up->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHWWv[imass]/BRHWWv[iref]*scaleFactor[1]);   
        out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Down->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHWWv[imass]/BRHWWv[iref]*scaleFactor[1]); 
      } else if(type[icard] == "ossf"){
        out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Up->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHTTv[imass]/BRHTTv[iref]*scaleFactor[0]);   
        out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Down->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHTTv[imass]/BRHTTv[iref]*scaleFactor[0]); 
        out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Up->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHWWv[imass]/BRHWWv[iref]*scaleFactor[1]);   
        out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Down->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHWWv[imass]/BRHWWv[iref]*scaleFactor[1]); 
      }
      out_histo_WH_htt_CMS_MVALepEffBoundingUp->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHTTv[imass]/BRHTTv[iref]*scaleFactor[0]); 
      out_histo_WH_htt_CMS_MVALepEffBoundingDown->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHTTv[imass]/BRHTTv[iref]*scaleFactor[0]);       
      out_histo_WH_htt_CMS_MVALepResBoundingUp->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHTTv[imass]/BRHTTv[iref]*scaleFactor[0]);
      out_histo_WH_htt_CMS_MVALepResBoundingDown->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHTTv[imass]/BRHTTv[iref]*scaleFactor[0]);
      out_histo_WH_htt_CMS_MVAMETResBoundingUp->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHTTv[imass]/BRHTTv[iref]*scaleFactor[0]);
      out_histo_WH_htt_CMS_MVAMETResBoundingDown->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHTTv[imass]/BRHTTv[iref]*scaleFactor[0]);
      out_histo_WH_htt_CMS_MVAJESBoundingUp->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHTTv[imass]/BRHTTv[iref]*scaleFactor[0]);
      out_histo_WH_htt_CMS_MVAJESBoundingDown->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHTTv[imass]/BRHTTv[iref]*scaleFactor[0]);
      out_histo_WH_hww_CMS_MVALepEffBoundingUp->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHWWv[imass]/BRHWWv[iref]*scaleFactor[1]); 
      out_histo_WH_hww_CMS_MVALepEffBoundingDown->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHWWv[imass]/BRHWWv[iref]*scaleFactor[1]);       
      out_histo_WH_hww_CMS_MVALepResBoundingUp->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHWWv[imass]/BRHWWv[iref]*scaleFactor[1]); 
      out_histo_WH_hww_CMS_MVALepResBoundingDown->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHWWv[imass]/BRHWWv[iref]*scaleFactor[1]);       
      out_histo_WH_hww_CMS_MVAMETResBoundingUp->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHWWv[imass]/BRHWWv[iref]*scaleFactor[1]); 
      out_histo_WH_hww_CMS_MVAMETResBoundingDown->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHWWv[imass]/BRHWWv[iref]*scaleFactor[1]);       
      out_histo_WH_hww_CMS_MVAJESBoundingUp->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHWWv[imass]/BRHWWv[iref]*scaleFactor[1]);    
      out_histo_WH_hww_CMS_MVAJESBoundingDown->Scale((ZHv[imass]+WHv[imass]+ttHv[imass])/(ZHv[iref]+WHv[iref]+ttHv[iref])*BRHWWv[imass]/BRHWWv[iref]*scaleFactor[1]);  

      TString outshapename = newDir + TString("/") + shapenamesv[icard];
      TFile outshapefile(outshapename,"RECREATE");
      out_histo_Data	 ->Write();
      out_histo_WH_htt   ->Write();
      out_histo_WH_hww   ->Write();
      out_histo_WZ	 ->Write();
      out_histo_ZZ	 ->Write();
      out_histo_Wjets	 ->Write();
      out_histo_Wgamma   ->Write();
      out_histo_VVV	 ->Write();
      out_histo_WH_htt_SM->Write();
      out_histo_WH_hww_SM->Write();
      if       (type[icard] == "sssf"){
      out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Up->Write();
      out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l1Down->Write();
      out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Up->Write();
      out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l1Down->Write();
      out_histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Up->Write();
      out_histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l1Down->Write();
      out_histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Up->Write();
      out_histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l1Down->Write();
      out_histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Up->Write();
      out_histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l1Down->Write();
      out_histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Up->Write();
      out_histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l1Down->Write();
      out_histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Up->Write();
      out_histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l1Down->Write();
      out_histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Up->Write();
      out_histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l1Down->Write();
      out_histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Up->Write();
      out_histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l1Down->Write();
      } else if(type[icard] == "ossf"){
      out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Up->Write();
      out_histo_WH_htt_CMS_MVAWH_httStatBounding_7TeV_vh3l2Down->Write();
      out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Up->Write();
      out_histo_WH_hww_CMS_MVAWH_hwwStatBounding_7TeV_vh3l2Down->Write();
      out_histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Up->Write();
      out_histo_WZ_CMS_MVAWZStatBounding_7TeV_vh3l2Down->Write();
      out_histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Up->Write();
      out_histo_ZZ_CMS_MVAZZStatBounding_7TeV_vh3l2Down->Write();
      out_histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Up->Write();
      out_histo_Wjets_CMS_MVAWjetsStatBounding_7TeV_vh3l2Down->Write();
      out_histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Up->Write();
      out_histo_Wgamma_CMS_MVAWgammaStatBounding_7TeV_vh3l2Down->Write();
      out_histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Up->Write();
      out_histo_VVV_CMS_MVAVVVStatBounding_7TeV_vh3l2Down->Write();
      out_histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Up->Write();
      out_histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_7TeV_vh3l2Down->Write();
      out_histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Up->Write();
      out_histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_7TeV_vh3l2Down->Write();
      }      
      out_histo_WH_htt_CMS_MVALepEffBoundingUp->Write();
      out_histo_WH_htt_CMS_MVALepEffBoundingDown->Write();
      out_histo_WH_hww_CMS_MVALepEffBoundingUp->Write();
      out_histo_WH_hww_CMS_MVALepEffBoundingDown->Write();
      out_histo_WZ_CMS_MVALepEffBoundingUp->Write();
      out_histo_WZ_CMS_MVALepEffBoundingDown->Write();
      out_histo_ZZ_CMS_MVALepEffBoundingUp->Write();
      out_histo_ZZ_CMS_MVALepEffBoundingDown->Write();
      out_histo_Wgamma_CMS_MVALepEffBoundingUp->Write();
      out_histo_Wgamma_CMS_MVALepEffBoundingDown->Write();
      out_histo_VVV_CMS_MVALepEffBoundingUp->Write();
      out_histo_VVV_CMS_MVALepEffBoundingDown->Write();
      out_histo_WH_htt_SM_CMS_MVALepEffBoundingUp->Write();
      out_histo_WH_htt_SM_CMS_MVALepEffBoundingDown->Write();
      out_histo_WH_hww_SM_CMS_MVALepEffBoundingUp->Write();
      out_histo_WH_hww_SM_CMS_MVALepEffBoundingDown->Write();
      out_histo_WH_htt_CMS_MVALepResBoundingUp->Write();
      out_histo_WH_htt_CMS_MVALepResBoundingDown->Write();
      out_histo_WH_hww_CMS_MVALepResBoundingUp->Write();
      out_histo_WH_hww_CMS_MVALepResBoundingDown->Write();
      out_histo_WZ_CMS_MVALepResBoundingUp->Write();
      out_histo_WZ_CMS_MVALepResBoundingDown->Write();
      out_histo_ZZ_CMS_MVALepResBoundingUp->Write();
      out_histo_ZZ_CMS_MVALepResBoundingDown->Write();
      out_histo_Wgamma_CMS_MVALepResBoundingUp->Write();
      out_histo_Wgamma_CMS_MVALepResBoundingDown->Write();
      out_histo_VVV_CMS_MVALepResBoundingUp->Write();
      out_histo_VVV_CMS_MVALepResBoundingDown->Write();
      out_histo_WH_htt_SM_CMS_MVALepResBoundingUp->Write();
      out_histo_WH_htt_SM_CMS_MVALepResBoundingDown->Write();
      out_histo_WH_hww_SM_CMS_MVALepResBoundingUp->Write();
      out_histo_WH_hww_SM_CMS_MVALepResBoundingDown->Write();
      out_histo_WH_htt_CMS_MVAMETResBoundingUp->Write();
      out_histo_WH_htt_CMS_MVAMETResBoundingDown->Write();
      out_histo_WH_hww_CMS_MVAMETResBoundingUp->Write();
      out_histo_WH_hww_CMS_MVAMETResBoundingDown->Write();
      out_histo_WZ_CMS_MVAMETResBoundingUp->Write();
      out_histo_WZ_CMS_MVAMETResBoundingDown->Write();
      out_histo_ZZ_CMS_MVAMETResBoundingUp->Write();
      out_histo_ZZ_CMS_MVAMETResBoundingDown->Write();
      out_histo_Wgamma_CMS_MVAMETResBoundingUp->Write();
      out_histo_Wgamma_CMS_MVAMETResBoundingDown->Write();
      out_histo_VVV_CMS_MVAMETResBoundingUp->Write();
      out_histo_VVV_CMS_MVAMETResBoundingDown->Write();
      out_histo_WH_htt_SM_CMS_MVAMETResBoundingUp->Write();
      out_histo_WH_htt_SM_CMS_MVAMETResBoundingDown->Write();
      out_histo_WH_hww_SM_CMS_MVAMETResBoundingUp->Write();
      out_histo_WH_hww_SM_CMS_MVAMETResBoundingDown->Write();
      out_histo_WH_htt_CMS_MVAJESBoundingUp->Write();
      out_histo_WH_htt_CMS_MVAJESBoundingDown->Write();
      out_histo_WH_hww_CMS_MVAJESBoundingUp->Write();
      out_histo_WH_hww_CMS_MVAJESBoundingDown->Write();
      out_histo_WZ_CMS_MVAJESBoundingUp->Write();
      out_histo_WZ_CMS_MVAJESBoundingDown->Write();
      out_histo_ZZ_CMS_MVAJESBoundingUp->Write();
      out_histo_ZZ_CMS_MVAJESBoundingDown->Write();
      out_histo_Wgamma_CMS_MVAJESBoundingUp->Write();
      out_histo_Wgamma_CMS_MVAJESBoundingDown->Write();
      out_histo_VVV_CMS_MVAJESBoundingUp->Write();
      out_histo_VVV_CMS_MVAJESBoundingDown->Write();
      out_histo_WH_htt_SM_CMS_MVAJESBoundingUp->Write();
      out_histo_WH_htt_SM_CMS_MVAJESBoundingDown->Write();
      out_histo_WH_hww_SM_CMS_MVAJESBoundingUp->Write();
      out_histo_WH_hww_SM_CMS_MVAJESBoundingDown->Write();
      out_histo_Wjets_CMS_MVAWBoundingUp->Write();
      out_histo_Wjets_CMS_MVAWBoundingDown->Write();
      out_histo_WZ_CMS_MVAWZNLOBoundingUp->Write();
      out_histo_WZ_CMS_MVAWZNLOBoundingDown->Write();
      out_histo_ZZ_CMS_MVAZZNLOBoundingUp->Write();
      out_histo_ZZ_CMS_MVAZZNLOBoundingDown->Write();
      outshapefile.Close();      
      refshapefile.Close();
      cout << outshapename << " created!" << endl;
    }
  }
  cout << endl;
}
