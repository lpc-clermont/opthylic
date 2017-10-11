///////////////////////////////////////////////////////////////////////////////////
// This file is part the software OpTHyLiC
// Copyright © 2015 Laboratoire de Physique Corpusculaire de Clermont-Ferrand
// Developpers: David Calvet, Emmanuel Busato, Timothée Theveneaux-Pelzer
// Contact: opthylic@in2p3.fr
//
//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////
// Usage for interpreter mode:
//  in parent directory:
//  > make
//  > source setup.[c]sh
// then in examples directory:
//  > root -l load.C 'runSignificance.C("input1.dat")'
//  or
//  > root -l load.C 'runSignificance.C("input1.dat","input2.dat")'
//  ...
///////////////////////////////////////////////////////////
// Usage for compiled mode:
//  in parent directory:
//  > make
//  > source setup.[c]sh
// then in examples directory:
// > ./runSignificance.exe --files input1.dat input2.dat ...
///////////////////////////////////////////////////////////

#if defined EXECUTABLE || defined __CLING__

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

#include <TStopwatch.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TArrow.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TMath.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TString.h>

#include "OpTHyLiC.h"

using namespace std;
using namespace OTH;

#endif

void runSignificance(const std::string& file1, 
		     const std::string& file2="", 
		     const std::string& file3="",
		     const std::string& file4="", 
		     const std::string& file5="", 
		     const std::string& file6="", 
		     const std::string& file7="", 
		     const std::string& file8="") {
  TStopwatch w;
  w.Start();

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  // create OpTHyLiC instance
#if defined CPP11
  OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_mt19937); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_mt19937_64); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_minstd_rand); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_minstd_rand0); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_ranlux24_base); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_ranlux48_base); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_ranlux24); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_ranlux48); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_knuth_b); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::TR3); // using default TRandom3-based pseudo-random number generator, and automatic seed
//   OpTHyLiC oth(OTH::SystMclimit,OTH::StatNormal); // Provided by the C++11 standard library
#else
  OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatLogN); // using default TRandom3-based pseudo-random number generator, and automatic seed
#endif

  oth.addChannel("ch1",file1);
  if(file2!="")
    oth.addChannel("ch2",file2);
  if(file3!="")
    oth.addChannel("ch3",file3);
  if(file4!="")
    oth.addChannel("ch4",file4);
  if(file5!="")
    oth.addChannel("ch5",file5);
  if(file6!="")
    oth.addChannel("ch6",file6);
  if(file7!="")
    oth.addChannel("ch7",file7);
  if(file8!="")
    oth.addChannel("ch8",file8);

  oth.printSamples();
  
  const int Nexp=1e5;
  
  // Compute observed significance
  //
  // In order to compute the expected significance, replace OTH::SignifObserved by OTH::SignifExpectedP2sig, 
  // OTH::SignifExpectedP1sig, OTH::SignifExpectedMed, OTH::SignifExpectedM1sig or OTH::SignifExpectedM2sig 
  // depending on the kind of expected significance you want.
  //
  // By default, the signal strength is set to 1. If you want another value, just add a third argument to
  // oth.significance(...) with the desired value. 
  std::pair<double, double> s = oth.significance(OTH::SignifObserved,Nexp);
  const double p=s.first;
  const double z=s.second;

  cout << endl << "Results (cpu time=" << w.CpuTime()<< " sec, real time=" << w.RealTime() << " sec): " << endl;
  cout << " -> p=" << p << endl;
  cout << " -> z=" << z << endl;

  TH1F* hLLRb=(TH1F*)oth.getHistoLLRb();
  hLLRb->GetXaxis()->SetTitleSize(0.05);
  hLLRb->GetXaxis()->SetTitleOffset(0.88);
  hLLRb->GetXaxis()->SetTitle("q_{ #mu}");
  hLLRb->GetYaxis()->SetTitle("");
  hLLRb->SetLineColor(14);
  TH1F* hLLRsb=(TH1F*)oth.getHistoLLRsb();

  TCanvas *c1 = new TCanvas("c1", "c1",668,105,700,500);
  c1->SetLogy();
  hLLRb->Draw();
  hLLRsb->Draw("same");
  double qmuobs=oth.computeLLRdata();
  TArrow* arr = new TArrow(qmuobs,hLLRb->GetMaximum()/5.,qmuobs,0,0.02,"|>");
  arr->SetLineWidth(3);
  arr->SetLineColor(kRed);
  arr->SetFillColor(kRed);
  arr->Draw();
  TLatex latex1;
  latex1.SetTextSize(0.05);
  latex1.SetTextColor(kRed);
  latex1.DrawLatex(qmuobs,hLLRb->GetMaximum()/2.,"q_{ #mu}^{obs}");

  w.Stop();
}

#if defined EXECUTABLE
int main(int argc, char *argv[])
{
  const char* args[] = {"", "", "", "","", "", "", ""};
  std::vector<TString> vec(args, args + 8);
  int fileCounter=0;
  for (int i=1; i<argc; ++i) {
    std::string arg(argv[i]);
    if(arg=="--files") {
      int j=i+1;
      TString f(argv[j]);
      while(!f.Contains("--") && j<argc) {
	vec[fileCounter++]=f;
	++j; f=argv[j];
      }
      i=j-1;
    }
  }
  if(fileCounter==0) {
    cout << "ERROR! not input file specified" << endl;
    return -1;
  }
  cout << "Reading " << fileCounter << " files: " << vec[0]; 
  for(unsigned int i=1; i<vec.size(); ++i) {
    if(vec[i]!="") 
      cout << ", " << vec[i] ;
  }
  cout << endl;
  
  TApplication *theApp=new TApplication("myapp",&argc, argv);
  runSignificance(vec[0].Data(),vec[1].Data(),vec[2].Data(),vec[3].Data(),vec[4].Data(),vec[5].Data(),vec[6].Data(),vec[7].Data());
  cout<<endl<<"Press Ctrl+C or Clic on File->Exit ROOT in a Canvas to exit"<<endl;
  theApp->Run();
  theApp->Terminate();
  return 0;
}
#endif
