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

#include <iostream>
using namespace std;

#include "TH1.h"
#include "TMath.h"

#include "OTHBase.h"
#include "OTHAlgorithms.h"
using namespace OTH;

Algorithms::Algorithms() 
{}

Algorithms::~Algorithms()
{}

double Algorithms::computeCLs(TH1 *pLLRsb,const TH1 *pLLRb,const double llr)
{
  const int bin=pLLRsb->FindBin(llr);
  const int maxBin=pLLRsb->GetNbinsX()+1;
  double clsb=0,clb=0;
  for(int b=bin ; b<=maxBin ; ++b) {
    clsb+=pLLRsb->GetBinContent(b);
    clb+=pLLRb->GetBinContent(b);
  }
  if (clb>1e-5) return clsb/clb;
  return -1;
}

double Algorithms::sigStrengthExclusion(Base &clgen,const double mu0,const double mu0Step,
					const int nbExp,const int type,double &cls,const double confLevel,
					const bool extrapol)
{
  // coarse scan of mu
  double mu=mu0,muPrev=0,muStep=mu0Step,clsPrev=0;
  cout << "---> Searching for a reasonable mu interval, from " << mu << endl;
  //  searching for a non-zero value of CLs
  cls=0;
  double muFactor=10;
  enum {None,Up,Down};
  int direction=None;
  for(int i=0 ; cls<=0 || cls>0.9 ; ++i) {
    if (i>50) {
      cout << "######################################################" << endl
	   << "##### still no correct value for CLs, aborting ! #####" << endl
	   << "######################################################" << endl;
      return 0;
    }
    cls=clgen.generateForCLs(mu,nbExp,type);
    cout << "-> scanning first point: mu=" << mu << ", CLs=" << cls << endl;
    if (cls<=0) {
      if (direction==Up) muFactor/=2;
      mu/=muFactor;
      direction=Down;
    } else if (cls>0.9) {
      if (direction==Down) muFactor/=2;
      mu*=muFactor;
      direction=Up;
    }
  }
  
  //  searching for a second point
  if (cls<1-confLevel) muStep=1/muStep;    
  muPrev=mu;
  clsPrev=cls;
  mu*=muStep;
  
  cls=0;
  direction=None;
  for(int i=0 ; cls<=0 || cls>0.9 ; ++i) {
    if (i>50) {
      cout << "######################################################" << endl
	   << "##### still no correct value for CLs, aborting ! #####" << endl
	   << "######################################################" << endl;
      return 0;
    }
    cls=clgen.generateForCLs(mu,nbExp,type);
    cout << "-> scanning second point: mu=" << mu << ", CLs=" << cls << endl;
    if (cls<=0) {
      if (direction==Up) muFactor/=2;
      mu/=muFactor;
      direction=Down;
    } else if (cls>0.9) {
      if (direction==Down) muFactor/=2;
      mu*=muFactor;
      direction=Up;
    }
    if (mu==muPrev) mu=muPrev*1.1;
  }

  const double targCLs=1-confLevel;
  const double logTargCLs=TMath::Log(targCLs);
  const double minCLs=targCLs*0.95;
  const double maxCLs=targCLs*1.05;
  const double precMu=0.01;
  double muMin=muPrev,muMax=mu,clsMin=clsPrev,clsMax=cls;
  if (muMin>muMax) {
    muMin=mu;
    clsMin=cls;
    muMax=muPrev;
    clsMax=clsPrev;
  }
  double logClsMin=TMath::Log(clsMin);
  double logClsMax=TMath::Log(clsMax);

  if (extrapol) {
    // direct extrapolation of mu
    cout << "---> Extrapolating mu from references: " << muMin << ", " << muMax << endl;
    mu=muMin+(muMax-muMin)*(logTargCLs-logClsMin)/(logClsMax-logClsMin);
    if (mu<0) mu=0;
    cls=clgen.generateForCLs(mu,nbExp,type);
    cout << "---> Extrapolated mu=" << mu << ", CLs=" << cls << endl;

  } else {
    // finer scan of mu (dichotomy)
    cout << "---> Log-dichotomy search with mu references: " << muMin << ", " << muMax << endl;
    for(int i=0 ; (muMax-muMin)/muMin>precMu ; ++i) {
      if (i>50) {
	cout << "##################################################" << endl
	     << "##### no convergence of mu found, aborting ! #####" << endl
	     << "##################################################" << endl;
	return 0;
      }
      if (clsMax<0.00001||clsMax==clsMin) mu=muMin+(muMax-muMin)*(targCLs-clsMin)/(-clsMin);
      else {
	mu=muMin+(muMax-muMin)*(logTargCLs-logClsMin)/(logClsMax-logClsMin);
	if (mu<0) mu=0;
      }
      cls=clgen.generateForCLs(mu,nbExp,type);
      cout << "-> searching for mu=" << mu << ", CLs=" << cls << " (refs: " << muMin << ", " << muMax << ")" << endl;
      if (cls>minCLs && cls<maxCLs) {
	cout << "---> Close enough to " << targCLs << ", stopping" << endl;
	return mu;
      } else if (cls>targCLs) {
	if (mu>muMax) {
	  muMin=muMax;
	  clsMin=clsMax;
	  logClsMin=logClsMax;
	  muMax=mu;
	  clsMax=cls;
	  logClsMax=TMath::Log(cls);
	} else {
	  muMin=mu;
	  clsMin=cls;
	  logClsMin=TMath::Log(cls);
	}
      } else {
	if (mu<muMin) {
	  muMax=muMin;
	  clsMax=clsMin;
	  logClsMax=logClsMin;
	  muMin=mu;
	  clsMin=cls;
	  logClsMin=TMath::Log(cls);
	} else {
	  muMax=mu;
	  clsMax=cls;
	  logClsMax=TMath::Log(cls);
	}
      }
    }
    mu=(muMin+muMax)/2;
    cls=clgen.generateForCLs(mu,nbExp,type);
    cout << "---> Best mu=" << mu << " +- " << (muMax-muMin)/2 << ", CLs=" << cls << endl;
  }
  return mu;
}

double Algorithms::getCLsFromLLR(const int type,TH1 *pLLRsb,const TH1 *pLLRb)
{
  // compute quantiles
  const int nbQuant=5;
  if (type>=nbQuant) return 0;
  const double cdf[nbQuant]={0.0228,0.1587,0.5,0.8413,0.9772};
  double sum=0,sumPrev=0,varPrev=0;
  for(int b=0 ; b<=pLLRb->GetNbinsX()+1 ; ++b) {
    const double val=pLLRb->GetBinContent(b);
    if (val>0) {
      sum+=val;
      const double var=pLLRb->GetBinLowEdge(b);
      if (sum>cdf[type]) {
	if (0==sumPrev) return computeCLs(pLLRsb,pLLRb,var);
	else {
	  const double distance=(cdf[type]-sumPrev)/(sum-sumPrev);
	  return computeCLs(pLLRsb,pLLRb,varPrev+distance*(var-varPrev));
	}
      }
      varPrev=var;
      sumPrev=sum;
    }
  }
  return -1;
}
    
vector<double> Algorithms::getQuantiles(const TH1 *pHisto,const bool print)
{
  // compute quantiles
  const int nbQuant=5;
  const double cdf[nbQuant]={0.0228,0.1587,0.5,0.8413,0.9772};
  vector<double> varQ(nbQuant,0);
  double sum=0,sumPrev=0,varPrev=0;
  for(int b=0,q=0 ; b<=pHisto->GetNbinsX()+1 && q<nbQuant ; ++b) {
    const double val=pHisto->GetBinContent(b);
    if (val>0) {
      sum+=val;
      const double var=pHisto->GetBinLowEdge(b);
      for(; sum>cdf[q] && q<nbQuant ; ++q) {
	if (0==sumPrev) varQ[q]=var;
	else {
	  const double distance=(cdf[q]-sumPrev)/(sum-sumPrev);
	  varQ[q]=varPrev+distance*(var-varPrev);
	}
      }
      varPrev=var;
      sumPrev=sum;
    }
  }

  if(print) {
    cout << "======= Expected signal strengths ==============" << endl;
    const char *names[nbQuant]={"-2 sig","-1 sig","median","+1 sig","+2 sig"};
    for(int q=0 ; q<nbQuant ; ++q) {
      cout << "---> var(" << names[q] << ")=" << varQ[q] << endl;    
    }
  }

  return varQ;
}

