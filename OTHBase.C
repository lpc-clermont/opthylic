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
#include <stdexcept>
using namespace std;

#include "TGraph.h"
#include "TMath.h"
#include "TH1.h"

#include "OTHBase.h"
#include "OTHAlgorithms.h"
using namespace OTH;

Base::Base() :
  m_additiveSystComb(false),
  m_pExpMu(0),
  m_pCLs(0),
  m_pMuObs(0),
  m_pCLsMu(0),
  m_confLevel(0.95)
{}

Base::~Base()
{
  // do not delete histos and graphs, belong to ROOT
}

void Base::setConfLevel(const double cl)
{
  if (cl<1 && cl>0) m_confLevel = cl;
  else throw runtime_error("Wrong confidence level value provided !");
}

void Base::scanCLsVsMu(const double muMin,const double muMax,const int steps,const int nbExp,const int type)
{
  if (m_pCLsMu) {
    delete m_pCLsMu;
    m_pCLsMu=0;
  }

  if (muMin>=muMax || steps<2) return;
  m_pCLsMu=new TGraph(steps);
  m_pCLsMu->SetMarkerStyle(kCircle);

  const double muStep=(muMax-muMin)/static_cast<double>(steps-1);
  cout << "---> Interval for mu=[" << muMin << "," << muMax << "], " << steps << " steps" << endl;

  double mu=muMin;
  for(int step=0 ; step<steps ; ++step) {
    const double cls=generateForCLs(mu,nbExp,type);
    cout << "-> scanning for mu=" << mu << ", CLs=" << cls << endl;
    m_pCLsMu->SetPoint(step,mu,cls);
    mu+=muStep;
  }
}

pair<double,double> Base::significance(const SignifType type,const int nbExp,const double mu) 
{
  setSigStrength(mu);
  generateDistrLLR(nbExp);
  double pval=0.;
  if(SignifObserved==type) {
    pval=pValueData();
  }
  else {
    TH1 *hLLRsb = getHistoLLRsb();
    if(!hLLRsb) {
      throw runtime_error("hLLRsb not found !");
    }
    vector<double> quantiles=Algorithms::getQuantiles(hLLRsb,false);
    double quantile=0;
    if(type<=SignifExpectedM2sig) {
      quantile=quantiles[type];
    }
    else {
      throw runtime_error("unknown significance type !");
    }
    TH1 *hLLRb = getHistoLLRb();
    if(!hLLRb) {
      throw runtime_error("hLLRb not found !");
    }
    pval=hLLRb->Integral(0,hLLRb->FindBin(quantile));
  }
  return make_pair(pval,sqrt(2.)*TMath::ErfInverse(1-2*pval));
}
