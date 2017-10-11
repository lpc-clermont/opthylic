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

#include <stdexcept>
#include <iostream>
using namespace std;

#include "TH1.h"
#include "TMath.h"

#include "OTHRdmGenerator.h"

#include "OTHSystematics.h"
using namespace OTH;

Systematics::Systematics(RdmGenerator* rdmGen, const SystType systInterpExtrapStyle) :
  m_variations(),
  m_pRdmGen(rdmGen),
  m_pH(0),
  m_names(),
  m_table(),
  m_pSF(0)
{
  m_pH=new TH1I("hSystSig","Systematics;Sigmas;Entries",240,-6,6);
  
  if(systInterpExtrapStyle==SystMclimit) m_pSF = &Systematics::getScaleFactorMCLimit;
  else if(systInterpExtrapStyle==SystLinear) m_pSF = &Systematics::getScaleFactorLinear;
  else if(systInterpExtrapStyle==SystExpo) m_pSF = &Systematics::getScaleFactorExpo;
  else if(systInterpExtrapStyle==SystPolyexpo) m_pSF = &Systematics::getScaleFactorPolyExpo;
  else {
    cerr << "OpTHyLiC Error ! Unknown systematic uncertainty style " 
	 << systInterpExtrapStyle << " !" << endl;
    throw runtime_error("Unknown systematic uncertainty style !");
  }
}

Systematics::~Systematics()
{
  // do not delete m_pH, belongs to ROOT
}

unsigned int Systematics::add(const string &name)
{
  map<string,unsigned int>::const_iterator it=m_table.find(name);
  if (it!=m_table.end()) return it->second;

  m_variations.push_back(0);
  m_names.push_back(name);
  if (m_names.size()!=m_variations.size()) throw runtime_error("Sizes of systematics names and variations differ !");
  const unsigned int index=m_variations.size()-1;
  m_table[name]=index;
  return index;
}

void Systematics::variate()
{
  for(unsigned int i=0 ; i<m_variations.size() ; ++i) {
    // find a variation in sigmas, within +-5
    double var=0;
    do {
      var=m_pRdmGen->gaus(0,1);
    } while (var<-5 || var>5);
    m_variations[i]=var;
    m_pH->Fill(var);
  }
}

double Systematics::getVariation(const string &name) const
{
  map<string,unsigned int>::const_iterator it=m_table.find(name);
  if (it!=m_table.end()) return m_variations[it->second];
  throw runtime_error("Unknown systematics name !");  
}

double Systematics::getVariation(const unsigned int index) const
{
  if (index<m_variations.size()) return m_variations[index];
  throw runtime_error("Unknown systematics index !");
}

string Systematics::getName(const unsigned int index) const
{
  if (index<m_names.size()) return m_names[index];
  throw runtime_error("Unknown systematics index !");
}

void Systematics::print() const
{
  cout << "======= List of systematics =============" << endl;
  for(unsigned int i=0 ; i<m_names.size() ; ++i) {
    cout << "-> '" << m_names[i] << "': var=" << m_variations[i]
	 << " sigma" << endl;
  }
}


double Systematics::getScaleFactor(const unsigned int index,
				  const double low,const double high) const
{
  return (this->*m_pSF)(index, low, high);
}

double Systematics::getScaleFactorMCLimit(const unsigned int index,
					  const double low,const double high) const

{
  if (index<m_variations.size()) {
    const double var=m_variations[index];
    double sig=-low;
    if (var>0) sig=high;
    const double quadMatch=var*(high-low)/2 + var*var*(high+low)/2;
    const double rf=1/(1+3*TMath::Abs(var));
    const double bridge=var*sig*(1-rf) + rf*quadMatch;
    double lnB=0;
    if (bridge<0) lnB=TMath::Exp(bridge);
    else lnB=bridge+1;
    return lnB;
  }
  throw runtime_error("Unknown systematics index !");
}

double Systematics::getScaleFactorLinear(const unsigned int index,
					 const double low,const double high) const

{
  if (index<m_variations.size()) {
    const double var=m_variations[index];
    double sf=1+var*high;
    if(var<0) sf=1-var*low;
    if(sf<0) sf=0;
    return sf;
  }
  throw runtime_error("Unknown systematics index !");
}

double Systematics::getScaleFactorExpo(const unsigned int index,
				       const double low,const double high) const

{
  if (index<m_variations.size()) {
    double var=m_variations[index];
    double sf=0;
    if((var>=0 && high>-1)||(var<0 && low>-1)) { // exponential interp/extrap
      double sig=high+1;
      if(var<0) {
	sig=low+1;
	var*=-1;
      }
      sf=TMath::Power(sig,var);
    }
    else { // linear
      sf=getScaleFactorLinear(index, low, high);
    }
    return sf;
  }
  throw runtime_error("Unknown systematics index !");
}

double Systematics::getScaleFactorPolyExpo(const unsigned int index,
					 const double low,const double high) const

{
  if (index<m_variations.size()) { // polynomial interpolation
    double var=m_variations[index];
    double sf=0;
    if(var>-1&&var<1) {
      double pow_up       = 1+high;
      double pow_down     = 1+low;
      double pow_up_log   = (1+high) <= 0. ? 0. : pow_up*TMath::Log(1+high);
      double pow_down_log = (1+low) <= 0. ? 0. : -pow_down*TMath::Log(1+low);
      double pow_up_log2  = (1+high) <= 0. ? 0. : pow_up_log*TMath::Log(1+high);
      double pow_down_log2= (1+low) <= 0. ? 0. : pow_down_log*TMath::Log(1+low);
      
      double S0 = (pow_up+pow_down)/2;
      double A0 = (pow_up-pow_down)/2;
      double S1 = (pow_up_log+pow_down_log)/2;
      double A1 = (pow_up_log-pow_down_log)/2;
      double S2 = (pow_up_log2+pow_down_log2)/2;
      double A2 = (pow_up_log2-pow_down_log2)/2;
      double a = 1/8.*(      15*A0 -  7*S1 + A2);
      double b = 1/8.*(-24 + 24*S0 -  9*A1 + S2);
      double c = 1/4.*(    -  5*A0 +  5*S1 - A2);
      double d = 1/4.*( 12 - 12*S0 +  7*A1 - S2);
      double e = 1/8.*(    +  3*A0 -  3*S1 + A2);
      double f = 1/8.*( -8 +  8*S0 -  5*A1 + S2);
      
      double var2 = var*var ;
      double var3 = var2*var ;
      sf = 1 + a*var + b*var2 + c*var3 + d*var2*var2 + e*var3*var2 + f*var3*var3;
    }
    else { // exponential extrapolation
      sf=getScaleFactorExpo(index, low, high);
    }
    return sf;
  }
  throw runtime_error("Unknown systematics index !");
}
