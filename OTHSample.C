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
#include <sstream>
#include <iomanip>
using namespace std;

#include "TH1.h"
#include "TMath.h"

#include "OTHAlgorithms.h"
#include "OTHSample.h"
using namespace OTH;

Sample::Sample() :
  m_name(),
  m_nameLaTeX(),
  m_yield(),
  m_systLow(0),
  m_systHigh(0),
  m_systs(),
  m_pHyield(0)
{}

Sample::Sample(const string &name,const std::string &nameLaTeX,
	       const double nominal,const double stat) :
  m_name(name),
  m_nameLaTeX(nameLaTeX),
  m_yield(nominal,stat),
  m_systLow(0),
  m_systHigh(0),
  m_systs(),
  m_pHyield(0)
{}

void Sample::addSyst(const std::string &name,const unsigned int id,
		     const double low,const double high)
{
  if (0==low && 0==high) return;
  m_systs.push_back(SingleSyst(name,id,low,high));
  m_systs.back().createDistr(m_name);

  if (low<0 && low<high) {
    m_systLow=-TMath::Sqrt(m_systLow*m_systLow+low*low);
  } else if (high<0) {
    m_systLow=-TMath::Sqrt(m_systLow*m_systLow+high*high);
  }

  if (high>0 && high>low) {
    m_systHigh=TMath::Sqrt(m_systHigh*m_systHigh+high*high);
  } else if (low>0) {
    m_systHigh=TMath::Sqrt(m_systHigh*m_systHigh+low*low);
  }

  m_yield.setSystLow(m_systLow*m_yield.yield());
  m_yield.setSystHigh(m_systHigh*m_yield.yield());
}

string Sample::getLaTeXSystFromId(const unsigned int id,const int precision) const
{
  for(unsigned int s=0 ; s<m_systs.size() ; ++s) {
    if (m_systs[s].getId()==id) return getLaTeXSyst(getSystLow(s),getSystHigh(s),precision);
  }
  return "---";
}

string Sample::getLaTeXTotalSyst(const int precision) const
{
  return getLaTeXSyst(m_systLow,m_systHigh,precision);
}
    
TH1 *Sample::getSystDistr(const string &systName) const
{
  for(unsigned int i=0 ; i<m_systs.size() ; ++i) {
    if (m_systs[i].getName()==systName) return m_systs[i].getDistr();
  }
  throw runtime_error("Unknown systematics name !");
}

void Sample::print() const
{
  cout << "-> sample '" << m_name << "': yield = " << m_yield << endl;
  for(unsigned int i=0 ; i<m_systs.size() ; ++i) {
    m_systs[i].print();
  }
  cout << " -- total syst: " << m_systHigh*100 << "% "
       << m_systLow*100 << "% (syst)" << endl;
}

string Sample::getLaTeXSyst(const double low,const double high,const int precision)
{
  ostringstream str;
  str << fixed << setprecision(precision);

  if (TMath::Abs(low)==TMath::Abs(high)) {
    if (low<0 && high>0) {
      str << "$\\pm " << high*100 << "$";
      return str.str();
    } else if (low>0 && high<0) {
      str << "$\\mp " << low*100 << "$";
      return str.str();
    }
  }

  str << "$_{";
  if (low<0) str << low*100;
  else str << "+" << low*100;
  str << "}^{";
  if (high<0) str << high*100;
  else str << "+" << high*100;
  str << "}$";
  return str.str();
}

void Sample::createYieldHisto()
{
  if (m_pHyield) delete m_pHyield;
  string name=m_name;
  name+="_yield";
  m_pHyield=new TH1F(name.c_str(),";Events;Probability",1000,0,5*m_yield.yield());
}

void Sample::fillYieldHisto(const double yield) const
{
  if (m_pHyield) m_pHyield->Fill(yield);
}

YieldWithUncert Sample::getGeneratedYield() const
{
  YieldWithUncert yield;
  if (m_pHyield) {
    vector<double> quant=Algorithms::getQuantiles(m_pHyield,false);
    yield.setYield(quant[2]);
    yield.setSystHigh(quant[3]-quant[2]);
    yield.setSystLow(quant[1]-quant[2]);
  }
  return yield;
}

