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

#include <sstream>
#include <iomanip>
using namespace std;

#include "TMath.h"

#include "OTHYieldWithUncert.h"
using namespace OTH;

YieldWithUncert::YieldWithUncert() : 
  m_yield(0),
  m_stat(0),m_stat2(0),
  m_systLow(0),m_systLow2(0),
  m_systHigh(0),m_systHigh2(0)
{}

YieldWithUncert::YieldWithUncert(const double yld,const double stat) : 
  m_yield(yld),
  m_stat(stat),m_stat2(stat*stat),
  m_systLow(0),m_systLow2(0),
  m_systHigh(0),m_systHigh2(0)
{}

void YieldWithUncert::set(const double yld,const double stat,const double syst)
{
  m_yield=yld;
  setStat(stat);
  setSymSyst(syst);
}

void YieldWithUncert::setStat(const double stat)
{
  m_stat=stat;
  m_stat2=stat*stat;
}

void YieldWithUncert::setSymSyst(const double syst)
{
  setSystLow(syst);
  setSystHigh(syst);
}

void YieldWithUncert::setSystLow(const double syst)
{
  m_systLow=TMath::Abs(syst);
  m_systLow2=syst*syst;
}

void YieldWithUncert::setSystHigh(const double syst)
{
  m_systHigh=TMath::Abs(syst);
  m_systHigh2=syst*syst;
}

void YieldWithUncert::add(const double yld,const double stat,const double syst)
{
  m_yield+=yld;
  addStat(stat);
  addSymSyst(syst);
}

void YieldWithUncert::addStat(const double stat)
{
  m_stat2+=stat*stat;
  m_stat=TMath::Sqrt(m_stat2);
}

void YieldWithUncert::addSymSyst(const double syst)
{
  addSystLow(syst);
  addSystHigh(syst);
}

void YieldWithUncert::addSystLow(const double syst)
{
  m_systLow2+=syst*syst;
  m_systLow=TMath::Sqrt(m_systLow2);
}

void YieldWithUncert::addSystHigh(const double syst)
{
  m_systHigh2+=syst*syst;
  m_systHigh=TMath::Sqrt(m_systHigh2);
}

YieldWithUncert &YieldWithUncert::operator+=(const YieldWithUncert &y)
{
  m_yield+=y.m_yield;
  m_stat2+=y.m_stat2;
  m_stat=TMath::Sqrt(m_stat2);
  m_systLow2+=y.m_systLow2;
  m_systLow=TMath::Sqrt(m_systLow2);
  m_systHigh2+=y.m_systHigh2;
  m_systHigh=TMath::Sqrt(m_systHigh2);
  return *this;
}

void YieldWithUncert::rescale(const double factor)
{
  m_yield*=factor;
  m_stat*=factor;
  m_stat2=m_stat*m_stat;
  m_systLow*=factor;
  m_systLow2=m_systLow*m_systLow;
  m_systHigh*=factor;
  m_systHigh2=m_systHigh*m_systHigh;
}
    
double YieldWithUncert::totalUncertLow() const
{
  return TMath::Sqrt(m_stat2+m_systLow2);
}

double YieldWithUncert::totalUncertHigh() const
{
  return TMath::Sqrt(m_stat2+m_systHigh2);
}

string YieldWithUncert::getLaTeX(const int precision,const bool showStat) const
{
  double max=(showStat?m_stat:0);
  if (m_systLow>max) max=m_systLow;
  if (m_systHigh>max) max=m_systHigh;

  int prec=precision;
  if (-1==prec) {
    if (max<10 && max>0) {
      prec=1;
      for(double val=1 ; val>1e-10 ; val/=10,++prec) {
	if (max>val) break;
      }
    } else prec=0;
  }

  ostringstream str;
  str << fixed << setprecision(prec);
  str << "$" << m_yield;
  if (showStat) str << " \\pm " << m_stat;
  if (m_systLow!=0 || m_systHigh!=0) {
    ostringstream strLow,strHigh;
    strLow << fixed << setprecision(prec) << m_systLow;
    strHigh << fixed << setprecision(prec) << m_systHigh;
    if (strLow.str()==strHigh.str()) str << " \\pm " << strLow.str();
    else str << " ^{+" << strHigh.str() << "}_{-" << strLow.str() << "}";
  }
  str << "$";
  return str.str();
}
      
YieldWithUncert operator+(const YieldWithUncert &y1,const YieldWithUncert &y2)
{
  YieldWithUncert result(y1);
  return result+=y2;
}

ostream &operator<<(ostream &out,const YieldWithUncert &y)
{
  out << y.yield() << " +- " << y.statUncert() << " (stat)";
  if (y.systUncertLow()!=0 || y.systUncertHigh()!=0) {
    if (y.systUncertLow()==y.systUncertHigh()) out << " +- " << y.systUncertLow();
    else out << " +" << y.systUncertHigh() << " -" << y.systUncertLow();
    out << " (syst)";
  }
  return out;
}


