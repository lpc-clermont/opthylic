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

#ifndef OTH_YIELDWITHUNCERT_H
#define OTH_YIELDWITHUNCERT_H

#include <iostream>

namespace OTH {

  class YieldWithUncert {

  public:

    YieldWithUncert();
    YieldWithUncert(const double yld,const double stat);

    void set(const double yld,const double stat=0,const double syst=0);
    inline void setYield(const double yld) {m_yield=yld;}
    void setStat(const double stat);
    void setSymSyst(const double syst);
    void setSystLow(const double syst);
    void setSystHigh(const double syst);
    
    void add(const double yld,const double stat=0,const double syst=0);
    void addStat(const double stat);
    void addSymSyst(const double syst);
    void addSystLow(const double syst);
    void addSystHigh(const double syst);
    YieldWithUncert &operator+=(const YieldWithUncert &y);
    
    void rescale(const double factor);
    
    inline double yield() const {return m_yield;}
    inline double statUncert() const {return m_stat;}
    inline double systUncertLow() const {return m_systLow;}
    inline double systUncertHigh() const {return m_systHigh;}
    double totalUncertLow() const;
    double totalUncertHigh() const;
    
    std::string getLaTeX(const int precision=-1,const bool showStat=true) const;
    
  private:
    double m_yield;
    double m_stat,m_stat2;
    double m_systLow,m_systLow2;
    double m_systHigh,m_systHigh2;
  };
  
}

OTH::YieldWithUncert operator+(const OTH::YieldWithUncert &y1,const OTH::YieldWithUncert &y2);
std::ostream &operator<<(std::ostream &out,const OTH::YieldWithUncert &y);


#endif // OTH_YIELDWITHUNCERT_H
