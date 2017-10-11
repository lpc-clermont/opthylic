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

#ifndef OTH_SAMPLE_H
#define OTH_SAMPLE_H

#include <deque>

#include "OTHYieldWithUncert.h"
#include "OTHSingleSyst.h"

class TH1;

namespace OTH {

  class Sample {
    
  public:
    
    Sample();
    Sample(const std::string &name,const std::string &nameLaTeX,
	   const double nominal,const double stat);

    void addSyst(const std::string &name,const unsigned int id,
		 const double low,const double high);

    inline void setNameLaTeX(const std::string name) {m_nameLaTeX=name;}
    
    inline std::string getName() const {return m_name;}
    inline std::string getNameLaTeX() const {return m_nameLaTeX;}
    inline YieldWithUncert getYield() const {return m_yield;}
    inline double getNominal() const {return m_yield.yield();}
    inline double getStat() const {return m_yield.statUncert();}
    inline double getSystLow() const {return m_systLow;}
    inline double getSystHigh() const {return m_systHigh;}

    inline unsigned int getSystSize() const {return m_systs.size();}
    inline double getSystLow(const unsigned int i) const {return m_systs[i].getLow();}
    inline double getSystHigh(const unsigned int i) const {return m_systs[i].getHigh();}
    inline unsigned int getSystId(const unsigned int i) const {return m_systs[i].getId();}

    std::string getLaTeXSystFromId(const unsigned int id,const int precision) const;
    std::string getLaTeXTotalSyst(const int precision) const;

    inline void fillSystDistr(const unsigned int i,const double value) const {m_systs[i].fillDistr(value);}
    
    TH1 *getSystDistr(const std::string &systName) const;
    void print() const;

    void createYieldHisto();
    void fillYieldHisto(const double yield) const;
    TH1 *getYieldHisto() const {return m_pHyield;}
    YieldWithUncert getGeneratedYield() const;
    
  private:
    static std::string getLaTeXSyst(const double low,const double high,const int precision);

    std::string m_name,m_nameLaTeX;
    YieldWithUncert m_yield;
    double m_systLow,m_systHigh;
    std::deque<SingleSyst> m_systs;
    TH1 *m_pHyield;
  };

}

#endif // OTH_SAMPLE_H
