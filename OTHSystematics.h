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

#ifndef OTH_SYSTEMATICS_H
#define OTH_SYSTEMATICS_H

#include <string>
#include <deque>
#include <map>

class TH1;

#include "OTHTypes.h"

namespace OTH {

  class RdmGenerator;

  class Systematics {
    
  public:
    
    Systematics(RdmGenerator* rdmGen, const SystType systInterpExtrapStyle);
    virtual ~Systematics();

    unsigned int add(const std::string &name);
    virtual void variate();
    
    double getScaleFactor(const unsigned int index,
			  const double low,const double high) const;
				  
    double getScaleFactorMCLimit(const unsigned int index,
				 const double low,const double high) const;
    double getScaleFactorLinear(const unsigned int index,
				const double low,const double high) const;
    double getScaleFactorExpo(const unsigned int index,
			      const double low,const double high) const;
    double getScaleFactorPolyExpo(const unsigned int index,
				  const double low,const double high) const;

    double getVariation(const std::string &name) const;
    double getVariation(const unsigned int index) const;

    inline unsigned int getSize() const {return m_names.size();}
    std::string getName(const unsigned int index) const;
    TH1 *getDistr() const {return m_pH;}
    void print() const;
    
  protected:
    std::deque<double> m_variations;
    RdmGenerator* m_pRdmGen;
    TH1 *m_pH;

  private:
    Systematics();
    Systematics(const Systematics&);
    Systematics &operator=(const Systematics&);

    std::deque<std::string> m_names;
    std::map<std::string,unsigned int> m_table;
    
    // this pointer-to-function will point to one of the getScaleFactorXXX functions above
    double (Systematics::*m_pSF) (const unsigned int index,
				  const double low,const double high) const;
  };

}

#endif // OTH_SYSTEMATICS_H
