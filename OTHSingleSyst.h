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

#ifndef OTH_SINGLESYST_H
#define OTH_SINGLESYST_H

#include <string>

class TH1;

namespace OTH {

  class SingleSyst {
    
  public:
    
    SingleSyst();

    SingleSyst(const std::string &name,const unsigned int id,
	       const double low,const double high);

    virtual ~SingleSyst();

    // not to be used once distribution has been produced
    SingleSyst(const SingleSyst &syst);
    SingleSyst &operator=(const SingleSyst &syst);

    void createDistr(const std::string &smplName);
    void fillDistr(const double value) const;
    
    inline std::string getName() const {return m_name;}
    inline unsigned int getId() const {return m_id;}
    inline double getLow() const {return m_low;}
    inline double getHigh() const {return m_high;}
    
    void print() const;
    inline TH1 *getDistr() const {return m_pH;}
    
  private:
    std::string m_name;
    unsigned int m_id;
    double m_low,m_high;
    TH1 *m_pH;
  };

}

#endif // OTH_SINGLESYST_H
