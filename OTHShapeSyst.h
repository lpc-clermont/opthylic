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

#ifndef OTH_SHAPESYST_H
#define OTH_SHAPESYST_H

#include <vector>
class TH1;

namespace OTH {
  
  class ShapeSyst {
    
  public:
    
    ShapeSyst(const std::string& name, TH1* histoUp, TH1* histoDown);

    virtual ~ShapeSyst();

    inline std::string getName() const {return m_name;}

    int getNbins() const;
    float getBinContentUp(const int i) const;
    float getBinContentDown(const int i) const;

    void print() const;

  private:
    ShapeSyst();
    ShapeSyst(const ShapeSyst&);
    ShapeSyst &operator=(const ShapeSyst&);

    std::string m_name;
    TH1* m_histoUp;
    TH1* m_histoDown;
  };
}

#endif // OTH_SHAPESYST_H
