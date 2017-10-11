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

#ifndef OTH_SHAPE_H
#define OTH_SHAPE_H

#include <iostream>
#include <deque>

#include "OTHShapeSyst.h"
#include "OTHSingleSyst.h"

namespace OTH {
  
  class Shape {
    
  public:
    
    Shape(const std::string &name,TH1* histoNominal);

    virtual ~Shape();

    inline std::string getName() const {return m_name;}
    inline int getNbins() const {return m_histoNominal->GetNbinsX();}
    float getBinContent(const int i) const;
    float getBinError(const int i) const;
    float getBinCenter(const int i) const;

    inline void setNameLaTeX(const std::string name) {m_nameLaTeX=name;}
    inline std::string getNameLaTeX() const {return m_nameLaTeX;}

    void addShapeSystematic(ShapeSyst* shapeSyst);
    void addGlobalSystematic(SingleSyst* globalSyst);
    inline std::deque<ShapeSyst*> getShapeSystematics() const {return m_shapeSyst;}
    inline std::deque<SingleSyst*> getGlobalSystematics() const {return m_globalSyst;}

    void print(const bool printStatUncert=true) const;
    void writeInputFile(std::ostream &out,const int iBin) const;

  private:
    Shape();
    Shape(const Shape&);
    Shape &operator=(const Shape&);

    std::string m_name;
    std::string m_nameLaTeX;
    TH1* m_histoNominal;
    std::deque<ShapeSyst*> m_shapeSyst;
    std::deque<SingleSyst*> m_globalSyst;
  };
}

#endif // OTH_SHAPE_H
