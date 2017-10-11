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

#ifndef OTH_PDFGENERATOR_H
#define OTH_PDFGENERATOR_H

#include "OTHTypes.h"

namespace OTH {

  class RdmGenerator;
  
  class PdfGenerator {
    
  public:
    PdfGenerator(RdmGenerator* rdmGen, const StatType statSampling);
    virtual ~PdfGenerator();
    
    double draw(const double mean, const double sigma);
    
    double drawNormal(const double mean, const double sigma);
    double drawLogN(const double mean, const double sigma);
    double drawGammaHyper(const double mean, const double sigma);
    double drawGammaUni(const double mean, const double sigma);
    double drawGammaJeffreys(const double mean, const double sigma);
    
    double drawGamma(const double mean, const double sigma, const float shapeParameterShift);
    
    int poisson(const double expected);

  protected:
    RdmGenerator* m_pRdmGen;

  private:
    PdfGenerator();
    PdfGenerator(const PdfGenerator&);
    PdfGenerator &operator=(const PdfGenerator&);
    
    // this pointer-to-function will point to one of the drawXXX functions above
    double (PdfGenerator::*m_pDraw) (const double mean, const double sigma);
  };

}

#endif // OTH_PDFGENERATOR_H
