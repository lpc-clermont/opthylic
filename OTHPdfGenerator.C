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

#include <iostream>
#include <stdexcept>
using namespace std;

#include "OTHRdmGenerator.h"

#include "OTHPdfGenerator.h"

using namespace OTH;

PdfGenerator::PdfGenerator(RdmGenerator* rdmGen, const StatType statSampling) :
  m_pRdmGen(rdmGen),
  m_pDraw(0)
{  
  if(statSampling==StatNormal) m_pDraw = &PdfGenerator::drawNormal;
  else if(statSampling==StatLogN) m_pDraw = &PdfGenerator::drawLogN;
  else if(statSampling==StatGammaHyper) m_pDraw = &PdfGenerator::drawGammaHyper;
  else if(statSampling==StatGammaUni) m_pDraw = &PdfGenerator::drawGammaUni;
  else if(statSampling==StatGammaJeffreys) m_pDraw = &PdfGenerator::drawGammaJeffreys;
  else {
    cerr << "OpTHyLiC Error ! Unknown sampling method "
	 << statSampling << " !" << endl;
    throw runtime_error("Unknown sampling method !");
  }
}

PdfGenerator::~PdfGenerator()
{}

int PdfGenerator::poisson(const double expected)
{
  return m_pRdmGen->poisson(expected);
}


double PdfGenerator::draw(const double mean, const double sigma)
{
  return (this->*m_pDraw)(mean, sigma);
}

double PdfGenerator::drawNormal(const double mean, const double sigma)
{
  double rand=0;
  do {
    rand=m_pRdmGen->gaus(mean,sigma);
  } while (rand<0);
  return rand;
}

double PdfGenerator::drawLogN(const double mean, const double sigma)
{
  double rand=0;
  if (0!=mean) {
    rand=m_pRdmGen->logNormal(mean,sigma);
  }
  else {
    rand = drawNormal(mean, sigma);
  }
  return rand;
  
}

double PdfGenerator::drawGammaHyper(const double mean, const double sigma)
{
  return drawGamma(mean, sigma, 0.);
}

double PdfGenerator::drawGammaUni(const double mean, const double sigma)
{
  return drawGamma(mean, sigma, 1.);
}

double PdfGenerator::drawGammaJeffreys(const double mean, const double sigma)
{
  return drawGamma(mean, sigma, 0.5);
}

double PdfGenerator::drawGamma(const double mean, const double sigma, const float shapeParameterShift)
{
  double rand=0;

  if (0!=mean) {
    rand = m_pRdmGen->gamma(mean, sigma, shapeParameterShift);
  }
  else {
    rand = drawNormal(mean, sigma);
  }
  return rand;
}
