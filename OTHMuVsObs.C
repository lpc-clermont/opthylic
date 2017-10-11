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

#include "OTHMuVsObs.h"
using namespace OTH;

MuVsObs::MuVsObs() :
  m_obsMin(-1),
  m_obsMax(-1),
  m_muMin(0),
  m_muMax(0)
{}

MuVsObs::MuVsObs(const MuVsObs &muObs) :
  m_obsMin(muObs.m_obsMin),
  m_obsMax(muObs.m_obsMax),
  m_muMin(muObs.m_muMin),
  m_muMax(muObs.m_muMax)
{}

MuVsObs::~MuVsObs()
{}

MuVsObs &MuVsObs::operator=(const MuVsObs &muObs)
{
  if (this!=&muObs) {
    m_obsMin=muObs.m_obsMin;
    m_obsMax=muObs.m_obsMax;
    m_muMin=muObs.m_muMin;
    m_muMax=muObs.m_muMax;
  }
  return *this;
}

void MuVsObs::add(const int obs,const double mu)
{
  if (m_obsMin<0 || m_obsMax<0) {
    m_obsMin=obs;
    m_obsMax=obs;
    m_muMin=mu;
    m_muMax=mu;
  } else {
    if (obs<m_obsMin) {
      m_obsMin=obs;
      m_muMin=mu;
    } else if (obs>m_obsMax) {
      m_obsMax=obs;
      m_muMax=mu;
    }
  }
}

void MuVsObs::reset()
{
  m_obsMin=-1;
  m_obsMax=-1;
  m_muMin=0;
  m_muMax=0;
} 

double MuVsObs::interpolateMu(const int obs) const
{
  if (m_obsMin!=m_obsMax) {
    const double distance=static_cast<double>(obs-m_obsMin)/static_cast<double>(m_obsMax-m_obsMin);
    double mu=m_muMin+distance*(m_muMax-m_muMin);
    if (mu<0) mu=0.001;
    return mu;
  }
  return -1;
}

    
