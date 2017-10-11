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
using namespace std;

#include "OTHObserved.h"
using namespace OTH;

Observed::Observed() :
  m_events()
{}

Observed::Observed(const unsigned int size) :
  m_events(size,0)
{}

Observed::Observed(const Observed &obs) :
  m_events(obs.m_events)
{}

Observed::~Observed()
{}

bool Observed::operator==(const Observed &obs) const
{
  for(unsigned int i=0 ; i<m_events.size() ; ++i) {
    if (m_events[i]!=obs.m_events[i]) return false;
  }
  return true;
}

bool Observed::operator<(const Observed &obs) const
{
  for(unsigned int i=0 ; i<m_events.size() ; ++i) {
    if (m_events[i]<obs.m_events[i]) return true;
    if (m_events[i]>obs.m_events[i]) return false;
  }
  return false;
}

void Observed::print() const
{
  for(unsigned int i=0 ; i<m_events.size() ; ++i) {
    cout << m_events[i] << " ";
  }
}
