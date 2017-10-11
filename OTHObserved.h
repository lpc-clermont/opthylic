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

#ifndef OTH_OBSERVED_H
#define OTH_OBSERVED_H

#include <vector>

namespace OTH {

  class Observed {

  public:

    Observed();
    
    Observed(const unsigned int size);
    
    Observed(const Observed &obs);

    ~Observed();

    inline void resize(const unsigned int size) {m_events.resize(size,0);}
    inline void set(const unsigned int i,const int obs) {m_events[i]=obs;}

    bool operator==(const Observed &obs) const;
    bool operator<(const Observed &obs) const;

    void print() const;
    
  private:
    Observed &operator=(const Observed&);

    std::vector<int> m_events;
    
  };

}

#endif // OTH_OBSERVED_H
