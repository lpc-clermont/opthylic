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

#ifndef OTH_MUVSOBS_H
#define OTH_MUVSOBS_H

namespace OTH {

  class MuVsObs {

  public:

    MuVsObs();
    
    MuVsObs(const MuVsObs &muObs);

    ~MuVsObs();

    MuVsObs &operator=(const MuVsObs &muObs);

    void add(const int obs,const double mu);
    void reset();

    double interpolateMu(const int obs) const;
    
  private:

    int m_obsMin,m_obsMax;
    double m_muMin,m_muMax;
  };

}

#endif // OTH_MUVSOBS_H
