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

#ifndef OTH_ALGORITHMS_H
#define OTH_ALGORITHMS_H

#include <vector>

class TH1;

namespace OTH {

  class Base;

  class Algorithms {

  public:

    Algorithms();
    
    ~Algorithms();

    static double computeCLs(TH1 *pLLRsb,const TH1 *pLLRb,const double llr);
    
    static double sigStrengthExclusion(Base &clgen,const double mu0,const double mu0Step,
				       const int nbExp,const int type,double &cls,const double confLevel,
				       const bool extrapol=false);

    static double getCLsFromLLR(const int type,TH1 *pLLRsb,const TH1 *pLLRb);
    
    static std::vector<double> getQuantiles(const TH1 *pExpMu,const bool print=true);
    
  };

}

#endif // OTH_ALGORITHMS_H
