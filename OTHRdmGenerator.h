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

#ifndef OTH_RDMGENERATOR_H
#define OTH_RDMGENERATOR_H

#if defined CPP11
#include <random>
#endif

#include "TRandom3.h"


namespace OTH {
  /// Generic random number generator
  class RdmGenerator {
  public:
    RdmGenerator(const int seed=0);
    virtual ~RdmGenerator();
    int getInitSeed() const;
    virtual int poisson(const double expected)=0; // poisson distribution
    virtual double gaus(const double mean=0., const double sigma=1.)=0; // normal distribution
    virtual double logNormal(const double mean, const double sigma)=0; // lognormal distribution
    virtual double gamma(const double mean, const double sigma, const float shapeParameterShift)=0; // gamma distribution
    virtual double uniform()=0; // uniform distribution
  protected:
    int m_seed;//the original seed is kept
  };
  
  /// Random number generator using TRandom3
  class RdmGenerator_TR3 : public RdmGenerator {
  public:
    RdmGenerator_TR3(const int seed=0);
    virtual ~RdmGenerator_TR3();
    int poisson(const double expected); // poisson distribution
    double gaus(const double mean=0., const double sigma=1.); // normal distribution
    double logNormal(const double mean, const double sigma); // lognormal distribution
    double gamma(const double mean, const double sigma, const float shapeParameterShift); // gamma distribution
    double uniform(); // uniform distribution
  private:
    TRandom3 m_rdm;
  };
  
#if defined CPP11
  /// Template class for random number generator using pseudo-random number engines in the C++11 std library
  template <class T> class RdmGenerator_STD : public RdmGenerator {
  public:
    RdmGenerator_STD(const int seed=0);
    virtual ~RdmGenerator_STD();
    int poisson(const double expected); // poisson distribution
    double gaus(const double mean=0., const double sigma=1.); // normal distribution
    double logNormal(const double mean, const double sigma); // lognormal distribution
    double gamma(const double mean, const double sigma, const float shapeParameterShift); // gamma distribution
    double uniform(); // uniform distribution
  private:
    // pseudo-random number engine
    T m_engine;
  };
#endif

}

#endif // OTH_RDMGENERATOR_H
