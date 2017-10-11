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
using namespace OTH;

#include "TMath.h"

/// Generic random number generator
RdmGenerator::RdmGenerator(const int seed) :
  m_seed(seed)
{}

RdmGenerator::~RdmGenerator()
{}

int RdmGenerator::getInitSeed() const
{
  return m_seed;
}


/// Random number generator using TRandom3
RdmGenerator_TR3::RdmGenerator_TR3(const int seed) : 
  RdmGenerator(seed),
  m_rdm(seed)
{}

RdmGenerator_TR3::~RdmGenerator_TR3()
{}

int RdmGenerator_TR3::poisson(const double expected)
{
  return m_rdm.Poisson(expected);
}

double RdmGenerator_TR3::uniform()
{
  return m_rdm.Uniform();
}

double RdmGenerator_TR3::gaus(const double mean, const double sigma)
{
  return m_rdm.Gaus(mean, sigma);
}

double RdmGenerator_TR3::logNormal(const double mean, const double sigma)
{
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // We use the following definition :
  //   f(x;mu,sig)=1/(x*sqrt(2pi*sig^2)) TMath::Exp(-(ln x - mu)^2/(2sig^2))
  //
  // mu and sig are the expectation and standard deviation of the normal r.v. equal to ln x.
  //
  // Relation between mu, sigma and expectation (E), standard deviation (SD) of logN r.v. is as follows :
  //   mu=ln(E^2/sqrt(E^2+SD^2))
  //   sig=sqrt(ln(1+SD^2/E^2))
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  double rand=gaus(TMath::Log(mean*mean/TMath::Sqrt(mean*mean+sigma*sigma)),TMath::Sqrt(TMath::Log(1+sigma*sigma/(mean*mean))));
  rand=TMath::Exp(rand);
  return rand;
}

double RdmGenerator_TR3::gamma(const double mean, const double sigma, const float shapeParameterShift)
{
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // It uses the following definition
  //  f(x;beta,gamma) = \frac{x^{\gamma-1}\cdot\exp^{-x/\beta}{\Gamma(\gamma)\cdot\beta^{\gamma}}
  //  Relation between beta,gamma and expectation(E),standard deviation(SD) is as follows :
  //    beta = SD^2/E
  //    gamma = (E/SD)^2
  //
  // Ref:  Marsaglia, G. and Tsang, W. W.
  //       A Simple Method for Generating Gamma Variables
  //       ACM Transactions on Mathematical Software, Vol. 26, No. 3, September 2000
  ///////////////////////////////////////////////////////////////////////////////////////////////
  double rand = 0.;
  double gammavar=mean*mean/(sigma*sigma)+shapeParameterShift;
  double beta=sigma*sigma/mean;
  while(1) {
    double d=gammavar-1./3.;
    if(d<=0) {
      cerr << "OpTHyLiC Error ! can't generate gamma random number (change constraint type to OTH::StatLogN or OTH::StatNormal) -> quitting" << endl;
      throw runtime_error("gamma random number not supported");
    }
    double c=1./TMath::Sqrt(9.*d);
    double xgen=0;
    double v=0;
    while(v<=0.) {
      xgen=gaus();
      v=1.+c*xgen;
    }
    v = v*v*v; 
    double u=uniform();
    if(u<1.-.0331*xgen*xgen*xgen*xgen) {
      if(d*v*beta<1e4 && d*v*beta>0) {
	rand=d*v*beta;
	break;
      } 
    } 
    if(TMath::Log(u)<0.5*xgen*xgen+d*(1.-v+TMath::Log(v))) { 
      if(d*v*beta<1e4 && d*v*beta>0) {
	rand=d*v*beta ;
	break;
      }
    }
  }
  return rand;
}


#if defined CPP11
/// Template class for random number generator using pseudo-random number engines in the C++11 std library
template <class T>
RdmGenerator_STD<T>::RdmGenerator_STD(const int seed) : RdmGenerator(seed)
{
  if (m_seed==0) {// if seed is 0, generate one with entropy-based random device
    std::random_device rd;
    unsigned int sd = rd();
    m_engine.seed(sd);
  }
  else {
    m_engine.seed(seed);
  }
}

template <class T>
RdmGenerator_STD<T>::~RdmGenerator_STD()
{}

template <class T>
int RdmGenerator_STD<T>::poisson(const double expected)
{
  std::poisson_distribution<int> distr(expected);
  int r = distr(m_engine);
  return r;
}

template <class T>
double RdmGenerator_STD<T>::gaus(const double mean, const double sigma)
{
  std::normal_distribution<double> distr(mean,sigma);
  double r = distr(m_engine);
  return r;
}

template <class T>
double RdmGenerator_STD<T>::logNormal(const double mean, const double sigma)
{
  std::lognormal_distribution<double> distr(TMath::Log(mean*mean/TMath::Sqrt(mean*mean+sigma*sigma)),TMath::Sqrt(TMath::Log(1+sigma*sigma/(mean*mean))));
  double r = distr(m_engine);
  return r;
}

template <class T>
double RdmGenerator_STD<T>::gamma(const double mean, const double sigma, const float shapeParameterShift)
{
  double alpha=mean*mean/(sigma*sigma)+shapeParameterShift;
  double beta=sigma*sigma/mean;
  std::gamma_distribution<double> distr(alpha,beta);
  double r = distr(m_engine);
  return r;
}

template <class T>
double RdmGenerator_STD<T>::uniform()
{
  std::uniform_real_distribution<double> distr(0.,1.);
  double r = distr(m_engine);
  return r;
}

#if !defined __CLING__
// explicit template class instanciations*
template class OTH::RdmGenerator_STD<std::minstd_rand>;
template class OTH::RdmGenerator_STD<std::minstd_rand0>;
template class OTH::RdmGenerator_STD<std::mt19937>;
template class OTH::RdmGenerator_STD<std::mt19937_64>;
template class OTH::RdmGenerator_STD<std::ranlux24_base>;
template class OTH::RdmGenerator_STD<std::ranlux48_base>;
template class OTH::RdmGenerator_STD<std::ranlux24>;
template class OTH::RdmGenerator_STD<std::ranlux48>;
template class OTH::RdmGenerator_STD<std::knuth_b>;
#endif
#endif
