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

#ifndef OTH_BASE_H
#define OTH_BASE_H

class TH1;
class TGraph;

#include "OTHTypes.h"

namespace OTH {

  class Base {

  public:

    Base();
    
    virtual ~Base();
    
    virtual void setConfLevel(const double cl); // set confidence level of computed limits
    inline double getConfLevel() const {return m_confLevel;}; // get confidence level of computed limits
    
    // setting of signal strength to be used for all computations
    virtual void setSigStrength(const double mu) =0;

    // generation of nbExp pseudo-experiments to compute the LLR distributions
    // must be called before trying to compute any CLs or p-value
    virtual void generateDistrLLR(const int nbExp) =0;

    // generation of nbExp pseudo-experiments with the signal strength mu
    // computation of the CLs for the observed events in data
    // (calls setSigStrength, generateDistrLLR and computeCLsData)
    // purely virtual function declared here for use in OTH::Algorithms
    virtual double generateForCLs(const double mu,const int nbExp,const int type) =0;
    
    // methods called for observed and expected (median, -+1 sigma, +-2 sigma) limit computation
    // type of limit is from the above enum
    virtual double sigStrengthExclusion(const LimitType type,const int nbExp,double &cls,
					const double muHint=1,const OTH::MethType method=OTH::MethDichotomy)=0;
    
    // methods called for observed and expected (median, -+1 sigma, +-2 sigma) significance computation
    std::pair<double,double> significance(const SignifType type,const int nbExp,const double mu=1);

    // computation of the p-value for the observed number of events
    // the LLR distributions must have been generated before
    virtual double pValueData() const =0;

    // retrieve LLR histos
    virtual TH1 *getHistoLLRsb() const =0;
    virtual TH1 *getHistoLLRb() const =0;

    // get distribution of expected mu
    inline TH1 *getDistrExpMu() const {return m_pExpMu;}

    // get expected mu versus number of observed events
    inline TGraph *getExpMuVsObs() const {return m_pMuObs;}

    // get distribution of computed CLs
    inline TH1 *getDistrCLs() const {return m_pCLs;}

    void scanCLsVsMu(const double muMin,const double muMax,const int steps,const int nbExp,const int type);

    // get CLS as a function of mu
    inline TGraph *getCLsVsMu() const {return m_pCLsMu;}
    
  protected:    
    
    bool m_additiveSystComb; // true if combination type for systematics is additive
    
    TH1 *m_pExpMu,*m_pCLs; // expected mu, computed CLs
    TGraph *m_pMuObs,*m_pCLsMu; // mu_up vs obs, CLs vs mu_up
    
    double m_confLevel; // confidence level of computed limits
    
  private:
    Base(const Base&);
    Base &operator=(const Base&);
  };
}

#endif // OTH_BASE_H
