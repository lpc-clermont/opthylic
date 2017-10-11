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

#ifndef OTH_CHANNEL_H
#define OTH_CHANNEL_H

#include <map>

class TH1;

#include "OTHSample.h"
#include "OTHAlgorithms.h"
#include "OTHMuVsObs.h"
#include "OTHBase.h"

namespace OTH {
  class RdmGenerator;
  class Systematics;
  class PdfGenerator;
  
  class Channel: public Base {

  public:

    Channel(const std::string &name,Systematics &syste,PdfGenerator &statSampling);
    
    virtual ~Channel();
    
    std::string getName() const {return m_name;}
    std::string getNameLaTeX() const {return m_nameLaTeX;}
    const std::deque<Sample> &getBkgSamples() const {return m_bgSamples;}
    const Sample &getSigSample() const {return m_sigSample;}
    YieldWithUncert getYieldBkg() const;

    // setting of samples yields and uncertainties
    void addSamples(const std::string &fileName);
    
    unsigned int addBkgSample(const std::string &name,const double nominal,const double stat);
    void addBkgSystematics(const unsigned int iSample,
			   const std::string &systName,
			   const double up,const double down);
    
    void setSigSample(const std::string &name,const double nominal,const double stat);
    void addSigSystematics(const std::string &systName,
			   const double up,const double down);
    void setYieldData(const int obs) {m_yieldData=obs;}
    void setYieldDataToBkg() {m_yieldData=static_cast<int>(m_yieldBg);}
    int getYieldData() const {return m_yieldData;}
    void saveYieldData() {m_yieldSaved=m_yieldData;}
    void restoreYieldData() {m_yieldData=m_yieldSaved;}
    
    // setting of signal strength to be used for all computations
    virtual void setSigStrength(const double mu);
    double getSigStrength() {return m_sigStrength;}
    
    // computation of the LLR value for the given number of observed events
    double computeLLR(const int obs) const;
    double computeLLRdata() const {return computeLLR(m_yieldData);}
    void initDistrLLR(double &llrMin,double &llrMax);
    void generateSinglePseudoExp(double &llrB,double &llrSB);
    
    // generation of nbExp pseudo-experiments to compute the LLR distributions
    // must be called before trying to compute any CLs or p-value
    virtual void generateDistrLLR(const int nbExp);
    // normalisation of distributions after generateDistrLLR
    void endDistrLLR(const int nbExp);
    
    // computation of the p-value for the given number of observed events
    // the LLR distributions must have been generated before
    double pValue(const int obs) const;

    // computation of the p-value for the observed number of events
    // the LLR distributions must have been generated before
    virtual double pValueData() const;
    
    // computation of the CLs value for the given number of observed events
    // the LLR distributions must have been generated before
    double computeCLs(const int obs) const;
    
    // find maxmimal number of observed events for a 95% CL exclusion
    // the LLR distributions must have been generated before
    int findObsExclusion() const;
    
    // generation of nbExp pseudo-experiments with the signal strength mu
    // computation of the CLs for the given number of observed events
    // (calls setSigStrength, generateDistrLLR and computeCLs)
    virtual double generateForCLs(const double mu,const int nbExp,const int type);

    // generation of yield distribution
    void generateDistrYield(const int nbExp);
    YieldWithUncert getGeneratedYieldBkg() const;
    
    // methods called for observed and expected (median, -+1 sigma, +-2 sigma) limit computation
    virtual double sigStrengthExclusion(const LimitType type,const int nbExp,double &cls,
					const double muHint=1,const OTH::MethType method=MethDichotomy);

    // computation of the distribution of the observed signal strengths if no signal exists
    // computation of the quantiles of this distribution (returns the median)
    double expectedSigStrengthExclusion(const int nbMu,const int nbExp);
    
    int generateSinglePseudoData(const double mu=0);
   
    // combination type for systematics
    inline void setCombinationType(const bool additive) {m_additiveSystComb=additive;}

    // print samples
    void printSamples() const;
    
    // get histograms of systematics distributions
    TH1 *getSigSystDistr(const std::string &systName) const;
    TH1 *getBkgSystDistr(const std::string &bkgName,const std::string &systName) const;
    
    // get histograms of yield distributions
    TH1 *getSigYieldDistr() const {return m_sigSample.getYieldHisto();}
    TH1 *getBkgYieldDistr(const std::string &bkgName) const;
    
    enum {hDistrBg,hDistrSB,hLLRb,hLLRsb,hYieldBg,nbHistos};
    // histograms (see enum above)
    TH1* getHisto(const int i) const;
    virtual TH1 *getHistoLLRsb() const;
    virtual TH1 *getHistoLLRb() const;

  private:
    Channel();
    Channel(const Channel&);
    Channel &operator=(const Channel&);
    
    int generateSinglePseudoExpBg(double &expected) const;
    double generateSingleSample(const OTH::Sample &sample,const double mu=1) const;

    std::string m_name,m_nameLaTeX; // channel name
    
    Systematics &m_syste; // list of systematic uncertainties
    PdfGenerator &m_statSampling; // implementation of statistical uncertainty variation
    std::deque<Sample> m_bgSamples; // background samples
    Sample m_sigSample; // signal sample
    double m_sigStrength; // signal strength (scale factor of signal)
    int m_yieldData,m_yieldSaved; // number of observed events in data
    
    // internal numbers
    double m_yieldBg,m_yieldSB; // expected yields in b or mu*s+b
    mutable double m_cacheLog,m_cacheYieldS; // cache for LLR computation
    mutable bool m_cached; // cache flag
    
    // distributions
    std::vector<TH1*> m_pHs; // main histos
    std::vector<std::string> m_hNames; // histo names;

    // these are for a single channel only limit
    std::map<int,double> m_muObs; // values of mu_95 for given observed events
    MuVsObs m_muVsObs; // interpolation
  };

}

#endif // OTH_CHANNEL_H
